#include <iostream>
#include <stdlib.h>
#include <hip/hip_runtime.h>
#include <hipblas/hipblas.h>

#include "gpuSymMatMult_hip.hpp"


#define IDX2C(i,j,ld)  (((j)*(ld))+(i))

gpuSymMatMult_hip::gpuSymMatMult_hip()
{
    m_isLoaded = false;
    m_rows = m_cols = 0;
    m_handle = nullptr;
    m_A = nullptr;
    m_x = nullptr;
    m_y = nullptr;
    m_z = nullptr;
}

gpuSymMatMult_hip::~gpuSymMatMult_hip()
{
    if (m_A) hipFree(m_A);
    if (m_x) hipFree(m_x);
    if (m_y) hipFree(m_y);
    if (m_z) hipFree(m_z);
    hipblasDestroy((hipblasHandle_t)m_handle);
}

int gpuSymMatMult_hip::set_matrix(int rank, size_t g_col_start, size_t n_rows, size_t n_cols, const float *A)
{
    hipError_t hipStat;
    hipblasStatus_t hipblasStat;
    size_t gpu_free, gpu_total;
    int device_id = -1;
    int ret = EXIT_SUCCESS;

    m_global_col_start = g_col_start;
    m_rows = n_rows;
    m_cols = n_cols;
    std::cout << "[" << rank << "] "
              << "A in HIP = " << A << std::endl;

    hipStat = hipSetDevice(rank);
    if (hipStat != hipSuccess) {
        std::cerr << "[" << rank << "] "
                  << "Error: hipSetDevice() has failed: " << hipStat << std::endl;
        ret = EXIT_FAILURE;
        goto out;
    }

    hipStat = hipMemGetInfo(&gpu_free, &gpu_total);
    if (hipStat != hipSuccess) {
        std::cerr << "[" << rank << "] "
                  << "Error: hipMemGetInfo() has failed: " << hipStat << std::endl;
        ret = EXIT_FAILURE;
        goto out;
    }

    std::cout << "[" << rank << "] "
              << "GPU memory: " << gpu_free << " / " << gpu_total
              << " requesting: " << m_rows*m_cols*sizeof(*m_A) << std::endl;

    hipStat = hipMalloc((void**)&m_A, m_rows*m_cols*sizeof(*m_A));
    if (hipStat != hipSuccess) {
        std::cerr << "[" << rank << "] "
                  << "Error: hipMalloc() for A has failed: " << hipStat << std::endl;
        ret = EXIT_FAILURE;
        goto out;
    }
    hipStat = hipMalloc((void**)&m_x, m_rows*sizeof(*m_x));
    if (hipStat != hipSuccess) {
        std::cerr << "[" << rank << "] "
                  << "Error: hipMalloc() for x has failed: " << hipStat << std::endl;
        ret = EXIT_FAILURE;
        goto out;
    }
    hipStat = hipMalloc((void**)&m_y, m_cols*sizeof(*m_y));
    if (hipStat != hipSuccess) {
        std::cerr << "[" << rank << "] "
                  << "Error: hipMalloc() for y has failed: " << hipStat << std::endl;
        ret = EXIT_FAILURE;
        goto out;
    }
    hipStat = hipMalloc((void**)&m_z, m_rows*sizeof(*m_z));
    if (hipStat != hipSuccess) {
        std::cerr << "[" << rank << "] "
                  << "Error: hipMalloc() for z has failed: " << hipStat << std::endl;
        ret = EXIT_FAILURE;
        goto out;
    }
    hipblasStat = hipblasCreate((hipblasHandle_t *)&m_handle);
    if (hipblasStat != HIPBLAS_STATUS_SUCCESS) {
        std::cerr << "[" << rank << "] "
                  << "Error: hipblasCreate() has failed: " << hipblasStat << std::endl;
        ret = EXIT_FAILURE;
        goto out;
    }
    hipblasStat = hipblasSetMatrix(m_rows, m_cols, sizeof(*m_A), A, m_rows, m_A, m_rows);
    if (hipblasStat != HIPBLAS_STATUS_SUCCESS) {
        std::cerr << "[" << rank << "] "
                  << "Error: hipblasSetMatrix() has failed: " << hipblasStat << std::endl;
        ret = EXIT_FAILURE;
        goto out;
    }

    m_isLoaded = true;
out:
    if (ret == EXIT_FAILURE) {
        m_rows = m_cols = 0;
        m_global_col_start = 0;
        if (m_A) hipFree(m_A);
        if (m_x) hipFree(m_x);
        if (m_y) hipFree(m_y);
        if (m_z) hipFree(m_z);
        if (m_handle) hipblasDestroy((hipblasHandle_t)m_handle);

        m_isLoaded = false;
        m_A = nullptr;
        m_x = nullptr;
        m_y = nullptr;
        m_z = nullptr;
        m_handle = nullptr;
    }
    return ret;
}

// Perform BB^T*x in two steps where B is a submatrix from A formed by B = A_{:,start_col:(end_col-1)}.
//   i)  y = B^T*x
//   ii) z = By
int gpuSymMatMult_hip::sym_sgemv_range(int rank, size_t g_start_col, size_t g_end_col, size_t n_elem, const float *x, float *ret)
{
    hipError_t hipStat;
    hipblasStatus_t hipblasStat;
    size_t start_col, end_col, n_subcols;
    const float *B;
    float alpha = 1.0, beta = 0.0;

    if (n_elem != m_rows) {
        std::cerr << "[" << rank << "] "
                  << "Error: size mismatch: "
                  << "m_rows: " << m_rows << " n_elem: " << n_elem << std::endl;
        return EXIT_FAILURE;
    }
    if (g_start_col >= g_end_col) {
        std::cerr << "[" << rank << "] "
                  << "Error: invalid range = ["
                  << g_start_col << ", " << g_end_col << ")" << std::endl;
        return EXIT_FAILURE;
    }

    memset(ret, 0, n_elem*sizeof(*ret));

    // if there is no overlap, return immediately.
    if (g_end_col <= m_global_col_start || g_start_col >= m_global_col_start+m_cols) {
        return EXIT_SUCCESS;
    }

    start_col = std::max(g_start_col, m_global_col_start) - m_global_col_start;
    end_col = std::min(g_end_col, m_global_col_start+m_cols) - m_global_col_start;

    hipStat = hipMemcpy(m_x, x, n_elem*sizeof(*x), hipMemcpyHostToDevice);
    if (hipStat != hipSuccess) {
        std::cerr << "[" << rank << "] "
                  << "Error: hipMemcpy (host to Device) failed: " << hipStat
                  << " m_x = " << m_x << " x = " << x << std::endl;
        return EXIT_FAILURE;
    }

    n_subcols = end_col - start_col;
    B = m_A + m_rows*start_col;

    hipblasStat = hipblasSgemv((hipblasHandle_t)m_handle, HIPBLAS_OP_T, m_rows, n_subcols, &alpha, B, m_rows, m_x, 1, &beta, m_y, 1);
    if (hipblasStat != HIPBLAS_STATUS_SUCCESS) {
        std::cerr << "[" << rank << "] "
                  << "Error: hipblasSgemv failed (HIPBLAS_OP_T): " << hipblasStat << std::endl;
        return EXIT_FAILURE;
    }
    hipblasStat = hipblasSgemv((hipblasHandle_t)m_handle, HIPBLAS_OP_N, m_rows, n_subcols, &alpha, B, m_rows, m_y, 1, &beta, m_z, 1);
    if (hipblasStat != HIPBLAS_STATUS_SUCCESS) {
        std::cerr << "[" << rank << "] "
                  << "Error: hipblasSgemv failed (HIPBLAS_OP_N): " << hipblasStat << std::endl;
        return EXIT_FAILURE;
    }
    hipStat = hipMemcpy(ret, m_z, n_elem*sizeof(*ret), hipMemcpyDeviceToHost);
    if (hipStat != hipSuccess) {
        std::cerr << "[" << rank << "] "
                  << "Error: hipMemcpy (device to host) failed: " << hipStat << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

// Perform AA^T*x in two steps.
//   i)  y = A^T*x
//   ii) z = Ay
int gpuSymMatMult_hip::sym_sgemv(int rank, size_t n_elem, const float *x, float *ret)
{
    hipError_t hipStat;
    hipblasStatus_t hipblasStat;
    float alpha = 1.0, beta = 0.0;

    if (n_elem != m_rows) {
        std::cerr << "[" << rank << "] "
                  << "Error: size mismatch: "
                  << "m_rows: " << m_rows << " n_elem: " << n_elem << std::endl;
        return EXIT_FAILURE;
    }

    // Initialize the result buffer to zero
    memset(ret, 0, n_elem * sizeof(*ret));

    // Copy the input vector to device memory
    hipStat = hipMemcpy(m_x, x, n_elem * sizeof(*x), hipMemcpyHostToDevice);
    if (hipStat != hipSuccess) {
        std::cerr << "[" << rank << "] "
                  << "Error: hipMemcpy (host to Device) failed: " << hipStat
                  << " m_x = " << m_x << " x = " << x << std::endl;
        return EXIT_FAILURE;
    }

    // First stage of the matrix-vector multiplication: y = A^T * x
    hipblasStat = hipblasSgemv((hipblasHandle_t)m_handle, HIPBLAS_OP_T, m_rows, m_cols, &alpha, m_A, m_rows, m_x, 1, &beta, m_y, 1);
    if (hipblasStat != HIPBLAS_STATUS_SUCCESS) {
        std::cerr << "[" << rank << "] "
                  << "Error: hipblasSgemv failed (HIPBLAS_OP_T): " << hipblasStat << std::endl;
        return EXIT_FAILURE;
    }

    // Second stage of the matrix-vector multiplication: z = A * y
    hipblasStat = hipblasSgemv((hipblasHandle_t)m_handle, HIPBLAS_OP_N, m_rows, m_cols, &alpha, m_A, m_rows, m_y, 1, &beta, m_z, 1);
    if (hipblasStat != HIPBLAS_STATUS_SUCCESS) {
        std::cerr << "[" << rank << "] "
                  << "Error: hipblasSgemv failed (HIPBLAS_OP_N): " << hipblasStat << std::endl;
        return EXIT_FAILURE;
    }

    // Copy the result vector back to host memory
    hipStat = hipMemcpy(ret, m_z, n_elem * sizeof(*ret), hipMemcpyDeviceToHost);
    if (hipStat != hipSuccess) {
        std::cerr << "[" << rank << "] "
                  << "Error: hipMemcpy (device to host) failed: " << hipStat << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
