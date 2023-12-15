# gwPheWAS-Summit

SAIGE is a tool that allows users to run GWAS for a particular cohort. 
Our goal is to run SAIGE on the Oakridge Summit infrastructure to take advantage of the GPU power available on the PowerPC9 computers.
SAIGE is run on standard CPU infrastrucutre and does not take advantage of GPUs. We have modifed the SAIGE code for step1 (SAIGE_fit_glmm), to use GPUs instead of CPUs. In order to do so, we first calculate the GRM and then make the iterative calculations.
The modified code is available in this repo. For step2, the SAIGE developer (Wei) provided us with an improved version for step 2 which is faster due to how the inputs are read and stored.

In the following sections, we will outline how SAIGE was run on Summit once the code was optimized to run on GPUs.

## Files needed to run SAIGE on Summit
1. HARE files
2. PCs files
3. Filtered Genotype BED files
4. Imputed BGEN files with their indexes
5. Phenotype files
6. Covariate files which include Age, Sex and any other covariates needed

## Steps to run SAIGE on Summit
1. Load the modules

```
module load gcc/10.2.0
module load bzip2
module load openblas
module load hdf5
module load python/2.7.15-anaconda2-5.3.0
module load cmake
module load mercurial
module load openblas
module load  cuda/11.1.1
R_LIB=/gpfs/alpine/proj-shared/med112/task0101113/tools/R-2/R-4.0.3/library/
export PATH=$PATH:/gpfs/alpine/proj-shared/med112/task0101113/tools/R-2/R-4.0.3/bin
```

2. Get list of phenotype ids to run through SAIGE

`/gpfs/alpine/med112/proj-shared/data/phenotype_data/processed/phenotype_list/final_selection_phenotype_variable_metadata_qc.csv`

3. Create the submit files

```
for i in `cut -f2 /gpfs/alpine/med112/proj-shared/data/phenotype_data/processed/phenotype_list/final_selection_phenotype_variable_metadata_qc.csv | tail -n +2`; do echo $i; if [ ! -d /ccs/home/arodriguez/med112/task0101113/output/FullRun/$i ]; then python ./get_phetype_v1.py --phecode $i; fi; done
```

4. Create the task list for each population group for step1 and step2

```
Rscript merge_step1_rds_v1.r ASN
Rscript merge_step1_rds_v1.r AFR
Rscript merge_step1_rds_v1.r EUR
Rscript merge_step1_rds_v1.r HIS
Rscript ./merge_step2_rds_v1.r ASN
Rscript ./merge_step2_rds_v1.r AFR
Rscript ./merge_step2_rds_v1.r EUR
Rscript ./merge_step2_rds_v1.r HIS
```

5. Submit step1  task list to a large pool of nodes and wait for completion

```
$ more submit.step1.all.AFR.nvlm.txt
#!/bin/bash
#BSUB -nnodes 1300
#BSUB -W 24:00
#BSUB -q batch
#BSUB -P MED112
#BSUB -o all.s1.run3.AFR.stdout
#BSUB -e all.s1.run3.AFR.stderr
#BSUB -J all.s1.AFR
#BSUB -alloc_flags "nvme"

module load gcc/10.2.0
module load cuda/11.1.1
module load openblas
module load bzip2
module load hdf5
jsrun -n1300 -a1 -g6 -r1 --smpiargs="-disable_gpu_hooks" -E NVBLAS_CONFIG_FILE=~/gemm/nvblas.conf -E LD_PRELOAD=$OLCF_CUDA_ROOT/lib64/libnvblas.so /ccs/home/arodriguez/med112/task0101113/tools/R-2/R-4.0.3/bin/Rscript /gpfs/alpine/proj-shared/med112/task0101113/batch/pre-gwas/scripts_for_run/submit_step1_with_nvme.r /gpfs/alpine/proj-shared/med112/task0101113/output/FullRun/SUBMIT/complete_step1.tasks.AFR.rds all AFR


```

```
bsub submit.step1.all.AFR.nvlm.txt
```

7. QC check on results for step1 and resubmit if necessary
8. Submit step2 task list to a large pool of nodes and wait for completion

```
bsub submit.step2.all.AFR.nvlm.txt
```

10. QC check on results for step2 and resubmit if necessary

## Post-SAIGE Meta analysis
