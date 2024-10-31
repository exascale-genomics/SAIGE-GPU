# Installation

## Install from source on OLCF Summit

### Dependencies

The dependencies listed below are modules used to install SAIGE from source code in the [ORNL OLCF Summit HPC](https://docs.olcf.ornl.gov/systems/summit_user_guide.html).

```
module load cuda/11.0.2
module load python/2.7.15-anaconda2-5.3.0
module load r/4.0.5
module load cmake
module load openblas
module load spectrum-mpi/10.4.0.3-20210112
module load libxml2
module load gcc/11.1.0
```

There are several R libraries needed to be installed:

```
mkdir saige
cd saige

R_LIB=/path/to/your/R_lib
Rscript -e "install.packages(c(\
    'R.utils', \
    'Rcpp', \
    'RcppParallel', \
    'data.table', \
    'RcppEigen', \
    'Matrix', \
    'BH', \
    'optparse', \
    'SPAtest',
    'SKAT', \
    'RcppArmadillo', \
    'qlcMatrix', \
    'RhpcBLASctl'), \
  lib='/path/to/your/R_lib', \
  repos='https://cran.rstudio.com')"

```

### Install and quick test on HPC from Source Code

You can install SAIGE from source code. We have merged SAIGE-GPU with the latest version of [SAIGE](https://saigegit.github.io/SAIGE-doc/). This GitHub contains the latest merged version and can be installed following these instructions. You will need the above mentioned R libraries as they are not contained in the branch.

```
export PATH=${PATH}:${HOME}/.local/summit/anaconda2/5.3.0/2.7/bin
git clone https://github.com/exascale-genomics/SAIGE-GPU.git
cd SAIGE-GPU/src
R_LIB=/path/to/your/R_lib
R CMD INSTALL SAIGE --library=$R_LIB
```

You can submit a quick test job on a GPU machine. For the test you only need 1 GPU available of minimum 16 Gigabytes of memory.

The following file was used to submit to the HPC:

```
#!/bin/bash
#BSUB -nnodes 1
#BSUB -W 4:00
#BSUB -q batch-hm
#BSUB -P MED112
#BSUB -o s1.synt_data.gpu.haoyu.cate_var.stdout
#BSUB -e s1.synt_data.gpu.haoyu.cate_var.stderr
#BSUB -J s1.synt_data
#BSUB -alloc_flags nvme

module load python/2.7.15-anaconda2-5.3.0
module load r/4.0.5
module load gcc/9.1.0
module load cmake
module load openblas
module load cuda/11.0.2
module load spectrum-mpi/10.4.0.3-20210112

path_to_saige="/gpfs/alpine/proj-shared/med112/task0101113/tools/saige_20220326/SAIGE-DOE/extdata/"
TMPDIR="/path/to/output"

jsrun -n6 -a1 -c1 -g1 Rscript $path_to_saige/step1_fitNULLGLMM.R \
   --plinkFile=$path_to_saige/100k_arrays \
   --phenoFile=$path_to_saige/phenotypes.tsv \
   --invNormalize=FALSE \
   --phenoCol=case \
   --covarColList=sex_at_birth_Male,age,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
   --qCovarColList=sex_at_birth_Male \
   --sampleIDColinphenoFile=IID \
   --traitType=binary \
   --outputPrefix=$TMPDIR/GPU_step1_output \
   --minMAFforGRM 0.01 \
   --LOCO FALSE \
   --IsOverwriteVarianceRatioFile=TRUE \
   --nThreads=1; gsutil -m cp $TMPDIR/GPU_step1_output* $OUT_DIR
```

## Install via Conda on ALCF Polaris

### Prerequisites

1. **Access to Polaris**: Ensure you have access to the Polaris system.
2. **Conda**: Conda should be available on Polaris.

```bash
module use /soft/modulefiles/
module load conda/2024-04-29
```

### Installation Steps

#### Step 1: Create directory and clone the SAIGE-DOE GitHub repository

Create and activate a new Conda environment for SAIGE:

```bash
mkdir SAIGE-GPU
cd SAIGE-GPU
git clone https://github.com/exascale-genomics/SAIGE-DOE.git
cd SAIGE-DOE/
git pull origin SAIGE-step1GPU-step2Wei-openmpi
```

#### Step 2: Create a Conda Environment

Create and activate a new Conda environment for SAIGE. We will be using the existing YML file from the repository.
In this example I am naming my environment `RSAIGE_GPU_V2`, you can replace with your desired name. 
In addition, due to the size of the conda environment and packages, I will be installing the environment in a different mount where I have more space. You can create your environment either in your home directory if you have enough space, or your project space.

```bash
conda env create  --file=./conda_env/environment-RSAIGE.yml -p /grand/projects/GeomicVar/rodriguez/conda_envs/RSAIGE_GPU_V2
conda activate /grand/projects/GeomicVar/rodriguez/conda_envs/RSAIGE_GPU_V2
```

#### Step 3: Install Dependencies
SAIGE requires several packages and libraries. Install these with the following commands:


```bash
pip3 install cget click
conda install cuda -c nvidia/label/cuda-11.4.3
```

You need to install openMPI version `4.1.5`. This is difficult to perform within Conda, so we will install separately, but then include it in our Conda environment:

```bash
cd /grand/projects/GeomicVar/rodriguez/conda_envs/pkgs
wget https://download.open-mpi.org/release/open-mpi/v4.1/openmpi-4.1.5.tar.gz tar -xzf openmpi-4.1.5.tar.gz
cd openmpi-4.1.5
./configure --prefix=/grand/projects/GeomicVar/rodriguez/conda_envs/RSAIGE_GPU_V2/opt/openmpi
make -j4
make install
export PATH=/grand/projects/GeomicVar/rodriguez/conda_envs/RSAIGE_GPU_V2/opt/openmpi/bin:$PATH
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/grand/projects/GeomicVar/rodriguez/conda_envs/RSAIGE_GPU_V2/opt/openmpi/lib

mkdir -p $CONDA_PREFIX/etc/conda/activate.d
mkdir -p $CONDA_PREFIX/etc/conda/deactivate.d
echo 'export PATH=$CONDA_PREFIX/opt/openmpi/bin:$PATH' > $CONDA_PREFIX/etc/conda/activate.d/env_vars.sh
echo 'export LD_LIBRARY_PATH=$CONDA_PREFIX/opt/openmpi/lib:$LD_LIBRARY_PATH' >> $CONDA_PREFIX/etc/conda/activate.d/env_vars.sh
```

Install other required libraries, such as pbdMPI, savvy, superlu:

```bash
Rscript -e 'install.packages("pbdMPI", repos=c("https://cloud.r-project.org"))'
conda install -c conda-forge -c bioconda savvy
conda install conda-forge::superlu
```

#### Step 4: Compile SAIGE
To compile SAIGE-GPU, first clean any previous builds and then run make:

```bash
cd ~/SAIGE-GPU
R CMD INSTALL SAIGE-DOE
```

If you encounter linking errors, ensure that the PKG_LIBS line in the Makevars file correctly references the MPI library.

#### Step 5: Verify Installation
Check if the installation was successful by running the following commands where the output should be the list of parameter options:

```bash
path_to_saige=~/SAIGE-GPU_3/SAIGE-DOE
Rscript $path_to_saige/extdata/step1_fitNULLGLMM.R --help
```

If the help information is displayed for each command, the installation is complete.

You can also run a test with the provided test input files. You can replace `mpirun -n 4` in the command below with the appropriate number of GPUs you have available.

```bash
# ask for a node
qsub -A geomicVar -I -l select=1 -l walltime=1:00:00 -l filesystems=home:eagle -q debug

# once the node is provided
module use /soft/modulefiles/
module load conda
conda activate /grand/projects/GeomicVar/rodriguez/conda_envs/RSAIGE_GPU_V2

path_to_saige=~/SAIGE-GPU_3/SAIGE-DOE
mpirun -n 4 Rscript $path_to_saige/extdata/step1_fitNULLGLMM.R \
--plinkFile=$path_to_saige/extdata/input/plinkforGRM_1000samples_10kMarkers \
--phenoFile=$path_to_saige/extdata/input/pheno_1000samples.txt \
--invNormalize=FALSE \
--phenoCol=y \
--covarColList=x1,x2 \
--sampleIDColinphenoFile=IID \
--traitType=binary \
--outputPrefix=./GPU_step1_output \
--minMAFforGRM 0.01 \
--LOCO  \
 --IsOverwriteVarianceRatioFile=TRUE \
--nThreads=1
```

You should see a succesful run where all GPUs are used. The log should provide the IDs of the GPUs used.


## Run on HPC Systems Using Singularity Container

It is not necessary to build SAIGE-GPU from source code. Instead you can pull the latest Docker container and run a small example from the SAIGE package itself:

```
singularity build saige-doe.sif docker://tnnandi/saige-doe:2

mpirun -n <number of gpus>
   singularityexec --bind /SAIGE_container/SAIGE-DOE/extdata saige_1.1.9.sif /SAIGE_container/SAIGE-DOE/extdata/step1_fitNULLGLMM.R \
   --plinkFile=$ARRAYS_DIR/100k_arrays \
   --phenoFile=$INPUT_DIR/phenotypes.tsv \
   --invNormalize=FALSE \
   --phenoCol=case \
   --covarColList=sex_at_birth_Male,age,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
   --qCovarColList=sex_at_birth_Male \
   --sampleIDColinphenoFile=IID \
   --traitType=binary \
   --outputPrefix=$TMPDIR/GPU_step1_output \
   --minMAFforGRM 0.01 \
   --LOCO FALSE \
   --IsOverwriteVarianceRatioFile=TRUE \
   --nThreads=1; gsutil -m cp $TMPDIR/GPU_step1_output* $OUT_DIR/
```



## Install on Google Cloud Platform

I will place instructions on how to run on GCP shortly...
