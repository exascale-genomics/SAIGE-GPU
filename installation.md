# Installation

## Dependencies

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
R
>install.packages(c())

```

## Install and quick test on HPC from Source Code

```
git clone git://path.to.git
cd SAIGE-GPU
export PATH=${PATH}:${HOME}/.local/summit/anaconda2/5.3.0/2.7/bin
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

## Run on HPC Systems Using Singularity Container

It is not necessary to build SAIGE-GPU from source code. Instead you can pull the latest Docker container and run a small example from the SAIGE package itself:

```
singularity build saige-doe.sif docker://tnnandi/saige-doe:2

singularity exec --bind /SAIGE_container/SAIGE-DOE/extdata saige_1.1.9.sif /SAIGE_container/SAIGE-DOE/extdata/step1_fitNULLGLMM.R \
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
