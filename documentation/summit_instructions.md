---
layout: page
title: Install on HPC
---

We will be installing SAIGE-GPU from the GitHub repository.

## Installation on OLCF ORNL Summit

The OLCF ORNL [Summit](https://docs.olcf.ornl.gov/systems/summit_user_guide.html) HPC was used to develop the SAIGE-GPU software.

```
module load cuda/11.0.2
module load python/2.7.15-anaconda2-5.3.0
module load r/4.0.5
module load cmake
module load openblas
module load spectrum-mpi/10.4.0.3-20210112
module load libxml2
module load gcc/11.1.0

export PATH=${PATH}:${HOME}/.local/summit/anaconda2/5.3.0/2.7/bin
R_LIB=/ccs/home/arodriguez/med112/task0101113/tools/saige_20220326/SAIGE-DOE/extdata/R_lib
R CMD INSTALL SAIGE --library=$R_LIB
```

## Test on Summit

We used a synthetic genotype dataset representing the AFR population based from the 1000 Genome project.
The synthetic data used had 150,000 individuals. The genotype file has 100,000 variants present.
This data was built following the instructions [here](https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/COXHAP).

The data was downloaded and pruned to be able to build the BED and BGEN files needed to run SAIGE-GPU.

We used 6 GPUs to distribute the full GRM evenly. Step 1 completed within 8 minutes.
Similarly, we submitted the same dataset through the CPU version of SAIGE and the job completed in 47 minutes.

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

#module load libxml2
#module load gcc/11.1.0

module load python/2.7.15-anaconda2-5.3.0
module load r/4.0.5
module load gcc/9.1.0
module load cmake
module load openblas
module load cuda/11.0.2
module load spectrum-mpi/10.4.0.3-20210112

#path_to_saige="/ccs/home/arodriguez/med112/task0101113/tools/saige_20230621/SAIGE/extdata/"
path_to_saige="/gpfs/alpine/proj-shared/med112/task0101113/tools/saige_20220326/SAIGE-DOE/extdata/"

jsrun -n6 -a1 -c1 -g1 Rscript $path_to_saige/step1_fitNULLGLMM.R --plinkFile=/gpfs/alpine/proj-shared/med112/task0101113/synthetic_data/Haoyu/AFR_merge_chrs.100k.forCate_vr --phenoFile=/gpfs/alpine/proj-shared/med112/task0101113/synthetic_data/Haoyu/AFR_pheno_rho_3_GA_5 --phenoCol=pheno_1 --sampleIDColinphenoFile=IID --traitType=quantitative --outputPrefix=/gpfs/alpine/proj-shared/med112/task0101113/synthetic_data/Haoyu/output/synth_quant.AFR.100k.cate_var --nThreads=1 --IsOverwriteVarianceRatioFile=TRUE --isCovariateOffset=FALSE --minMAFforGRM 0.01 --maxiterPCG 500 --maxiter 25 --LOCO TRUE --invNormalize=TRUE > /gpfs/alpine/proj-shared/med112/task0101113/synthetic_data/Haoyu/output/synth_quant.AFR.cate_var.100k.log
```
