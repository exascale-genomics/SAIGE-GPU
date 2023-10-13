---
layout: page
title: Install on HPC
---

We will be installing SAIGE-GPU from the GitHub repository.

## Installation on OLCF ORNL Summit

The [OLCF ORNL Summit HPC](https://docs.olcf.ornl.gov/systems/summit_user_guide.html) was used to develope the SAIGE-GPU software.

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
