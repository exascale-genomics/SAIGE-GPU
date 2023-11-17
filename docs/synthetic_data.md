# Synthetic Data Test

We used a synthetic genotype dataset representing the AFR population based from the 1000 Genome project.
The synthetic data used had 150,000 individuals. The genotype file has 100,000 variants present.
This data was built following the instructions [here](https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/COXHAP).

The data was downloaded and pruned to be able to build the BED and BGEN files needed to run SAIGE-GPU.

We used 6 GPUs to distribute the full GRM evenly. Step 1 completed within 8 minutes.
Similarly, we submitted the same dataset through the CPU version of SAIGE and the job completed in 47 minutes.
