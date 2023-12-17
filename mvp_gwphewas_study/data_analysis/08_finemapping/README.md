# Fine-Mapping Analysis Scripts

This directory contains all the scripts used to complete the fine-mapping analysis in Verma et al. The scripts were run on the Oak Ridge Leadership Computing Facility's Andes supercomputer. The scripts are broken down into the following five steps:

## 1) Locus Defining ##
These scripts define the loci used for fine-mapping by tiling the genome and merging together adjacent tiles with significant variants. The key shell script in this directory, **define_sig_loci.sh**, will submit the locus-defining process and merge the loci together into a comma-delimited list upon completion.

## 2) Matrix Making ##
These scripts create the LD matrices specific to each population group, trait, and locus which are needed to accurately fine-map the locus-trait pairs in each group. The key shell script in this directory, **make_matrices.sh**, will submit the matrix-making process, ultimately outputting a matrix and map file for each locus.

## 3) Fine-Mapping ##
These scripts execute the SuSiE-based fine-mapping of each population-trait-locus combination using the in-sample LD matrices. The key shell script in this directory, **map_all_loci.sh**, will submit the fine-mapping process, ultimately outputting a ".rds", R-object, file for each successful mapping. During testing, single-trait CAFEH was also considered in lieu SuSiE, but as the results were nearly identical, SuSiE was used due to ease of processing the results. The corresponding CAFEH mapping scripts are available in this directory as are the plotting scripts needed to plot the results of individual loci using both methods.  

## 4) Signal-Merging ##
These scripts combine the results of all the individual SuSiE output ".rds" files into two tab-delimited files by using the Jaccard-similarity method to merge the signals across population groups. The key shell script in this directory, **merge_all_loci.sh**, will submit the merging process. 

## 5) Synthesis and Analysis ##
These scripts are the downstream analysis scripts that analyze and plot all the summary-level results. These scripts created the fine-mapping figures referenced in Verma et al. The figures and the scripts that generated them are as follows:

Fig. 3 - **make_summary_plots.R**  
fig. S3 - Sankey plot is from **make_sankey_plot.R** and the bar plot is from **make_summary_plots.R**. The flow diagram was manually created.  
fig. S4 - **make_summary_plots.R**  
fig. S5 - **make_summary_plots.R**  
fig. S6 - **make_summary_plots.R**  
fig. S7 - **make_summary_plots.R**  

For questions about how individual scripts function, please see individual file headers and help messages.
