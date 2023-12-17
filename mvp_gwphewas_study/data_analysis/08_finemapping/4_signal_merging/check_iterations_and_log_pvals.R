################################################################################

#check_iterations_and_log_pvals.R

#This script is designed to open up a SuSiE results RDS file check and report 
#number of iterations run and recalculate any negative log gwas/residual 
#p-values.  

################################################################################

# 0) Load Needed Libraries ====

# 1) Read in Needed Arguments from Command ====

#Need to read in the following: 
  # 1) File of SuSiE fine-mapping results

args <- commandArgs(trailingOnly = TRUE);
if (length(args) == 1) {
  #Print confirmation output
  print("Necessary Parameters Input")
  #Set variables from the arguments in this case
  susie_results <- args[1]
} else {
  print("Usage: %> Rscript check_iterations_and_log_pvals.R susie_results");
  quit(save="no");
}

# 2) Read in file ====

#Read in rds file
fitted_rss1 <- readRDS(susie_results)

#Add a catch that will keep us from checking the empty rds objects
if (is.null(fitted_rss1)) {
  quit(save="no")
}

#Check for an infinite p-value 
if(max(fitted_rss1$neglogp_gwas) == Inf || suppressWarnings(max(fitted_rss1$neglogp_absolute)) == Inf || suppressWarnings(max(fitted_rss1$neglogp_precede)) == Inf){
  print(paste0("WARNING: Need to recalculate negative log p-values for ", susie_results))
}

#Report out the number of iterations
print(paste0("NOTE: The number of iterations ", susie_results, " ran is ", fitted_rss1$niter))



