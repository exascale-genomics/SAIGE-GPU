################################################################################

#finemap_SuSiE.R

#This script is designed to execute SuSiE fine-mapping given a fixed region 
#when called by finemap_SuSiE.py after temp file generation. It is set up for 
#in-sample fine mapping only! The code saves its output as a .rds file that can 
#later be used for plotting or for downstream analysis.

################################################################################

# 0) Load Needed Libraries ====

library(coloc)
library(stringr)

# 1) Read in Needed Arguments from Command ====

#Need to read in the following: 
  # 1) file_prefix (TRAIT.CHR.START.END.ANC)
  # 2) random seed
  # 3) confidence level
  # 4) summary statistics processed by MungeSumStats
  # 5) list of snps from .map file
  # 6) ld matrix file in .ld or similar format
  # 7) file of case/control counts
  # 8) output file folder

args <- commandArgs(trailingOnly = TRUE);
if (length(args) == 8) {
  #Print confirmation statement
  print("Necessary Parameters Input")
  #Set variables from the arguments in this case
  file_prefix <- args[1]
  random_seed <- as.integer(args[2])
  confidence <- as.numeric(args[3])
  summary_stats_file <- args[4]
  ld_ref_snps <- args[5]
  ld_ref_matrix <- args[6]
  count_loc <- args[7]
  out_folder <- args[8]
} else {
  print("Usage: %> Rscript finemap_SuSiE.R file_prefix random_seed confidence_level summary_stats.txt snp_list.map ld_matrix.ld count_loc output_folder");
  quit(save="no");
}

# 2) Prep for executing SuSiE ====

#Define out_loc
if(substr(out_folder, nchar(out_folder), nchar(out_folder)) != "/"){
  out_loc <- paste0(out_folder, "/")
} else {
  out_loc <- out_folder
}
if (file.exists(out_loc) == FALSE) {
  print("ERROR: Invalid output folder given.")
  quit(save="no");
}
rds_out_loc <- paste0(out_loc, file_prefix, ".rds")

#Identify trait and ancestry
trait <- unlist(strsplit(file_prefix, split = "[.]"))[1]
ancestry <- unlist(strsplit(file_prefix, split = "[.]"))[length(unlist(strsplit(file_prefix, split = "[.]")))]
#Read in trait registry and filter for our trait
trait_register <- read.csv(count_loc, header = TRUE, sep = "\t", fill = TRUE)
trait_register <- trait_register[which(trait_register$Trait == trait),]
#Determine Trait Type
trait_type <- trait_register[1,"Trait_Type"]
#Verify that we have at least some individuals for the ancestry
if(is.na(trait_register[1,paste0("num_samples.", ancestry)])){
  print(paste0("ERROR: Attempted to map a trait with no samples for trait/locus/anc: ", file_prefix))
  quit(save="no");
}

#Read in summary_statistics and snps list from map file
summary_stats_raw <- read.table(summary_stats_file, header = TRUE)
filt_snps <- read.table(ld_ref_snps, header = FALSE)
snps_list <- as.character(filt_snps[,2])
#Filter summary_statistics for the SNPs we care about
summary_stats_filt <- summary_stats_raw[which(summary_stats_raw$ID %in% snps_list),]

#Check for cc and make lists to execute coloc
#Add sample_size and proportion cases fields on for any cases and controls
if(trait_type == "binary") {
  #Check that we have both cases and controls in this case
  if(is.na(trait_register[1,paste0("num_cases.", ancestry)]) || is.na(trait_register[1,paste0("num_controls.", ancestry)])){
    print(paste0("ERROR: Attempted to map a binary trait with either no cases or no controls for trait/locus/anc: ", file_prefix))
    quit(save="no");
  } else {
  #Add prop_cases
  summary_stats_filt <- cbind.data.frame(summary_stats_filt, prop_cases=(trait_register[1,paste0("num_cases.", ancestry)]/(trait_register[1,paste0("num_cases.", ancestry)] + trait_register[1,paste0("num_controls.", ancestry)])))
  }
}

#These summary statistics are ordered in positional order, the same order that the SNPs
#from the MAP file are in.  This makes sense because I used the summary statistics SNPs
#to filter the plink files

#Read in the LD-matrix
ld_ref_raw <- read.table(ld_ref_matrix, header = FALSE)
ld_ref_raw <- as.matrix(ld_ref_raw)
colnames(ld_ref_raw) <- snps_list
rownames(ld_ref_raw) <- snps_list

#Check allele frequencies. Increase any 0% up to to a minimum value
#Get minimum values of the allele frequency columns
min_FRQ <- min(min(summary_stats_filt[which(summary_stats_filt$maf > 0), "maf"]), 1/(2*max(summary_stats_filt$num_samples)))
#Over write values
summary_stats_filt$maf <- ifelse(summary_stats_filt$maf < min_FRQ, min_FRQ, summary_stats_filt$maf)
#Also overwrite max values
summary_stats_filt$maf <- ifelse(summary_stats_filt$maf > 1-min_FRQ, 1-min_FRQ, summary_stats_filt$maf)

#Set variables for SuSiE
z <- summary_stats_filt$z
n <- max(summary_stats_filt$num_samples)
R <- ld_ref_raw
beta <- summary_stats_filt$beta
varbeta <- (summary_stats_filt$sebeta)^2

#Check for cc and make lists to create coloc datasets
print("Creating coloc datasets")
if(trait_type == "cc") {
  coloc_data <- list(z=z, snp=summary_stats_filt$ID, position=summary_stats_filt$pos37,
                     type="cc", N=n, MAF=summary_stats_filt$maf, s=summary_stats_filt$prop_cases, LD=R, beta=beta, varbeta=varbeta) 
} else {
  coloc_data <- list(z=z, snp=summary_stats_filt$ID, position=summary_stats_filt$pos37,
                     type="quant", N=n, MAF=summary_stats_filt$maf, LD=R, beta=beta, varbeta=varbeta)
}
#Check dataset
if (is.null(check_dataset(coloc_data))) {
  print("Dataset check passed.")
} else{
  print("Dataset check failed.")
  check_dataset(coloc_data)
}

# 3) Set error checking function ====

check <- function(expression){
  
  withCallingHandlers(expression,
                      
                      warning = function(w){
                      },
                      error = function(e){
                        message(paste0("ERROR: ", file_prefix, " is unmappable due to negative estimated residual variance."))
                        #Write an empty output file and quit
                        fitted_rss1 = NULL
                        saveRDS(fitted_rss1, file = rds_out_loc)
                        quit(save="no");
                      },
                      finally = {
                      })
}

# 4) Execute SuSiE ====

#Execute SuSiE
set.seed(random_seed) #Nothing up until here has been random
fitted_rss1 = check(runsusie(coloc_data, n=n, min_abs_corr=0, coverage=confidence, estimate_residual_variance = TRUE, L = 5, maxit=10000, repeat_until_convergence=FALSE))

#Verify convergence
if (fitted_rss1$converged == FALSE) {
  print(paste0("ERROR: ", file_prefix, " did not converge in < 10,000 iterations. Outputting empty rds file."))
  #Write an empty output file and quit
  fitted_rss1 = NULL
  saveRDS(fitted_rss1, file = rds_out_loc)
  quit(save="no");
}

## Calculate residuals ##
#Get credible sets
cred_sets=names(fitted_rss1$sets$cs)
cred_sets=vapply(cred_sets, str_replace, character(1), pattern="L", replacement="")
cred_sets=vapply(cred_sets, as.numeric, 0)
cred_sets=sort(cred_sets)

#Calculate residuals
for (i in cred_sets) {
  cred_set <- names(fitted_rss1$sets$cs[[paste0("L", i)]])
  #Calculate preceding residuals
  if (i > 2) {
    prediction = colSums(fitted_rss1$mu[seq(1,i - 1, 1),] * fitted_rss1$alpha[seq(1,i - 1, 1),])
    test = matrix(sqrt(n - 1) * z) - ((n - 1) * R %*% prediction)
    test2 = (test)/sqrt(fitted_rss1$sigma2 * (n-1) * diag(R))
  } else if(i == 2) {
    prediction = fitted_rss1$mu[i-1,] * fitted_rss1$alpha[i-1,]
    test = matrix(sqrt(n - 1) * z) - (n - 1) * R %*% prediction
    test2 = (test)/sqrt(fitted_rss1$sigma2 * (n-1) * diag(R))
  } else{
    test = matrix(sqrt(n - 1) * z)
    test2 = z
  }
  neglogp_precede = (-log(2)-pt(-abs(test2), df = n - 1, lower.tail = TRUE, log.p = TRUE))/log(10)
  #Calculate absolute residuals
  a = fitted_rss1$mu[-i,]
  b = fitted_rss1$alpha[-i,]
  if (is.null(ncol(a)) && is.null(ncol(b))) {
    prediction = colSums(t(a) * t(b))
  } else{
    prediction = colSums(a * b)
  }
  test = matrix(sqrt(n - 1) * z) - ((n - 1) * R %*% prediction)
  test2 = (test)/sqrt(fitted_rss1$sigma2 * (n-1) * diag(R))
  neglogp_absolute = (-log(2)-pt(-abs(test2), df = n - 1, lower.tail = TRUE, log.p = TRUE))/log(10)
  #Check if i == 1 to process info to data frames correctly
  if(i == 1){
    neglogp_precede_df = as.data.frame(neglogp_precede)
    neglogp_absolute_df = as.data.frame(neglogp_absolute)
  } else{
    neglogp_precede_df = cbind.data.frame(neglogp_precede_df, neglogp_precede)
    neglogp_absolute_df = cbind.data.frame(neglogp_absolute_df, neglogp_absolute)
  }
}

#Check if the residual dfs were actually populated (They won't be if no credible sets were identifiable)
if (exists("neglogp_precede_df") && exists("neglogp_absolute_df")) {
  #Set column names
  colnames(neglogp_precede_df) <- cred_sets
  colnames(neglogp_absolute_df) <- cred_sets
} else {
  print(paste0("WARNING: No credible sets detected for ", file_prefix))
  neglogp_precede_df <- NULL
  neglogp_absolute_df <- NULL
}

#Append residuals to SuSiE object
fitted_rss1$neglogp_precede <- neglogp_precede_df
fitted_rss1$neglogp_absolute <- neglogp_absolute_df

#Append chr, bp, and neglogp_gwas to SuSiE object
fitted_rss1$chr <- median(summary_stats_filt$chrom)
fitted_rss1$bp <- summary_stats_filt$pos37
fitted_rss1$bp38 <- summary_stats_filt$pos38
fitted_rss1$neglogp_gwas <- summary_stats_filt$X.logP

#Save susie_object
saveRDS(fitted_rss1, file = rds_out_loc)

