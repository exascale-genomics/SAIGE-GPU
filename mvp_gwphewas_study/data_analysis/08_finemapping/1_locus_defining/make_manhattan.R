################################################################################

#make_manhattan.R

#The purpose of this script is to make a Manhattan plot for either a given region
#or genome-wide.

################################################################################

# 0) Read in Arguments from command ====

#Need to read in inp_dir and out_dir
args <- commandArgs(trailingOnly = TRUE);
if (length(args) == 2) {
  
  #Print confirmation output
  print("Necessary Parameters Input")
  
  #Threshold for Expression (Because of the conversion to TPM the expression levels should be independent of the sequencing depth)
  inp_file <- as.character(args[1])
  out_file <- as.character(args[2])
  chromo <- FALSE
  rstart <- FALSE
  rend <- FALSE
  
} else if (length(args) == 5) {
  
  #Print confirmation output
  print("Necessary Parameters Input")
  
  #Threshold for Expression (Because of the conversion to TPM the expression levels should be independent of the sequencing depth)
  inp_file <- as.character(args[1])
  out_file <- as.character(args[2])
  chromo <- as.numeric(args[3])
  rstart <- as.numeric(args[4])
  rend <- as.numeric(args[5])
  
}else {
  print("Usage: %> Rscript make_manhattan.R inp_file out_file (chromosome) (region_start) (region_end)");
  quit(save="no");
}


# 1) Load libraries ====

library(qqman)

# 2) FUNCTION make_qq(input_file, out_file) ====

make_qq <- function(input_file, out_file, chromo=FALSE, rstart=FALSE, rend=FALSE){
    
    #Read in results file
    raw_results <- read.table(input_file, sep="\t", header=TRUE, stringsAsFactors=FALSE)

    #Substring the chromosome column
    raw_results$chrom <- as.integer(raw_results$chrom)

    #Set minimum p_value=1E-250
    raw_results$pval <- ifelse(raw_results$pval < 5e-324, 5e-324, as.numeric(raw_results$pval))    
    
    #Filter file
    raw_results <- raw_results[c("ID","chrom","pos","pval")]
    raw_results$pos <- as.numeric(raw_results$pos)
    colnames(raw_results) <- c("SNP","CHR","BP","P")
    
    #Open file
    jpeg(out_file, width = 11.5, height = 8, units = 'in', res = 500)
    #Make manhattan plot
    #Check if a specific chromosome/range was given
    attach(raw_results)
    if (rstart != FALSE) {
      manhattan(subset(raw_results, CHR == chromo & BP >= rstart & BP <= rend), suggestiveline = F, annotatePval = 5.0e-8, xlim = c(rstart, rend))
    } else if (chromo != FALSE) {
      manhattan(subset(raw_results, CHR == chromo), suggestiveline = F, annotatePval = 5.0e-8)
    } else {
      manhattan(raw_results, suggestiveline = F, chrlabs = c(1:22, "X"), annotatePval = 5.0e-8)
    }
    dev.off()
    detach(raw_results)
    
}

# 3) Set Call function ====

make_qq(inp_file, out_file, chromo, rstart, rend)
