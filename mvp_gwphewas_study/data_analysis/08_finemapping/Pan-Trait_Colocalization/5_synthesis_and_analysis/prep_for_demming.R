################################################################################

#prep_for_demming.R

#The purpose of this script is to run ancillary analyses in preparation for the
#Demming regression analysis of variant effect sizes.

#This code runs with the R/4.2 module

################################################################################

# 0) Call libraries and set variables and directories ====

#Call libraries
library(tidyverse)

#Set directory
work_dir <- "/project/voight_GWAS/mconery/downstream_analyses/"

# 1) Read in locus file, filter for high-pip signals, and analyze number of signals per trait ====

#Define ancestries
ancestries <- c("AFR", "AMR", "EAS", "EUR")

#Read in file of loci
loci_raw <- read.table(paste0(work_dir, "final_output.grch37.high_thresh.broad.csv"), sep = ",", header = FALSE, stringsAsFactors = FALSE)
colnames(loci_raw) <- c("chr", "start", "end", "length", "traits")

#Make a list of all the trait-locus pairs
trait_locus_list <- str_split(loci_raw$traits, "/")
names(trait_locus_list) <- paste0("chr",loci_raw$chr,".",loci_raw$start,".",loci_raw$end)
add_locus_name <- function(i, trait_locus_list){
  temp <- trait_locus_list[[i]]
  temp_name <- names(trait_locus_list)[i]
  return(paste0(temp, ".", temp_name))
}
temp <- lapply(seq(1, length(trait_locus_list)), FUN = add_locus_name, trait_locus_list=trait_locus_list)
trait_locus_list <- unlist(temp)
remove(temp)

#Convert the list to a dataframe
trait_locus_df <- t(data.frame(str_split(trait_locus_list, "[.]")))
rownames(trait_locus_df) <- trait_locus_list
colnames(trait_locus_df) <- c("trait", "chr", "start", "end")
trait_locus_df <- as.data.frame(trait_locus_df)

#Recast variable type
trait_locus_df$start <- as.numeric(trait_locus_df$start)
trait_locus_df$end <- as.numeric(trait_locus_df$end)

#Read in necessary file
merged_append <- read.table(file = paste0(work_dir, "master.signals.appended.txt"), sep = "\t", header=TRUE)

#Add on a column of binary info for whether the signal was found in one or more ancestries
merged_append <- cbind.data.frame(merged_append, suggestive_evidence=ifelse((is.na(merged_append$AFR.suggestive_variants) == FALSE & merged_append$AFR.suggestive_variants != "Mapped") | 
                                                                              (is.na(merged_append$AMR.suggestive_variants) == FALSE & merged_append$AMR.suggestive_variants != "Mapped") | 
                                                                              (is.na(merged_append$EAS.suggestive_variants) == FALSE & merged_append$EAS.suggestive_variants != "Mapped")| 
                                                                              (is.na(merged_append$EUR.suggestive_variants) == FALSE & merged_append$EUR.suggestive_variants != "Mapped"), "Suggestively Associated\nin Unmapped Population", NA))
#Bind on a column that classifies how to color the signal
merged_append <- cbind.data.frame(merged_append, coloring=ifelse(is.na(merged_append$suggestive_evidence), ifelse(merged_append$max_power_unmapped_ancestries < 0.8, "Underpowered to Detect Suggestive\nAssociation in Unmapped Population", NA), merged_append$suggestive_evidence))
#Add on binary columns for ancestries 
merged_append <- cbind.data.frame(merged_append,
                                  AFR=ifelse(is.na(merged_append$AFR.num_variants), FALSE, TRUE),
                                  AMR=ifelse(is.na(merged_append$AMR.num_variants), FALSE, TRUE),
                                  EAS=ifelse(is.na(merged_append$EAS.num_variants), FALSE, TRUE),
                                  EUR=ifelse(is.na(merged_append$EUR.num_variants), FALSE, TRUE))
merged_append$AFR <- as.logical(merged_append$AFR)
merged_append$AMR <- as.logical(merged_append$AMR)
merged_append$EAS <- as.logical(merged_append$EAS)
merged_append$EUR <- as.logical(merged_append$EUR)

#Make a function that takes 2 ancestries and outputs a comparison box plot
calc_anc_variant_count <- function(merged_append, ANC1, ANC2, work_dir){
  #Filter down for signals found in both ancestries
  merged_ANC1_ANC2 <- merged_append[which(merged_append[,ANC1] == 1 & merged_append[,ANC2] == 1),]
  temp <- merged_ANC1_ANC2[which(merged_ANC1_ANC2[,paste0(ANC1,".num_variants")] == 1 & merged_ANC1_ANC2[,paste0(ANC2,".num_variants")] == 1),]
  temp2 <- table(temp$trait)
  jpeg(paste0(work_dir,ANC1,".",ANC2,".high_pip.shared_signals.beta_compare.jpeg"))
  print(hist(temp2, xlab="Number of Signals per Trait", ylab="Number of Traits", main=NULL, breaks=seq(0,max(temp2),1)))
  dev.off()
  return(c(nrow(merged_ANC1_ANC2), nrow(temp)))
}
#Call function
calc_anc_variant_count(merged_append = merged_append, ANC1 = "AFR", ANC2 = "AMR", work_dir = work_dir)
calc_anc_variant_count(merged_append = merged_append, ANC1 = "AFR", ANC2 = "EUR", work_dir = work_dir)
calc_anc_variant_count(merged_append = merged_append, ANC1 = "AMR", ANC2 = "EUR", work_dir = work_dir)

# 2) Make the files necessary for the Demming Regressions of the High-PiP Variants ====

#Filter for the high-confidence signals in each of the 3 2-way ancestry comparisons
get_comparison_filtered_df <- function(merged_append, ANC1, ANC2){
  #Filter down for signals found in both ancestries
  merged_ANC1_ANC2 <- merged_append[which(merged_append[,ANC1] == 1 & merged_append[,ANC2] == 1),]
  temp <- merged_ANC1_ANC2[which(merged_ANC1_ANC2[,paste0(ANC1,".num_variants")] == 1 & merged_ANC1_ANC2[,paste0(ANC2,".num_variants")] == 1),]
  temp <- temp[,c("trait", "locus", "signal", "merged.best_variant")]
  temp <- temp %>% mutate(trait_variant_combo = paste(trait, merged.best_variant, sep = "."))
  return(temp)
}
AFR_AMR_shared_df <- get_comparison_filtered_df(merged_append, "AFR", "AMR")
AFR_EUR_shared_df <- get_comparison_filtered_df(merged_append, "AFR", "EUR")
AMR_EUR_shared_df <- get_comparison_filtered_df(merged_append, "AMR", "EUR")

#Read in variant-level results file
variant_raw <- read.csv(paste0(work_dir, "master.variants.supplement.txt"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
#Append on variant-trait combo id to each row
variant_raw <- cbind.data.frame(variant_raw, trait_variant_combo=paste0(variant_raw$phenotype, '.', variant_raw$mvp_id))

#Make a function to merge the variant level info onto the signal info and write the resulting table to file
merge_signal_variant_info <- function(ANC_shared_df, variant_raw, ANC1, ANC2, out_dir=work_dir){
  temp_1 <- variant_raw %>% filter(ancestry == ANC1)
  temp_2 <- variant_raw %>% filter(ancestry == ANC2)
  temp_merge <- inner_join(temp_1, temp_2, by = "trait_variant_combo")
  temp <- inner_join(ANC_shared_df, temp_merge, by = "trait_variant_combo")
  colnames(temp) <- str_replace_all(colnames(temp), "[.]x", ".ANC1")
  colnames(temp) <- str_replace_all(colnames(temp), "[.]y", ".ANC2")
  write.table(temp, file = paste0(out_dir, ANC1, ".", ANC2, ".demming_input_table.tsv"), quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
}
merge_signal_variant_info(AFR_AMR_shared_df, variant_raw, "AFR", "AMR")
merge_signal_variant_info(AFR_EUR_shared_df, variant_raw, "AFR", "EUR")
merge_signal_variant_info(AMR_EUR_shared_df, variant_raw, "AMR", "EUR")
