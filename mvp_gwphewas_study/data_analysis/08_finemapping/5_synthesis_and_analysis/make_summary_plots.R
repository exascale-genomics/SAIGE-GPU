################################################################################

#make_summary_plots.R

#The purpose of this script is to run some summary level calculations and make 
#summary level plots for the results of our fine-mapping experiment. These 
#include the following:

# 1) Plots of effect sizes vs allele frequencies for the high-confidence set
# 2) A barplot of coding enrichment 
# 3) A donut plot of the number of variants per credible set
# 4) A scatter plot showing sample sizes vs the number of signals per phenotype.
# 5) An upset plot showing signal sharing across ancestries
# 6) A box plot showing the difference in counts for AFR/EUR shared signals

################################################################################

# 0) Call libraries and set directories ====

#Call libraries
library(plyr)
library(ggpubr)
library(ggplot2)
library(stringr)
library(pbapply)
library(ComplexUpset)
library(data.table)
library(tidyverse)
library(PairedData)
library(ggallin)
library(pbapply)
library(ggrepel)
library(scales)
library(ggpmisc)
library(gridExtra)
library(pwr)
library(VennDiagram)

#Set directory
work_dir <- "/project/voight_GWAS/mconery/downstream_analyses/"
#work_dir <- "C:/Users/mitch/Documents/UPenn/Voight_Lab/MVP_Finemapping/Figure_Building_and_Synthesis/"

# 1) Read in locus file and make dataframe ====

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

# 2) Read in merged results file ====

#Read if file of merged signals
merged_raw <- read.csv(paste0(work_dir, "master.signals.txt"), sep = "\t", header = TRUE)

#Get variant ids for all variants in any credible sets
all_variants <- unique(unlist(str_split(merged_raw$merged.variant_ids, pattern = ",")))

#Extract traits for all variants so that Anurag/Jenny can use them to get mafs/meta p-values/betas
extract_traits <- function(merged_row){
  trait <- merged_row[1]
  num_traits <- length(unlist(str_split(merged_row[13], pattern = ",")))
  return(rep(trait, num_traits))
} 
all_variants_traits <- cbind.data.frame(trait=unlist(apply(merged_raw, MARGIN = 1, FUN = extract_traits)),
                                        variant=unlist(str_split(merged_raw$merged.variant_ids, pattern = ",")))
all_variants_traits <- unique(all_variants_traits)
#Make pip dictionaries
get_ancestry_pips <- function(merged_raw, ANC){
  ANC_high_pips <- merged_raw[which(merged_raw[,paste0(ANC,".max_overall_pip")] >= 0.95), c(paste0(ANC,".best_variant"), "trait", paste0(ANC,".max_overall_pip"))]
  ANC_high_pips <- ANC_high_pips[order(ANC_high_pips[,paste0(ANC,".max_overall_pip")]),]
  ANC_high_pips <- cbind.data.frame(ANC_high_pips, locus_variant=paste0(ANC_high_pips[,"trait"], ".", ANC_high_pips[,paste0(ANC,".best_variant")]))
  ANC_high_pips <- ANC_high_pips[which(duplicated(ANC_high_pips$locus_variant) == FALSE),]
  row.names(ANC_high_pips) <- ANC_high_pips$locus_variant
  return(ANC_high_pips)
}
AFR_high_pips <- get_ancestry_pips(merged_raw, "AFR")
AMR_high_pips <- get_ancestry_pips(merged_raw, "AMR")
EAS_high_pips <- get_ancestry_pips(merged_raw, "EAS")
EUR_high_pips <- get_ancestry_pips(merged_raw, "EUR")
#Bind high_pip info onto the all_variants_traits array
all_variants_traits <- cbind.data.frame(all_variants_traits, 
                                        AFR_high_pip=AFR_high_pips[paste0(all_variants_traits$trait, ".", all_variants_traits$variant),"AFR.max_overall_pip"],
                                        AMR_high_pip=AMR_high_pips[paste0(all_variants_traits$trait, ".", all_variants_traits$variant),"AMR.max_overall_pip"],
                                        EAS_high_pip=EAS_high_pips[paste0(all_variants_traits$trait, ".", all_variants_traits$variant),"EAS.max_overall_pip"],
                                        EUR_high_pip=EUR_high_pips[paste0(all_variants_traits$trait, ".", all_variants_traits$variant),"EUR.max_overall_pip"])
#Write the file
write.table(all_variants_traits, paste0(work_dir, "all_fine-mapped_variants.txt"), quote = FALSE, row.names = FALSE, col.names = TRUE)

# 3) Calculate Number of Double Counted Signals and Overlapping loci ====

#Run convert_to_mappable_list.py from the "1_locus_defining" folder of the repo to 
#split the list of loci up into a single line for each locus trait pair
#Read in the results
loci_sep <- read.table(paste0(work_dir, "final_output.grch37.high_thresh.broad.trait_separated.csv"), sep = ",", header = FALSE, stringsAsFactors = FALSE)
colnames(loci_sep) <- c("chr", "start", "end", "length", "trait")

#Filter the loci_raw file for overlapping loci
find_overlapping <- function(locus_row, locus_sep){
  loc_chromo = as.numeric(locus_row["chr"])
  loc_start = as.numeric(locus_row["start"])
  loc_end = as.numeric(locus_row["end"])
  trait = locus_row["trait"]
  locus_filt = locus_sep[which(locus_sep$chr == loc_chromo & locus_sep$trait == trait &
                                 ((locus_sep$start < loc_start & locus_sep$end > loc_start) |
                                    (locus_sep$start < loc_end & locus_sep$end > loc_end))),]
  return(ifelse(nrow(locus_filt) > 0, TRUE, FALSE))
}
loci_overlap <- cbind.data.frame(loci_sep, overlap=apply(loci_sep, MARGIN = 1, FUN = find_overlapping, locus_sep = loci_sep))
loci_overlap <- loci_overlap[which(loci_overlap$overlap == TRUE),]
loci_overlap <- loci_overlap[order(loci_overlap$trait, loci_overlap$chr, loci_overlap$start),]
#Count loci-trait pairs that overlap
nrow(loci_overlap[which(loci_overlap$overlap == TRUE),]) # It's 1,583 locus-trait pairs (6% of all locus-trait pairs)

#Filter down the merged_raw dataframe for signals found in the overlapping loci
overlapping_loci_traits <- paste0(loci_overlap$trait, ".chr", loci_overlap$chr, ".", loci_overlap$start, ".", loci_overlap$end)
merged_potential_overlap <- merged_raw[which(paste0(merged_raw$trait, ".", merged_raw$locus) %in% overlapping_loci_traits),]
#Create function that will check each signal in the overlapping loci against all others in all ancestries
check_signal_overlap <- function(merged_row, merged_potential_overlap, ancestry_names=ancestries){
  merged_trait = as.character(merged_row["trait"])
  merged_chr = as.character(merged_row["chr"])
  merged_locus = as.character(merged_row["locus"])
  merged_filt <- merged_potential_overlap[which(merged_potential_overlap$trait == merged_trait & 
                                                  merged_potential_overlap$chr == merged_chr & 
                                                  merged_potential_overlap$locus != merged_locus),] #This is the list of loci that aren't the same as the merged row but are for same trait on same chromosome
  merged_row_variants <- list(unlist(str_split(merged_row["AFR.variant_ids"], ",")), unlist(str_split(merged_row["AMR.variant_ids"], ",")),
                              unlist(str_split(merged_row["EAS.variant_ids"], ",")), unlist(str_split(merged_row["EUR.variant_ids"], ",")))
  names(merged_row_variants) <- ancestry_names
  variant_intersects <- lapply(ancestries, FUN = check_anc_signal_overlap, merged_row_variants=merged_row_variants, merged_filt=merged_filt)
  return(paste0(unique(variant_intersects[!(is.na(variant_intersects))]), collapse = ","))
}
check_anc_signal_overlap <- function(ANC, merged_row_variants, merged_filt){
  anc_variants <- unlist(merged_row_variants[ANC])
  if (all(is.na(anc_variants))) {
    return(NA)
  }
  check_variants <- str_split(merged_filt[,paste0(ANC, ".variant_ids")], ",")
  merged_intersect <- merged_filt[which(ifelse(lapply(lapply(check_variants, FUN = intersect, y=anc_variants),length) > 0, TRUE, FALSE)),]
  if (nrow(merged_intersect) > 0) {
    return(paste0(merged_intersect$trait, ".", merged_intersect$locus, ".", merged_intersect$signal, collapse = ","))
  } else {
    return(NA)
  }
}
merged_overlap <- cbind.data.frame(merged_potential_overlap, overlap=pbapply(merged_potential_overlap, FUN = check_signal_overlap, MARGIN = 1, merged_potential_overlap=merged_potential_overlap))

#Calculate the number of signals that are duplicated due to locus overlap
nrow(merged_overlap[which(merged_overlap$overlap != ""),]) # 231 signals that are due to duplication
sum(unlist(lapply(str_split(merged_overlap[which(merged_overlap$overlap != ""),"overlap"], ","), length))!=1) #All are pairs except for 1
nrow(merged_overlap[which(merged_overlap$overlap != ""),])/nrow(merged_raw) #0.4% of signals are affected by the locus overlap problem
#Calculate lengths of merged loci if we had to rerun
extract_element <- function(overlap_trait_locus, element_num){str_split(overlap_trait_locus,"[.]")[[1]][element_num]}
starts_1 <- as.numeric(vapply(str_replace(merged_overlap$overlap[merged_overlap$overlap != ""], ",", "[.]"), FUN = extract_element, element_num = 3, FUN.VALUE = character(1)))
ends_1 <- as.numeric(vapply(str_replace(merged_overlap$overlap[merged_overlap$overlap != ""], ",", "[.]"), FUN = extract_element, element_num = 4, FUN.VALUE = character(1)))
starts_2 <- as.numeric(merged_overlap[which(merged_overlap$overlap != ""),"start"])
ends_2 <- as.numeric(merged_overlap[which(merged_overlap$overlap != ""),"end"])
get_list_function <- function(element_num, list_function, list_1, list_2){return(list_function(list_1[element_num], list_2[element_num]))}
max(vapply(seq(1, length(ends_1), 1), FUN = get_list_function, FUN.VALUE = 0, list_function=max, list_1=ends_1, list_2=ends_2) - vapply(seq(1, length(starts_1), 1), FUN = get_list_function, FUN.VALUE = 0, list_function=min, list_1=starts_1, list_2=starts_2) + 1) 
#Largest merged locus would be 3.25Mb

# 4) Read in the list of all significant variants with their betas and p-values and process accordingly (Commented out as the work is already completed and it's intensive to run) ====

##Read in the RData file
load(paste0(work_dir, "all_fine-mapped_variants.withPIP.withsummarystats.withR2_Panel.RData"))
#
##Create unique identifiers for the rows in final
rownames(final) <- paste0(final$ID, ".", final$Trait)

##Set suffixes for ancestry specific columns in final
ancestry_column_suffixes <- c(".num_variants", ".max_overall_pip", ".best_variant", 
                              ".best_variant_cs_pip", ".best_variant_ancestry.specific_neglogp", ".variant_ids")
#Make a function that extracts ancestry specific column info
extract_ancestry_columns <- function(ancestry, merged_row, ancestry_column_suffixes, final, trait){
  temp <- merged_row[paste0(ancestry, ancestry_column_suffixes)]
  temp_beta <- final[paste0(temp[paste0(ancestry, ".best_variant")],".",trait),paste0("beta.",ancestry)]
  names(temp_beta) <- paste0(ancestry, ".best_variant_beta")
  return(c(temp[1:3], temp_beta, temp[4:6]))
}
##Make a function that checks for possible suggestively significant variants
check_suggestive <- function(ancestry, final, mapped_ancestries, output_row, suggest_sig=1e-3){
  #Check if we can exit immediately
  if (mapped_ancestries[ancestry]) {
    return("Mapped")
  } else {
    merged_variants <- unlist(str_split(output_row["merged.variant_ids"],","))
    trait <- output_row["trait"]
    meta_betas <- final[paste0(merged_variants, ".", trait),"beta.META"]
    anc_betas_pvals <- final[paste0(merged_variants, ".", trait),c(paste0("beta.",ancestry),paste0("pval.",ancestry))]
    matching_anc_pvals <- anc_betas_pvals[ifelse(meta_betas*anc_betas_pvals[,1] > 0, TRUE, FALSE),2]
    suggest_sig_index <- ifelse(is.na(matching_anc_pvals), FALSE, matching_anc_pvals<=suggest_sig)
    if(sum(suggest_sig_index)>0){
      return(paste0(paste0(merged_variants[suggest_sig_index], "(", anc_betas_pvals[which(suggest_sig_index),2], ")"), collapse = ","))
    }
    else{
      return(NA)
    }
  }
}
#Make a function that takes in a list of merged variants for a signal as well as their ancestry specific breakdown and returns the minimum absolute beta for each
check_element <- function(character_vector, variant){return(variant %in% character_vector)}
extract_element <- function(inp_list, element_num){return(inp_list[element_num])}
get_min_abs_beta <- function(variant, abs_check_betas, check_ancestries_variants){
  relevant_ancestries <- names(which(unlist(lapply(check_ancestries_variants, FUN = check_element, variant=variant))))
  return(min(abs_check_betas[variant, paste0("beta.", relevant_ancestries)]))
}
str_split_unlist <- function(string, pattern){return(unlist(unlist(str_split(string, pattern = pattern))))} #Splits a string and unlists the output to a vector
#Make a function that runs a power test for a given variant in all unmapped ancestries
check_variant_power_all_unmapped <- function(variant, unmapped_sample_sizes, min_abs_check_betas, check_ses, unmapped_ancestries){
  ancestry_sample_sizes <- as.data.frame(unmapped_sample_sizes[variant,])
  colnames(ancestry_sample_sizes) <- colnames(unmapped_sample_sizes)
  rownames(ancestry_sample_sizes) <- variant
  ancestry_ses <- as.data.frame(check_ses[variant,])
  colnames(ancestry_ses) <- colnames(check_ses)
  rownames(ancestry_ses) <- variant
  abs_beta <- min_abs_check_betas[variant]
  powers <- vapply(names(ancestry_sample_sizes), FUN = check_variant_power_one_unmapped, FUN.VALUE = 0, abs_beta=abs_beta, ancestry_sample_sizes=ancestry_sample_sizes, ancestry_ses=ancestry_ses)
  return(max(powers, na.rm = TRUE))
}
check_variant_power_one_unmapped <- function(ancestry, abs_beta, ancestry_sample_sizes, ancestry_ses, sig_level = 1e-3){
  reject_point <- qt(p = sig_level/2, lower.tail = FALSE, df = as.numeric(ancestry_sample_sizes[,ancestry]) - 2) * as.numeric(ancestry_ses[,ancestry])
  power = pt((reject_point - abs_beta)/as.numeric(ancestry_ses[,ancestry]), df = as.numeric(ancestry_sample_sizes[,ancestry]) - 2, lower.tail = FALSE)
  power = power + pt((-reject_point - abs_beta)/as.numeric(ancestry_ses[,ancestry]), df = as.numeric(ancestry_sample_sizes[,ancestry]) - 2, lower.tail = TRUE)
}
#Make a function that checks for lack of power
check_power <- function(output_row, final, mapped_ancestries){
  #Check if we can exit immediately
  if (sum(mapped_ancestries) == length(mapped_ancestries)) {
    return(NA)
  } else {
    #Get basic information
    check_ancestries <- names(mapped_ancestries)[mapped_ancestries] #Get ancestries from which to check the betas
    unmapped_ancestries <- names(mapped_ancestries)[!(mapped_ancestries)]
    merged_variants <- unlist(str_split(output_row["merged.variant_ids"],","))
    trait <- output_row["trait"]
    #Get breakdown of variants in credible set for each mapped ancestries
    check_ancestries_variants <- unlist(output_row[paste0(check_ancestries,".variant_ids")])
    names(check_ancestries_variants) <- check_ancestries
    check_ancestries_variants <- lapply(check_ancestries_variants, FUN = str_split_unlist, pattern = ",")
    #Get absolute betas for each variant across the ancestries in which the variant appeared in the credible set
    abs_check_betas <- abs(as.data.frame(final[paste0(merged_variants, ".", trait),paste0("beta.", check_ancestries)]))
    rownames(abs_check_betas) <- merged_variants
    colnames(abs_check_betas) <- paste0("beta.", check_ancestries)
    #Get paired standard errors to calculate Cohen's Ds
    check_ses <- as.data.frame(final[paste0(merged_variants, ".", trait),paste0("se.", unmapped_ancestries)])
    rownames(check_ses) <- merged_variants
    colnames(check_ses) <- paste0(unmapped_ancestries)
    #Get sample sizes for unmapped_ancestries
    unmapped_sample_sizes <- as.data.frame(final[paste0(merged_variants,".",trait),paste0("N.", unmapped_ancestries)])
    rownames(unmapped_sample_sizes) <- paste0(merged_variants)
    colnames(unmapped_sample_sizes) <- paste0(unmapped_ancestries)
    #Get minimum absolute betas for each variant in the ancestries in which the variant was mapped
    min_abs_check_betas <- vapply(merged_variants, get_min_abs_beta, FUN.VALUE = 0, abs_check_betas=abs_check_betas, check_ancestries_variants=check_ancestries_variants)
    #Calculate maximum power for any variants in the list
    max_power_merged_variants <- vapply(merged_variants, FUN = check_variant_power_all_unmapped, FUN.VALUE = 0, unmapped_sample_sizes=unmapped_sample_sizes, min_abs_check_betas=min_abs_check_betas, check_ses=check_ses, unmapped_ancestries=unmapped_ancestries)
    return(max(max_power_merged_variants))
  }
}
#Make a function that adds data to the merged dataframe for a single row
add_data_to_merged_row <- function(merged_row, final){
  #Identify best ancestry
  best_ancestry <- ancestries[which.max(c(merged_row["AFR.max_overall_pip"], merged_row["AMR.max_overall_pip"], merged_row["EAS.max_overall_pip"], merged_row["EUR.max_overall_pip"]))]
  #Begin creating the output_row
  output_row <- merged_row[1:11]
  #Extract the meta-analyzed p-values and betas
  output_row <- c(output_row,final[paste0(merged_row["merged.best_variant"],".",merged_row["trait"]),c("pval.META","beta.META")])
  names(output_row)[(length(output_row)-1):length(output_row)] <- c("merged.best_variant.pval.META", "merged.best_variant.beta.META")
  output_row <- c(output_row, merged_row[12:13],best_ancestry=best_ancestry)
  #Append on ancestry specific info
  output_row <- c(output_row,
                  unlist(lapply(ancestries, FUN = extract_ancestry_columns, ancestry_column_suffixes=ancestry_column_suffixes, merged_row=merged_row, final=final, trait=merged_row["trait"])))
  #Identify ancestries in which signal is found
  mapped_ancestries <- as.vector(is.na(merged_row[paste0(ancestries, ".best_variant")])==FALSE)
  names(mapped_ancestries) <- ancestries
  #Check for suggestive association in unmapped ancestries
  temp <- unlist(lapply(ancestries, FUN = check_suggestive, final=final, mapped_ancestries=mapped_ancestries, output_row=output_row))
  names(temp) <- paste0(ancestries, ".suggestive_variants")
  output_row <- c(output_row, temp)
  #Check for lack of power as well
  output_row <- c(output_row, max_power_unmapped_ancestries=as.numeric(check_power(output_row, final, mapped_ancestries)))
  #Return the output_row
  return(unlist(output_row))
}
#Apply the outer function to make a new version of the merged file
merged_append <- as.data.frame(t(pbapply(merged_raw, FUN = add_data_to_merged_row, MARGIN = 1, final=final)))
#Write to file
write.table(merged_append, file = paste0(work_dir, "master.signals.appended.txt"), sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

# 5) Assess High-PIP Variants Betas vs MAFs and calculate statistics about high-confidence  ====

#Read in the RData file
load(paste0(work_dir, "all_fine-mapped_variants.withPIP.withsummarystats.withR2_Panel.RData"))
#Create unique identifiers for the rows in final
rownames(final) <- paste0(final$ID, ".", final$Trait)

#Prepare to merge dataframe of pips and dataframe of variant data (i.e. final)
rownames(all_variants_traits) <- paste0(all_variants_traits$variant, ".", all_variants_traits$trait)
#Merge dataframes
complete_high_pip_df <- all_variants_traits[which(rowSums(!(is.na(all_variants_traits[,c("AFR_high_pip", "AMR_high_pip", "EAS_high_pip", "EUR_high_pip")]))) > 0),]
complete_high_pip_df <- cbind.data.frame(complete_high_pip_df, final[rownames(complete_high_pip_df),]) #The info in final is from the GWAS summary statistics that Jenny procured

### Write list of high-pip variants to text file in format needed for vcf manipulation ###
temp <- str_split(complete_high_pip_df[which(!(duplicated(complete_high_pip_df$variant))),"variant"], pattern = ":", simplify = TRUE)
temp <- cbind.data.frame(temp[,c(1,2)])
write.table(temp, file = paste0(work_dir, "high_pip_targets.for_vcf.txt"), sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

#Read in the trait dictionary
trait_raw <- read.csv(paste0(work_dir, "MVP_R4.1000G_AGR.GIA.DataDictionary.revised.022023.txt"), header = TRUE, sep = "\t")
rownames(trait_raw) <- trait_raw$Trait
#Tack on description column to the complete_high_pip_df
complete_high_pip_df <- cbind.data.frame(complete_high_pip_df, trait_description=trait_raw[complete_high_pip_df$trait,"DescriptionShort"])

#Read in list of known/novel variant/trait combos
known_novel_raw <- read.table(paste0(work_dir, "all.pheno.finemap.known-novel.05152023.txt"), header = TRUE, sep = "\t")
known_novel_raw <- cbind.data.frame(known_novel_raw, unique_id=paste0(known_novel_raw$SNP, ".", known_novel_raw$phenotype))
known_novel_raw <- known_novel_raw[which(!(duplicated(known_novel_raw$unique_id))),]
rownames(known_novel_raw) <- known_novel_raw$unique_id
#Recode known_novel_raw column
known_novel_raw$Novel_label <- ifelse(known_novel_raw$Known_Association=="True", "Known Association",
                                ifelse(known_novel_raw$Novel_Association_Known_Signal=="True", "Unreported Assoc. Known Signal", 
                                       ifelse(known_novel_raw$Novel_Signal=="True", "Previously Unreported Signal", NA)))
#Tack novelty data onto complete_high_pip_df (data frame of all high-pip variants for all ancestries)
complete_high_pip_df <- cbind.data.frame(complete_high_pip_df, Novel=known_novel_raw[rownames(complete_high_pip_df,),"Novel_label"])

#Make a dataframe that splits each point with multiple ancestries-worth of high pips out onto separate lines
extract_single_ancestry_high_pips <- function(ancestry, complete_high_pip_df){
  if (ancestry != "META") {
    temp <- complete_high_pip_df[which(!(is.na(complete_high_pip_df[,paste0(ancestry,"_high_pip")]))),
                                 c("trait", "variant", "rsid", "ttype", "trait_description", "Novel", paste0(c("pval.", "eaf.", "beta.", "se.", "N."), ancestry))]
  } else {
    temp <- complete_high_pip_df[,c("trait", "variant", "rsid", "ttype", "trait_description", "Novel", paste0(c("pval.", "eaf.", "beta.", "se.", "N."), ancestry))]
  }
  colnames(temp) <- c("trait", "variant", "rsid", "ttype", "trait_description", "Novel", paste0(c("pval.", "eaf.", "beta.", "se.", "N."), "ANC"))
  temp <- cbind.data.frame(temp, ANC=rep(ancestry, nrow(temp)))
  return(temp)
}
#Call Function
separated_high_pip_df <- rbind.fill(lapply(ancestries, FUN = extract_single_ancestry_high_pips, complete_high_pip_df=complete_high_pip_df))

#Write the dataframe to file for Anurag
write.table(separated_high_pip_df, file = paste0(work_dir, "beta_v_maf.data_frame.txt"), sep = "\t", col.names = TRUE, quote = FALSE, row.names = FALSE)

#Calculate contributions to genetic variance
separated_high_pip_df <- cbind.data.frame(separated_high_pip_df, 
                                          heri_cont=2*separated_high_pip_df$beta.ANC^2*separated_high_pip_df$eaf.ANC*(1-separated_high_pip_df$eaf.ANC))
#Create labels for highest contributors
separated_high_pip_df <- separated_high_pip_df[order(rank(separated_high_pip_df$heri_cont, na.last = FALSE), decreasing = TRUE),]
outlier_labels = vector()
for (i in 1:nrow(separated_high_pip_df)) {
  clean_current_label <- str_replace(str_replace(str_replace(str_replace(paste0(separated_high_pip_df[i, c("trait_description","rsid","ANC")], collapse = "\n"), ", maximum", replacement = ""), ", mean", replacement = ""), ", minimum", replacement = ""), ", inv-norm transformed", replacement = "")
  if (!(clean_current_label %in% outlier_labels)) {
    outlier_labels <- append(outlier_labels, clean_current_label)
  } else {
    outlier_labels <- append(outlier_labels, NA)
  }
}
for(ancestry in ancestries){
  outlier_labels <- str_replace(outlier_labels, paste0("\n",ancestry), "")
}
separated_high_pip_df <- cbind.data.frame(separated_high_pip_df, outlier_labels)

#Read in derived vs ancestral allele info and clean
ancestral_raw_list <- lapply(paste0(work_dir, "high_pip_filtered_vcfs/", list.files(paste0(work_dir, "high_pip_filtered_vcfs"))), read.table, header=FALSE, stringsAsFactors = FALSE)
ancestral_raw <- as.data.frame(rbind.fill.matrix(ancestral_raw_list))
colnames(ancestral_raw) <- c("chr", "pos", "rsid", "ref", "alt", "junk1", "junk2", "AA")
ancestral_raw$chr <- as.numeric(ancestral_raw$chr)
ancestral_raw$pos <- as.numeric(ancestral_raw$pos)
#Fix the ancestral allele column by removing rows missing the info
ancestral_raw <- ancestral_raw[which(substr(ancestral_raw$AA, 1, 3) == "AA="),]
ancestral_raw$AA <- str_replace(ancestral_raw$AA, pattern = "AA=", replacement = "")
#Set rownames as rsids and use them to convert
rownames(ancestral_raw) <- ancestral_raw$rsid
#Tack on ancestral info to complete_high_pip_df and separated_high_pip_df
complete_high_pip_df <- cbind.data.frame(complete_high_pip_df, AA=ancestral_raw[complete_high_pip_df$rsid,"AA"])
separated_high_pip_df <- cbind.data.frame(separated_high_pip_df, AA=ancestral_raw[separated_high_pip_df$rsid,"AA"])

#Make dataframes for plotting
finemap_signal_summ_q <- separated_high_pip_df[which(separated_high_pip_df$ttype=="quantitative"),]
finemap_signal_summ_q <- finemap_signal_summ_q %>%
  mutate(eaf.ANC.up=ifelse(eaf.ANC>0.5, 1-eaf.ANC, eaf.ANC)) %>%
  mutate(beta.ANC.up=ifelse(eaf.ANC>0.5, -(beta.ANC), beta.ANC))
finemap_signal_summ_b <- separated_high_pip_df[which(separated_high_pip_df$ttype=="binary"),]
finemap_signal_summ_b <- finemap_signal_summ_b %>%
  mutate(eaf.ANC.up=ifelse(eaf.ANC>0.5, 1-eaf.ANC, eaf.ANC)) %>%
  mutate(or.ANC.up=ifelse(eaf.ANC>0.5, exp(-(beta.ANC)), exp(beta.ANC)))
#Recode novelty columns as factors
finemap_signal_summ_q$Novel <- factor(finemap_signal_summ_q$Novel, levels = c("Previously Unreported Signal", "Unreported Assoc. Known Signal", "Known Association"))
finemap_signal_summ_b$Novel <- factor(finemap_signal_summ_b$Novel, levels = c("Previously Unreported Signal", "Unreported Assoc. Known Signal", "Known Association"))
#Get minimum and maximum values for y-axis on plots
quantitative_beta_limit <- max(abs(finemap_signal_summ_q$beta.ANC.up), na.rm = TRUE)
binary_or_limit <- max(exp(abs(finemap_signal_summ_b$beta.ANC)), na.rm = TRUE)

#Make data tables for annotating the figures
Novel_Map <- c("Previously Unreported\nSignal", "Known\nAssociation", "Unreported Assoc.\nKnown Signal")
names(Novel_Map) <- c("Previously Unreported Signal", "Known Association", "Unreported Assoc. Known Signal")
quantitative_annotate_table <- as.data.frame(table(finemap_signal_summ_q[,c("ANC","Novel")]))
quantitative_annotate_table$Novel <- Novel_Map[quantitative_annotate_table$Novel]
binary_annotate_table <- as.data.frame(table(finemap_signal_summ_b[,c("ANC", "Novel")]))
binary_annotate_table$Novel <- Novel_Map[binary_annotate_table$Novel]
for (ANC in unique(quantitative_annotate_table$ANC)) {
  quantitative_annotate_table[which(quantitative_annotate_table$ANC == ANC),"Freq"] <-  quantitative_annotate_table[which(quantitative_annotate_table$ANC == ANC),"Freq"]/sum(quantitative_annotate_table[which(quantitative_annotate_table$ANC == ANC),"Freq"])
}
for (ANC in unique(binary_annotate_table$ANC)) {
  binary_annotate_table[which(binary_annotate_table$ANC == ANC),"Freq"] <-  binary_annotate_table[which(binary_annotate_table$ANC == ANC),"Freq"]/sum(binary_annotate_table[which(binary_annotate_table$ANC == ANC),"Freq"])
}

#Make inset plots for minor allele plots
quantitative_inset <- ggplot(quantitative_annotate_table, aes(y=Freq, x=ANC, alpha=Novel)) + 
  geom_bar(aes(fill=ANC), stat="identity") + 
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0%", "25%", "50%", "75%", "100%"), limits = c(0,1)) + 
  scale_fill_manual(values = c("#56B893", "#F87850", "#7B8DBF", "#DF71B6")) + 
  scale_alpha_manual(values = c(1, 0.6, 0.3, 0.1)) + 
  ylab(label = "Prior Report") + 
  theme(axis.title.x = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA), legend.title = element_blank(), 
        axis.line = element_line(colour = "black"), legend.position = "none", 
        plot.background = element_rect(fill = "transparent",colour = NA),
        axis.title.y = element_text(size = 14), axis.text.y = element_text(size = 14), axis.text.x = element_text(size = 14))
data.tb.quantitative <- tibble(x = 0.5, y = -quantitative_beta_limit, plot = list(quantitative_inset))
binary_inset <- ggplot(binary_annotate_table, aes(y=Freq, x=ANC, alpha=Novel)) + 
  geom_bar(aes(fill=ANC), stat="identity") + 
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0%", "25%", "50%", "75%", "100%"), limits = c(0,1)) + 
  scale_fill_manual(values = c("#56B893", "#F87850", "#7B8DBF", "#DF71B6")) + 
  scale_alpha_manual(values = c(1, 0.6, 0.3, 0.1)) + 
  ylab(label = "Prior Report") + 
  theme(axis.title.x = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA), legend.title = element_blank(), 
        axis.line = element_line(colour = "black"), legend.position = "none", 
        plot.background = element_rect(fill = "transparent",colour = NA), 
        axis.title.y = element_text(size = 14), axis.text.y = element_text(size = 14), axis.text.x = element_text(size = 14))
data.tb.binary <- tibble(x = 0.5, y = 1/binary_or_limit, plot = list(binary_inset))

#Make legend plot
legend_df <- cbind.data.frame(x = unlist(lapply(seq(1,10,3), FUN = rep, times = 3)), y = rep(seq(1.2,1,-0.1), 4), 
                              ANC=unlist(lapply(ancestries, FUN = rep, times = 3)), Novel=rep(levels(finemap_signal_summ_b$Novel), 4))
legend_labels_df <- cbind.data.frame(x=c(rep(0, 3), seq(1,10,3)), y=c(seq(1.2,1,-0.1), rep(1.33, 4)), 
                                     label=c(levels(finemap_signal_summ_b$Novel), ancestries), hjust = c(rep(1,3), rep(0.5,4)),
                                     angle=c(rep(0,3), rep(45,4)))
legend_df$Novel <- factor(legend_df$Novel, levels = levels(finemap_signal_summ_b$Novel))
legend_plot <- ggplot(data = legend_df) + 
  geom_point(mapping = aes(x = x, y = y, color = ANC, alpha = Novel, shape = Novel), size = 4) + 
  geom_text(data = legend_labels_df, mapping = aes(x = x, y = y, label = label, hjust = hjust, angle=angle), size = 4) + 
  xlim(-21, 10.1) + 
  ylim(0.95,1.4) + 
  scale_color_manual(values = c("#56B893", "#F87850", "#7B8DBF", "#DF71B6")) + 
  scale_alpha_manual(values = c(1, 0.6, 0.3, 0.1)) + 
  theme_void() + 
  theme(legend.title=element_blank(), axis.title=element_text(size=16), legend.position = "none", 
        legend.text = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_blank(), plot.background = element_blank(), axis.text = element_blank(), 
        axis.title.x = element_blank(), axis.title.y = element_blank(), axis.ticks = element_blank()) 
data.tb.legend.minor <- tibble(x = 0.50, y = quantitative_beta_limit, plot = list(legend_plot))
data.tb.legend.derived <- tibble(x = 1, y = quantitative_beta_limit, plot = list(legend_plot))

#Make minor allele plots for quantitative trait and binary trait
quant_minor_plot <- ggplot() +
  geom_point(data = finemap_signal_summ_q[which(finemap_signal_summ_q$ANC == "EUR"),], aes(x=eaf.ANC.up, y=beta.ANC.up, color=ANC, shape=Novel, alpha=Novel)) + 
  geom_point(data = finemap_signal_summ_q[which(finemap_signal_summ_q$ANC != "EUR"),], aes(x=eaf.ANC.up, y=beta.ANC.up, color=ANC, shape=Novel, alpha=Novel)) +
  scale_y_continuous(limits = c(-quantitative_beta_limit, quantitative_beta_limit), breaks = seq(-3,3,1), labels = format(seq(-3,3,1), width = 5, justify = "right")) + 
  xlim(c(0, 0.5)) + ylab("Beta of Minor Allele") + 
  geom_plot(data = data.tb.quantitative, aes(x, y, label = plot), vp.width = 0.6) +
  geom_plot(data = data.tb.legend.minor, aes(x, y, label = plot), vp.width = 0.8, vp.height = 0.20) + 
  theme_minimal() +
  scale_color_manual(values = c("#56B893", "#F87850", "#7B8DBF", "#DF71B6")) + 
  scale_alpha_manual(values = c(1, 0.6, 0.3, 0.1)) + 
  theme(legend.title=element_blank(), axis.title=element_text(size=16), legend.position = "none", 
        legend.text = element_text(size=12), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"), legend.box="vertical", axis.title.x = element_blank(), 
        axis.text = element_text(size = 18), axis.title.y = element_text(size = 18)) + 
  guides(colour = guide_legend(override.aes = list(size=6)), shape = guide_legend(override.aes = list(size=6))) 
binary_minor_plot <- ggplot() +
  geom_point(data = finemap_signal_summ_b[which(finemap_signal_summ_b$ANC == "EUR"),], aes(x=eaf.ANC.up, y=or.ANC.up, color=ANC, shape=Novel, alpha=Novel)) + 
  geom_point(data = finemap_signal_summ_b[which(finemap_signal_summ_b$ANC != "EUR"),], aes(x=eaf.ANC.up, y=or.ANC.up, color=ANC, shape=Novel, alpha=Novel)) +
  scale_y_continuous(limits = c(1/binary_or_limit, binary_or_limit), 
                     trans = "log", breaks = 10^(seq(-3,3,1)), labels = format(c("0.001","0.01", "0.1", "1", "10", "100", "1,000"), width = 5, justify = "right")) + 
  xlab("MAF") + xlim(c(0, 0.5)) + 
  ylab("Odds Ratio of Minor Allele") +  
  geom_plot(data = data.tb.binary, aes(x, y, label = plot), vp.width = 0.6) + 
  theme_minimal() +
  scale_color_manual(values = c("#56B893", "#F87850", "#7B8DBF", "#DF71B6")) + 
  scale_alpha_manual(values = c(1, 0.6, 0.3, 0.1)) + 
  theme(legend.title=element_blank(), axis.title=element_text(size=16), legend.position = "none", 
        legend.text = element_text(size=12), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"), legend.box="vertical", axis.text = element_text(size = 16),
        axis.title.x = element_text(size = 18), axis.title.y = element_text(size = 18)) + 
  guides(colour = guide_legend(override.aes = list(size=6)), shape = guide_legend(override.aes = list(size=6)))

#Make full combined minor allele plot
jpeg(paste0(work_dir, "minor_allele.effect_size_v_maf.scatter.jpeg"), width = 6, height = 11.5, units = 'in', res = 500)
ggarrange(quant_minor_plot, binary_minor_plot, ncol = 1, nrow = 2, common.legend = FALSE, align = "v", heights = c(4.25, 4.25))
dev.off()

### Make dataframes for plotting the ancestral allele frequency/effect sizes ###
extract_element <- function(overlap_trait_locus, element_num){str_split(overlap_trait_locus,":")[[1]][element_num]}
separated_high_pip_ea_df <- cbind.data.frame(separated_high_pip_df, ea=unlist(lapply(separated_high_pip_df$variant, extract_element, element_num=3)), nea=unlist(lapply(separated_high_pip_df$variant, extract_element, element_num=4)))
finemap_signal_derived_q <- separated_high_pip_ea_df[which(separated_high_pip_df$AA %in% c("A", "C", "G", "T")),] %>%
  filter(ttype=="quantitative") %>%
  mutate(derived_af=ifelse(ea == AA, 1-eaf.ANC, eaf.ANC)) %>%
  mutate(derived_beta=ifelse(ea == AA, -beta.ANC, beta.ANC))
finemap_signal_derived_b <- separated_high_pip_ea_df[which(separated_high_pip_df$AA %in% c("A", "C", "G", "T")),]  %>%
  filter(ttype=="binary") %>%
  mutate(derived_af=ifelse(ea == AA, 1-eaf.ANC, eaf.ANC)) %>%
  mutate(derived_or=ifelse(ea == AA, exp(-(beta.ANC)), exp(beta.ANC)))
#Recode novelty columns as factors
finemap_signal_derived_q$Novel <- factor(finemap_signal_derived_q$Novel, levels = c("Previously Unreported Signal", "Unreported Assoc. Known Signal", "Known Association"))
finemap_signal_derived_b$Novel <- factor(finemap_signal_derived_b$Novel, levels = c("Previously Unreported Signal", "Unreported Assoc. Known Signal", "Known Association"))

#Make data tables for annotating the figures
quantitative_derived_annotate_table <- as.data.frame(table(finemap_signal_derived_q[,c("ANC","Novel")]))
quantitative_derived_annotate_table$Novel <- Novel_Map[quantitative_derived_annotate_table$Novel]
binary_derived_annotate_table <- as.data.frame(table(finemap_signal_derived_b[,c("ANC", "Novel")]))
binary_derived_annotate_table$Novel <- Novel_Map[binary_derived_annotate_table$Novel]
for (ANC in unique(quantitative_derived_annotate_table$ANC)) {
  quantitative_derived_annotate_table[which(quantitative_derived_annotate_table$ANC == ANC),"Freq"] <-  quantitative_derived_annotate_table[which(quantitative_derived_annotate_table$ANC == ANC),"Freq"]/sum(quantitative_derived_annotate_table[which(quantitative_derived_annotate_table$ANC == ANC),"Freq"])
}
for (ANC in unique(binary_derived_annotate_table$ANC)) {
  binary_derived_annotate_table[which(binary_derived_annotate_table$ANC == ANC),"Freq"] <-  binary_derived_annotate_table[which(binary_derived_annotate_table$ANC == ANC),"Freq"]/sum(binary_derived_annotate_table[which(binary_derived_annotate_table$ANC == ANC),"Freq"])
}

#Make inset plots for derived allele plots
quantitative_derived_inset <- ggplot(quantitative_derived_annotate_table, aes(y=Freq, x=ANC, alpha=Novel)) + 
  geom_bar(aes(fill=ANC), stat="identity") + 
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0%", "25%", "50%", "75%", "100%"), limits = c(0,1)) + 
  scale_fill_manual(values = c("#56B893", "#F87850", "#7B8DBF", "#DF71B6")) + 
  scale_alpha_manual(values = c(1, 0.6, 0.3, 0.1)) + 
  ylab(label = "Prior Report") + 
  theme(axis.title.x = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA), legend.title = element_blank(), 
        axis.line = element_line(colour = "black"), legend.position = "none", 
        plot.background = element_rect(fill = "transparent",colour = NA),
        axis.title.y = element_text(size = 14), axis.text.y = element_text(size = 14), axis.text.x = element_text(size = 14))
data.tb.quantitative_derived <- tibble(x = 1, y = -quantitative_beta_limit, plot = list(quantitative_derived_inset))
binary_derived_inset <- ggplot(binary_derived_annotate_table, aes(y=Freq, x=ANC, alpha=Novel)) + 
  geom_bar(aes(fill=ANC), stat="identity") + 
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0%", "25%", "50%", "75%", "100%"), limits = c(0,1)) + 
  scale_fill_manual(values = c("#56B893", "#F87850", "#7B8DBF", "#DF71B6")) + 
  scale_alpha_manual(values = c(1, 0.6, 0.3, 0.1)) + 
  ylab(label = "Prior Report") + 
  theme(axis.title.x = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA), legend.title = element_blank(), 
        axis.line = element_line(colour = "black"), legend.position = "none", 
        plot.background = element_rect(fill = "transparent",colour = NA),
        axis.title.y = element_text(size = 14), axis.text.y = element_text(size = 14), axis.text.x = element_text(size = 14))
data.tb.binary_derived <- tibble(x = 1, y = 1/binary_or_limit, plot = list(binary_derived_inset))

#Make quantitative trait plot with ancestral alleles
quant_derived_plot <- ggplot() +
  geom_point(data = finemap_signal_derived_q[which(finemap_signal_derived_q$ANC == "EUR"),], aes(x=derived_af, y=derived_beta, color=ANC, shape=Novel, alpha=Novel)) + 
  geom_point(data = finemap_signal_derived_q[which(finemap_signal_derived_q$ANC != "EUR"),], aes(x=derived_af, y=derived_beta, color=ANC, shape=Novel, alpha=Novel)) +
  xlab("Derived AF") + xlim(c(0,1)) + 
  ylab("Beta of Derived Allele") + ylim(-quantitative_beta_limit, quantitative_beta_limit) + 
  geom_plot(data = data.tb.quantitative_derived, aes(x, y, label = plot), vp.width = 0.6) +
  geom_plot(data = data.tb.legend.derived, aes(x, y, label = plot), vp.width = 0.8, vp.height = 0.20) + 
  theme_minimal() +
  scale_color_manual(values = c("#56B893", "#F87850", "#7B8DBF", "#DF71B6")) + 
  scale_alpha_manual(values = c(1, 0.6, 0.3, 0.1)) + 
  theme(legend.title=element_blank(), axis.title=element_text(size=16), legend.position = "none", 
        legend.text = element_text(size=12), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"), legend.box="vertical", axis.title.x = element_blank(), 
        axis.text = element_text(size = 18), axis.title.y = element_text(size = 18)) + 
  guides(colour = guide_legend(override.aes = list(size=6)), shape = guide_legend(override.aes = list(size=6)))
#Make binary trait plot with ancestral alleles
binary_derived_plot <- ggplot() +
  geom_point(data = finemap_signal_derived_b[which(finemap_signal_derived_b$ANC == "EUR"),], aes(x=derived_af, y=derived_or, color=ANC, shape=Novel, alpha=Novel)) + 
  geom_point(data = finemap_signal_derived_b[which(finemap_signal_derived_b$ANC != "EUR"),], aes(x=derived_af, y=derived_or, color=ANC, shape=Novel, alpha=Novel)) +
  scale_y_continuous(limits = c(1/binary_or_limit, binary_or_limit), 
                     trans = "log", breaks = 10^(seq(-3,3,1)), labels = c("0.001","0.01", "0.1", "1", "10", "100", "1,000")) + 
  xlab("Derived AF") + xlim(c(0,1)) + 
  ylab("Odds Ratio of Derived Allele") +  
  geom_plot(data = data.tb.binary_derived, aes(x, y, label = plot), vp.width = 0.6) +
  theme_minimal() +
  scale_color_manual(values = c("#56B893", "#F87850", "#7B8DBF", "#DF71B6")) + 
  scale_alpha_manual(values = c(1, 0.6, 0.3, 0.1)) + 
  theme(legend.title=element_blank(), axis.title=element_text(size=16), legend.position = "none", 
        legend.text = element_text(size=12), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"), legend.box="vertical", axis.text = element_text(size = 16),
        axis.title.x = element_text(size = 18), axis.title.y = element_text(size = 18)) + 
  guides(colour = guide_legend(override.aes = list(size=6)), shape = guide_legend(override.aes = list(size=6)))

#Make full combined derived allele plot
jpeg(paste0(work_dir, "derived_allele.effect_size_v_maf.scatter.jpeg"), width = 6, height = 11.5, units = 'in', res = 500)
ggarrange(quant_derived_plot, binary_derived_plot, ncol = 1, nrow = 2, common.legend = FALSE, align = "v", heights = c(4.25, 4.25))
dev.off()

### Make plot for lead SNPS ###
#Read in lead snps file
leads_raw <- read.csv(paste0(work_dir, "Table-S3-all.pheno.meta.gia.lead.snps.w_summary.csv"), sep = ",", header = TRUE)
#Modify novelty data
leads_raw$Novel <- ifelse(leads_raw$Locus_Novelty == "Known Association", "Known Association", 
                          ifelse(leads_raw$Locus_Novelty == "Novel Association Known Signal", "Unreported Assoc. Known Signal", "Previously Unreported Signal"))

#Append on trait type information
leads_raw$ttype <- trait_raw[leads_raw$Trait,"Trait_Type"]

#Create plotting dataframes
leads_signal_summ_q <- leads_raw %>%
  filter(ttype=="quantitative") %>%
  mutate(AF.up=MAF) %>%
  mutate(Beta.up=ifelse(EAF>0.5, -(Beta), Beta))
leads_signal_summ_b <- leads_raw %>%
  filter(ttype=="binary") %>%
  mutate(AF.up=MAF) %>%
  mutate(or.up=ifelse(EAF>0.5, exp(-(Beta)), exp(Beta)))
#Recode novelty columns as factors
leads_signal_summ_q$Novel <- factor(leads_signal_summ_q$Novel, levels = c("Previously Unreported Signal", "Unreported Assoc. Known Signal", "Known Association"))
leads_signal_summ_b$Novel <- factor(leads_signal_summ_b$Novel, levels = c("Previously Unreported Signal", "Unreported Assoc. Known Signal", "Known Association"))

#Create inset plot dataframes
quantitative_annotate_lead_table <- as.data.frame(table(leads_signal_summ_q[,c("Novel")]))
colnames(quantitative_annotate_lead_table) <- c("Novel", "Freq")
binary_annotate_lead_table <- as.data.frame(table(leads_signal_summ_b[,c("Novel")]))
colnames(binary_annotate_lead_table) <- c("Novel", "Freq")
quantitative_annotate_lead_table[,"Freq"] <-  quantitative_annotate_lead_table[,"Freq"]/sum(quantitative_annotate_lead_table[,"Freq"])
binary_annotate_lead_table[,"Freq"] <-  binary_annotate_lead_table[,"Freq"]/sum(binary_annotate_lead_table[,"Freq"])

#Make inset plots for minor allele plots
quantitative_leads_inset <- ggplot(quantitative_annotate_lead_table, aes(x=1, y=Freq, alpha=Novel), fill="#36454F") + 
  geom_bar(stat="identity") + 
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0%", "25%", "50%", "75%", "100%"), limits = c(0,1)) +
  scale_x_discrete(breaks = c(1), labels = c("META")) + 
  scale_alpha_manual(values = c(1, 0.6, 0.3, 0.1)) + 
  ylab(label = "Prior Report") + 
  theme(axis.title.x = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = "transparent",colour = NA), legend.title = element_blank(), 
        axis.line = element_line(colour = "black"), legend.position = "none", 
        plot.background = element_rect(fill = "transparent",colour = NA),
        axis.title.y = element_text(size = 14), axis.text.y = element_text(size = 14), axis.text.x = element_text(size = 14))
data.tb.quantitative_leads <- tibble(x = 0.5, y = -quantitative_beta_limit, plot = list(quantitative_leads_inset))
binary_leads_inset <- ggplot(binary_annotate_lead_table, aes(x=1, y=Freq, alpha=Novel), fill="#36454F") + 
  geom_bar(stat="identity") + 
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0%", "25%", "50%", "75%", "100%"), limits = c(0,1)) + 
  scale_x_discrete(breaks = c(1), labels = c("META")) + 
  scale_alpha_manual(values = c(1, 0.6, 0.3, 0.1)) + 
  ylab(label = "Prior Report") + 
  theme(axis.title.x = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA), legend.title = element_blank(), 
        axis.line = element_line(colour = "black"), legend.position = "none", 
        plot.background = element_rect(fill = "transparent",colour = NA),
        axis.title.y = element_text(size = 14), axis.text.y = element_text(size = 14), axis.text.x = element_text(size = 14))
data.tb.binary_leads <- tibble(x = 0.5, y = 1/binary_or_limit, plot = list(binary_leads_inset))

#Make legend plot
legend_leads_df <- cbind.data.frame(x = rep(1, 3), y = seq(1.2,1,-0.1), 
                              Novel=levels(finemap_signal_summ_b$Novel))
legend_labels_leads_df <- cbind.data.frame(x=c(rep(0, 3), 1), y=c(seq(1.2,1,-0.1), 1.33), 
                                     label=c(levels(finemap_signal_summ_b$Novel), "META"), hjust = c(rep(1,3), 0.5),
                                     angle=c(rep(0,3), rep(45,1)))
legend_leads_df$Novel <- factor(legend_leads_df$Novel, levels = levels(leads_signal_summ_b$Novel))
legend_leads_plot <- ggplot(data = legend_leads_df) + 
  geom_point(mapping = aes(x = x, y = y, alpha = Novel, shape = Novel), size = 4, color="#36454F") + 
  geom_text(data = legend_labels_leads_df, mapping = aes(x = x, y = y, label = label, hjust = hjust, angle=angle), size = 4) + 
  xlim(-15, 1.1) + 
  ylim(0.95,1.4) + 
  scale_alpha_manual(values = c(1, 0.6, 0.3, 0.1)) + 
  theme_void() + 
  theme(legend.title=element_blank(), axis.title=element_text(size=16), legend.position = "none", 
        legend.text = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_blank(), plot.background = element_blank(), axis.text = element_blank(), 
        axis.title.x = element_blank(), axis.title.y = element_blank(), axis.ticks = element_blank()) 
data.tb.legend.leads.minor <- tibble(x = 0.50, y = quantitative_beta_limit, plot = list(legend_leads_plot))

#Make lead allele plots for quantitative trait and binary trait
quant_leads_plot <- ggplot() +
  geom_point(data = leads_signal_summ_q, aes(x=AF.up, y=Beta.up, shape=Novel, alpha=Novel), color="#36454F") + 
  scale_y_continuous(limits = c(-quantitative_beta_limit, quantitative_beta_limit), breaks = seq(-3,3,1), labels = format(seq(-3,3,1), width = 5, justify = "right")) + 
  xlim(c(0, 0.5)) + ylab("Beta of Minor Allele") + 
  geom_plot(data = data.tb.quantitative_leads, aes(x, y, label = plot)) +
  geom_plot(data = data.tb.legend.leads.minor, aes(x, y, label = plot), vp.width = 0.8, vp.height = 0.20) + 
  theme_minimal() +
  scale_alpha_manual(values = c(1, 0.6, 0.3, 0.1)) + 
  theme(legend.title=element_blank(), axis.title=element_text(size=16), legend.position = "none", 
        legend.text = element_text(size=12), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"), legend.box="vertical", axis.title.x = element_blank(),
        axis.text = element_text(size = 18), axis.title.y = element_text(size = 18)) + 
  guides(colour = guide_legend(override.aes = list(size=6)), shape = guide_legend(override.aes = list(size=6))) 
binary_leads_plot <- ggplot() +
  geom_point(data = leads_signal_summ_b, aes(x=AF.up, y=or.up, shape=Novel, alpha=Novel), color="#36454F") + 
  scale_y_continuous(limits = c(1/binary_or_limit, binary_or_limit), 
                     trans = "log", breaks = 10^(seq(-3,3,1)), labels = format(c("0.001","0.01", "0.1", "1", "10", "100", "1,000"), width = 5, justify = "right")) + 
  xlab("MAF") + xlim(c(0, 0.5)) + 
  ylab("Odds Ratio of Minor Allele") +  
  geom_plot(data = data.tb.binary_leads, aes(x, y, label = plot)) + 
  theme_minimal() +
  scale_alpha_manual(values = c(1, 0.6, 0.3, 0.1)) + 
  theme(legend.title=element_blank(), axis.title=element_text(size=16), legend.position = "none", 
        legend.text = element_text(size=12), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"), legend.box="vertical", axis.text = element_text(size = 16),
        axis.title.x = element_text(size = 18), axis.title.y = element_text(size = 18)) + 
  guides(colour = guide_legend(override.aes = list(size=6)), shape = guide_legend(override.aes = list(size=6)))

#Make full combined minor allele plot
jpeg(paste0(work_dir, "minor_allele.effect_size_v_maf.leads.scatter.jpeg"), width = 6, height = 11.5, units = 'in', res = 500)
ggarrange(quant_leads_plot, binary_leads_plot, ncol = 1, nrow = 2, common.legend = FALSE, align = "v", heights = c(4.25, 4.25))
dev.off()

# 6) Calculate random Stats ====

#Get counts of ancestral alleles
high_pip_AA_table <- table(complete_high_pip_df$AA)
high_pip_AA_table
high_pip_AA_table/sum(high_pip_AA_table)
#Get counts of novel associations
high_pip_novel_table <- table(complete_high_pip_df$Novel)
high_pip_novel_table
high_pip_novel_table/sum(high_pip_novel_table)

#Make list of results for variants with high-derived allele frequencies (>0.9)
temp <- finemap_signal_derived_q[which(finemap_signal_derived_q$derived_af > 0.9),] 
temp_2 <- finemap_signal_derived_b[which(finemap_signal_derived_b$derived_af > 0.9),] 
colnames(temp) <- str_replace(colnames(temp), "derived_beta", "derived_effect_size")
colnames(temp_2) <- str_replace(colnames(temp_2), "derived_or", "derived_effect_size")
high_frequency_derived_df <- rbind.data.frame(temp, temp_2)
high_frequency_derived_df <- subset(high_frequency_derived_df, select = -c(outlier_labels))
write.table(high_frequency_derived_df, paste0(work_dir, "high_confidence.high_derived_allele_frequency.variants.tsv"), row.names = FALSE, col.names = TRUE, quote = FALSE, sep="\t")

#Calculate percentages of variants with increasing and decreasing effects
nrow(finemap_signal_summ_q[which(finemap_signal_summ_q$beta.ANC.up > 0),])/nrow(finemap_signal_summ_q) # 0.4976001
nrow(finemap_signal_derived_q[which(finemap_signal_derived_q$derived_beta > 0),])/nrow(finemap_signal_derived_q) # 0.4964292
nrow(leads_signal_summ_q[which(leads_signal_summ_q$Beta.up > 0),])/nrow(leads_signal_summ_q) #0.4741175
nrow(finemap_signal_summ_b[which(finemap_signal_summ_b$or.ANC.up > 1),])/nrow(finemap_signal_summ_b) # 0.7165862
nrow(finemap_signal_derived_b[which(finemap_signal_derived_b$derived_or > 1),])/nrow(finemap_signal_derived_b) # 0.7128447
nrow(leads_signal_summ_b[which(leads_signal_summ_b$or.up > 1),])/nrow(leads_signal_summ_b) #0.6353746

#Calculate percentage of minor alleles that are derived alleles and vice versa
separated_high_pip_ma_df <- cbind.data.frame(separated_high_pip_ea_df, ma = ifelse(separated_high_pip_ea_df$eaf.ANC < 0.5, separated_high_pip_ea_df$ea, separated_high_pip_ea_df$nea))
separated_high_pip_ma_df <- separated_high_pip_ma_df[which(separated_high_pip_ma_df$AA %in% c("A","C","G","T")),]
sum(separated_high_pip_ma_df[which(!(is.na(separated_high_pip_ma_df$AA))),"AA"] != separated_high_pip_ma_df[which(!(is.na(separated_high_pip_ma_df$AA))),"ma"])/nrow(separated_high_pip_ma_df)

#Calculate the number of high confidence (PIP > 0.95) fine-mapped variants, present only in 1 population
length(unique(complete_high_pip_df[which(rowSums(!(is.na(complete_high_pip_df[,c("AFR_high_pip", "AMR_high_pip", "EAS_high_pip", "EUR_high_pip")]))) == 1),"variant"])) #6,299 high-pip variants mapped in a single ancestry
length(unique(complete_high_pip_df[which(rowSums(!(is.na(complete_high_pip_df[,c("AFR_high_pip", "AMR_high_pip", "EAS_high_pip", "EUR_high_pip")]))) == 1 & !(is.na(complete_high_pip_df[,c("AFR_high_pip")]))),"variant"]))
length(unique(complete_high_pip_df[which(rowSums(!(is.na(complete_high_pip_df[,c("AFR_high_pip", "AMR_high_pip", "EAS_high_pip", "EUR_high_pip")]))) == 1 & !(is.na(complete_high_pip_df[,c("AMR_high_pip")]))),"variant"]))
length(unique(complete_high_pip_df[which(rowSums(!(is.na(complete_high_pip_df[,c("AFR_high_pip", "AMR_high_pip", "EAS_high_pip", "EUR_high_pip")]))) == 1 & !(is.na(complete_high_pip_df[,c("EAS_high_pip")]))),"variant"]))
length(unique(complete_high_pip_df[which(rowSums(!(is.na(complete_high_pip_df[,c("AFR_high_pip", "AMR_high_pip", "EAS_high_pip", "EUR_high_pip")]))) == 1 & !(is.na(complete_high_pip_df[,c("EUR_high_pip")]))),"variant"]))
nrow(complete_high_pip_df[which(rowSums(!(is.na(complete_high_pip_df[,c("AFR_high_pip", "AMR_high_pip", "EAS_high_pip", "EUR_high_pip")]))) == 1),]) #14,423 high-pip variant phenotype pairs mapped in a single ancestry
nrow(complete_high_pip_df[which(rowSums(!(is.na(complete_high_pip_df[,c("AFR_high_pip", "AMR_high_pip", "EAS_high_pip", "EUR_high_pip")]))) == 1 & !(is.na(complete_high_pip_df[,c("AFR_high_pip")]))),])
nrow(complete_high_pip_df[which(rowSums(!(is.na(complete_high_pip_df[,c("AFR_high_pip", "AMR_high_pip", "EAS_high_pip", "EUR_high_pip")]))) == 1 & !(is.na(complete_high_pip_df[,c("AMR_high_pip")]))),])
nrow(complete_high_pip_df[which(rowSums(!(is.na(complete_high_pip_df[,c("AFR_high_pip", "AMR_high_pip", "EAS_high_pip", "EUR_high_pip")]))) == 1 & !(is.na(complete_high_pip_df[,c("EAS_high_pip")]))),])
nrow(complete_high_pip_df[which(rowSums(!(is.na(complete_high_pip_df[,c("AFR_high_pip", "AMR_high_pip", "EAS_high_pip", "EUR_high_pip")]))) == 1 & !(is.na(complete_high_pip_df[,c("EUR_high_pip")]))),])

### Identify high-confidence variant-phenotype pairs that were fine-mapped in only one ancestry and monomorphic in others ###
#Write the list of high-confidence variants to file (Jenny needs this so we can identify which are monomorphic across ancestries)
write.table(unique(complete_high_pip_df$variant), file = paste0(work_dir, "high-confidence_variants.txt"), sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
#Read in the list of allele frequencies and imputation scores (from Jenny)
raw_imp_af_df <- read.table(paste0(work_dir, "high-confidence_variants.withR2_AF.txt"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
rownames(raw_imp_af_df) <- raw_imp_af_df$MVPID
#Bind on info about whether the GWAS had enough samples to be run from trait_raw onto complete_high_pip_df
af_high_pip_df <- cbind.data.frame(complete_high_pip_df,
                                   AFR_trait_run=ifelse(is.na(trait_raw[complete_high_pip_df$trait, "dbgap_cat.AFR"]), FALSE, TRUE),
                                   AMR_trait_run=ifelse(is.na(trait_raw[complete_high_pip_df$trait, "dbgap_cat.AMR"]), FALSE, TRUE),
                                   EAS_trait_run=ifelse(is.na(trait_raw[complete_high_pip_df$trait, "dbgap_cat.EAS"]), FALSE, TRUE),
                                   EUR_trait_run=ifelse(is.na(trait_raw[complete_high_pip_df$trait, "dbgap_cat.EUR"]), FALSE, TRUE))
#Bind on imputation R2 scores
af_high_pip_df <- cbind.data.frame(af_high_pip_df, 
                                   R2_AFR=raw_imp_af_df[af_high_pip_df$variant,"R2_AFR"],
                                   R2_AMR=raw_imp_af_df[af_high_pip_df$variant,"R2_AMR"],
                                   R2_EAS=raw_imp_af_df[af_high_pip_df$variant,"R2_EAS"],
                                   R2_EUR=raw_imp_af_df[af_high_pip_df$variant,"R2_EUR"])
#Identify variants that were not tested due to allele frequency issues
af_high_pip_df <- cbind.data.frame(af_high_pip_df,
                                   mono_AFR=ifelse(af_high_pip_df$AFR_trait_run == TRUE & af_high_pip_df$R2_AFR > 0.3 & is.na(af_high_pip_df$pval.AFR), TRUE, FALSE),
                                   mono_AMR=ifelse(af_high_pip_df$AMR_trait_run == TRUE & af_high_pip_df$R2_AMR > 0.3 & is.na(af_high_pip_df$pval.AMR), TRUE, FALSE),
                                   mono_EAS=ifelse(af_high_pip_df$EAS_trait_run == TRUE & af_high_pip_df$R2_EAS > 0.3 & is.na(af_high_pip_df$pval.EAS), TRUE, FALSE),
                                   mono_EUR=ifelse(af_high_pip_df$EUR_trait_run == TRUE & af_high_pip_df$R2_EUR > 0.3 & is.na(af_high_pip_df$pval.EUR), TRUE, FALSE))
#Identify variants that were fine-mapped in only one ancestry and monomorphic in others 
check_mono_anc <- function(anc, af_high_pip_df, ancestries, min_anc_af=0){
  if (length(anc) == 1) {
    temp <- af_high_pip_df[which(rowSums(is.na(af_high_pip_df[,c("AFR_high_pip", "AMR_high_pip", "EAS_high_pip", "EUR_high_pip")])) == 3 & 
                                   !(is.na(af_high_pip_df[,c(paste0(anc,"_high_pip"))])) &
                                   rowSums(af_high_pip_df[,paste0("mono_", ancestries[ancestries != anc])]) == 3 &
                                   (af_high_pip_df[,paste0("eaf.", anc)] > min_anc_af & af_high_pip_df[,paste0("eaf.", anc)] < 1 - min_anc_af)),]
  } else if (length(anc) == 2) {
    temp <- af_high_pip_df[which(rowSums(is.na(af_high_pip_df[,c("AFR_high_pip", "AMR_high_pip", "EAS_high_pip", "EUR_high_pip")])) == 4 - length(anc) & 
                                   rowSums(!(is.na(af_high_pip_df[,c(paste0(anc,"_high_pip"))]))) == length(anc) &
                                   rowSums(af_high_pip_df[,paste0("mono_", ancestries[!(ancestries %in% anc)])]) == 4 - length(anc) &
                                   (rowSums(af_high_pip_df[,paste0("eaf.", anc)] > min_anc_af) == length(anc) & 
                                      rowSums(af_high_pip_df[,paste0("eaf.", anc)] < 1 - min_anc_af) == length(anc))),]
  } else {
    temp <- af_high_pip_df[which(rowSums(is.na(af_high_pip_df[,c("AFR_high_pip", "AMR_high_pip", "EAS_high_pip", "EUR_high_pip")])) == 4 - length(anc) & 
                                   rowSums(!(is.na(af_high_pip_df[,c(paste0(anc,"_high_pip"))]))) == length(anc) &
                                   af_high_pip_df[,paste0("mono_", ancestries[!(ancestries %in% anc)])] == 4 - length(anc) &
                                   (rowSums(af_high_pip_df[,paste0("eaf.", anc)] > min_anc_af) == length(anc) & 
                                      rowSums(af_high_pip_df[,paste0("eaf.", anc)] < 1 - min_anc_af) == length(anc))),]
  }
  return(temp)
}
#Call function
mono_list_all <- lapply(ancestries, check_mono_anc, af_high_pip_df=af_high_pip_df, ancestries=ancestries)
names(mono_list_all) <- ancestries
mono_list_common <- lapply(ancestries, check_mono_anc, af_high_pip_df=af_high_pip_df, ancestries=ancestries, min_anc_af=0.05)
names(mono_list_common) <- ancestries

# 7) Read in variant information file and use to convert to rsids (Commented out as the work is already completed and it's intensive to run) ====

### Get ID to rsid conversion ###
##Read in file version
#annotation_raw <- read.table(gzfile(paste0(work_dir, "Annotation_gwPheWAS_dbGaP_b38.GIA.txt.gz")), header = TRUE, stringsAsFactors = FALSE)
##Create lists from the annotations to transform rsids to mvp ids
#mvp_to_rsid <- as.list(annotation_raw$RSID)
#names(mvp_to_rsid) <- annotation_raw$MVPID
#rsid_to_mvp <- as.list(annotation_raw$MVPID)
#names(rsid_to_mvp) <- annotation_raw$RSID
#rsid_to_grch38 <- as.list(paste0(annotation_raw$CHROM_b38,":",annotation_raw$POS_b38))
#names(rsid_to_grch38) <- annotation_raw$RSID
#
##Convert all variants to rsids
#finemapped_rsids <- unlist(mvp_to_rsid[all_variants])
#names(finemapped_rsids) <- NULL
##Remove the  full annotation raw file
#remove(annotation_raw)
##Save down RDSs
#saveRDS(finemapped_rsids, paste0(work_dir, "finemapped_rsids.rds"))
#saveRDS(mvp_to_rsid, paste0(work_dir, "mvp_to_rsid.rds"))
#saveRDS(rsid_to_mvp, paste0(work_dir, "rsid_to_mvp.rds"))
#saveRDS(rsid_to_grch38, paste0(work_dir, "rsid_to_grch38.rds"))
#
##Read in list of all significant variants across 4 ancestries and meta analysis
#all_assoc_list <- lapply(list.files(path=work_dir, pattern=".GIA.GWS.txt", all.files=TRUE, full.names=TRUE), read.csv, header = TRUE, sep="\t")
#all_assoc_raw <- rbind.fill(all_assoc_list)
##Create vector of all the significant variants wince we'll need to get vep annotations for those
#assoc_rsids <- unique(all_assoc_raw$SNP_ID)
#saveRDS(assoc_rsids, paste0(work_dir, "assoc_rsids.rds"))
#
#### Save down a list of the variants with PIP>0.95 in any ancestry ###
##Get 0.95% rsids
#high_pip_rsids <- mvp_to_rsid[unique(c(merged_raw[which(merged_raw$AFR.max_overall_pip >= 0.95),"AFR.best_variant"],
#                                       merged_raw[which(merged_raw$AMR.max_overall_pip >= 0.95),"AMR.best_variant"],
#                                       merged_raw[which(merged_raw$EAS.max_overall_pip >= 0.95),"EAS.best_variant"],
#                                       merged_raw[which(merged_raw$EUR.max_overall_pip >= 0.95),"EUR.best_variant"]))]
#high_pip_rsids <- unlist(high_pip_rsids)
#saveRDS(high_pip_rsids, paste0(work_dir, "high_pip_rsids.rds"))    

### Save down a list of the variants in credible sets with 1 variants ###
#one_rsids <- mvp_to_rsid[unique(merged_raw[which(merged_raw$merged.num_variants == 1),"merged.best_variant"])]
#one_rsids <- unlist(one_rsids)
#saveRDS(one_rsids, paste0(work_dir, "one_rsids.rds")) 

### Save down a list of the variants in credible sets with 2-5 variants ###
#two_to_five_rsids <- mvp_to_rsid[unique(c(merged_raw[which(merged_raw$merged.num_variants <= 5 & merged_raw$merged.num_variants >= 2),"AFR.best_variant"],
#                                       merged_raw[which(merged_raw$merged.num_variants <= 5 & merged_raw$merged.num_variants >= 2),"AMR.best_variant"],
#                                       merged_raw[which(merged_raw$merged.num_variants <= 5 & merged_raw$merged.num_variants >= 2),"EAS.best_variant"],
#                                       merged_raw[which(merged_raw$merged.num_variants <= 5 & merged_raw$merged.num_variants >= 2),"EUR.best_variant"]))]
#two_to_five_rsids <- unlist(two_to_five_rsids)
#saveRDS(two_to_five_rsids, paste0(work_dir, "two_to_five_rsids.rds")) 

### Save down a list of the variants in credible sets with 6-10 variants ###
#six_to_ten_rsids <- mvp_to_rsid[unique(c(merged_raw[which(merged_raw$merged.num_variants <= 10 & merged_raw$merged.num_variants >= 6),"AFR.best_variant"],
#                                       merged_raw[which(merged_raw$merged.num_variants <= 10 & merged_raw$merged.num_variants >= 6),"AMR.best_variant"],
#                                       merged_raw[which(merged_raw$merged.num_variants <= 10 & merged_raw$merged.num_variants >= 6),"EAS.best_variant"],
#                                       merged_raw[which(merged_raw$merged.num_variants <= 10 & merged_raw$merged.num_variants >= 6),"EUR.best_variant"]))]
#six_to_ten_rsids <- unlist(six_to_ten_rsids)
#saveRDS(six_to_ten_rsids, paste0(work_dir, "six_to_ten_rsids.rds")) 

### Get a list of the study-wide significant rsids too ###
###Extract list of rsids that are significant at study-wide threshold in META GWAS (5th GWAS)
#study_wide_rsids <- all_assoc_list[[5]][which(all_assoc_list[[5]][,"pval"] <= 4.6e-11),"SNP_ID"]
#saveRDS(study_wide_rsids, paste0(work_dir, "study_wide_rsids.rds"))

# 8) Use Ensembl Online portal to get Variant Feature Annotations from list of rsids ====

### Read in objects (Using objects generated above; run in lieu of block above after the objects are generated) ###
finemapped_rsids <- readRDS(paste0(work_dir, "finemapped_rsids.rds"))
assoc_rsids <- readRDS(paste0(work_dir, "assoc_rsids.rds"))
high_pip_rsids <- readRDS(paste0(work_dir, "high_pip_rsids.rds"))
one_rsids <- readRDS(paste0(work_dir, "one_rsids.rds"))
two_to_five_rsids <- readRDS(paste0(work_dir, "two_to_five_rsids.rds"))
six_to_ten_rsids <- readRDS(paste0(work_dir, "six_to_ten_rsids.rds"))
study_wide_rsids <- readRDS(paste0(work_dir, "study_wide_rsids.rds"))
mvp_to_rsid <- readRDS(paste0(work_dir, "mvp_to_rsid.rds"))

#Make list of all rsids
all_rsids <- unique(c(finemapped_rsids, assoc_rsids))

### Commented out the next block since it's already complete ###
##Export list of all fine-mapped rsids in chunks of 100K
#for (i in seq(0,length(all_rsids)%/%100000,1)) {
#  write.table(all_rsids[seq(1+100000*i,100000*i+100000,1)], file = paste0(work_dir, "vep_variant_ids/all_rsids.chunk", i, ".txt"), quote = FALSE, col.names = FALSE, row.names = FALSE)
#}
##Now upload results to Online portal and then download them

#Read in the VEP results
vep_list <- lapply(list.files(path=paste0(work_dir, "vep_results/"), pattern=".txt", all.files=TRUE, full.names=TRUE), read.table, header = FALSE)
#Combine VEP results
vep_raw <- rbind.fill(vep_list)
#Retain only relevant columns and add names
vep_raw <- vep_raw[,c(1,2,4)]
colnames(vep_raw) <- c("variant", "location", "consequence")

#Read in the vep categories file
vep_categories <- read.table(paste0(work_dir, "vep_categories.prioritized.txt"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
rownames(vep_categories) <- vep_categories$VEP.Annotation

#Group the variant annotations
vep_grouping_map <- vep_categories$Grouped.Annotation
names(vep_grouping_map) <- vep_categories$VEP.Annotation

#Bind grouped categories onto vep_raw
vep_raw <- cbind.data.frame(vep_raw, group_annotate=vep_grouping_map[vep_raw$consequence])
#Remove variant rows that were lacking grouped annotations (They didn't have any consequences either. They look like garbage rows)
vep_raw <- vep_raw[which(!(is.na(vep_raw$group_annotate))),]

#Read in file that converts rsids to mvp ids
rsid_to_mvp <- readRDS(paste0(work_dir, "rsid_to_mvp.rds"))
#Convert rsids in vep 
temp <- rsid_to_mvp[vep_raw$variant]
temp <- unlist(temp)
vep_raw <- cbind.data.frame(vep_raw, mvpid=temp)
remove(temp)
#Filter down the vep dataframe so that each MVPid has 1 and only 1 row (We'll retain the patched annotations )
vep_filt <- vep_raw[order(vep_raw$location),c("variant", "mvpid", "consequence", "location", "group_annotate")]
vep_filt <- vep_filt[!duplicated(vep_filt$mvpid, fromLast=TRUE),c("variant", "mvpid", "consequence", "group_annotate")] #From last will ensure that patches are retained
rownames(vep_filt) <- vep_filt$variant

#Make list of annotation priorities
vep_priority <- vep_categories$Priority
names(vep_priority) <- vep_categories$VEP.Annotation

# 9) Export lists of high-confidence variants for heterogeneity analysis ====

#Make a full coding data frame to get top "top" coding variants and write to rds object
full_variant_df <- cbind.data.frame(final, vep_cat=vep_filt[final$rsid,"consequence"], vep_group=vep_filt[final$rsid,"group_annotate"])
full_coding_df <- full_variant_df[which(vep_categories[full_variant_df$vep_cat,"Coding"] == 1),]
full_coding_df <- cbind.data.frame(full_coding_df, novelty=known_novel_raw[rownames(full_coding_df),"Novel_label"])
saveRDS(full_coding_df, file = paste0(work_dir, "full_coding_df.rds"))

#Calculate summary-level statistics for high-confidence fine-mapped variants in multiple ancestries
nrow(full_variant_df[which(rowSums(!(is.na(full_variant_df[,c("AFR_high_pip", "AMR_high_pip", "EAS_high_pip", "EUR_high_pip")]))) >= 2),]) #622 variant-phenotype pairs fine-mapped with PIP > 0.95 in 2+ ancestries
nrow(full_variant_df[which(rowSums(!(is.na(full_variant_df[,c("AFR_high_pip", "AMR_high_pip", "EAS_high_pip", "EUR_high_pip")]))) >= 2 & #241 coding pairs fine-mapped with PIP > 0.95 in 2+ ancestries
                             vep_categories[full_variant_df$vep_cat,"Coding"] == 1),])

#Export lists of high-confidence variants for heterogeneity analysis.
write.table(full_variant_df[which(rowSums(!(is.na(full_variant_df[,c("AFR_high_pip", "AMR_high_pip", "EAS_high_pip", "EUR_high_pip")]))) >= 2),], 
            file = paste0(work_dir, "high_confidence.shared_ancestry.full_data.txt"), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(full_variant_df[which(rowSums(!(is.na(full_variant_df[,c("AFR_high_pip", "AMR_high_pip", "EAS_high_pip", "EUR_high_pip")]))) >= 2 & 
                                    vep_categories[full_variant_df$vep_cat,"Coding"] == 1),], 
            file = paste0(work_dir, "high_confidence.coding.shared_ancestry.full_data.txt"), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

# 10) Make upset plot for fine-mapped variants by ancestry intersections ====

#Read in previously written merged_append file
merged_append <- read.table(file=paste0(work_dir, "master.signals.appended.txt"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)

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

#Make an upset plot of signals by ancestry
jpeg(paste0(work_dir, "cs.upset_plot.jpeg"), width = 8, height = 8, units = 'in', res = 500)
upset(
  merged_append,
  ancestries, 
  base_annotations=list(
    'Intersection size'=intersection_size(
      mapping=aes(fill=coloring), 
      text_mapping = aes(label = ifelse(!!get_size_mode('exclusive_intersection') <= 1000, !!get_size_mode('exclusive_intersection'), "")
      )
    ) + 
      scale_fill_manual(values=c("Suggestively Associated\nin Unmapped Population"='#E41A1C', "Underpowered to Detect Suggestive\nAssociation in Unmapped Population"='#0072BB')) + 
      ylab('Number of Credible Sets') + 
      scale_y_continuous(breaks = seq(10000,50000,10000), labels = paste0(seq(1,5,1), "0,000")) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.title=element_blank(),
            legend.justification=c(1,1), legend.position=c(1,1), legend.text = element_text(size=12),
            legend.key.height= unit(0.8, 'cm'), legend.key.width= unit(0.8, 'cm'), 
            axis.title.y = element_text(size=16), axis.text.y=element_text(size = 12)
            )
  ), 
  themes=upset_modify_themes(
    list('intersections_matrix'=theme(text=element_text(size=16)))
    ),
  width_ratio=0.15, 
  height_ratio=0.6, #Adjusts ratio of intersection matrix height to Intersection bar plot
  name = 'Population Intersection',
  set_sizes=(
    upset_set_size()
    + theme(axis.text.x=element_text(angle=90, size = 12)) + 
      scale_y_continuous(breaks = seq(20000,50000,20000), labels = paste0(seq(2,5,2), "0K"), trans = "reverse")
  )
)
dev.off()

#Calculate single-ancestry statistics
nrow(merged_append[which(merged_append$AFR==TRUE & merged_append$AMR==FALSE & merged_append$EAS==FALSE & merged_append$EUR==FALSE),])
nrow(merged_append[which(merged_append$AFR==FALSE & merged_append$AMR==TRUE & merged_append$EAS==FALSE & merged_append$EUR==FALSE),])
nrow(merged_append[which(merged_append$AFR==FALSE & merged_append$AMR==FALSE & merged_append$EAS==TRUE & merged_append$EUR==FALSE),])
nrow(merged_append[which(merged_append$AFR==FALSE & merged_append$AMR==FALSE & merged_append$EAS==FALSE & merged_append$EUR==TRUE),])
#Calculate portions of single-ancestry bars that are suggestive or underpowered 
nrow(merged_append[which(merged_append$AFR==TRUE & merged_append$AMR==FALSE & merged_append$EAS==FALSE & merged_append$EUR==FALSE & !(is.na(merged_append$coloring))),])
nrow(merged_append[which(merged_append$AFR==FALSE & merged_append$AMR==TRUE & merged_append$EAS==FALSE & merged_append$EUR==FALSE & !(is.na(merged_append$coloring))),])
nrow(merged_append[which(merged_append$AFR==FALSE & merged_append$AMR==FALSE & merged_append$EAS==TRUE & merged_append$EUR==FALSE & !(is.na(merged_append$coloring))),])
nrow(merged_append[which(merged_append$AFR==FALSE & merged_append$AMR==FALSE & merged_append$EAS==FALSE & merged_append$EUR==TRUE & !(is.na(merged_append$coloring))),])

# 11) Make a bar plot of VEP Annotations ====

#Make a function that gets the highest priority vep annotation for a group of variants
get_best_vep_for_cs <- function(variants, vep_priority, vep_grouping_map){
  variants <- unlist(strsplit(variants, split = ","))
  return(names(vep_priority[min(vep_priority[vep_grouping_map[vep_filt[variants,"consequence"]]])]))
}
#Make a function that returns the grouped annotations 
get_grouped_vep <- function(rsid_vector, vep_filt){
  temp <- vep_filt[rsid_vector,"group_annotate"]
  return(replace(temp, is.na(temp), "Non-Gene Unknown"))
}
#Call function that returns the grouped annotations
high_pip_veps <- as.data.frame(table(get_grouped_vep(rsid_vector = high_pip_rsids, vep_filt = vep_filt)), row.names = 1)
one_veps <- as.data.frame(table(get_grouped_vep(rsid_vector = one_rsids, vep_filt = vep_filt)), row.names = 1)
two_to_five_veps <- as.data.frame(table(get_grouped_vep(rsid_vector = two_to_five_rsids, vep_filt = vep_filt)), row.names = 1)
six_to_ten_veps <- as.data.frame(table(get_grouped_vep(rsid_vector = six_to_ten_rsids, vep_filt = vep_filt)), row.names = 1)
study_wide_veps <- as.data.frame(table(get_grouped_vep(rsid_vector = study_wide_rsids, vep_filt = vep_filt)), row.names = 1)

#Bind a dataframe together
bar_plot_df <- rbind.data.frame(data.frame(vep_category=rownames(study_wide_veps), rsid_type=rep("Study-Wide\nSignificant", length(study_wide_veps)), value=study_wide_veps$Freq),
                                data.frame(vep_category=rownames(six_to_ten_veps), rsid_type=rep("In CS with\n6-10 Variants", length(six_to_ten_veps)), value=six_to_ten_veps$Freq),
                                data.frame(vep_category=rownames(two_to_five_veps), rsid_type=rep("In CS with\n2-5 Variants", length(two_to_five_veps)), value=two_to_five_veps$Freq),
                                data.frame(vep_category=rownames(one_veps), rsid_type=rep("In CS with\n1 Variant", length(one_veps)), value=one_veps$Freq))
#Change factor levels
bar_plot_df$rsid_type <- factor(bar_plot_df$rsid_type, levels = c("Study-Wide\nSignificant", "In CS with\n6-10 Variants", "In CS with\n2-5 Variants", "In CS with\n1 Variant"))
included_vep_categories = vep_categories[which(vep_categories$Grouped.Annotation %in% unique(bar_plot_df$vep_category)),]
included_vep_categories = unique(included_vep_categories[order(included_vep_categories$Priority),"Grouped.Annotation"])
bar_plot_df$vep_category <- factor(bar_plot_df$vep_category, levels = included_vep_categories)

### Run Fisher's Exact Tests ###
#Extract mapping of grouped categories to coding binary
temp <- vep_categories[which(duplicated(vep_categories$Grouped.Annotation) == FALSE),c("Grouped.Annotation", "Coding")]
coding_map <- temp[,"Coding"]
names(coding_map) <- temp$Grouped.Annotation
remove(temp)
#Bind coding info on to vep variables
high_pip_veps <- cbind.data.frame(high_pip_veps, coding=coding_map[rownames(high_pip_veps)])
one_veps <- cbind.data.frame(one_veps, coding=coding_map[rownames(one_veps)])
two_to_five_veps <- cbind.data.frame(two_to_five_veps, coding=coding_map[rownames(two_to_five_veps)])
six_to_ten_veps <- cbind.data.frame(six_to_ten_veps, coding=coding_map[rownames(six_to_ten_veps)])
study_wide_veps <- cbind.data.frame(study_wide_veps, coding=coding_map[rownames(study_wide_veps)])
#Create function that run the fisher exact test
vep_fisher_coding <- function(first_veps, second_veps){
  first_veps_coding <- sum(first_veps$coding*first_veps$Freq)
  second_veps_coding <- sum(second_veps$coding*second_veps$Freq)
  first_veps_non_coding <- sum((1-first_veps$coding)*first_veps$Freq)
  second_veps_non_coding <- sum((1-second_veps$coding)*second_veps$Freq)
  vep_matrix <- matrix(data = c(first_veps_coding, first_veps_non_coding, second_veps_coding, second_veps_non_coding), nrow = 2, ncol = 2, byrow = TRUE)
  return(fisher.test(vep_matrix))
}
vep_fisher_annotation <- function(first_veps, second_veps, annotation_type){
  first_veps_annotate <- ifelse(annotation_type %in% rownames(first_veps), first_veps[annotation_type,"Freq"], 0)
  second_veps_annotate <- ifelse(annotation_type %in% rownames(second_veps), second_veps[annotation_type,"Freq"], 0)
  first_veps_non_annotate <- sum(first_veps$Freq) - first_veps_annotate
  second_veps_non_annotate <- sum(second_veps$Freq) - second_veps_annotate
  vep_matrix <- matrix(data = c(first_veps_annotate, first_veps_non_annotate, second_veps_annotate, second_veps_non_annotate), nrow = 2, ncol = 2, byrow = TRUE)
  return(fisher.test(vep_matrix))
}
#Execute test for all comparisons
coding_list <- list()
annotate_list <- list()
confidence_buckets <- list(high_pip_veps, one_veps, two_to_five_veps, six_to_ten_veps, study_wide_veps)
names(confidence_buckets) <- c("high_pip_veps", "one_veps", "two_to_five_veps", "six_to_ten_veps", "study_wide_veps")
for (i in 1:(length(confidence_buckets)-1)) {
  for (j in (i+1):length(confidence_buckets)) {
    coding_list[[length(coding_list) + 1]] <- vep_fisher_coding(confidence_buckets[[i]], confidence_buckets[[j]])
    names(coding_list)[length(coding_list)] <- paste0(names(confidence_buckets)[i],":",names(confidence_buckets)[j])
    #Get union of vep annotations
    union_annotations <- union(row.names(confidence_buckets[[i]]), row.names(confidence_buckets[[j]]))
    #Loop over union of vep_annotations and fill annotate_list
    for (annotation in union_annotations) {
      annotate_list[[length(annotate_list) + 1]] <- vep_fisher_annotation(confidence_buckets[[i]], confidence_buckets[[j]], annotation)
      names(annotate_list)[length(annotate_list)] <- paste0(names(confidence_buckets)[i],":",names(confidence_buckets)[j],":",annotation)
    }
  }
}

#Create function that can extract p-values from test lists
extract_pval <- function(test_result){return(test_result$p.value)}
extract_estimate <- function(test_result){return(test_result$estimate)}
create_fisher_results_df <- function(fisher_list){
  fisher_pvals <- unlist(lapply(fisher_list, extract_pval))
  fisher_estimates <- unlist(lapply(fisher_list, extract_estimate))
  temp <- str_split(names(fisher_pvals), ":")
  temp <- as.data.frame(matrix(unlist(temp), ncol = length(temp[[1]]), byrow = TRUE))
  rownames(fisher_pvals) <- NULL
  rownames(fisher_estimates) <- NULL
  fisher_pvals <- cbind.data.frame(temp, fisher_estimates, fisher_pvals)
  if (ncol(fisher_pvals) == 4) {
    colnames(fisher_pvals) <- c("category_1", "category_2", "OR", "Fisher_pval")
  } else if (ncol(fisher_pvals) == 5){
    colnames(fisher_pvals) <- c("category_1", "category_2", "annotation_type", "OR", "Fisher_pval")
  } else {
    print("ERROR: Unexpected Input")
  }
  rownames(fisher_pvals) <- NULL
  fisher_pvals <- cbind.data.frame(fisher_pvals, Bonf_pval=p.adjust(fisher_pvals$Fisher_pval, method = "bonferroni"))
  return(fisher_pvals)
}
coding_fisher_df <- create_fisher_results_df(coding_list)
annotation_fisher_df <- create_fisher_results_df(annotate_list)
#Save RDS objects
saveRDS(coding_fisher_df, paste0(work_dir, "coding_fisher_df.rds"))
saveRDS(annotation_fisher_df, paste0(work_dir, "annotation_fisher_df.rds"))

#Calculate ymins for the coding boxes
coding_ymins <- vector()
for (i in 1:length(levels(bar_plot_df$rsid_type))) {
  subset_bar_plot_df <- bar_plot_df[which(bar_plot_df$rsid_type == levels(bar_plot_df$rsid_type)[i]),]
  box_raw_min <- 1-sum(coding_map[subset_bar_plot_df$vep_category]*subset_bar_plot_df$value)/sum(subset_bar_plot_df$value)
  box_pad_min <- box_raw_min
  coding_ymins = append(coding_ymins, box_pad_min)
}

#Rescale the bar_plot_df
for (category in c("Study-Wide\nSignificant", "In CS with\n6-10 Variants", "In CS with\n2-5 Variants", "In CS with\n1 Variant")) {
  bar_plot_df[which(bar_plot_df$rsid_type == category),"value"] <-  bar_plot_df[which(bar_plot_df$rsid_type == category),"value"]/sum(bar_plot_df[which(bar_plot_df$rsid_type == category),"value"])
}

#Adjust coding_fisher_df to use as p-values for barplot
plot_fisher_df <- coding_fisher_df[c(5,7,8,10), c("category_1", "category_2", "Fisher_pval", "Bonf_pval")] #Remove categories not on the plot manually by row number
plot_fisher_df <- cbind.data.frame(c(1.05, 1.15, 1.05, 1.05), plot_fisher_df)
colnames(plot_fisher_df) <- c("y.position", "group1", "group2", "p", "p.adj")
plot_fisher_df$p <- format(plot_fisher_df$p, scientific = TRUE, digits = 2)
plot_fisher_df$p.adj <- format(plot_fisher_df$p.adj, scientific = TRUE, digits = 2)
bar_plot_label_map <- c("Study-Wide\nSignificant", "In CS with\n6-10 Variants", "In CS with\n2-5 Variants", "In CS with\n1 Variant")
names(bar_plot_label_map) <- c("study_wide_veps", "six_to_ten_veps", "two_to_five_veps", "one_veps")
plot_fisher_df$group1 <- bar_plot_label_map[plot_fisher_df$group1] 
plot_fisher_df$group2 <- bar_plot_label_map[plot_fisher_df$group2]

#Make plot
jpeg(paste0(work_dir, "vep_comparison.bar.jpeg"), width = 8.5, height = 8, units = 'in', res = 500)
ggplot(bar_plot_df, aes(y=value, x=rsid_type)) + 
  geom_bar(aes(fill=vep_category), stat="identity") + 
  annotate("rect", xmin = 0.55, xmax = 1.45, ymin = coding_ymins[1] - 0.006, ymax = 1.001, fill = "#00000000", color = "black", linetype = "solid") +  
  annotate("rect", xmin = 1.55, xmax = 2.45, ymin = coding_ymins[2] - 0.007, ymax = 1.001, fill = "#00000000", color = "black", linetype = "solid") + 
  annotate("rect", xmin = 2.55, xmax = 3.45, ymin = coding_ymins[3] - 0.008, ymax = 1.001, fill = "#00000000", color = "black", linetype = "solid") + 
  annotate("rect", xmin = 3.55, xmax = 4.45, ymin = coding_ymins[4] - 0.012, ymax = 1.001, fill = "#00000000", color = "black", linetype = "solid") + 
  stat_pvalue_manual(plot_fisher_df, label = "p", bracket.shorten = -0.1, label.size = 6) + 
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0%", "25%", "50%", "75%", "100%"), limits = c(0,1.2)) + 
  guides(fill=guide_legend(ncol = 3, byrow = TRUE)) + 
  theme(axis.title.y = element_blank(), axis.title.x = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), legend.title = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "bottom", legend.text = element_text(size = 14), axis.text = element_text(size = 16)) 
dev.off()

# 12) Make Donut plot of variant counts ====

#Discretize the variant counts
discretize_counts <- function(number){
  if (is.na(number)) {
    return("N/A")
  }else if (number == 1) {
    return("1")
  } else if(number >= 2 && number <= 5){
    return("2-5")
  } else if(number >= 6 && number <= 10){
    return("6-10")
  } else if(number >= 11 && number <= 20){
    return("11-20")
  } else if(number >= 21 && number <= 50){
    return("21-50")
  } else if(number >50){
    return(">50")
  } else{
    return("error")
  }
}
merged_append <- cbind.data.frame(merged_append,
                                  merged_discretized=factor(vapply(merged_append$merged.num_variants, discretize_counts, FUN.VALUE = character(1)), levels=c(">50", "21-50", "11-20", "6-10", "2-5", "1")),
                                  AFR_discretized=factor(vapply(merged_append$AFR.num_variants, discretize_counts, FUN.VALUE = character(1)), levels=c("N/A", ">50", "21-50", "11-20", "6-10", "2-5", "1")),
                                  AMR_discretized=factor(vapply(merged_append$AMR.num_variants, discretize_counts, FUN.VALUE = character(1)), levels=c("N/A", ">50", "21-50", "11-20", "6-10", "2-5", "1")),
                                  EAS_discretized=factor(vapply(merged_append$EAS.num_variants, discretize_counts, FUN.VALUE = character(1)), levels=c("N/A", ">50", "21-50", "11-20", "6-10", "2-5", "1")),
                                  EUR_discretized=factor(vapply(merged_append$EUR.num_variants, discretize_counts, FUN.VALUE = character(1)), levels=c("N/A", ">50", "21-50", "11-20", "6-10", "2-5", "1")))

### Make Donut Plot ###
#Create a function for making the plots
make_labeled_donut_plot <- function(annotations, plot_out_loc, plot_label){
  #Calculate fractions of VEP annotations for best variant case
  data <- as.data.frame(table(annotations))
  colnames(data) <- c("class", "count")
  # Compute percentages
  data <- data %>% 
    arrange(desc(class)) %>%
    mutate(prop = count / sum(count)) %>%
    mutate(lab.ypos = cumsum(prop) - 0.5*prop)
  #Make a donut plot of VEP annotations
  temp_plot <- ggplot(data, aes(x = 2, y = prop, fill = class)) +
    geom_bar(stat = "identity", color = "black") +
    coord_polar(theta = "y", start = 0) + # Try to remove that to understand how the chart is built initially
    geom_text(aes(y = lab.ypos, label=class), size = 10) + 
    annotate(geom = 'text', x = 0.5, y = 0, label = plot_label, size = 10) + 
    xlim(.5, 2.5) + # Try to remove that to see how to make a pie chart
    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y = element_blank(),
          axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.x = element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
          legend.position="none")
  jpeg(plot_out_loc, width = 8, height = 8, units = 'in', res = 500)
  print(temp_plot)
  dev.off()
}
##Call function
make_labeled_donut_plot(merged_append$merged_discretized, paste0(work_dir, "variant_count.donut.jpeg"), "Number of\nVariants per\nCredible Set")

# 13) Make scatter plot of Trait Categories ====

#Read in the trait dictionary
trait_raw <- read.csv(paste0(work_dir, "MVP_R4.1000G_AGR.GIA.DataDictionary.revised.022023.txt"), header = TRUE, sep = "\t")
rownames(trait_raw) <- trait_raw$Trait
trait_raw <- cbind.data.frame(trait_raw, num_samples.META=(ifelse(is.na(trait_raw$num_samples.AFR), 0, trait_raw$num_samples.AFR) + 
                                                             ifelse(is.na(trait_raw$num_samples.AMR), 0, trait_raw$num_samples.AMR) + 
                                                             ifelse(is.na(trait_raw$num_samples.EAS), 0, trait_raw$num_samples.EAS) + 
                                                             ifelse(is.na(trait_raw$num_samples.EUR), 0, trait_raw$num_samples.EUR)))

#Make a new data set for the bar plot
label_data <- as.data.frame(table(merged_raw$trait))
colnames(label_data) <- c("Phenotype", "Signals")
label_data$Phenotype <- as.character(label_data$Phenotype)
label_data <- cbind.data.frame(label_data, Category=trait_raw[label_data$Phenotype,"Category"], 
                               Trait_Type=trait_raw[label_data$Phenotype,"Trait_Type"], 
                               num_samples.META=trait_raw[label_data$Phenotype,"num_samples.META"],
                               num_cases.META=trait_raw[label_data$Phenotype,"num_cases.META"],
                               num_controls.META=trait_raw[label_data$Phenotype,"num_controls.META"])
#Calculate effective sample sizes
label_data <- cbind.data.frame(label_data, 
                               eff_samples.META=ifelse(label_data$Trait_Type == "binary", 
                                                       2/((1/label_data$num_cases.META)+(1/label_data$num_controls.META)), 
                                                       label_data$num_samples.META))
#Reorder trait categories
label_data$Category <- factor(label_data$Category, levels = c('PheCodes', 'Labs', 'Vitals', 'Baseline_Survey', 'Lifestyle_Survey'))

#Set color palette
cbPalette <- c("#009E73", "#0072B2", "#56B4E9", "#F0E442", "#E69F00", "#D55E00", "#CC79A7")

#Make scatter_plot
jpeg(paste0(work_dir, "trait_category.scatter.jpeg"), width = 8, height = 8, units = 'in', res = 500)
ggplot(label_data, aes(y=Signals, x=eff_samples.META)) + 
  geom_point(aes(color=Category)) + 
  scale_y_continuous(breaks = c(1,10,100,1000,10000), minor_breaks = rep(1:9, 4)*(10^rep(0:3, each = 9)), limits = c(0, max(label_data$Signals)), labels = c(1,10,100,"1,000","10,000"), trans = pseudolog10_trans) + 
  scale_x_continuous(breaks = c(1000,10000,100000,1000000), limits = c(1000, max(label_data$eff_samples.META)), labels = c("   1,000","   10,000","   100,000","   1,000,000"), trans = pseudolog10_trans) + 
  xlab("(Effective) Meta-Analysis Sample Size") + ylab("Number of Mapped Signals") + 
  guides(colour = guide_legend(override.aes = list(size=6))) +
  scale_color_manual(labels=c('Phecodes', 'Labs', 'Vitals', 'Baseline Survey', 'Lifestyle Survey'), values=cbPalette[1:5]) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), legend.title = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(color="black", size=16, angle=90, hjust = 1), axis.title=element_text(size=24), 
        legend.position = "bottom", axis.text.y = element_text(color="black", size=16), 
        legend.text = element_text(size=13), legend.key=element_blank()) 
dev.off()

#Calculate absolute numbers of signals in each Trait Category
table(trait_raw[merged_raw[,"trait"],"Category"])

# 14) Compare Number of Variants in EUR/AFR/AMR Shared Signals ====

#Make a function that takes 2 ancestries and outputs a comparison box plot
compare_variant_count <- function(merged_append, ANC1, ANC2, work_dir){
  #Filter down for signals found in both ancestries
  merged_ANC1_ANC2 <- merged_append[which(merged_append[,ANC1] == 1 & merged_append[,ANC2] == 1),]
  merged_ANC1_ANC2 <- cbind.data.frame('Number.Variants'=c(merged_ANC1_ANC2[,paste0(ANC1,".num_variants")], merged_ANC1_ANC2[,paste0(ANC2,".num_variants")]),
                                       Ancestry=c(rep(ANC1, nrow(merged_ANC1_ANC2)), rep(ANC2, nrow(merged_ANC1_ANC2))), 
                                       Unique.ID=c(paste0(merged_ANC1_ANC2[,"trait"],".",merged_ANC1_ANC2[,"locus"],".",merged_ANC1_ANC2[,"signal"]), paste0(merged_ANC1_ANC2[,"trait"],".",merged_ANC1_ANC2[,"locus"])))
  merged_ANC1_ANC2 <- merged_ANC1_ANC2[order(merged_ANC1_ANC2$Unique.ID),]
  ##Make comparison plot (Commented out since it doesn't illustrate the results of the test well)
  #temp<-ggboxplot(merged_ANC1_ANC2, x = "Ancestry", y = "Number.Variants",
  #                color = "Ancestry", palette = "jco", xlab = "Ancestry", ylab = "Number of Variants") + 
  #  stat_compare_means(comparisons = list(c(ANC1,ANC2)), paired = TRUE, method = "wilcox.test", label = "p.signif") + 
  #  scale_y_continuous(trans=pseudolog10_trans)
  #jpeg(paste0(work_dir, ANC1, ".", ANC2, ".variant_count_comparison.box.jpeg"), width = 8, height = 8, units = 'in', res = 500)
  #print(temp)
  #dev.off()
  #Calculate Differences
  ANC1_subset <- subset(merged_ANC1_ANC2,  Ancestry == ANC1, Number.Variants, drop = TRUE)
  ANC2_subset <- subset(merged_ANC1_ANC2,  Ancestry == ANC2, Number.Variants, drop = TRUE)
  differences <- cbind.data.frame(Unique.ID=merged_ANC1_ANC2[which(merged_ANC1_ANC2$Ancestry == ANC1),"Unique.ID"], ANC1_min_ANC2=ANC1_subset - ANC2_subset)
  res <- wilcox.test(differences$ANC1_min_ANC2)
  #Create annotation df
  annotations <- data.frame(xpos = c(-Inf), ypos =  c(Inf), annotateText = c(paste0("Wilcoxon*','~~italic(p)==", str_replace(format(res$p.value, scientific = TRUE, digits = 3), pattern = "e", replacement = "~~x~~10^"))), hjustvar = c(0), vjustvar = c(1))
  #Make difference plot
  temp <- ggplot(data = differences, mapping = aes(y = ANC1_min_ANC2)) + geom_boxplot() + scale_y_continuous(trans=pseudolog10_trans) + 
    ylab(paste0("Difference in Number of Variants (", ANC1, " - ", ANC2, ")")) + 
    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.x = element_blank(), axis.line.y = element_line(color = "black"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
          axis.text.y = element_text(size = 12), axis.title.y = element_text(size = 15)) + 
    geom_text(data=annotations,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText), parse = TRUE, size = 5.5)
  jpeg(paste0(work_dir, ANC1, ".", ANC2, ".variant_count_difference.box.jpeg"), width = 8, height = 8, units = 'in', res = 500)
  print(temp)
  dev.off()
}
#Call function
compare_variant_count(merged_append = merged_append, ANC1 = "AFR", ANC2 = "EUR", work_dir = work_dir)
compare_variant_count(merged_append = merged_append, ANC1 = "AFR", ANC2 = "AMR", work_dir = work_dir)
compare_variant_count(merged_append = merged_append, ANC1 = "AMR", ANC2 = "EUR", work_dir = work_dir)

# 15) Export Summary-Level Supplemental Table ====

##Read in appended master table and reload table of all fine-mapped gwas results (called "final")
merged_append <- read.table(paste0(work_dir, "master.signals.appended.txt"), stringsAsFactors = FALSE, header = TRUE)
load(paste0(work_dir, "all_fine-mapped_variants.withPIP.withsummarystats.withR2_Panel.RData"))
rownames(final) <- paste0(final$ID,".",final$Trait)

#Remove irrelevant columns from the table before writing
merged_export <- subset(merged_append, select = -c(AFR.suggestive_variants, AMR.suggestive_variants, EAS.suggestive_variants, EUR.suggestive_variants, max_power_unmapped_ancestries, component_signals, merged.best_variant_ancestry.specific_neglogp))

#Convert p-value columns 
merged_export$AFR.best_variant_ancestry.specific_neglogp <- final[paste0(merged_export$AFR.best_variant,".",merged_export$trait),c("pval.AFR")]
merged_export$AMR.best_variant_ancestry.specific_neglogp <- final[paste0(merged_export$AMR.best_variant,".",merged_export$trait),c("pval.AMR")]
merged_export$EAS.best_variant_ancestry.specific_neglogp <- final[paste0(merged_export$EAS.best_variant,".",merged_export$trait),c("pval.EAS")]
merged_export$EUR.best_variant_ancestry.specific_neglogp <- final[paste0(merged_export$EUR.best_variant,".",merged_export$trait),c("pval.EUR")]
colnames(merged_export) <- stringr::str_replace(colnames(merged_export), "best_variant_ancestry.specific_neglogp", "best_variant.pval.ANC")
colnames(merged_export) <- vapply(colnames(merged_export), FUN = gsub, FUN.VALUE = character(1), pattern="best_variant_ancestry.specific_neglogp", replacement="best_variant.pval.ANC")
colnames(merged_export) <- vapply(colnames(merged_export), FUN = gsub, FUN.VALUE = character(1), pattern="best_variant_beta", replacement="best_variant.beta.ANC")

#Append on trait category and description info
merged_export <- cbind.data.frame(merged_export, trait_raw[merged_export$trait,c("Category", "Description")])
colnames(merged_export) <- vapply(colnames(merged_export), FUN = gsub, FUN.VALUE = character(1), pattern="Category", replacement="category")
colnames(merged_export) <- vapply(colnames(merged_export), FUN = gsub, FUN.VALUE = character(1), pattern="Description", replacement="description")
colnames(merged_export) <- vapply(colnames(merged_export), FUN = gsub, FUN.VALUE = character(1), pattern="trait", replacement="phenotype")

#Reorder the columns
merged_export <- merged_export[,
                                 c("phenotype","category","description","chr","start","end","locus","signal",
                                 "merged.num_variants","merged.max_overall_pip","merged.best_variant",
                                 "merged.best_variant_cs_pip","merged.best_variant.beta.META","merged.best_variant.pval.META","merged.variant_ids",
                                 "best_ancestry",
                                 "AFR.num_variants","AFR.max_overall_pip","AFR.best_variant",
                                 "AFR.best_variant_cs_pip","AFR.best_variant.beta.ANC","AFR.best_variant.pval.ANC","AFR.variant_ids",
                                 "AMR.num_variants","AMR.max_overall_pip","AMR.best_variant",
                                 "AMR.best_variant_cs_pip","AMR.best_variant.beta.ANC","AMR.best_variant.pval.ANC","AMR.variant_ids",
                                 "EAS.num_variants","EAS.max_overall_pip","EAS.best_variant",
                                 "EAS.best_variant_cs_pip","EAS.best_variant.beta.ANC","EAS.best_variant.pval.ANC","EAS.variant_ids",
                                 "EUR.num_variants","EUR.max_overall_pip","EUR.best_variant",
                                 "EUR.best_variant_cs_pip","EUR.best_variant.beta.ANC","EUR.best_variant.pval.ANC","EUR.variant_ids")
                               ]

##Write to file
write.table(merged_export, file = paste0(work_dir, "master.signals.supplement.txt"), col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

# 16) Make plots for non-coding variation ====

#Make non-coding dataframe (with all variants including coding ones)
non_coding_df <- rbind.data.frame(cbind.data.frame(category= rep("In CS with\n1 Variant", length(one_rsids)), rsid=one_rsids, vep_annotation=vep_filt[one_rsids,"group_annotate"]),
                                  cbind.data.frame(category= rep("In CS with\n2-5 Variants", length(two_to_five_rsids)), rsid=two_to_five_rsids, vep_annotation=vep_filt[two_to_five_rsids,"group_annotate"]),
                                  cbind.data.frame(category= rep("In CS with\n6-10 Variants", length(six_to_ten_rsids)), rsid=six_to_ten_rsids, vep_annotation=vep_filt[six_to_ten_rsids,"group_annotate"]),
                                  cbind.data.frame(category= rep("Study-Wide\nSignificant", length(study_wide_rsids)), rsid=study_wide_rsids, vep_annotation=vep_filt[study_wide_rsids,"group_annotate"]))
#Filter out coding variants
non_coding_df <- non_coding_df[which(coding_map[non_coding_df$vep_annotation] == 0),]

#Add on GRCh38 positions
rsid_to_grch38 <- readRDS(file=paste0(work_dir, "rsid_to_grch38.rds"))
non_coding_df <- cbind.data.frame(non_coding_df, grch38=unlist(rsid_to_grch38[non_coding_df$rsid]))

#Write out a bed file format file of SNPs in GRCh38 positions
non_coding_bed <- cbind.data.frame(rbind.fill.matrix(lapply(str_split(non_coding_df$grch38,":"), t)), non_coding_df$rsid)
non_coding_bed[,1] <- paste0("chr",non_coding_bed[,1])
non_coding_bed[,2] <- as.integer(non_coding_bed[,2])
non_coding_bed[,4] <- non_coding_bed[,3]
non_coding_bed[,3] <- non_coding_bed[,2] + 1
non_coding_bed[,5] <- rownames(non_coding_bed)
write.table(non_coding_bed, paste0(work_dir, "non_coding_bed.bed"), col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

#Read in RegulomeDB Data
regulomedb_raw <- read.table(file="/project/voight_datasets_01/RegulomeDB/ENCFF250UJY.tsv", header = TRUE)
#Filter for non-coding variants in the regulomedb index and calculate loss across the categories
non_coding_df_filt <- non_coding_df[which(non_coding_df$rsid %in% regulomedb_raw$rsid),]
table(non_coding_df_filt$category)/table(non_coding_df$category)
#Append on regulomedb ranking and probability score
rownames(regulomedb_raw) <- regulomedb_raw$rsid
non_coding_df_filt <- cbind.data.frame(non_coding_df_filt, regulomedb_raw[non_coding_df_filt$rsid,c("ranking", "probability_score")])
saveRDS(non_coding_df_filt, paste0(work_dir, "non_coding_df_filt.rds"))
#non_coding_df_filt <- readRDS(paste0(work_dir, "non_coding_df_filt.rds"))

#Check allele frequencies for categories
allele_freq_raw <- fread(file = paste0(work_dir, "MVP_R4.100G_AGR.AF_R2_Overall.txt.gz"), header = TRUE)
allele_freq_raw <- allele_freq_raw[,c("MVPID", "AF")]
rownames(non_coding_df_filt) <- NULL
non_coding_df_filt$MVPID <- unlist(rsid_to_mvp[non_coding_df_filt$rsid])
non_coding_df_append <- merge(x=non_coding_df_filt, y=allele_freq_raw, by="MVPID", all.x=TRUE, all.y=FALSE)
non_coding_df_append$maf <- ifelse(non_coding_df_append[,"AF"] < 0.5, non_coding_df_append[,"AF"], 1 - non_coding_df_append[,"AF"])
non_coding_df_append <- non_coding_df_append[!(duplicated(paste0(non_coding_df_append$rsid, non_coding_df_append$category))),]
non_coding_df_append$category <- factor(non_coding_df_append$category, levels = c("Study-Wide\nSignificant", "In CS with\n6-10 Variants", "In CS with\n2-5 Variants", "In CS with\n1 Variant"))
saveRDS(non_coding_df_append, paste0(work_dir, "non_coding_df_append.rds"))

#Make functions to generate wilcoxon results
run_wilcox_tests <- function(test_field, non_coding_df_append){
  temp_list <- list()
  for (i in 1:(length(levels(non_coding_df_append$category))-1)) {
    category_i <- levels(non_coding_df_append$category)[i]
    category_j <- levels(non_coding_df_append$category)[i+1]
    temp <- wilcox.test(non_coding_df_append[which(non_coding_df_append$category == category_i),test_field], non_coding_df_append[which(non_coding_df_append$category == category_j),test_field])
    temp_list[[length(temp_list) + 1]] <- t(c(category_i, category_j , temp$p.value))
  }
  temp <- as.data.frame(rbind.fill.matrix(temp_list))
  colnames(temp) <- c("group1", "group2", "raw_p")
  temp$p <- format(p.adjust(temp$raw_p, method = "BH"), scientific = TRUE, digits = 2)
  return(temp)
}
run_wilcox_tests_base_compare <- function(test_field, non_coding_df_append){
  temp_list <- list()
  for (i in 2:(length(levels(non_coding_df_append$category)))) {
    category_i <- levels(non_coding_df_append$category)[1]
    category_j <- levels(non_coding_df_append$category)[i]
    temp <- wilcox.test(non_coding_df_append[which(non_coding_df_append$category == category_i),test_field], non_coding_df_append[which(non_coding_df_append$category == category_j),test_field])
    temp_list[[length(temp_list) + 1]] <- t(c(category_i, category_j , temp$p.value))
  }
  temp <- as.data.frame(rbind.fill.matrix(temp_list))
  colnames(temp) <- c("group1", "group2", "raw_p")
  temp$p <- format(p.adjust(temp$raw_p, method = "BH"), scientific = TRUE, digits = 2)
  return(temp)
}

#Run wilcoxon tests for the allele frequencies
maf_wilcox_df <- run_wilcox_tests("maf", non_coding_df_append = non_coding_df_append)
maf_wilcox_df$y.position <- rep(0.58, nrow(maf_wilcox_df))
#Make allele frequencies violin plot
jpeg(paste0(work_dir, "non-coding.maf.violin.jpeg"), width = 8.5, height = 8, units = 'in', res = 500)
ggplot(non_coding_df_append[which(is.na(non_coding_df_append$maf) == FALSE),], aes(x=category, y=maf)) + 
  geom_violin(trim=FALSE) + stat_summary(fun="median", geom="crossbar", width=0.2, aes(color = "black")) + 
  stat_summary(fun="mean", geom="crossbar", width=0.2, aes(color = "red")) + 
#  stat_pvalue_manual(maf_wilcox_df, label = "p", bracket.shorten = 0.1) + 
  scale_y_continuous(breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5), limits = c(0,0.6), name = "Minor Allele Frequency") + 
  scale_colour_manual(values = c("red", "black"), labels=c("mean", "median")) + 
  theme(axis.title.x = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), legend.title = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "bottom", axis.title.y = element_text(size = 15), legend.text = element_text(size = 15), 
        axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 12)) 
dev.off()

#Run wilcoxon tests for the probability scores
prob_score_wilcox_df <- run_wilcox_tests_base_compare("probability_score", non_coding_df_append = non_coding_df_append)
prob_score_wilcox_df$y.position <- c(1.05, 1.15, 1.25)
#Make probability score violin plot
jpeg(paste0(work_dir, "non-coding.prob_score.violin.jpeg"), width = 8.5, height = 8, units = 'in', res = 500)
ggplot(non_coding_df_append[which(is.na(non_coding_df_append$probability_score) == FALSE),], aes(x=category, y=probability_score)) + 
  geom_violin(trim=FALSE) + stat_summary(fun="median", geom="crossbar", width=0.2, aes(color = "black")) + 
  stat_summary(fun="mean", geom="crossbar", width=0.2, aes(color = "red")) + 
  stat_pvalue_manual(prob_score_wilcox_df, label = "p", bracket.shorten = 0.1) + 
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0,1.3)) + 
  scale_colour_manual(values = c("red", "black"), labels=c("mean", "median")) + 
  ylab("RegulomeDB Functional Probability Score") + 
  theme(axis.title.x = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), legend.title = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "bottom", axis.title.y = element_text(size = 15),  legend.text = element_text(size = 15),
        axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 12)) 
dev.off()

#Calculate median values of probability scores
median(non_coding_df_append[which(non_coding_df_append$category=="Study-Wide\nSignificant"),"probability_score"])
median(non_coding_df_append[which(non_coding_df_append$category=="In CS with\n6-10 Variants"),"probability_score"])
median(non_coding_df_append[which(non_coding_df_append$category=="In CS with\n2-5 Variants"),"probability_score"])
median(non_coding_df_append[which(non_coding_df_append$category=="In CS with\n1 Variant"),"probability_score"])

mean(non_coding_df_append[which(non_coding_df_append$category=="Study-Wide\nSignificant"),"probability_score"])
mean(non_coding_df_append[which(non_coding_df_append$category=="In CS with\n6-10 Variants"),"probability_score"])
mean(non_coding_df_append[which(non_coding_df_append$category=="In CS with\n2-5 Variants"),"probability_score"])
mean(non_coding_df_append[which(non_coding_df_append$category=="In CS with\n1 Variant"),"probability_score"])

# 17) Calculate the number of non-EUR high-confidence signals not found in EUR or found with less precision ====

#Identify the AFR/AMR/EAS only credible sets
non_eur_high_pip_df <- af_high_pip_df[which(is.na(af_high_pip_df$EUR_high_pip)),]
#Create a function that takes a row of the non_eur_high_pip_df and checks it for matches against EUR cred sets in the merged_raw_df
check_non_eur_high_pip <- function(non_eur_high_pip_row, merged_raw){
  temp <- merged_raw[which(merged_raw$trait == paste0(non_eur_high_pip_row["trait"]) & merged_raw$chr == paste0("chr",non_eur_high_pip_row["chrom"]) & !(is.na(merged_raw$EUR.num_variants))),]
  temp2 <- str_split(temp$EUR.variant_ids, pattern = ",")
  temp3 <- unlist(lapply(temp2, FUN=check_array, check_element=non_eur_high_pip_row["variant"]))
  return(ifelse(sum(temp3) > 0, temp[temp3,"EUR.num_variants"], NA))
}
check_array <- function(check_array, check_element){
  return(ifelse(check_element %in% check_array, TRUE, FALSE))
}
#Call the function
non_eur_high_pip_df <- cbind.data.frame(non_eur_high_pip_df, in_eur_cs=pbapply(non_eur_high_pip_df, FUN = check_non_eur_high_pip, MARGIN = 1, merged_raw=merged_raw))
#Append on RegulomeDB annotations
non_eur_high_pip_df <- cbind.data.frame(non_eur_high_pip_df, non_coding_df_append[as.character(non_eur_high_pip_df$ID),c("ranking", "probability_score")])
#Save the RDS object
saveRDS(non_eur_high_pip_df, file = paste0(work_dir, "non_eur_high_pip_df.rds"))

#Calculate the needed statistics for discovery gain
nrow(non_eur_high_pip_df[which(is.na(non_eur_high_pip_df$in_eur_cs)),]) #2,069 high-confidence phenotype-variant associations that were found only by using non-EUR ancestry
length(unique(non_eur_high_pip_df[which(is.na(non_eur_high_pip_df$in_eur_cs)),"variant"])) #974 unique high-confidence variants found only by using non-EUR ancestry
length(unique(non_eur_high_pip_df[which(is.na(non_eur_high_pip_df$in_eur_cs)),"trait"])) #271 phenotypes were in these

#Calculate the needed statistics for precision gain
nrow(non_eur_high_pip_df[which(!(is.na(non_eur_high_pip_df$in_eur_cs))),]) #124 phenotype-variant associations were fine-mapped with high-confidence only by using non-EUR ancestry
median(non_eur_high_pip_df[which(!(is.na(non_eur_high_pip_df$in_eur_cs))),"in_eur_cs"]) #These had a median size of 3 variants in the EUR ancestry
mean(non_eur_high_pip_df[which(!(is.na(non_eur_high_pip_df$in_eur_cs))),"in_eur_cs"]) #Average size of 7.3
length(unique(non_eur_high_pip_df[which(!(is.na(non_eur_high_pip_df$in_eur_cs))),"variant"])) #45 unique variants
length(unique(non_eur_high_pip_df[which(!(is.na(non_eur_high_pip_df$in_eur_cs))),"trait"])) #83 unique traits

non_eur_high_pip_df[which(!(is.na(non_eur_high_pip_df$in_eur_cs)) & non_eur_high_pip_df$in_eur_cs > 100),] #Identify the most extreme improvement in fine-mapping precision in non-EUR ancestry

#Calculate minor allele frequencies
check_min_1_minus <- function(value){return(min(value, 1-value))}
non_eur_high_pip_df <- cbind.data.frame(non_eur_high_pip_df, 
                                        maf.AFR=vapply(non_eur_high_pip_df$eaf.AFR, FUN = check_min_1_minus, FUN.VALUE = 0),
                                        maf.AMR=vapply(non_eur_high_pip_df$eaf.AMR, FUN = check_min_1_minus, FUN.VALUE = 0),
                                        maf.EAS=vapply(non_eur_high_pip_df$eaf.EAS, FUN = check_min_1_minus, FUN.VALUE = 0),
                                        maf.EUR=vapply(non_eur_high_pip_df$eaf.EUR, FUN = check_min_1_minus, FUN.VALUE = 0))

#Prep file for export
non_eur_high_pip_df$ID <- as.character(non_eur_high_pip_df$ID)
colnames(non_eur_high_pip_df) <- str_replace(colnames(non_eur_high_pip_df), "trait", "phenotype")
colnames(non_eur_high_pip_df) <- str_replace(colnames(non_eur_high_pip_df), "R2_", "R2.")
colnames(non_eur_high_pip_df) <- str_replace(colnames(non_eur_high_pip_df), "probability_score", "RegulomeDB score")
colnames(non_eur_high_pip_df) <- str_replace(colnames(non_eur_high_pip_df), "ranking", "RegulomeDB ranking")

non_eur_high_pip_export <- cbind.data.frame(non_eur_high_pip_df, rsid=unlist(mvp_to_rsid[non_eur_high_pip_df$ID]))
non_eur_high_pip_export <- cbind.data.frame(non_eur_high_pip_export, vep_annotation=vep_filt[non_eur_high_pip_export$rsid,"consequence"])
non_eur_high_pip_export <- cbind.data.frame(non_eur_high_pip_export, Trait_Type=trait_raw[non_eur_high_pip_export$phenotype,"Category"])
non_eur_high_pip_export <- non_eur_high_pip_export[,c("phenotype", "ID", "rsid", "chrom", "pos_b38", "Trait_Type", "Category", "phenotype_description", "vep_annotation",
                                                      "RegulomeDB score", "RegulomeDB ranking",
                                                      "AFR_high_pip", "AMR_high_pip", "EAS_high_pip","Novel", "in_eur_cs", 
                                                      "maf.AFR", "maf.AMR", "maf.EAS", "maf.EUR", 
                                                      "R2.AFR", "R2.AMR", "R2.EAS", "R2.EUR",
                                                      "AFR_phenotype_run", "AMR_phenotype_run", "EAS_phenotype_run", "EUR_phenotype_run")]
write.table(non_eur_high_pip_export, file = paste0(work_dir, "non_eur_high_pip.supplement.tsv"), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

# 18) Make a Bar Plot of Precision by Trait Category ====

#Append discretized variant counts onto each signal
merged_append <- cbind.data.frame(merged_append,
                                  merged_discretized=factor(vapply(merged_append$merged.num_variants, discretize_counts, FUN.VALUE = character(1)), levels=c(">50", "21-50", "11-20", "6-10", "2-5", "1")),
                                  AFR_discretized=factor(vapply(merged_append$AFR.num_variants, discretize_counts, FUN.VALUE = character(1)), levels=c("N/A", ">50", "21-50", "11-20", "6-10", "2-5", "1")),
                                  AMR_discretized=factor(vapply(merged_append$AMR.num_variants, discretize_counts, FUN.VALUE = character(1)), levels=c("N/A", ">50", "21-50", "11-20", "6-10", "2-5", "1")),
                                  EAS_discretized=factor(vapply(merged_append$EAS.num_variants, discretize_counts, FUN.VALUE = character(1)), levels=c("N/A", ">50", "21-50", "11-20", "6-10", "2-5", "1")),
                                  EUR_discretized=factor(vapply(merged_append$EUR.num_variants, discretize_counts, FUN.VALUE = character(1)), levels=c("N/A", ">50", "21-50", "11-20", "6-10", "2-5", "1")))
#Append on category info
merged_append <- cbind.data.frame(merged_append, Category=trait_raw[merged_append$trait,"Category"])

#Make category bar plot df
category_bar_plot_df <- table(merged_append$merged_discretized, merged_append$Category)
temp <- cbind.data.frame(Category=rep(colnames(category_bar_plot_df)[1], nrow(category_bar_plot_df)), 
                         num_variants=rownames(category_bar_plot_df), values=category_bar_plot_df[,1]/sum(category_bar_plot_df[,1]))
rownames(temp) <- NULL
for (n in 2:ncol(category_bar_plot_df)) {
  temp <- rbind.data.frame(temp,
                           cbind.data.frame(Category=rep(colnames(category_bar_plot_df)[n], nrow(category_bar_plot_df)), 
                                            num_variants=rownames(category_bar_plot_df), values=category_bar_plot_df[,n]/sum(category_bar_plot_df[,n])))
  rownames(temp) <- NULL
}
category_bar_plot_df <- temp

#Factorize the num_variants column
category_bar_plot_df$num_variants <- factor(category_bar_plot_df$num_variants, levels=c(">50", "21-50", "11-20", "6-10", "2-5", "1"))

#Make barplot
jpeg(paste0(work_dir, "precision_by_category.bar.jpeg"), width = 8.5, height = 8, units = 'in', res = 500)
ggplot(category_bar_plot_df, aes(y=values, x=Category)) + 
  geom_bar(aes(fill=num_variants), stat="identity") + 
  guides(fill = guide_legend(title = "Variants\nper Signal")) + 
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0%", "25%", "50%", "75%", "100%"), limits = c(0,1)) + 
  scale_x_discrete(breaks = c("Baseline_Survey", "Labs", "Lifestyle_Survey", "PheCodes", "Vitals"), labels=c("Baseline Survey", "Labs", "Lifestyle Survey", "PheCodes", "Vitals")) + 
  theme(axis.title.y = element_blank(), axis.title.x = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "bottom", legend.title = element_text(size = 12), legend.text = element_text(size = 12),
        axis.text = element_text(size = 12)) 
dev.off()

# 19) Check signal Sharing within Precision Buckets ====

#Read in merged append file
merged_append <- read.table(file=paste0(work_dir, "master.signals.appended.txt"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
#Add on a column of binary info for whether the signal was found in one or more ancestries
merged_append <- cbind.data.frame(merged_append, suggestive_evidence=ifelse((is.na(merged_append$AFR.suggestive_variants) == FALSE & merged_append$AFR.suggestive_variants != "Mapped") | 
                                                                              (is.na(merged_append$AMR.suggestive_variants) == FALSE & merged_append$AMR.suggestive_variants != "Mapped") | 
                                                                              (is.na(merged_append$EAS.suggestive_variants) == FALSE & merged_append$EAS.suggestive_variants != "Mapped")| 
                                                                              (is.na(merged_append$EUR.suggestive_variants) == FALSE & merged_append$EUR.suggestive_variants != "Mapped"), "Suggestively Associated\nin Unmapped Ancestry", NA))
#Rebind on a column that classifies how to color the signal
merged_append <- cbind.data.frame(merged_append, coloring=ifelse(is.na(merged_append$suggestive_evidence), ifelse(merged_append$max_power_unmapped_ancestries < 0.8, "Underpowered to Detect Suggestive\nAssociation in Unmapped Ancestry", NA), merged_append$suggestive_evidence))
#Readd on binary columns for ancestries 
merged_append <- cbind.data.frame(merged_append,
                                  AFR=ifelse(is.na(merged_append$AFR.num_variants), FALSE, TRUE),
                                  AMR=ifelse(is.na(merged_append$AMR.num_variants), FALSE, TRUE),
                                  EAS=ifelse(is.na(merged_append$EAS.num_variants), FALSE, TRUE),
                                  EUR=ifelse(is.na(merged_append$EUR.num_variants), FALSE, TRUE))
merged_append$AFR <- as.logical(merged_append$AFR)
merged_append$AMR <- as.logical(merged_append$AMR)
merged_append$EAS <- as.logical(merged_append$EAS)
merged_append$EUR <- as.logical(merged_append$EUR)
#Append discretized variant counts onto each signal
merged_append <- cbind.data.frame(merged_append,
                                  merged_discretized=factor(vapply(merged_append$merged.num_variants, discretize_counts, FUN.VALUE = character(1)), levels=c(">50", "21-50", "11-20", "6-10", "2-5", "1")),
                                  AFR_discretized=factor(vapply(merged_append$AFR.num_variants, discretize_counts, FUN.VALUE = character(1)), levels=c("N/A", ">50", "21-50", "11-20", "6-10", "2-5", "1")),
                                  AMR_discretized=factor(vapply(merged_append$AMR.num_variants, discretize_counts, FUN.VALUE = character(1)), levels=c("N/A", ">50", "21-50", "11-20", "6-10", "2-5", "1")),
                                  EAS_discretized=factor(vapply(merged_append$EAS.num_variants, discretize_counts, FUN.VALUE = character(1)), levels=c("N/A", ">50", "21-50", "11-20", "6-10", "2-5", "1")),
                                  EUR_discretized=factor(vapply(merged_append$EUR.num_variants, discretize_counts, FUN.VALUE = character(1)), levels=c("N/A", ">50", "21-50", "11-20", "6-10", "2-5", "1")))

#Make precision bar plot df
merged_append$num_ancestries = rowSums(merged_append[,ancestries]) 
precision_bar_plot_df <- table(merged_append$merged_discretized, merged_append$num_ancestries)
temp <- cbind.data.frame(num_ancestries=rep(colnames(precision_bar_plot_df)[1], nrow(precision_bar_plot_df)), 
                         num_variants=rownames(precision_bar_plot_df), values=precision_bar_plot_df[,1]/rowSums(precision_bar_plot_df))
rownames(temp) <- NULL
for (n in 2:ncol(precision_bar_plot_df)) {
  temp <- rbind.data.frame(temp,
                           cbind.data.frame(num_ancestries=rep(colnames(precision_bar_plot_df)[n], nrow(precision_bar_plot_df)), 
                                            num_variants=rownames(precision_bar_plot_df), values=precision_bar_plot_df[,n]/rowSums(precision_bar_plot_df)))
  rownames(temp) <- NULL
}
precision_bar_plot_df <- temp
#Factorize the num_variants column
precision_bar_plot_df$num_variants <- factor(precision_bar_plot_df$num_variants, levels=c("1","2-5","6-10","11-20","21-50",">50"))
#Make a stacked bar-plot by the discretized variant count categories
jpeg(paste0(work_dir, "sharing_by_precision.bar.jpeg"), width = 8.5, height = 8, units = 'in', res = 500)
ggplot(precision_bar_plot_df, aes(y=values, x=num_variants)) + 
  geom_bar(aes(fill=num_ancestries), stat="identity") + 
  guides(fill = guide_legend(title = "Ancestries\nper Signal")) + 
  scale_fill_manual(values = c("steelblue1", "steelblue2", "steelblue3", "steelblue4")) + 
  scale_x_discrete(breaks = c("1","2-5","6-10","11-20","21-50",">50")) +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0%", "25%", "50%", "75%", "100%"), limits = c(0,1)) + 
  theme(axis.title.y = element_blank(), axis.title.x = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "bottom", legend.title = element_text(size = 12), legend.text = element_text(size = 12),
        axis.text = element_text(size = 12)) 
dev.off()

### Do the same thing as above, but clump together the upper categories ###
#Rediscretize the data set
merged_append$discretized_two <- factor(ifelse(merged_append$merged_discretized %in% c("1","2-5","6-10"), as.character(merged_append$merged_discretized), ">10"),
                                        levels = c("1","2-5","6-10",">10"))
merged_append$single_ancestry <- ifelse(merged_append$num_ancestries == 1, TRUE, FALSE)
#Remake the precision bar plot dataframe with the high-count categories regrouped
precision_grouped_bar_plot_df <- table(merged_append$discretized_two, merged_append$num_ancestries)
temp <- cbind.data.frame(num_ancestries=rep(colnames(precision_grouped_bar_plot_df)[1], nrow(precision_grouped_bar_plot_df)), 
                         num_variants=rownames(precision_grouped_bar_plot_df), values=precision_grouped_bar_plot_df[,1]/rowSums(precision_grouped_bar_plot_df))
rownames(temp) <- NULL
for (n in 2:ncol(precision_grouped_bar_plot_df)) {
  temp <- rbind.data.frame(temp,
                           cbind.data.frame(num_ancestries=rep(colnames(precision_grouped_bar_plot_df)[n], nrow(precision_grouped_bar_plot_df)), 
                                            num_variants=rownames(precision_grouped_bar_plot_df), values=precision_grouped_bar_plot_df[,n]/rowSums(precision_grouped_bar_plot_df)))
  rownames(temp) <- NULL
}
precision_grouped_bar_plot_df <- temp
precision_grouped_bar_plot_df$num_variants <- factor(precision_grouped_bar_plot_df$num_variants, levels = c("1","2-5","6-10",">10"))

#Run Fisher tests for enrichment in the single-ancestry shared category
sharing_fisher_results <- list()
temp <- table(merged_append[,c("discretized_two","single_ancestry")])
for (i in 1:(length(levels(merged_append$discretized_two))-1)) {
  sharing_fisher_results <- append(sharing_fisher_results, fisher.test(temp[levels(merged_append$discretized_two)[c(i,i+1)],])$p.value)
  names(sharing_fisher_results)[length(sharing_fisher_results)] <- paste0(levels(merged_append$discretized_two)[c(i,i+1)], collapse = "_")
}
#Create fisher dataframe
precision_plot_fisher_df <- data.frame(rbind.fill.matrix(lapply(str_split(names(sharing_fisher_results), "_"), matrix, nrow=1, ncol=2)))
colnames(precision_plot_fisher_df) <- c("group1", "group2")
precision_plot_fisher_df <- cbind.data.frame(precision_plot_fisher_df, p=format(unlist(sharing_fisher_results), scientific = TRUE, digits = 2), y.position=rep(1.05, nrow(precision_plot_fisher_df)))

#Make a stacked bar-plot by the discretized variant count categories
jpeg(paste0(work_dir, "sharing_by_precision.grouped.bar.jpeg"), width = 8.5, height = 8, units = 'in', res = 500)
ggplot(precision_grouped_bar_plot_df, aes(y=values, x=num_variants)) + 
  geom_bar(aes(fill=num_ancestries), stat="identity") + 
  annotate("rect", xmin = 0.55, xmax = 1.45, ymin = 1-precision_grouped_bar_plot_df[which(precision_grouped_bar_plot_df$num_variants == "1" & precision_grouped_bar_plot_df$num_ancestries=="1"),"values"], ymax = 0, fill = "#00000000", color = "black", linetype = "solid") +  
  annotate("rect", xmin = 1.55, xmax = 2.45, ymin = 1-precision_grouped_bar_plot_df[which(precision_grouped_bar_plot_df$num_variants == "2-5" & precision_grouped_bar_plot_df$num_ancestries=="1"),"values"], ymax = 0, fill = "#00000000", color = "black", linetype = "solid") + 
  annotate("rect", xmin = 2.55, xmax = 3.45, ymin = 1-precision_grouped_bar_plot_df[which(precision_grouped_bar_plot_df$num_variants == "6-10" & precision_grouped_bar_plot_df$num_ancestries=="1"),"values"], ymax = 0, fill = "#00000000", color = "black", linetype = "solid") + 
  annotate("rect", xmin = 3.55, xmax = 4.45, ymin = 1-precision_grouped_bar_plot_df[which(precision_grouped_bar_plot_df$num_variants == ">10" & precision_grouped_bar_plot_df$num_ancestries=="1"),"values"], ymax = 0, fill = "#00000000", color = "black", linetype = "solid") + 
  guides(fill = guide_legend(title = "Ancestries\nper Signal")) + 
  scale_fill_manual(values = c("steelblue1", "steelblue2", "steelblue3", "steelblue4")) + 
  scale_x_discrete(breaks = c("1","2-5","6-10",">10")) + 
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0%", "25%", "50%", "75%", "100%"), limits = c(0,1.2)) + 
  stat_pvalue_manual(precision_plot_fisher_df, label = "p", bracket.shorten = 0.1) + 
  theme(axis.title.y = element_blank(), axis.title.x = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "bottom", legend.title = element_text(size = 12), legend.text = element_text(size = 12),
        axis.text = element_text(size = 12)) 
dev.off()

# 20) Export Variant-Level Supplemental Table ====

#Read in variant-level data
variant_level_raw <- read.table(paste0(work_dir, "master.variants.txt"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)

#Set rownames in order to extract needed data from "final" and "vep_filt"
rownames(final) <- paste(final$Trait, final$ID, sep = ".")
rownames(vep_filt) <- vep_filt$mvpid
#Extract trait, vep, and rsid data
variant_level_appended <- cbind.data.frame(variant_level_raw, final[paste(variant_level_raw$phenotype, variant_level_raw$mvp_id, sep = "."), c("rsid")],
                                           trait_raw[variant_level_raw$phenotype,c("Category", "Description")], 
                                           vep_filt[variant_level_raw$mvp_id, "consequence"])
colnames(variant_level_appended)[17:20] <- c("rsid", "category", "description", "vep_annotation")
variant_level_appended$rsid <- as.character(variant_level_appended$rsid)
#Extract GWAS result data from final
extract_variant_gwas <- function(variant_level_row, final){
  unique_id <- paste0(variant_level_row[c("phenotype","mvp_id")],collapse = ".")
  ancestry <- variant_level_row[c("ancestry")]
  temp <- final[unique_id,paste0(c("pval.", "eaf.", "beta.", "se.", "N."),ancestry)]
  names(temp) <- paste0(c("pval.", "eaf.", "beta.", "se.", "N."),"ANC") 
  return(temp)
}
temp <- pbapply(variant_level_appended, FUN = extract_variant_gwas, MARGIN = 1, final=final)
variant_level_appended <- cbind.data.frame(variant_level_appended, rbind.fill.matrix(temp))

#Tack on the novelty info specifically for the high-pip variants
variant_level_export <- cbind.data.frame(variant_level_appended, 
                                         High_PIP_Novelty=ifelse(variant_level_appended$overall_pip >= 0.95, complete_high_pip_df[paste0(variant_level_appended$mvp_id, ".", variant_level_appended$phenotype),"Novel"], NA))
variant_level_export$High_PIP_Novelty <- ifelse(is.na(variant_level_export$High_PIP_Novelty), "-", 
                                                ifelse(variant_level_export$High_PIP_Novelty == "Novel Signal", "Novel SNP to GWAS Catalog", 
                                                       ifelse(variant_level_export$High_PIP_Novelty=="Novel Association Known Signal", "Novel Association", variant_level_export$High_PIP_Novelty)))

#Reorganize columns
variant_level_export <- variant_level_export[,c("phenotype", "category", "description",
                                                "mvp_id","rsid","bp","bp38","vep_annotation",
                                                "chr","locus","merged_signal","ancestry","anc_signal",
                                                "eaf.ANC","beta.ANC","se.ANC","pval.ANC","N.ANC",
                                                "overall_pip","cs_pip","mu","mu2","cs_log_bayes_factor", "High_PIP_Novelty")]
rownames(variant_level_export) <- NULL

#Reset merged signal numbers
check_pheno_locus <- function(pheno_locus, variant_level_export){
  temp2 <- variant_level_export[which(paste0(variant_level_export$phenotype, ".", variant_level_export$locus) == pheno_locus),]
  bad_signals <- unique(temp2$merged_signal)
  signal_map <- seq(1, length(bad_signals), 1)
  names(signal_map) <- as.character(bad_signals)
  temp2$merged_signal <- signal_map[as.character(temp2$merged_signal)]
  return(temp2)
}
variant_level_export_list <- lapply(unique(paste0(variant_level_export$phenotype, ".", variant_level_export$locus)), FUN = check_pheno_locus, variant_level_export=variant_level_export)
variant_level_export_clean <- rbind.fill(variant_level_export_list)

#Export the file
write.table(variant_level_export_clean, file = paste0(work_dir, "master.variants.supplement.txt"), col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

# 21) Append Fine-Mapping Results onto the Trait Overview Supplementary Table (Commented out as it's now implemented elsewhere) ====

##Read in file of lead SNPs and associations summarised by trait
#trait_overview_raw <- read.csv(paste0(work_dir, "trait_associations_summary_S4.tsv"), sep = "\t", header = TRUE)
#
##Append on number of signals and variants in credible sets
#count_signals <- function(trait_overview_row, merged_raw){
#  trait <- trait_overview_row["Trait"]
#  temp <- merged_raw[which(merged_raw$trait == trait),]
#  num_signals <- nrow(temp)
#  num_variants <- length(unique(unlist(str_split(temp$merged.variant_ids, ","))))
#  return(c(num_signals, num_variants))
#}
#temp <- pbapply(trait_overview_raw, MARGIN = 1, FUN = count_signals, merged_raw=merged_raw)
#temp <- rbind.fill.matrix(temp)
#temp <- t(temp)
#colnames(temp) <- c("Fine-Mapped.Signals", "Variants.in.95%.CS")
#trait_overview_append <- cbind.data.frame(trait_overview_raw, temp)
#
##Append on number of high-PIP variants with novelty info
#count_high_pip_variants <- function(trait_overview_row, complete_high_pip_df){
#  trait <- trait_overview_row["Trait"]
#  temp <- complete_high_pip_df[which(complete_high_pip_df$trait == trait),]
#  num_high_pip_variants <- nrow(temp)
#  num_known_variants <- sum(temp$Novel == "Known Association")
#  num_novel_variants <- sum(temp$Novel %in% c("Novel Association Known Signal", "Novel Signal"))
#  num_novel_variants_unknown_gwas <- sum(temp$Novel %in% c("Novel Signal"))
#  return(c(num_high_pip_variants, num_known_variants, num_novel_variants,num_novel_variants_unknown_gwas))
#}
#temp <- pbapply(trait_overview_append, MARGIN = 1, FUN = count_high_pip_variants, complete_high_pip_df=complete_high_pip_df)
#temp  <- rbind.fill.matrix(temp)
#temp <- t(temp)
#colnames(temp) <- c("High-PiP.Variants", "Known.High-PIP.Variants.in.GWAS.Catalog", "Novel.High-PIP.Variant.Associations", "Novel.High-PIP.Variants.not.found.in.GWAS.Catalog")
#trait_overview_append <- cbind.data.frame(trait_overview_append, temp)

#Write to file
#write.table(trait_overview_append, file = paste0(work_dir, "trait_associations_summary_S4_w_fine-map.tsv"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


# 22) Make a scatter plot of credible set sizes by effective sample sizes ====

#Read in the trait dictionary
trait_raw <- read.csv(paste0(work_dir, "MVP_R4.1000G_AGR.GIA.DataDictionary.revised.022023.txt"), header = TRUE, sep = "\t")
rownames(trait_raw) <- trait_raw$Trait

#Make an ancestry separated data frame of signals
anc_sep_num_var_df <- cbind.data.frame(rbind.data.frame(merged_raw[,c("trait", "locus", "signal")],merged_raw[,c("trait", "locus", "signal")],
                                       merged_raw[,c("trait", "locus", "signal")],merged_raw[,c("trait", "locus", "signal")]),
                                       ANC=(c(rep("AFR",nrow(merged_raw)), rep("AMR",nrow(merged_raw)), rep("EAS",nrow(merged_raw)), rep("EUR",nrow(merged_raw)))),
                                       ANC.num_variants=(c(merged_raw[,c("AFR.num_variants")],merged_raw[,c("AMR.num_variants")],
                                                                        merged_raw[,c("EAS.num_variants")],merged_raw[,c("EUR.num_variants")])))
#Clear out missing data
anc_sep_num_var_df <- anc_sep_num_var_df[which(is.na(anc_sep_num_var_df$ANC.num_variants) == FALSE),]
#Append on sample sizes
get_eff_sample_size <- function(anc_sep_num_var_row, trait_raw){
  ANC <- as.character(anc_sep_num_var_row["ANC"])
  trait <- as.character(anc_sep_num_var_row["trait"])
  trait_type <- trait_raw[trait,"Trait_Type"]
  trait_sizes <- trait_raw[trait,c(paste0("num_samples.", ANC), paste0("num_cases.", ANC), paste0("num_controls.", ANC))]
  eff_size <- ifelse(trait_type == "binary", 2/((1/trait_sizes[2])+(1/trait_sizes[3])), ifelse(trait_type == "quantitative", trait_sizes[1], NA))
  return(eff_size)
}
temp <- unlist(pbapply(anc_sep_num_var_df, FUN = get_eff_sample_size, MARGIN = 1, trait_raw=trait_raw))
anc_sep_num_var_df <- cbind.data.frame(anc_sep_num_var_df, eff_sample_size=temp)

#Factorize ancestry
anc_sep_num_var_df$ANC <- factor(anc_sep_num_var_df$ANC, levels = ancestries)

#Make Plot
jpeg(paste0(work_dir, "cs_size_v_sample_size.scatter.jpeg"), width = 8, height = 8, units = 'in', res = 500)
ggplot() +
  geom_point(data = anc_sep_num_var_df[which(anc_sep_num_var_df$ANC == "EUR"),], aes(x=eff_sample_size, y=ANC.num_variants, color=ANC)) + 
  geom_point(data = anc_sep_num_var_df[which(anc_sep_num_var_df$ANC != "EUR"),], aes(x=eff_sample_size, y=ANC.num_variants, color=ANC)) + 
  scale_x_continuous(breaks = c(1000,10000,100000,1000000), limits = c(min(anc_sep_num_var_df$eff_sample_size), max(anc_sep_num_var_df$eff_sample_size)), 
                     labels = c("1,000","10,000","100,000","1,000,000"), trans = pseudolog10_trans) + 
  ylab("Credible Set Size") + xlab("Effective Sample Size") + 
  theme_minimal() +
  scale_color_manual(values = c("#56B893", "#F87850", "#7B8DBF", "#DF71B6")) + 
  theme(legend.title=element_blank(), axis.title=element_text(size=18), legend.position = "bottom",
        legend.text = element_text(size=18), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"), legend.box="vertical", axis.title.x = element_text(size = 18), 
        axis.text = element_text(size = 18), axis.title.y = element_text(size = 18)) + 
  guides(colour = guide_legend(override.aes = list(size=6)), shape = guide_legend(override.aes = list(size=6))) 
dev.off()

#Run a mean cs by trait version as well / start with binding data onto trait dataframe
get_mean_trait_cs_sizes <- function(trait, anc_sep_num_var_df){
  temp <- anc_sep_num_var_df[which(anc_sep_num_var_df$trait == trait),]
  temp2 <- c(mean(temp[which(temp$ANC == "AFR"),"ANC.num_variants"]),mean(temp[which(temp$ANC == "AMR"),"ANC.num_variants"]),
             mean(temp[which(temp$ANC == "EAS"),"ANC.num_variants"]), mean(temp[which(temp$ANC == "EUR"),"ANC.num_variants"]))
  return(temp2)
}
temp <- vapply(rownames(trait_raw), FUN = get_mean_trait_cs_sizes, FUN.VALUE = numeric(4), anc_sep_num_var_df=anc_sep_num_var_df)
temp <- t(temp)
colnames(temp) <- paste0("mean_cs_size.", ancestries)
trait_append <- cbind.data.frame(trait_raw, temp)
#Calculate effective sample sizes on trait dataframe
trait_append <- cbind.data.frame(trait_append, 
                                 num_eff.AFR=ifelse(trait_append$Trait_Type=="binary", 2/((1/trait_append$num_controls.AFR)+(1/trait_append$num_cases.AFR)), trait_append$num_samples.AFR),
                                 num_eff.AMR=ifelse(trait_append$Trait_Type=="binary", 2/((1/trait_append$num_controls.AMR)+(1/trait_append$num_cases.AMR)), trait_append$num_samples.AMR),
                                 num_eff.EAS=ifelse(trait_append$Trait_Type=="binary", 2/((1/trait_append$num_controls.EAS)+(1/trait_append$num_cases.EAS)), trait_append$num_samples.EAS),
                                 num_eff.EUR=ifelse(trait_append$Trait_Type=="binary", 2/((1/trait_append$num_controls.EUR)+(1/trait_append$num_cases.EUR)), trait_append$num_samples.EUR))
#Convert trait df to a plottable format
trait_plot_df <- cbind.data.frame(trait=c(rownames(trait_append), rownames(trait_append), rownames(trait_append), rownames(trait_append)),
                                  ANC=c(rep("AFR", nrow(trait_append)), rep("AMR", nrow(trait_append)), rep("EAS", nrow(trait_append)), rep("EUR", nrow(trait_append))),
                                  num_eff=c(trait_append$num_eff.AFR, trait_append$num_eff.AMR, trait_append$num_eff.EAS, trait_append$num_eff.EUR),
                                  mean_cs_size=c(trait_append$mean_cs_size.AFR, trait_append$mean_cs_size.AMR, trait_append$mean_cs_size.EAS, trait_append$mean_cs_size.EUR))
#Remove Missing rows
trait_plot_df <- trait_plot_df[which(!(is.na(trait_plot_df$num_eff) | is.na(trait_plot_df$mean_cs_size))),]

#Make Plot
jpeg(paste0(work_dir, "mean_cs_size_v_sample_size.scatter.jpeg"), width = 8, height = 8, units = 'in', res = 500)
ggplot() +
  geom_point(data = trait_plot_df[which(trait_plot_df$ANC == "EUR"),], aes(x=num_eff, y=mean_cs_size, color=ANC)) + 
  geom_point(data = trait_plot_df[which(trait_plot_df$ANC != "EUR"),], aes(x=num_eff, y=mean_cs_size, color=ANC)) + 
  scale_x_continuous(breaks = c(1000,10000,100000,1000000), limits = c(min(trait_plot_df$num_eff), max(trait_plot_df$num_eff)), 
                     labels = c("1,000","10,000","100,000","1,000,000"), trans = pseudolog10_trans) + 
  scale_y_continuous(breaks = c(1,10,100,1000), limits = c(min(trait_plot_df$mean_cs_size), max(max(trait_plot_df$mean_cs_size), 1000)), 
                     labels = c("1","10","100","1,000"), trans = pseudolog10_trans) + 
  ylab("Mean Credible Set Size") + xlab("Effective Sample Size") + 
  theme_minimal() +
  scale_color_manual(values = c("#56B893", "#F87850", "#7B8DBF", "#DF71B6")) + 
  theme(legend.title=element_blank(), axis.title=element_text(size=18), legend.position = "bottom",
        legend.text = element_text(size=18), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"), legend.box="vertical", axis.title.x = element_text(size = 18), 
        axis.text = element_text(size = 18), axis.title.y = element_text(size = 18)) + 
  guides(colour = guide_legend(override.aes = list(size=6)), shape = guide_legend(override.aes = list(size=6))) 
dev.off()
