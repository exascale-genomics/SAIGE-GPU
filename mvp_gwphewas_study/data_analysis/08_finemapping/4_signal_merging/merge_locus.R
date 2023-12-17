################################################################################

#merge_locus.R

#This script is designed to merge the credible sets for a locus across 
#ancestries. It takes in 4 RDS files corresponding to AFR, AMR, EAS, and EUR 
#ancestries. It merges the credible sets using Jaccard indices and outputs a 
#text file of results and a tree diagram showing the relationship between 
#credible sets at the locus.

################################################################################

# 0) Load Needed Libraries ====

library(igraph)
library(stringr)
library(plyr)

# 1) Read in Needed Arguments from Command ====

#Need to read in the following (Either input all or no optional parameters): 
  # 1) AFR SuSiE RDS File
  # 2) AMR SuSiE RDS File
  # 3) EAS SuSiE RDS File
  # 4) EUR SuSiE RDS File
  # 5) output file folder
  # 6) significance threshold OPTIONAL (default is 5e-8)
  # 7) significance type OPTIONAL (gwas, absolute_residual, or preceding_residual / default is gwas)
  # 8) purity threshold OPTIONAL (default is 0.1)
  # 9) cutting height OPTIONAL (default is 0.9)
  # 10) rds naming file convention OPTIONAL (default is TRAIT.CHR.START.END.ANC.rds)
  # 11) random seed OPTIONAL (default is 5)

args <- commandArgs(trailingOnly = TRUE);
if (length(args) == 5){
  #Print confirmation statement
  print("Necessary Parameters Input")
  #Set variables from the arguments in this case
  afr_rds_file <- args[1]
  amr_rds_file <- args[2]
  eas_rds_file <- args[3]
  eur_rds_file <- args[4]
  out_folder <- args[5]
  significance_thresh <- 5e-8
  significance_type <- "gwas"
  purity_thresh <- 0.1
  cut_height <- 0.9
  name_conv <- "PHENO.chr.start.end.ANC.rds"
  random_seed <- 5
} else if (length(args) == 11) {
  #Print confirmation statement
  print("Necessary and Optional Parameters Input")
  #Set variables from the arguments in this case
  afr_rds_file <- args[1]
  amr_rds_file <- args[2]
  eas_rds_file <- args[3]
  eur_rds_file <- args[4]
  out_folder <- args[5]
  significance_thresh <- as.numeric(args[6])
  significance_type <- args[7]
  purity_thresh <- as.numeric(args[8])
  cut_height<- as.numeric(args[9])
  name_conv <- args[10]
  random_seed <- as.numeric(args[11])
} else {
  print("ERROR: Usage: %> Rscript merge_locus.R afr_rds_file amr_rds_file eas_rds_file eur_rds_file out_folder [significance_thresh] [significance_type] [purity_thresh] [cut_height] [name_conv] [random_seed]");
  quit(save="no");
}

# 2) Verify Inputs ====

#Set ancestries
inp_ancestries <- c("AFR", "AMR", "EAS", "EUR")
inp_ancestry_files <- c(afr_rds_file, amr_rds_file, eas_rds_file, eur_rds_file)
names(inp_ancestry_files) <- inp_ancestries

#Remove any ancestries/files that were determined to have missing rds files and were given NULL locations
ancestries <- vector()
for (ancestry in inp_ancestries) {
  if (inp_ancestry_files[ancestry] != "NA" && !(is.na(inp_ancestry_files[ancestry]))) {
    ancestries <- append(ancestries, ancestry)
  }
}
ancestry_files <- inp_ancestry_files[ancestries]
names(ancestry_files) <- ancestries

#Verify that all rds files have same location and fit the naming convention (Also set the out file name)
name_split = strsplit(name_conv, "[.]")
check_rds_file_name <- function(anc_rds_file, name_split){
  anc_base = strsplit(basename(anc_rds_file), "[.]")[[1]]
  #Swap the ancestry from the base name list and return it
  anc_base[which(name_split[[1]] == "ANC")] = "MERGE"
  return(paste0(anc_base, collapse = "."))
}
out_file <- vapply(ancestry_files, check_rds_file_name, FUN.VALUE = character(1), name_split=name_split)
if (length(unique(out_file)) != 1) {
  print("ERROR: RDS file names either do not match each other or the naming convention.")
  quit(save = "no")
} else {
  out_file = out_file[1]
}
#Get name of trait locus
trait_locus_name <- str_replace(string = out_file,pattern = ".MERGE.rds", replacement = "")
trait <- str_split(trait_locus_name, "[.]")[[1]][1]
chrom <- str_split(trait_locus_name, "[.]")[[1]][2]
start_pos <- str_split(trait_locus_name, "[.]")[[1]][3]
end_pos <- str_split(trait_locus_name, "[.]")[[1]][4]

#Verify that RDS file locations exist
miss_rds_ind = FALSE
for (ancestry in names(ancestry_files)) {
  if (file.exists(ancestry_files[ancestry]) == FALSE) {
    print(paste0("ERROR: Invalid ", ancestry," RDS file location given."))
    miss_rds_ind = TRUE
  }
}
if (miss_rds_ind == TRUE) {
  quit(save="no");
}

#Verify that valid significance_type was given
if (!(significance_type %in% c("gwas", "absolute_residual", "preceding_residual"))) {
  print("ERROR: Invalid significance type given.")
  quit(save="no");
}

#Define out_loc
if(substr(out_folder, nchar(out_folder), nchar(out_folder)) != "/"){
  out_loc <- paste0(out_folder, "/")
} else {
  out_loc <- out_folder
}
if (dir.exists(out_loc) == FALSE) {
  print("ERROR: Invalid output folder given.")
  quit(save="no");
}
df_out_loc <- paste0(out_loc, out_file)
df_out_loc <- str_replace(df_out_loc, pattern = ".rds", replacement = ".txt")
plot_out_loc <- str_replace(df_out_loc, pattern = ".txt", replacement = ".jpg")

# 3) Read in and merge RDS Files ====

#Read in files
raw_rds <- lapply(ancestry_files, readRDS)
names(raw_rds) = ancestries

#Remove any ancestries with empty files
finish_ancestries <- vector()
for (ancestry in ancestries) {
  if (is.null(raw_rds[[ancestry]]) == FALSE) {
    finish_ancestries <- append(finish_ancestries, ancestry)
  }
}
#Quit if all files are empty
if (length(finish_ancestries) == 0) {
  print(paste0("WARNING: All ancestries failed to converge for ", trait_locus_name))
  write.table(NULL, file = df_out_loc, quote = FALSE, row.names = FALSE, col.names = FALSE) #Write an empty file for accounting purposes
  quit(save = "no");
}

### Verify that each of the 4 ancestries appears to be complete ###
bad_write_ind = FALSE
for (ancestry in finish_ancestries) {
  variant_count = ncol(raw_rds[[ancestry]][["alpha"]])
  if (ncol(raw_rds[[ancestry]][["mu"]]) != variant_count || ncol(raw_rds[[ancestry]][["mu2"]]) != variant_count || ncol(raw_rds[[ancestry]][["lbf_variable"]]) != variant_count || length(raw_rds[[ancestry]][["pi"]]) != variant_count || length(raw_rds[[ancestry]][["XtXr"]]) != variant_count || length(raw_rds[[ancestry]][["pip"]]) != variant_count || length(raw_rds[[ancestry]][["bp"]]) != variant_count || length(raw_rds[[ancestry]][["neglogp_gwas"]]) != variant_count || length(raw_rds[[ancestry]][["bp38"]]) != variant_count) {
    print(paste0("ERROR: ", ancestry, " results appear to have written with errors for ", trait_locus_name))
    bad_write_ind = TRUE
  }
}
if (bad_write_ind == TRUE) {
  quit(save = "no");
}

### Verify that none of the 4 ancestries hit the maximum value threshold for p-values and absolute residual p-values ###
min_p_ind = FALSE
for (ancestry in finish_ancestries) {
  if (suppressWarnings(max(raw_rds[[ancestry]][["neglogp_gwas"]])) == Inf || suppressWarnings(max(raw_rds[[ancestry]][["neglogp_absolute"]])) == Inf){
    print(paste0("WARNING: Minimum p-value threshold breached for ", trait_locus_name, ".", ancestry))
    min_p_ind = TRUE
  }
}
if (min_p_ind == TRUE) {
  print(paste0("ERROR: Removed ", trait_locus_name, " from temporary consideration due to minimum p-value threshold"))
  quit(save="no")
}

### Create a table of all credible set combinations from the 4 ancestries ###
#Create functions that take pairwise minimums and maximums of elements in lists and another that uses them to calculate jaccard distances
take_min <- function(index, cs_pip, cs2_pip){
  return(min(cs_pip[index], cs2_pip[index]))
}
take_max <- function(index, cs_pip, cs2_pip){
  return(max(cs_pip[index], cs2_pip[index]))
}
calc_jaccard_dist <- function(cs2_pip, cs_pip, SNPs_intersect){
  return(1 - (sum(vapply(SNPs_intersect, take_min, FUN.VALUE = 0, cs_pip=cs_pip, cs2_pip=cs2_pip))/sum(vapply(SNPs_intersect, take_max, FUN.VALUE = 0, cs_pip=cs_pip, cs2_pip=cs2_pip))))
}

#Create empty vectors to hold data
jaccard_dists <- vector()
signals_1 <- vector()
signals_2 <-vector()
#Loop over the ancestries twice
for (ancestry in finish_ancestries) {
  for (ancestry2 in finish_ancestries) {
    #Next loop over the credible sets
    for (cs in raw_rds[[ancestry]][["sets"]][["cs_index"]]) {
      #Verify that there is at least one cs
      if (is.null(cs)) {
        print(paste0("WARNING: No Signals detected for ", trait_locus_name, ".", ancestry))
        next
      }
      for (cs2 in raw_rds[[ancestry2]][["sets"]][["cs_index"]]) {
        #Verify that there is at least one cs
        if (is.null(cs2)) {
          print(paste0("WARNING: No Signals detected for ", trait_locus_name, ".", ancestry2))
          next
        }
        #Extract PIPs
        cs_pip = raw_rds[[ancestry]][["alpha"]][cs,]
        cs2_pip = raw_rds[[ancestry2]][["alpha"]][cs2,]
        #Get names of SNPs
        SNPs = names(cs_pip)
        SNPs2 = names(cs2_pip)
        SNPs_intersect = intersect(SNPs, SNPs2)
        #Calculate jaccard
        jaccard_dist = calc_jaccard_dist(cs2_pip, cs_pip, SNPs_intersect)
        #Append values to vectors
        jaccard_dists <- append(jaccard_dists, jaccard_dist)
        signals_1 <- append(signals_1, paste(ancestry, cs, sep = "."))
        signals_2 <- append(signals_2, paste(ancestry2, cs2, sep = "."))
      }
    }
  }
}
jaccard_df <- cbind.data.frame(signal_1=signals_1, signal_2=signals_2, jaccard_dist=jaccard_dists)

### Build distance matrix and make tree ###
#First check if more than 1 credible set was found across all ancestries (Need at least two to do clustering)
if (nrow(jaccard_df) > 1) {
  #Make matrix
  jaccard_G <- graph.data.frame(jaccard_df,directed=FALSE)
  jaccard_A <- as_adjacency_matrix(jaccard_G,names=TRUE,sparse=FALSE,attr="jaccard_dist",type='lower')
  dist_mat <- as.dist(jaccard_A)
  #Cluster the signals
  set.seed(random_seed)
  hclust_signals <- hclust(dist_mat, method = 'complete')
  #Cut the tree
  cut_signals <- cutree(hclust_signals, h=cut_height)
  #Get unique merged cs
  unique_signals <- unique(cut_signals)
  
  #Output the merging tree
  jpeg(plot_out_loc, width = 8, height = 8, units = 'in', res = 500)
  plot(as.dendrogram(hclust_signals))
  if (min(hclust_signals$height) <= cut_height) {
    if(length(dist_mat) > 1){
      rect.hclust(hclust_signals , h = cut_height, border = 2:6) 
    }
  }
  abline(h = cut_height, col = 'red')
  dev.off()
} else {
  #Spoof the cut signal results as if more than one credible set was found
  unique_signals <- 1
  cut_signals <- 1
  names(cut_signals) <- jaccard_df[1,"signal_1"]
}

### Build functions for extracting data for analysis ###
#Make a function that extracts first element from a list item
extract_pos <- function(list_array, pos){return(list_array[pos])}
#Make a function that extracts variant-level info from the RDS object for each vetted signal
extract_variant_info <- function(vetted_signal, raw_rds, merged_signal){
  signal_split <- strsplit(vetted_signal, "[.]")
  ancestry <- signal_split[[1]][1]
  row_index <- as.numeric(signal_split[[1]][2])
  variant_indices <- raw_rds[[ancestry]][["sets"]][["cs"]][[paste0("L",signal_split[[1]][2])]]
  all_variants <- names(raw_rds[[ancestry]][["sets"]][["cs"]][[paste0("L",signal_split[[1]][2])]])
  overall_pips <- raw_rds[[ancestry]][["pip"]][all_variants]
  cs_pips <- raw_rds[[ancestry]][["alpha"]][row_index,all_variants]
  mus <- raw_rds[[ancestry]][["mu"]][row_index,all_variants]
  mu2s<- raw_rds[[ancestry]][["mu2"]][row_index,all_variants]
  lbf_variables<- raw_rds[[ancestry]][["lbf_variable"]][row_index,all_variants]
  bps<- raw_rds[[ancestry]][["bp"]][variant_indices]
  bp38s <- raw_rds[[ancestry]][["bp38"]][variant_indices]
  neglogp_gwas <- raw_rds[[ancestry]][["neglogp_gwas"]][variant_indices]
  neglogp_absolutes <- raw_rds[[ancestry]][["neglogp_absolute"]][all_variants,as.character(row_index)]
  temp <- cbind.data.frame(ancestry=rep(ancestry, length(all_variants)), anc_signal=rep(row_index, length(all_variants)), 
                           merged_signal=rep(merged_signal, length(all_variants)), mvp_id=all_variants, chr=rep(raw_rds[[ancestry]][["chr"]], length(all_variants)), 
                           bp=bps, bp38=bp38s, overall_pip=overall_pips, cs_pip=cs_pips, mu=mus, mu2=mu2s, cs_log_bayes_factor=lbf_variables, 
                           neglogp_gwas=neglogp_gwas, neglogp_absolute=neglogp_absolutes)
  return(temp)
}
#Make a function to extract info from the RDS for each vetted signal
extract_signal_info <- function(vetted_signal, variant_level_df){
  signal_split <- strsplit(vetted_signal, "[.]")
  temp_df <- variant_level_df[which(variant_level_df$ancestry == signal_split[[1]][1] & variant_level_df$anc_signal == signal_split[[1]][2]),]
  all_variants <- temp_df$mvp_id
  max_overall_pip <- max(temp_df$overall_pip)
  best_variant <- all_variants[which.max(temp_df$overall_pip)]
  best_cs_pip <- temp_df$cs_pip[which.max(temp_df$overall_pip)]
  best_neglogp <- temp_df$neglogp_gwas[which.max(temp_df$overall_pip)]
  return(c(vetted_signal, signal_split[[1]][1], length(all_variants), max_overall_pip, best_variant, best_cs_pip, best_neglogp, paste0(all_variants, collapse = ",")))
}
#Make a function that gets maximum values and unions of variant sets from a list of info on all cs merged into a master cs
get_output_info <- function(vetted_signals_info, ancestry){
  if (length(vetted_signals_info) > 1) {
    #Get maximum PIP values for the merged cred set
    max_merged_overall_pip_index <- which.max(lapply(vetted_signals_info, FUN = extract_pos, pos=4))
    best_merged_overall_variant <- lapply(vetted_signals_info, FUN = extract_pos, pos=5)[[max_merged_overall_pip_index]]
    max_merged_overall_pip <- lapply(vetted_signals_info, FUN = extract_pos, pos=4)[[max_merged_overall_pip_index]]
    max_merged_overall_cs_pip <- lapply(vetted_signals_info, FUN = extract_pos, pos=6)[[max_merged_overall_pip_index]]
    best_merged_overall_neglogp <-lapply(vetted_signals_info, FUN = extract_pos, pos=7)[[max_merged_overall_pip_index]]
    #Get variants in credible set
    merged_variants_union <- unique(unlist(lapply(lapply(vetted_signals_info, FUN = extract_pos, pos=8),strsplit, split=",")))
  } else {
    best_merged_overall_variant <- vetted_signals_info[[1]][5]
    max_merged_overall_pip <- vetted_signals_info[[1]][4]
    max_merged_overall_cs_pip <- vetted_signals_info[[1]][6]
    best_merged_overall_neglogp <-vetted_signals_info[[1]][7]
    #Get variants in credible set
    merged_variants_union <- unique(unlist(strsplit(vetted_signals_info[[1]][8], split = ",")))
  }
  #Return outputs
  temp_list <- c(length(merged_variants_union), max_merged_overall_pip, best_merged_overall_variant, max_merged_overall_cs_pip, best_merged_overall_neglogp, paste0(merged_variants_union, collapse = ","))
  names(temp_list) <- paste0(ancestry, c(".num_variants", ".max_overall_pip", ".best_variant", ".best_variant_cs_pip", ".best_variant_ancestry-specific_neglogp", ".variant_ids"))
  return(temp_list)
}

### Extract info from raw_rds for unique signals and output to tables ###
#Make lists to hold summary-level and variant-level results
out_list <- list()
variant_level_list <- list()
#Loop over unique signals and extract info for each merged signal
for (i in unique_signals){
  component_signals <- names(cut_signals[which(cut_signals == unique_signals[i])])
  if (length(component_signals) == 1) { #Add a two part check for no signals found
    if (is.na(component_signals)) {
      next
    }
  }
  #Vet component_signals for purity/significance filters
  component_signals_vetted <- vector()
  for (signal in component_signals){
    signal_split <- strsplit(signal, "[.]")
    signal_purity <- raw_rds[[signal_split[[1]][1]]][["sets"]][["purity"]][paste0("L",signal_split[[1]][2]),"min.abs.corr"]
    if (significance_type == "gwas") {
      signal_variants <- raw_rds[[signal_split[[1]][1]]][["sets"]][["cs"]][[paste0("L",signal_split[[1]][2])]]
      max_signal_sig <- max(raw_rds[[signal_split[[1]][1]]][["neglogp_gwas"]][signal_variants])
    } else if(significance_type == "absolute_residual"){
      signal_variants <- names(raw_rds[[signal_split[[1]][1]]][["sets"]][["cs"]][[paste0("L",signal_split[[1]][2])]])
      max_signal_sig <- max(raw_rds[[signal_split[[1]][1]]][["neglogp_absolute"]][signal_variants, signal_split[[1]][2]])
    } else {
      signal_variants <- names(raw_rds[[signal_split[[1]][1]]][["sets"]][["cs"]][[paste0("L",signal_split[[1]][2])]])
      max_signal_sig <- max(raw_rds[[signal_split[[1]][1]]][["neglogp_precede"]][signal_variants, signal_split[[1]][2]])
    }
    if (max_signal_sig >= -log10(significance_thresh) && signal_purity >= purity_thresh) {
      component_signals_vetted <- append(component_signals_vetted, signal)
    }
  }
  #Make sure at least one signal passed the filters
  if (length(component_signals_vetted) == 0){
    next #Skip to next signal
  }
  #Loop over the vetted signals and extract relevant variant-level info into a dataframe using an applied function
  complete_variant_info <- rbind.fill.matrix(lapply(component_signals_vetted, extract_variant_info, raw_rds=raw_rds, merged_signal=i))
  complete_variant_info <- cbind.data.frame(phenotype=rep(trait, nrow(complete_variant_info)), 
                                            locus=rep(paste(chrom, start_pos, end_pos, sep = "."), nrow(complete_variant_info)), complete_variant_info)
  #Summarize the variant level info now into signal level info
  vetted_signals_info <- lapply(component_signals_vetted, extract_signal_info, variant_level_df=complete_variant_info)
  names(vetted_signals_info) <- component_signals_vetted
  #Get output values for the merged cred set
  merged_values <- get_output_info(vetted_signals_info, "merged")
  #Count signals per ancestry
  component_signals_by_ancestry = table(factor(unlist(lapply(strsplit(component_signals_vetted, "[.]"), FUN = extract_pos, pos=1)), levels = ancestries))
  #Collect CS-level info
  temp_list <- c(trait, chrom, start_pos, end_pos, paste(chrom, start_pos, end_pos, sep = "."), length(out_list) + 1, max(as.numeric(component_signals_by_ancestry)))
  names(temp_list) <- c("trait", "chr", "start", "end", "locus", "signal", "component_signals")
  temp_list <- c(temp_list, merged_values)
  for (ancestry in inp_ancestries) {
    if (!(ancestry %in% names(component_signals_by_ancestry))) {
      #Append empty data to temp_list
      empty_list <- c(NA, NA, NA, NA, NA, NA)
      names(empty_list) <- paste0(ancestry, c(".num_variants", ".max_overall_pip", ".best_variant", ".best_variant_cs_pip", ".best_variant_ancestry-specific_neglogp", ".variant_ids"))
      temp_list <- c(temp_list, empty_list)
    } else if (component_signals_by_ancestry[ancestry] != 0) {
      ancestry_vetted_signals_info <- vetted_signals_info[lapply(vetted_signals_info, FUN = extract_pos, pos=2) == ancestry]
      temp_list <- c(temp_list, get_output_info(ancestry_vetted_signals_info, ancestry))
    } else {
      #Append empty data to temp_list
      empty_list <- c(NA, NA, NA, NA, NA, NA)
      names(empty_list) <- paste0(ancestry, c(".num_variants", ".max_overall_pip", ".best_variant", ".best_variant_cs_pip", ".best_variant_ancestry-specific_neglogp", ".variant_ids"))
      temp_list <- c(temp_list, empty_list)
    }
  }
  #Append outputs onto lists
  out_list[[length(out_list) + 1]] <- as.data.frame(t(temp_list))
  variant_level_list[[length(variant_level_list) + 1]] <- complete_variant_info
}

#Merge output lists into dataframes
out_summary_df <- rbind.fill(out_list)
out_variant_df <- rbind.fill(variant_level_list)
#Write out the files
write.table(out_summary_df, file=df_out_loc, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(out_variant_df, file=str_replace(df_out_loc, "MERGE", "MERGE_VARIANT"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)


