################################################################################

#nominate_genes.R

#The purpose of this script is to nominate genes at signals and provide some 
#ancillary analyses related to the nominations.

#This code runs with the R/4.2 module

################################################################################

# 0) Call libraries and set variables and directories ====

#Call libraries
library(org.Hs.eg.db)
library(plyr)
library(ggpubr)
library(ggplot2)
library(stringr)
library(pbapply)
library(data.table)
library(tidyverse)
library(PairedData)
library(ggrepel)
library(scales)
library(ggpmisc)
library(gridExtra)
library(pwr)
library(bedr)
library(parallel)
library(clusterProfiler)
library(biomartr)
library(org.Hs.eg.db)

#Set directory
work_dir <- "/project/voight_GWAS/mconery/downstream_analyses/"
gene_nom_dir <- "/project/voight_GWAS/mconery/gwphewas_gene_nominations/"
#Set other file locations
gene_ann_loc <- "/project/voight_datasets/annot/gencode.v19.annotation.gtf.gz"

#Set global variables
pip_thresh=0.01

# 1) Check Distribution of max-PIPs for each variant/trait combo ====

#Read in variant-level results file
variant_raw <- read.csv(paste0(work_dir, "master.variants.supplement.txt"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
#Append on variant-trait combo id to each row
variant_raw <- cbind.data.frame(variant_raw, trait_variant_combo=paste0(variant_raw$phenotype, '.', variant_raw$mvp_id))

#Make a dataframe where each row is a specific trait/variant combo
get_max_pip <- function(trait_variant_combo, variant_raw){
  return(max(variant_raw[variant_raw$trait_variant_combo == trait_variant_combo, "overall_pip"]))
}
cl <- makeCluster(detectCores())
trait_variant_combo_max_pips <- unlist(parLapply(cl = cl, X = unique(variant_raw$trait_variant_combo), 
                                                 fun = get_max_pip, variant_raw=variant_raw))
stopCluster(cl)

jpeg(paste0(gene_nom_dir,"trait_variant_max_pips.hist.jpeg"))
print(hist(trait_variant_combo_max_pips, xlab="Max PIP", ylab="Freq of Trait/Variant Combos", main=NULL, breaks=seq(0,1,0.01)))
dev.off()

# 2) Use Exons to Nominate Genes ====

#Count before filtering
length(unique(paste0(variant_raw$phenotype, ".", variant_raw$mvp_id))) #894,652 unique variant trait pairs
length(unique(variant_raw$phenotype)) #936 traits represented

#Filter for variant/trait combos with a minimum PIP of 0.01 
variant_filt <- variant_raw[variant_raw$overall_pip >= pip_thresh,]
length(unique(paste0(variant_filt$phenotype, ".", variant_filt$mvp_id))) #508,808 unique variant trait pairs
length(unique(variant_raw$phenotype)) #936 traits represented

#Convert table to bed file format to use bedtools format
variant_bed <- cbind.data.frame(chr=paste0("chr",variant_filt$chr),
                                start=variant_filt$bp,
                                end=variant_filt$bp+1,
                                variant_filt[,colnames(variant_filt)[-match("chr", colnames(variant_filt))]])

#Read in v19 gencode gene annotations
gene_ann_raw <- read.csv(gene_ann_loc, header = FALSE, sep = "\t")
colnames(gene_ann_raw) <- c("chr", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")
#Extract the exons
gene_ann_exon <- gene_ann_raw[which(gene_ann_raw$feature == "exon"),]
#Break up attribute column
max_zero <- function(value){return(max(value, 0))} #Function to prevent neagative start positions
temp <- str_locate(string = gene_ann_exon$attribute, "gene_name")[,"start"]
temp2 <- str_locate(string = gene_ann_exon$attribute, "gene_type")[,"start"]
gene_ann_exon_org <- cbind.data.frame(chr=gene_ann_exon[,c("chr")], 
                                      start = vapply(gene_ann_exon$start,
                                                     FUN = max_zero, FUN.VALUE = numeric(1))-1, 
                                      end = gene_ann_exon$end,
                                      score=gene_ann_exon[,c("score")], 
                                      name=substr(gene_ann_exon$attribute, temp + 10, str_locate(substr(gene_ann_exon$attribute, temp, nchar(gene_ann_exon$attribute)), ";") + temp - 2),
                                      type=substr(gene_ann_exon$attribute, temp2 + 10, str_locate(substr(gene_ann_exon$attribute, temp2, nchar(gene_ann_exon$attribute)), ";") + temp2 - 2),
                                      gene_ann_exon[,c("strand", "feature", "frame")])
remove(temp, temp2, gene_ann_exon)
#Remove exons where the start position equals the end position (These don't make any sense)
gene_ann_exon_org <- gene_ann_exon_org[gene_ann_exon_org$start != gene_ann_exon_org$end,]

#Sort bed files
variant_bed_sort <- bedr.sort.region(variant_bed, method = "lexicographical")
gene_ann_exon_sort <- bedr.sort.region(gene_ann_exon_org, method = "lexicographical")

#Intersect the variants and exons to identify coding variants
gencode_genes_intersect <- bedr(
  input = list(a = variant_bed_sort, b = gene_ann_exon_sort), 
  method = "intersect", 
  params = "-wa -wb -sorted"
)

#Count number of trait-variant pairs and traits after intersecting
length(unique(gencode_genes_intersect$phenotype)) #684 traits after gencode intersect
length(unique(paste0(gencode_genes_intersect$phenotype, ".", gencode_genes_intersect$mvp_id))) #40,783 variant trait pairs after intersect
length(unique(paste0(gencode_genes_intersect$phenotype, ".", gencode_genes_intersect$mvp_id, ".", gencode_genes_intersect$name))) #42,660 trait-variant-gene combos after intersect

#Read in the vep categories file that groups together variant annotations
vep_categories <- read.table(paste0(work_dir, "vep_categories.prioritized.txt"), header = TRUE, row.names = 1, sep = "\t")
#Append grouped annotations onto gencode_genes_intersect
gencode_genes_intersect <- cbind.data.frame(gencode_genes_intersect, 
                                            grouped_vep=vep_categories[gencode_genes_intersect$vep_annotation,"Grouped.Annotation"])
#Filter for genes that are implicated via missense variants, splice/start/stop Gain/Loss variants
coding_nominations <- gencode_genes_intersect[gencode_genes_intersect$grouped_vep %in% c("Splice/Start/Stop Gain/Loss", "Missense"),]
length(unique(coding_nominations$phenotype)) # 582 traits after removing synonymous/extraneous vep annotations
length(unique(paste0(coding_nominations$phenotype, ".", coding_nominations$mvp_id))) #9,915 variant trait pairs
length(unique(paste0(coding_nominations$phenotype, ".", coding_nominations$mvp_id, ".", coding_nominations$name))) #10,465 trait-variant-gene combos

#Filter for protein coding genes
coding_nominations <- coding_nominations[coding_nominations$type == "protein_coding",]
length(unique(coding_nominations$phenotype)) # 579 traits after removing non-protein coding transcripts
length(unique(paste0(coding_nominations$phenotype, ".", coding_nominations$mvp_id))) #9,616 variant trait pairs
length(unique(paste0(coding_nominations$phenotype, ".", coding_nominations$mvp_id, ".", coding_nominations$name))) #9,915 trait-variant-gene combos

#Get unique ids of all merged signals that should be excluded from the ABC nomination process
coding_nominations <- cbind.data.frame(coding_nominations, 
                                       trait_variant_gene=paste(coding_nominations$phenotype, coding_nominations$mvp_id, coding_nominations$name, sep="."),
                                       merged_signal_id=paste(coding_nominations$phenotype, coding_nominations$locus, coding_nominations$merged_signal, sep="."))

#Convert to a single row for each trait/variant/gene combo and make into a dataframe
get_trait_variant_genes <- function(trait_variant_gene, coding_nominations, ancestries=c("AFR","AMR","EAS","EUR")) {
  temp <- coding_nominations[coding_nominations$trait_variant_gene == trait_variant_gene,]
  trait <- temp$phenotype[1]
  coding_signal_gene <- temp$name[1]
  variant <- temp$mvp_id[1]
  rsid <- temp$rsid[1]
  vep_annot <- temp$vep_annotation[1]
  grouped_vep <- temp$grouped_vep[1]
  mapped_ancestries <- unique(temp$ancestry)
  max_overall_pip <- max(temp$overall_pip)
  get_anc_info <- function(ancestry, temp_df){
    temp2 <- temp_df[temp_df$ancestry == ancestry,][1,]
    return(c(temp2$overall_pip, temp2$eaf.ANC, temp2["beta.ANC"], temp2$mu))
  }
  anc_info <- unlist(lapply(ancestries, get_anc_info, temp_df=temp))
  names(anc_info) <- NULL
  merged_signal_ids <- unique(temp$merged_signal_id)
  return(c(trait, coding_signal_gene, variant, rsid, grouped_vep, vep_annot, 
           paste0(mapped_ancestries, collapse = ","), paste0(merged_signal_ids, collapse = ","),  
           max_overall_pip, anc_info))
}
cl <- makeCluster(detectCores())
  trait_variant_genes_list <- parLapply(cl = cl, X = unique(coding_nominations$trait_variant_gene), fun = get_trait_variant_genes, coding_nominations=coding_nominations)
stopCluster(cl)
coding_trait_variant_gene_df <- cbind.data.frame(unique(coding_nominations$trait_variant_gene), as.data.frame(matrix(unlist(trait_variant_genes_list), byrow = TRUE, ncol = 25)))
colnames(coding_trait_variant_gene_df) <- c("trait_variant_gene", "trait", "gene", "mvp_id", "rsid", "grouped_vep", "vep_annot", 
                                          "mapped_ancestries", "merged_signals", "max_overall_pip", 
                                          "AFR.overall_pip", "AFR.eaf", "AFR.beta", "AFR.susie_mu",
                                          "AMR.overall_pip", "AMR.eaf", "AMR.beta", "AMR.susie_mu",
                                          "EAS.overall_pip", "EAS.eaf", "EAS.beta", "EAS.susie_mu",
                                          "EUR.overall_pip", "EUR.eaf", "EUR.beta", "EUR.susie_mu")
#Save down dataframe
write.table(coding_trait_variant_gene_df, file=paste0(gene_nom_dir, "coding_nominations.trait_variant_gene_table.tsv"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
coding_trait_variant_gene_df <- read.table(file=paste0(gene_nom_dir, "coding_nominations.trait_variant_gene_table.tsv"), header=TRUE, sep="\t")

#Count statistics again as a check only since no additional filtering
length(unique(coding_trait_variant_gene_df$trait)) #579 traits again
length(unique(paste0(coding_trait_variant_gene_df$trait, ".", coding_trait_variant_gene_df$mvp_id))) #9,616 trait-variant pairs
nrow(coding_trait_variant_gene_df) #9,915 trait-variant-gene combos again

# 3) Use ABC Connections to Nominate Genes ====

#Read in the ABC Connections
abc_raw <- read.table(paste0(gene_nom_dir, "AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.txt.gz"), sep = "\t", header = TRUE)
#Bed tools sort the dataframe
abc_sort <- bedr.sort.region(abc_raw, method = "lexicographical")

#Make a bed file where each fine-mapped variant is represented once
unique_variant_bed_sort <- variant_bed_sort[duplicated(variant_bed_sort$mvp_id) == FALSE,c("chr", "start", "end", "mvp_id")]

#Intersect the variants and exons to identify coding variants
abc_intersect <- bedr(
  input = list(a = unique_variant_bed_sort, b = abc_sort), 
  method = "intersect", 
  params = "-wa -wb -sorted"
)
#Write to file
write.table(abc_intersect, file=paste0(gene_nom_dir, "abc_intersect.bed"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
abc_intersect <- read.table(paste0(gene_nom_dir, "abc_intersect.bed"), header=TRUE)

#Combine rows of abc_intersect together so that each variant-gene combination represents a single row
abc_intersect[,"variant_gene"] <- paste(abc_intersect$mvp_id, abc_intersect$TargetGene, sep = ".")
abc_variant_genes <- unique(abc_intersect$variant_gene)
get_abc_info <- function(abc_variant_gene, abc_intersect){
  temp <- abc_intersect[abc_intersect$variant_gene == abc_variant_gene,]
  cell_types <- temp$CellType
  num_cell_types <- length(cell_types)
  variant <- temp$mvp_id[1]
  gene <- temp$TargetGene[1]
  enhancer_type <- temp$class[1]
  implication_type <- ifelse(enhancer_type == "promoter" && temp$isSelfPromoter == "True", "abc_promoter_overlap", "abc_interaction")
  max_abc_score <- max(temp$ABC.Score)
  return(c(variant, gene, implication_type, enhancer_type, max_abc_score, num_cell_types, paste(cell_types, collapse = ",")))
}
cl <- makeCluster(detectCores())
abc_nominations_list <- parLapply(cl = cl, X = abc_variant_genes, fun = get_abc_info, abc_intersect=abc_intersect)
stopCluster(cl)
abc_nominations <- as.data.frame(matrix(unlist(abc_nominations_list), ncol = 7, byrow = TRUE))
colnames(abc_nominations) <- c("mvp_id", "implicated_gene", "implication_type", "enhancer_type", "max_abc_score", "num_cell_types", "cell_type")

#Write dataframe to file
write.table(abc_nominations, file=paste0(gene_nom_dir, "abc_nominations.tsv"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
abc_nominations <- read.table(file=paste0(gene_nom_dir, "abc_nominations.tsv"), sep = "\t", header = TRUE)
#Set rownames of the abc_nominations df that reflects unique identifiers of variants and genes
rownames(abc_nominations) <- paste0(abc_nominations$mvp_id,".",abc_nominations$implicated_gene)
#Filter down the list of variants for just those with nominations
variant_abc <- variant_filt[variant_filt$mvp_id %in% unique(abc_nominations$mvp_id),]
#Combined rows of variant_abc together into a single row for each variant/trait like in the coding data frame
combine_by_trait_variants <- function(trait_variant, variant_abc, vep_categories, ancestries=c("AFR","AMR","EAS","EUR")) {
  temp <- variant_abc[variant_abc$trait_variant_combo == trait_variant,]
  trait <- temp$phenotype[1]
  variant <- temp$mvp_id[1]
  rsid <- temp$rsid[1]
  vep_annot <- temp$vep_annotation[1]
  grouped_vep <- vep_categories[vep_annot,"Grouped.Annotation"]
  mapped_ancestries <- unique(temp$ancestry)
  max_overall_pip <- max(temp$overall_pip)
  get_anc_info <- function(ancestry, temp_df){
    temp2 <- temp_df[temp_df$ancestry == ancestry,][1,]
    return(c(temp2$overall_pip, temp2$eaf.ANC, temp2["beta.ANC"], temp2$mu))
  }
  anc_info <- unlist(lapply(ancestries, get_anc_info, temp_df=temp))
  names(anc_info) <- NULL
  merged_signal_ids <- unique(temp$merged_signal_id)
  return(c(trait, variant, rsid, grouped_vep, vep_annot, 
           paste0(mapped_ancestries, collapse = ","), paste0(merged_signal_ids, collapse = ","),  
           max_overall_pip, anc_info))
}
cl <- makeCluster(detectCores())
combined_variant_abc_list <- parLapply(cl = cl, X = unique(variant_abc$trait_variant_combo), fun = combine_by_trait_variants, variant_abc=variant_abc, vep_categories = vep_categories)
stopCluster(cl)
combined_variant_abc <- cbind.data.frame(unique(variant_abc$trait_variant_combo), 
                                         as.data.frame(matrix(unlist(combined_variant_abc_list), byrow = TRUE, ncol = 24)))
colnames(combined_variant_abc) <- c("trait_variant", "trait", "mvp_id", "rsid", "grouped_vep", "vep_annot", 
                                          "mapped_ancestries", "merged_signals", "max_overall_pip", 
                                          "AFR.overall_pip", "AFR.eaf", "AFR.beta", "AFR.susie_mu",
                                          "AMR.overall_pip", "AMR.eaf", "AMR.beta", "AMR.susie_mu",
                                          "EAS.overall_pip", "EAS.eaf", "EAS.beta", "EAS.susie_mu",
                                          "EUR.overall_pip", "EUR.eaf", "EUR.beta", "EUR.susie_mu")

#Tie the abc nominations to each variant 
tie_abc_noms_to_variants <- combined_variant_abc %>% inner_join(abc_nominations, by = "mvp_id", relationship = "many-to-many")
#Count traits, trait-variants, and trait-variant genes again
length(unique(tie_abc_noms_to_variants$trait)) #727 traits
length((unique(paste0(tie_abc_noms_to_variants$trait_variant)))) #82,353 trait-variant combos
length((unique(paste0(tie_abc_noms_to_variants$trait_variant, ".", tie_abc_noms_to_variants$implicated_gene)))) #486,334 trait-variant combos

#Homogenize the coding and abc dataframes and then merge them
tie_abc_noms_to_variants <- 
  tie_abc_noms_to_variants[,c("trait_variant", "trait", "implicated_gene", 
                              "mvp_id", "rsid", "grouped_vep", "vep_annot", 
                            "mapped_ancestries", "merged_signals", "max_overall_pip", 
                            "AFR.overall_pip", "AFR.eaf", "AFR.beta", "AFR.susie_mu",
                            "AMR.overall_pip", "AMR.eaf", "AMR.beta", "AMR.susie_mu",
                            "EAS.overall_pip", "EAS.eaf", "EAS.beta", "EAS.susie_mu",
                            "EUR.overall_pip", "EUR.eaf", "EUR.beta", "EUR.susie_mu",
                            "implication_type", "max_abc_score", "num_cell_types", "cell_type")]
tie_abc_noms_to_variants$trait_variant <- paste0(tie_abc_noms_to_variants$trait_variant, ".", tie_abc_noms_to_variants$implicated_gene)
colnames(tie_abc_noms_to_variants) <- c("trait_variant_gene", "trait", "gene", 
                                        "mvp_id", "rsid", "grouped_vep", "vep_annot", 
                                        "mapped_ancestries", "merged_signals", "max_overall_pip", 
                                        "AFR.overall_pip", "AFR.eaf", "AFR.beta", "AFR.susie_mu",
                                        "AMR.overall_pip", "AMR.eaf", "AMR.beta", "AMR.susie_mu",
                                        "EAS.overall_pip", "EAS.eaf", "EAS.beta", "EAS.susie_mu",
                                        "EUR.overall_pip", "EUR.eaf", "EUR.beta", "EUR.susie_mu",
                                        "implication_type", "max_abc_score", "num_cell_types", "cell_type")
coding_trait_variant_gene_df <- cbind.data.frame(coding_trait_variant_gene_df,
                                                 implication_type = rep("nonsynonymous_coding_variant", nrow(coding_trait_variant_gene_df)),
                                                 max_abc_score = rep(NA, nrow(coding_trait_variant_gene_df)),
                                                 num_cell_types = rep(NA, nrow(coding_trait_variant_gene_df)),
                                                 cell_type = rep(NA, nrow(coding_trait_variant_gene_df)))
merged_nominations_df <- rbind.data.frame(coding_trait_variant_gene_df, tie_abc_noms_to_variants)
#Write file to table
write.table(merged_nominations_df, file=paste0(gene_nom_dir, "merged_nominations.tsv"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# 4) Analyze the Gene Nominations ====

#Read in the merged nominations
merged_nominations_df <- read.table(file=paste0(gene_nom_dir, "merged_nominations.tsv"), sep = "\t", header = TRUE)

#Count the number of traits, trait-variants, and trait-variant-genes again
length(unique(merged_nominations_df$trait)) #786 traits
length((unique(paste0(merged_nominations_df$trait, ".", merged_nominations_df$mvp_id)))) #89,652 trait-variant combos
length((unique(paste0(merged_nominations_df$trait_variant_gene)))) #494,181 trait-variant combos

#Remove ABC nominations of nonsynonymous variants for the genes in which the variants are located
#Get list of trait_variant_genes that have both abc and coding variant implications
double_counted_combos <- merged_nominations_df %>% group_by(trait_variant_gene) %>% summarize(count_implications = n()) %>%
  filter(count_implications > 1) %>% select(trait_variant_gene) %>% as.data.frame
double_counted_combos <- double_counted_combos[,1]
merged_nominations_df <- merged_nominations_df %>% filter(!(trait_variant_gene %in% double_counted_combos) | 
                                                            implication_type == "nonsynonymous_coding_variant")

#Count the number of traits, trait-variants, and trait-variant-genes again
length(unique(merged_nominations_df$trait)) #786 traits
length((unique(paste0(merged_nominations_df$trait, ".", merged_nominations_df$mvp_id)))) #89,652 trait-variant combos
length((unique(paste0(merged_nominations_df$trait_variant_gene)))) #494,181 trait-variant combos
temp <- merged_nominations_df %>% filter(implication_type != "nonsynonymous_coding_variant")
length(unique(temp$trait)) #727 traits
length((unique(paste0(temp$trait, ".", temp$mvp_id)))) #81,608 trait-variant combos
length((unique(paste0(temp$trait_variant_gene)))) #484,266 trait-variant combos

#Extract the Gencode genes
gene_ann_gene <- gene_ann_raw[which(gene_ann_raw$feature == "gene"),]
#Split up the attribute column in the gencode genes list
temp <- str_locate(string = gene_ann_gene$attribute, "gene_name")[,"start"]
temp2 <- str_locate(string = gene_ann_gene$attribute, "gene_type")[,"start"]
gene_ann_gene_org <- paste0(substr(gene_ann_gene$attribute, temp + 10, str_locate(substr(gene_ann_gene$attribute, temp, nchar(gene_ann_gene$attribute)), ";") + temp - 2),"/",
                                      substr(gene_ann_gene$attribute, temp2 + 10, str_locate(substr(gene_ann_gene$attribute, temp2, nchar(gene_ann_gene$attribute)), ";") + temp2 - 2))
remove(temp, temp2, gene_ann_gene)
gene_ann_gene_org <- as.data.frame(matrix(unlist(str_split(unique(gene_ann_gene_org), "/")), byrow = TRUE, ncol = 2))
colnames(gene_ann_gene_org) <- c("gene", "gene_type")
#Remove duplicated gencode gene names with distinct gene_types (sort first so that transcript types with known functions will be at top)
gene_ann_gene_org <- rbind.data.frame(gene_ann_gene_org[gene_ann_gene_org$gene_type == "protein_coding",],
                                      gene_ann_gene_org[gene_ann_gene_org$gene_type == "miRNA",],
                                      gene_ann_gene_org[gene_ann_gene_org$gene_type == "snRNA",],
                                      gene_ann_gene_org[gene_ann_gene_org$gene_type == "snoRNA",],
                                      gene_ann_gene_org[gene_ann_gene_org$gene_type == "rRNA",],
                                      gene_ann_gene_org[gene_ann_gene_org$gene_type == "Mt_tRNA",],
                                      gene_ann_gene_org[gene_ann_gene_org$gene_type == "Mt_rRNA",],
                                      gene_ann_gene_org[!(gene_ann_gene_org$gene_type %in% c("protein_coding", "miRNA", "snRNA", "snoRNA", "rRNA", "Mt_tRNA", "Mt_rRNA")),])
gene_ann_gene_org <- gene_ann_gene_org[!(duplicated(gene_ann_gene_org$gene)),]
#Get gene types from gencode list
merged_nominations_df <- merged_nominations_df %>% inner_join(gene_ann_gene_org, by = "gene")

#Count the number of traits, trait-variants, and trait-variant-genes again
length(unique(merged_nominations_df$trait)) #785 traits
length((unique(paste0(merged_nominations_df$trait, ".", merged_nominations_df$mvp_id)))) #86,819 trait-variant combos
length((unique(paste0(merged_nominations_df$trait_variant_gene)))) #430,586 trait-variant combos
temp <- merged_nominations_df %>% filter(implication_type != "nonsynonymous_coding_variant")
length(unique(temp$trait)) #725 traits
length((unique(paste0(temp$trait, ".", temp$mvp_id)))) #78,657 trait-variant combos
length((unique(paste0(temp$trait_variant_gene)))) #420,671 trait-variant combos

#Make histogram of abc scores for abc nominations (Color by implication_type)
jpeg(paste0(gene_nom_dir,"abc_scores.hist.jpeg"))
merged_nominations_df %>% filter(implication_type != "nonsynonymous_coding_variant") %>%
ggplot(aes(max_abc_score, fill=implication_type)) + 
  geom_histogram(binwidth=0.01)
dev.off()
#Make same histogram but color by implication type
jpeg(paste0(gene_nom_dir,"abc_scores.color_by_gene_type.hist.jpeg"), width = 720, height=720)
merged_nominations_df %>% filter(implication_type != "nonsynonymous_coding_variant") %>%
  ggplot(aes(max_abc_score, fill=gene_type)) + 
  geom_histogram(binwidth=0.01)
dev.off()

#Make a bar plot again colored by protein coding proportion
bar_plot_df <- merged_nominations_df %>% filter(implication_type != "nonsynonymous_coding_variant") %>% 
  mutate(bin_name = floor(max_abc_score *100)/100 ) %>% 
  mutate(protein_coding = ifelse(gene_type == "protein_coding", "Protein Coding", "ncRNA")) %>% 
  group_by(protein_coding, bin_name) %>% summarize(bin_count = n())
bar_plot_df <- merged_nominations_df %>% filter(implication_type != "nonsynonymous_coding_variant") %>% 
  mutate(bin_name = floor(max_abc_score *100)/100 ) %>% 
  mutate(protein_coding = ifelse(gene_type == "protein_coding", "Protein Coding", "ncRNA")) %>%
  group_by(bin_name) %>% summarize(total_bin_count = n()) %>% inner_join(bar_plot_df, by="bin_name") %>%
  mutate(bin_pct = bin_count/total_bin_count) 
bar_plot_df$bin_name <- factor(as.character(bar_plot_df$bin_name), levels=as.character(seq(0,max(bar_plot_df$bin_name),0.01)))
jpeg(paste0(gene_nom_dir,"abc_scores.color_by_gene_type.barplot.jpeg"), width = 720, height=720)
  ggplot(bar_plot_df, aes(y=bin_pct, x=bin_name)) + geom_bar(aes(fill=protein_coding), stat="identity") + 
    xlab(label = "Maximum ABC Score") + ylab(label = "Proportion of Protein Coding Genes") + 
  theme(axis.title = element_text(size = 16),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), legend.title = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "bottom", legend.text = element_text(size = 14), axis.text.x = element_blank()) 
dev.off()

#Filter the merged_nominations for the protein-coding genes and only things with an abc score > 0.1
filtered_merged_nominations_df <- merged_nominations_df %>% filter(max_abc_score > 0.1 | is.na(max_abc_score))
temp <- filtered_merged_nominations_df %>% filter(implication_type != "nonsynonymous_coding_variant")
length(unique(temp$trait)) #529 traits
length((unique(paste0(temp$trait, ".", temp$mvp_id)))) #22,538 trait-variant combos
length((unique(paste0(temp$trait_variant_gene)))) #40,936 trait-variant combos
filtered_merged_nominations_df <- filtered_merged_nominations_df %>% filter(gene_type == "protein_coding") 
temp <- filtered_merged_nominations_df %>% filter(implication_type != "nonsynonymous_coding_variant")
length(unique(temp$trait)) #517 traits
length((unique(paste0(temp$trait, ".", temp$mvp_id)))) #21,798 trait-variant combos
length((unique(paste0(temp$trait_variant_gene)))) #38,141 trait-variant 

#Write table to file
write.table(filtered_merged_nominations_df, file=paste0(gene_nom_dir, "filtered_merged_nominations.tsv"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
filtered_merged_nominations_df <- read.table(paste0(gene_nom_dir, "filtered_merged_nominations.tsv"), header = TRUE, sep = "\t")

#Count the number of traits, trait-variants, and trait-variant-genes again
length(unique(filtered_merged_nominations_df$trait)) #668 traits
length((unique(paste0(filtered_merged_nominations_df$trait, ".", filtered_merged_nominations_df$mvp_id)))) #31,020 trait-variant combos
length((unique(paste0(filtered_merged_nominations_df$trait_variant_gene)))) #48,056 trait-variant combos
temp <- filtered_merged_nominations_df %>% mutate(trait_variant = paste0(trait, ".", mvp_id)) %>% 
  filter(implication_type == "nonsynonymous_coding_variant") %>% select(trait_variant) 
filtered_merged_nominations_df %>% mutate(trait_variant = paste0(trait, ".", mvp_id)) %>% 
  filter(implication_type != "nonsynonymous_coding_variant") %>% select(trait_variant) %>%
  intersect(temp) %>% nrow #394 trait-variant pairs that implicate different genes by coding and ABC 
filtered_merged_nominations_df %>% mutate(trait_variant = paste0(trait, ".", mvp_id)) %>%
  group_by(trait_variant) %>% summarize(unique_genes = length(unique(gene))) %>% 
  filter(unique_genes > 1) %>% nrow #8,039 variant-traits implicate multiple gene
filtered_merged_nominations_df %>% mutate(trait_gene = paste0(trait, ".", gene)) %>% 
  group_by(trait_gene) %>% summarize(unique_variants = length(unique(trait))) %>% 
  filter(unique_variants == 1) %>% nrow #0 trait-genes were implicated by multiple variants

#Make a histogram of the number of genes nominated per trait
jpeg(paste0(gene_nom_dir,"genes_per_trait.hist.jpeg"), width = 720, height=720)
filtered_merged_nominations_df %>% group_by(trait) %>% summarize(num_genes = n()) %>%
  ggplot(aes(num_genes)) + geom_histogram(binwidth = 5) + 
  xlab(label = "Number of Genes") + ylab(label = "Number of Traits") + 
  theme(axis.title = element_text(size = 16),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), legend.title = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "bottom", legend.text = element_text(size = 14), 
        axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14)) 
dev.off()

#Read in the file of traits
trait_raw <- read.csv(paste0(work_dir, "MVP_R4.1000G_AGR.GIA.DataDictionary.revised.022023.txt"), header = TRUE, sep = "\t")
rownames(trait_raw) <- trait_raw$Trait
trait_raw <- cbind.data.frame(trait_raw, num_samples.META=(ifelse(is.na(trait_raw$num_samples.AFR), 0, trait_raw$num_samples.AFR) + 
                                                             ifelse(is.na(trait_raw$num_samples.AMR), 0, trait_raw$num_samples.AMR) + 
                                                             ifelse(is.na(trait_raw$num_samples.EAS), 0, trait_raw$num_samples.EAS) + 
                                                             ifelse(is.na(trait_raw$num_samples.EUR), 0, trait_raw$num_samples.EUR)))
#Calculate effective sample sizes
trait_raw <- cbind.data.frame(trait_raw, 
                               eff_samples.META=ifelse(trait_raw$Trait_Type == "binary", 
                                                       2/((1/trait_raw$num_cases.META)+(1/trait_raw$num_controls.META)), 
                                                       trait_raw$num_samples.META))
colnames(trait_raw) <- str_replace(colnames(trait_raw),"Trait", "trait")
#Make a scatter plot of number of gene nominations vs sample size
cbPalette <- c("#009E73", "#0072B2", "#56B4E9", "#F0E442", "#E69F00", "#D55E00", "#CC79A7")
trait_append <- filtered_merged_nominations_df %>% group_by(trait) %>% summarize(num_genes = n()) %>%
  inner_join(trait_raw, by="trait")
trait_append$Category <- factor(trait_append$Category, levels = c("PheCodes", "Labs", "Vitals", "Baseline_Survey", "Lifestyle_Survey"))
jpeg(paste0(gene_nom_dir,"genes_per_trait_sample_size.scatter.jpeg"), width = 720, height=720)
  ggplot(trait_append, aes(x=eff_samples.META, y=num_genes)) + geom_point(aes(color = Category)) + 
  scale_y_continuous(breaks = c(1,10,100,1000,10000), minor_breaks = rep(1:9, 4)*(10^rep(0:3, each = 9)), limits = c(0, max(trait_append$num_genes)), labels = c(1,10,100,"1,000","10,000"), trans = "pseudo_log") + 
  scale_x_continuous(breaks = c(1000,10000,100000,1000000), limits = c(1000, max(trait_append$eff_samples.META)), labels = c("   1,000","   10,000","   100,000","   1,000,000"), trans = "pseudo_log") + 
  xlab("(Effective) Meta-Analysis Sample Size") + ylab("Number of Implicated Genes") + 
  guides(colour = guide_legend(override.aes = list(size=6))) +
  scale_color_manual(labels=c('Phecodes', 'Labs', 'Vitals', 'Baseline Survey', 'Lifestyle Survey'), values=cbPalette[1:5]) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), legend.title = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(color="black", size=16, angle=90, hjust = 1), axis.title=element_text(size=24), 
        legend.position = "bottom", axis.text.y = element_text(color="black", size=16), 
        legend.text = element_text(size=13), legend.key=element_blank()) 
dev.off()

#Read in signal-level results file
merged_append <- read.table(file = paste0(work_dir, "master.signals.appended.txt"), sep = "\t", header=TRUE)
#Append on number of signals info to trait dataframe
temp <- as.data.frame(table(merged_append$trait))
colnames(temp) <- c("trait", "num_signals")
trait_append <- trait_append %>% inner_join(temp, by = "trait")
plot_regression <- lm(log10(num_genes) ~ log10(num_signals), data = trait_append)
jpeg(paste0(gene_nom_dir,"genes_signals_per_trait.scatter.jpeg"), width = 720, height=720)
  ggplot(trait_append, aes(x=num_signals, y=num_genes)) + geom_point(aes(color = Category)) + 
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "darkgray", linewidth = 2) +  # Add 45-degree line
    geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "solid") +  # Add regression line
  scale_y_continuous(breaks = c(1,10,100,1000,10000), minor_breaks = rep(1:9, 4)*(10^rep(0:3, each = 9)), limits = c(0, max(trait_append$num_genes)), labels = c(1,10,100,"1,000","10,000"), trans = "pseudo_log") +
  scale_x_continuous(breaks = c(1,10,100,1000,10000), limits = c(0, max(trait_append$num_signals)), labels = c("     1","    10","   100"," 1,000","10,000"), trans = "pseudo_log") + 
  xlab("Number of Fine-Mapped Signals") + ylab("Number of Implicated Genes") + 
  guides(colour = guide_legend(override.aes = list(size=6))) +
  scale_color_manual(labels=c('Phecodes', 'Labs', 'Vitals', 'Baseline Survey', 'Lifestyle Survey'), values=cbPalette[1:5]) + 
  annotate("text", x = 0, y = max(trait_append$num_genes), 
           label = paste("y =", round(coef(plot_regression)[2], 2), "x +", round(coef(plot_regression)[1], 2)),
           hjust = 0, vjust = 1, size = 6, color = "black") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), legend.title = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(color="black", size=16, angle=90, hjust = 1), axis.title=element_text(size=24), 
        legend.position = "bottom", axis.text.y = element_text(color="black", size=16), 
        legend.text = element_text(size=13), legend.key=element_blank()) 
dev.off()

# 5) Remake the plots including the traits without signals/gene nominations ====

#Append on the number of genes column again but this  time retain traits without gene nominations
trait_append <- filtered_merged_nominations_df %>% group_by(trait) %>% summarize(num_genes = n()) %>%
  right_join(trait_raw, by="trait")
trait_append$num_genes <- ifelse(is.na(trait_append$num_genes), 0, trait_append$num_genes)
#Append on the number of signals again but this  time retain traits without signals
temp <- as.data.frame(table(merged_append$trait))
colnames(temp) <- c("trait", "num_signals")
trait_append <- trait_append %>% left_join(temp, by = "trait")
trait_append$num_signals <- ifelse(is.na(trait_append$num_signals), 0, trait_append$num_signals)
#Recode the Category column
trait_append$Category <- factor(trait_append$Category, levels = c("PheCodes", "Labs", "Vitals", "Baseline_Survey", "Lifestyle_Survey"))

#Remake the num genes vs. sample size
jpeg(paste0(gene_nom_dir,"genes_per_trait_sample_size.all_traits.scatter.jpeg"), width = 720, height=720)
ggplot(trait_append, aes(x=eff_samples.META, y=num_genes)) + geom_point(aes(color = Category)) + 
  scale_y_continuous(breaks = c(1,10,100,1000,10000), minor_breaks = rep(1:9, 4)*(10^rep(0:3, each = 9)), limits = c(0, max(trait_append$num_genes)), labels = c(1,10,100,"1,000","10,000"), trans = "pseudo_log") + 
  scale_x_continuous(breaks = c(1000,10000,100000,1000000), limits = c(1000, max(trait_append$eff_samples.META)), labels = c("   1,000","   10,000","   100,000","   1,000,000"), trans = "pseudo_log") + 
  xlab("(Effective) Meta-Analysis Sample Size") + ylab("Number of Implicated Genes") + 
  guides(colour = guide_legend(override.aes = list(size=6))) +
  scale_color_manual(labels=c('Phecodes', 'Labs', 'Vitals', 'Baseline Survey', 'Lifestyle Survey'), values=cbPalette[1:5]) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), legend.title = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(color="black", size=16, angle=90, hjust = 1), axis.title=element_text(size=24), 
        legend.position = "bottom", axis.text.y = element_text(color="black", size=16), 
        legend.text = element_text(size=13), legend.key=element_blank()) 
dev.off()

#Remake the num signals vs. sample size
jpeg(paste0(gene_nom_dir,"trait_category.all_traits.scatter.jpeg"), width = 720, height=720)
ggplot(trait_append, aes(x=eff_samples.META, y=num_signals)) + geom_point(aes(color = Category)) + 
  scale_y_continuous(breaks = c(1,10,100,1000,10000), minor_breaks = rep(1:9, 4)*(10^rep(0:3, each = 9)), limits = c(0, max(trait_append$num_genes)), labels = c(1,10,100,"1,000","10,000"), trans = "pseudo_log") + 
  scale_x_continuous(breaks = c(1000,10000,100000,1000000), limits = c(1000, max(trait_append$eff_samples.META)), labels = c("   1,000","   10,000","   100,000","   1,000,000"), trans = "pseudo_log") + 
  xlab("(Effective) Meta-Analysis Sample Size") + ylab("Number of Signals") + 
  guides(colour = guide_legend(override.aes = list(size=6))) +
  scale_color_manual(labels=c('Phecodes', 'Labs', 'Vitals', 'Baseline Survey', 'Lifestyle Survey'), values=cbPalette[1:5]) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), legend.title = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(color="black", size=16, angle=90, hjust = 1), axis.title=element_text(size=24), 
        legend.position = "bottom", axis.text.y = element_text(color="black", size=16), 
        legend.text = element_text(size=13), legend.key=element_blank()) 
dev.off()

#Remake the num genes vs. num signals plot
jpeg(paste0(gene_nom_dir,"genes_signals_per_trait.all_traits.scatter.jpeg"), width = 720, height=720)
ggplot(trait_append, aes(x=num_signals, y=num_genes)) + geom_point(aes(color = Category)) + 
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "darkgray", linewidth = 2) +  # Add 45-degree line
  scale_y_continuous(breaks = c(1,10,100,1000,10000), minor_breaks = rep(1:9, 4)*(10^rep(0:3, each = 9)), limits = c(0, max(trait_append$num_genes)), labels = c(1,10,100,"1,000","10,000"), trans = "pseudo_log") +
  scale_x_continuous(breaks = c(1,10,100,1000,10000), limits = c(0, max(trait_append$num_signals)), labels = c("     1","    10","   100"," 1,000","10,000"), trans = "pseudo_log") + 
  xlab("Number of Fine-Mapped Signals") + ylab("Number of Implicated Genes") + 
  guides(colour = guide_legend(override.aes = list(size=6))) +
  scale_color_manual(labels=c('Phecodes', 'Labs', 'Vitals', 'Baseline Survey', 'Lifestyle Survey'), values=cbPalette[1:5]) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), legend.title = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(color="black", size=16, angle=90, hjust = 1), axis.title=element_text(size=24), 
        legend.position = "bottom", axis.text.y = element_text(color="black", size=16), 
        legend.text = element_text(size=13), legend.key=element_blank()) 
dev.off()

# 6) Execute KEGG Enrichment ====

#Read in gene nominations
filtered_merged_nominations_df <- read.table(paste0(gene_nom_dir, "filtered_merged_nominations.tsv"), header = TRUE, sep = "\t")

#Collapse quantiatitve traits
collapsed_nominations_df <- filtered_merged_nominations_df %>% 
  mutate(trait_cut=str_replace_all(trait, "_M[[:alpha:]]*_INT", "")) %>% 
  mutate(trait_cut_variant_gene=paste(trait_cut, mvp_id, gene, sep=".")) %>% group_by(trait_cut_variant_gene) %>%
  summarize(trait=first(trait_cut), gene=first(gene), variant=first(mvp_id), implication_type=first(implication_type), 
            max_overall_pip=max(max_overall_pip)) %>% as.data.frame()
#Write to file
write.table(collapsed_nominations_df, file = paste0(gene_nom_dir, "collapsed_nominations.tsv"),row.names = FALSE, col.names = TRUE, sep = "\t", quote=FALSE)

#Get trait list
traits_with_noms <- unique(collapsed_nominations_df$trait) #668 traits with at least one gene nomination

#Make a function 
test_kegg_enrichments <- function(test_trait, collapsed_nominations_df){
  #Load needed libraries
  library(org.Hs.eg.db)
  library(clusterProfiler)
  #Extract genes for trait as hgnc codes
  hgnc_symbols <- unique(collapsed_nominations_df[collapsed_nominations_df$trait == test_trait,"gene"])
  #Map to uniprot ids
  uniprot_ids <- mapIds(org.Hs.eg.db, keys = hgnc_symbols, column = "UNIPROT", keytype = "SYMBOL")
  uniprot_ids <- sort(uniprot_ids, decreasing = TRUE)
  names(uniprot_ids) <- NULL
  uniprot_ids <- uniprot_ids[!(is.na(uniprot_ids))]
  # Run KEGG pathway enrichment analysis
  kk <- enrichKEGG(gene = uniprot_ids, organism = 'hsa', keyType = "uniprot", 
                   pvalueCutoff = 1, pAdjustMethod = "BH")
  #Make return table
  if (is.null(kk)) {
    return(NULL)
  } else {
    temp <- kk@result
    temp <- cbind.data.frame(temp, trait = rep(test_trait, nrow(temp))) 
    return(temp)
  }
}
#Call the function
cl <- makeCluster(detectCores())
traits_nominations_list <- parLapply(cl = cl, X = traits_with_noms, 
                                  fun = test_kegg_enrichments, collapsed_nominations_df=collapsed_nominations_df)
stopCluster(cl)
#Remove null entries
is.not.null <- function(item){return(!(is.null(item)))}
traits_nominations_list <- traits_nominations_list[unlist(lapply(traits_nominations_list, FUN = is.not.null))]
#Merge the list together
full_kegg_df <- as.data.frame(rbind.fill.matrix(traits_nominations_list))
#Readjust the p-values
full_kegg_df[,"p.adjust"] <- p.adjust(full_kegg_df[,"pvalue"], method = "BH")

#Filter for significant enrichments
sig_kegg_df <- full_kegg_df %>% filter(p.adjust < 0.05)

#Reorganize file columns for outputting
colnames(sig_kegg_df)[1] <- "Pathway_ID"
colnames(sig_kegg_df) <- str_replace(colnames(sig_kegg_df), "trait", "Trait")
sig_kegg_df <- sig_kegg_df[,c("Trait", "Pathway_ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "geneID")]
#Reconvert uniprot gene ids back to gene symbols
convert_uniprot_symbol <- function(geneids){
  uniprot_ids <- unlist(str_split(geneids, "/"))
  gene_symbols <- mapIds(org.Hs.eg.db, keys = uniprot_ids, column = "SYMBOL", keytype = "UNIPROT")
  return(paste0(gene_symbols, collapse = "/"))
}
sig_kegg_df[,"Gene_Symbols"] <- vapply(sig_kegg_df$geneID, FUN = convert_uniprot_symbol, FUN.VALUE = character(1))
#Sort table
sig_kegg_df <- sig_kegg_df[order(sig_kegg_df$Trait),]

#Write to file
write.table(sig_kegg_df, file = paste0(gene_nom_dir, "significant_kegg_enrichments.tsv"), sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

#Get list of traits with 2+ KEGG enrichments
traits_w_kegg <- sig_kegg_df %>% group_by(Trait) %>% summarize(count=n()) %>% filter(count > 1) %>% select(Trait) %>% as.data.frame()
traits_w_kegg <- traits_w_kegg[,"Trait"]
#Create output directory for cluego datasets if needed
cluego_dir = paste0(gene_nom_dir, "cluego_kegg_pathways/")
if (!file.exists(cluego_dir)) {
  dir.create(cluego_dir, recursive = TRUE)
}
#Write files into a format useable for ClueGo for each trait
write_cluego_kegg <- function(trait, sig_kegg_df, cluego_dir){
  cluego_whole_df <- sig_kegg_df %>% filter(Trait == trait)
  cluego_whole_df <- cluego_whole_df[,c("Pathway_ID", "p.adjust")]
  cluego_whole_df$Pathway_ID <- str_replace(cluego_whole_df$Pathway_ID, "hsa", "KEGG:")
  colnames(cluego_whole_df) <- c("GOID.PathwayID", "P value")
  write.table(cluego_whole_df, file = paste0(cluego_dir, trait, ".kegg_pathways.cluego.tsv"), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
}
lapply(traits_w_kegg, FUN = write_cluego_kegg, sig_kegg_df=sig_kegg_df, cluego_dir=cluego_dir)

# 7) Add Variant Novelty to Gene Nominations ====
  
#Read in gene nominations
filtered_merged_nominations_df <- read.table(paste0(gene_nom_dir, "filtered_merged_nominations.tsv"), header = TRUE, sep = "\t")

#Read in list of known/novel variant/trait combos
known_novel_raw <- read.table(paste0(work_dir, "all.pheno.finemap.known-novel.05152023.txt"), header = TRUE, sep = "\t")
known_novel_raw <- cbind.data.frame(known_novel_raw, unique_id=paste0(known_novel_raw$SNP, ".", known_novel_raw$phenotype))
known_novel_raw <- known_novel_raw[which(!(duplicated(known_novel_raw$unique_id))),]
rownames(known_novel_raw) <- known_novel_raw$unique_id
#Recode known_novel_raw column
known_novel_raw$Novel_label <- ifelse(known_novel_raw$Known_Association=="True", "Known Association",
                                      ifelse(known_novel_raw$Novel_Association_Known_Signal=="True", "Novel Association Known Signal", 
                                             ifelse(known_novel_raw$Novel_Signal=="True", "Novel Signal", NA)))

#Append on novelty information to gene nominations
filtered_merged_nominations_df[,"High_PIP_Novelty"] <- known_novel_raw[paste(filtered_merged_nominations_df$mvp_id, filtered_merged_nominations_df$trait, sep="."), "Novel_label"]
filtered_merged_nominations_df$High_PIP_Novelty <- ifelse(is.na(filtered_merged_nominations_df$High_PIP_Novelty), "Not High PIP", filtered_merged_nominations_df$High_PIP_Novelty)

#Write to file
write.table(filtered_merged_nominations_df, file=paste0(gene_nom_dir, "gene_nominations_supplement.tsv"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

#Read in list of all signals
trait_raw <- read.csv(paste0(work_dir, "MVP_R4.1000G_AGR.GIA.DataDictionary.revised.022023.txt"), header = TRUE, sep = "\t")
#Select trait and description columns
trait_desc <- trait_raw %>% select(Trait, Description)
colnames(trait_desc) <- c("trait", "description")
trait_desc <- trait_desc %>% filter(!duplicated(trait))

#Make an abridged table for reviewing variants
compressed_nominations <- filtered_merged_nominations_df %>% select(trait, mvp_id, gene, implication_type, High_PIP_Novelty) %>%
  inner_join(trait_desc, by = "trait") %>% group_by(mvp_id, gene) %>% 
  summarize(traits = paste(trait, collapse=","), descriptions = paste(description, collapse=","), 
            implication_type = paste0(unique(implication_type), collapse=","), High_PIP_Novelty = paste0(High_PIP_Novelty, collapse = ","))
split_comma_check_high_pip <- function(inp_string){
  temp <- (unlist(str_split(inp_string, pattern = ",")) != "Not High PIP")
  return(c(length(temp), sum(temp)))
}
temp <- matrix(unlist(vapply(compressed_nominations$High_PIP_Novelty, split_comma_check_high_pip, FUN.VALUE = integer(2))), ncol = 2, byrow = TRUE)
colnames(temp) <- c("num_traits", "num_High_PIP_traits")
compressed_nominations <- cbind.data.frame(compressed_nominations, temp)
#Write to file
write.table(compressed_nominations, file=paste0(gene_nom_dir, "compressed_nominations.tsv"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# 8) Analyze Results with Quantitative Traits collapsed ====

#Read in gene nominations
filtered_merged_nominations_df <- read.table(paste0(gene_nom_dir, "filtered_merged_nominations.tsv"), header = TRUE, sep = "\t")
#Collapse quantiatitve traits
collapsed_nominations_df <- filtered_merged_nominations_df %>% 
  mutate(trait_cut=str_replace_all(trait, "_M[[:alpha:]]*_INT", "")) %>% 
  mutate(trait_cut_variant_gene=paste(trait_cut, mvp_id, gene, sep=".")) %>% group_by(trait_cut_variant_gene) %>%
  dplyr::summarize(trait=first(trait_cut), gene=first(gene), variant=first(mvp_id), implication_type=first(implication_type))
#Calculate statistics
nrow(collapsed_nominations_df)
length(unique(collapsed_nominations_df$trait))
length(unique(collapsed_nominations_df$variant))
length(unique(collapsed_nominations_df$gene))

#Make donut plot of implication types
temp <- collapsed_nominations_df %>% group_by(implication_type) %>% dplyr::summarize(pct=scales::percent(n()/nrow(collapsed_nominations_df))) %>% as.data.frame()
rownames(temp) <- temp$implication_type
collapsed_nominations_df <- collapsed_nominations_df %>% 
  mutate(category_pct=temp[implication_type,"pct"]) %>%
  mutate(implication_plot=ifelse(implication_type=="nonsynonymous_coding_variant", "Non-Synonymous\nCoding", 
                                 ifelse(implication_type=="abc_promoter_overlap", "ABC\nPromoter", "ABC\nInteraction"))) %>%
  mutate(implication_plot=paste0(implication_plot, "\n", category_pct))
data <- as.data.frame(table(collapsed_nominations_df$implication_plot))
colnames(data) <- c("class", "count")
#Set plot label
plot_label <- "Implication\nTypes"
# Compute percentages
data <- data %>% 
  arrange(desc(class)) %>%
  mutate(prop = count / sum(count)) %>%
  mutate(lab.ypos = cumsum(prop) - 0.5*prop) %>%
  mutate(lebal=percent(round(prop,2)))
#Make a donut plot of VEP annotations
temp_plot <- ggplot(data, aes(x = 2, y = prop, fill = class)) +
  geom_bar(stat = "identity", color = "black") +
  coord_polar(theta = "y", start = 0) + # Try to remove that to understand how the chart is built initially
  geom_text(aes(y = lab.ypos, x=2,  label=lebal), size = 10) + 
  annotate(geom = 'text', x = 0.5, y = 0, label = plot_label, size = 10) + 
  xlim(.5, 2.5) + # Try to remove that to see how to make a pie chart
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y = element_blank(),
        axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.x = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
        legend.position="none")
jpeg(paste0(gene_nom_dir, "implication_types.collapsed.donut.jpeg"), width = 8, height = 8, units = 'in', res = 500)
print(temp_plot)
dev.off()

#Make histogram plotting dataframe
trait_filt <- trait_raw %>% mutate(num_samples.META=(ifelse(is.na(num_samples.AFR), 0, num_samples.AFR) + 
                                         ifelse(is.na(num_samples.AMR), 0, num_samples.AMR) + 
                                         ifelse(is.na(num_samples.EAS), 0, num_samples.EAS) + 
                                         ifelse(is.na(num_samples.EUR), 0, num_samples.EUR))) %>% 
  mutate(trait_cut=str_replace_all(Trait, "_M[[:alpha:]]*_INT", ""),
         eff_samples.META=ifelse(Trait_Type == "binary", 
                                 2/((1/num_cases.META)+(1/num_controls.META)), num_samples.META)) %>% 
  group_by(trait_cut) %>% summarize(eff_samples.META=mean(eff_samples.META), Category=first(Category))
collapsed_plot_df <- collapsed_nominations_df %>% group_by(trait) %>% 
  summarize(genes=length(unique(gene))) %>% inner_join(trait_filt, join_by(trait == trait_cut)) %>%
  select(trait, genes, Category, eff_samples.META)
collapsed_plot_df$Category <- factor(collapsed_plot_df$Category, levels = c("PheCodes", "Labs", "Vitals", "Baseline_Survey", "Lifestyle_Survey"))
#Calc trait gene combinations
collapsed_plot_df %>% select(genes) %>% sum
#Make histogram plot
cbPalette <- c("#009E73", "#0072B2", "#56B4E9", "#F0E442", "#E69F00", "#D55E00", "#CC79A7")
jpeg(paste0(gene_nom_dir,"genes_per_trait.collapsed_trait.hist.jpeg"))
  ggplot(data=collapsed_plot_df, aes(x=genes, fill=Category)) + geom_histogram() + 
  scale_y_continuous(breaks = c(1,10,100,1000), labels = c(1,10,100,"1,000"), trans = "pseudo_log") + 
  xlab("Nominated Genes Per Trait") + ylab("Number of Traits") + 
  scale_fill_manual(labels=c('Phecodes', 'Labs', 'Vitals', 'Baseline Survey', 'Lifestyle Survey'), values=cbPalette[1:5]) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), legend.title = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(color="black", size=16, angle=90, hjust = 1), axis.title=element_text(size=24), 
        legend.position = "bottom", axis.text.y = element_text(color="black", size=16), 
        legend.text = element_text(size=13), legend.key=element_blank()) 
dev.off()
#Make scatter plot
jpeg(paste0(gene_nom_dir,"genes_per_trait_sample_size.collapsed.scatter.jpeg"), width = 720, height=720)
ggplot(collapsed_plot_df, aes(x=eff_samples.META, y=genes)) + geom_point(aes(color = Category)) + 
  scale_y_continuous(breaks = c(1,10,100,1000,10000), minor_breaks = rep(1:9, 4)*(10^rep(0:3, each = 9)), limits = c(0, max(collapsed_plot_df$genes)), labels = c(1,10,100,"1,000","10,000"), trans = "pseudo_log") + 
  scale_x_continuous(breaks = c(1000,10000,100000,1000000), limits = c(1000, max(collapsed_plot_df$eff_samples.META)), labels = c("   1,000","   10,000","   100,000","   1,000,000"), trans = "pseudo_log") + 
  xlab("(Effective) Meta-Analysis Sample Size") + ylab("Number of Implicated Genes") + 
  guides(colour = guide_legend(override.aes = list(size=6))) +
  scale_color_manual(labels=c('Phecodes', 'Labs', 'Vitals', 'Baseline Survey', 'Lifestyle Survey'), values=cbPalette[1:5]) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), legend.title = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(color="black", size=16, angle=90, hjust = 1), axis.title=element_text(size=24), 
        legend.position = "bottom", axis.text.y = element_text(color="black", size=16), 
        legend.text = element_text(size=13), legend.key=element_blank()) 
dev.off()


# 9) Check for Allelic Series ====

#Read in gene nominations
filtered_merged_nominations_df <- read.table(paste0(gene_nom_dir, "filtered_merged_nominations.tsv"), header = TRUE, sep = "\t")

#Find any instances of the same gene being implicated by both coding and non-coding variation for a given trait
appended_nominations_df <- filtered_merged_nominations_df %>% mutate(trait_gene = paste0(trait, "/", gene)) 
trait_genes <- unique(appended_nominations_df$trait_gene)
check_coding_non_coding <- function(trait_gene_inp, appended_nominations_df){
  temp <- appended_nominations_df %>% filter(trait_gene == trait_gene_inp)
  check_cod <- ifelse(nrow(temp[temp$implication_type == "nonsynonymous_coding_variant",]) > 0, TRUE, FALSE)
  check_prom <- ifelse(nrow(temp[temp$implication_type == "abc_promoter_overlap",]) > 0, TRUE, FALSE)
  check_interact <- ifelse(nrow(temp[temp$implication_type == "abc_interaction",]) > 0, TRUE, FALSE)
  return(c(check_cod, check_prom, check_interact))
}
implication_type_check_list <- vapply(trait_genes, FUN = check_coding_non_coding, FUN.VALUE = logical(3), appended_nominations_df=appended_nominations_df)
implication_type_check_df <- cbind.data.frame(trait_gene=trait_genes,
                                              matrix(data = unlist(implication_type_check_list), ncol = 3, byrow = TRUE))
colnames(implication_type_check_df) <- c("trait_gene", "coding", "promoter", "interaction")

#Check how many rows have hits in coding and one of the other two columns
implication_type_check_df %>% filter(coding == TRUE & (promoter == TRUE | interaction == TRUE)) %>% nrow #2192
implication_type_check_df %>% filter(coding == TRUE & (interaction == TRUE)) %>% nrow #1339 
implication_type_check_df %>% filter(coding == TRUE & (promoter == TRUE)) %>% nrow #1321
implication_type_check_df %>% filter(coding == TRUE & (promoter == TRUE & interaction == TRUE)) %>% nrow #468

#Make an abridged dataframe for reviewing nominations
split_slash_take_one <- function(string_inp){return(str_split(string_inp, pattern = "/")[[1]][1])}
split_slash_take_two <- function(string_inp){return(str_split(string_inp, pattern = "/")[[1]][2])}
implication_type_check_df <- implication_type_check_df %>%  
  mutate(trait = vapply(trait_gene, FUN = split_slash_take_one, FUN.VALUE=character(1)), 
         gene = vapply(trait_gene, FUN = split_slash_take_two, FUN.VALUE=character(1))) %>% 
  inner_join(trait_desc, by = "trait") 
#Count number of unique genes
implication_type_check_df %>% filter(coding == TRUE & (promoter == TRUE & interaction == TRUE)) %>% 
  group_by(gene) %>% summarize(traits = paste(trait, collapse=","), descriptions = paste(description, collapse=",")) %>% 
  nrow() #209 unique genes are represented

#Append prevalency info to the implication data frame
trait_thin <- trait_raw %>% mutate(prevalency_AFR = ifelse(Trait_Type == "quantitative", "quantitative", num_cases.AFR/(num_cases.AFR + num_controls.AFR)),
                                   prevalency_AMR = ifelse(Trait_Type == "quantitative", "quantitative", num_cases.AMR/(num_cases.AMR + num_controls.AMR)),
                                   prevalency_EAS = ifelse(Trait_Type == "quantitative", "quantitative", num_cases.EAS/(num_cases.EAS + num_controls.EAS)),
                                   prevalency_EUR = ifelse(Trait_Type == "quantitative", "quantitative", num_cases.EUR/(num_cases.EUR + num_controls.EUR))) %>%
  select(Trait, prevalency_AFR, prevalency_AMR, prevalency_EAS, prevalency_EUR)
implication_type_check_df <- implication_type_check_df %>% inner_join(trait_thin, join_by(trait == Trait))
#Collapse the quantitative traits together
review_df <- implication_type_check_df %>%
  mutate(trait_cut=str_replace_all(trait, "_M[[:alpha:]]*_INT", ""), 
         description_cut=str_replace_all(description, " \\(.*transformed\\)", "")) %>% group_by(trait_cut, gene) %>%
  summarize(description=first(description_cut),
    coding_any=max(coding), promoter_any=max(promoter), interaction_any=max(interaction),
    prevalency_AFR=first(prevalency_AFR), prevalency_AMR=first(prevalency_AMR), prevalency_EAS=first(prevalency_EAS), prevalency_EUR=first(prevalency_EUR))
#Make a histogram of the trait-collapsed gene nomination dataframe
jpeg(paste0(gene_nom_dir,"traits_per_gene.collapsed_trait.hist.jpeg"))
  review_df %>% group_by(gene) %>% summarize(traits = n()) %>% 
    ggplot(aes(x=traits)) + geom_histogram(fill = "#E69F00") + 
      scale_y_continuous(breaks = c(1,10,100,1000), labels = c(1,10,100,"1,000"), trans = "pseudo_log") + 
            xlab("Traits per Nominated Gene") + ylab("Number of Genes") + 
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
            panel.background = element_blank(), legend.title = element_blank(), axis.line = element_line(colour = "black"),
            axis.text.x = element_text(color="black", size=16, angle=90, hjust = 1), axis.title=element_text(size=24), 
            legend.position = "bottom", axis.text.y = element_text(color="black", size=16), 
            legend.text = element_text(size=13), legend.key=element_blank()) 
dev.off()

#Write tables to file
write.table(implication_type_check_df, file = paste0(gene_nom_dir, "gene_trait_implications.tsv"), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(review_df, file = paste0(gene_nom_dir, "gene_trait_implications.collapsed.tsv"), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

#Count number of traits per gene
three_way_gene_pleiotropy_df <- review_df %>% filter(coding_any==1 & promoter_any==1 & interaction_any==1) %>%
  group_by(gene) %>% summarise(traits=n()) %>% arrange(desc(traits))
#Write to table
write.table(three_way_gene_pleiotropy_df, file = paste0(gene_nom_dir, "three_way_gene_pleiotropy.tsv"), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
