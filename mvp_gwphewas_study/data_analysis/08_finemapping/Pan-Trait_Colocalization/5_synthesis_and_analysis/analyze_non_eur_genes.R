################################################################################

#analyze_non_eur_genes.R

#The purpose of this script is to study the genes nominated by the ABC/coding
#analysis and identify any that are nominated by variants that are specific to 
#Non-EUR populations. 

#This code runs locally

################################################################################

# 0) Call libraries and set variables and directories ====

#Call libraries
library(tidyverse)
library(ggrepel)

#Set directory
gene_nom_dir <- "C:/Users/mitch/Documents/UPenn/Voight_Lab/MVP_Finemapping/gene_nominations/"

#Set trait file
trait_file <- "C:/Users/mitch/Documents/UPenn/Voight_Lab/MVP_Finemapping/Figure_Building_and_Synthesis/MVP_R4.1000G_AGR.GIA.DataDictionary.revised.022023.txt"

#Set a random seed
random_seed=5
set.seed(random_seed)

#Set a threshold on rarity
rarity_thresh=0.001

# 1) Read in gene nominations with Ancestry Info ====

#Read in file
nominations_w_anc <- read.csv(paste0(gene_nom_dir, "filtered_merged_nominations_w_anc_summary.tsv"), header = TRUE, sep = "\t")

#Add on a collapsed trait column
nominations_w_anc <- nominations_w_anc %>%  mutate(trait_cut=trait) %>%
  mutate(trait_cut_variant_gene=paste(trait_cut, mvp_id, gene, sep=".")) %>%
  mutate(trait_cut_gene=paste(trait_cut, gene, sep="."))

# 2) Calculate the number of genes and variant genes nominated via only non-EUR variation ====

#Filter for variants that are monomorphic in EUR study population or are ultra rare
non_eur_df <- nominations_w_anc %>% filter(is.na(maf_EUR) | maf_EUR < rarity_thresh)

#separate out monomorphic vs ultra rare and count for trait_cut version (collapsed traits)
non_eur_df %>% filter(is.na(maf_EUR)) %>% dplyr::select(trait_cut_variant_gene) %>% unique() %>% nrow #573 unique trait variant gene combos for monomorphic EUR variants
non_eur_df %>% filter(maf_EUR < rarity_thresh) %>% dplyr::select(trait_cut_variant_gene) %>% unique() %>% nrow #685 unique trait variant gene combos for ultra-rare EUR variants
non_eur_df %>% dplyr::select(trait_cut_variant_gene) %>% unique() %>% nrow # 1258 unique trait-variant-gene combos for both categories

#Check whether the trait-gene combo would have been found by any other variants
temp <- non_eur_df %>% dplyr::select(trait_cut_gene) %>% unique() %>% as.data.frame()
temp <- temp[,1]
eur_df <- nominations_w_anc %>% filter(is.na(maf_EUR)==FALSE & maf_EUR >= rarity_thresh)
temp2 <- ifelse(temp %in% eur_df$trait_cut_gene, FALSE, TRUE)
names(temp2) <- temp
non_eur_df <- non_eur_df %>% mutate(tg_non_eur_only=temp2[trait_cut_gene])
#Count the number of trait-gene combos that can only be found by non-EUR variants
non_eur_only_trait_genes <- non_eur_df %>% filter(tg_non_eur_only == TRUE) %>% dplyr::select(trait_cut_gene) %>% unique %>% nrow
total_non_eur_trait_genes <- non_eur_df %>% dplyr::select(trait_cut_gene) %>% unique %>% nrow
total_trait_genes <- nominations_w_anc %>% dplyr::select(trait_cut_gene) %>% unique %>% nrow
non_eur_only_trait_genes/total_non_eur_trait_genes #60.0% of trait-genes combos implicated via non-EUR variants were only implicated via non-EUR variants
non_eur_only_trait_genes/total_trait_genes #3.5% of all trait-gene combos were only implicated via non-EUR variants. 

#Make a bed file of non-EUR variants (I wanted to use this with the ABC viewer, but it appears their tool is broken.)
split_colon_take_x <- function(string_inp, x){return(str_split(string_inp, pattern = ":")[[1]][x])}
non_eur_bed <- non_eur_df %>% dplyr::select(mvp_id, rsid, max_overall_pip, pvalue_AFR, max_overall_pip) %>% 
  mutate(chromo=paste0("chr", split_colon_take_x(mvp_id, 1)), start=as.numeric(split_colon_take_x(mvp_id, 2))) %>%
  mutate(end=start+1) %>% dplyr::select(chromo, start, end, rsid, max_overall_pip, pvalue_AFR, max_overall_pip)
write.table(non_eur_bed, file = paste0(gene_nom_dir, "non_eur_variants.bed"), sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

# 3) Create a version of the output table to save as the supplemental table ====

supplement_df <- nominations_w_anc %>% 
  dplyr::select(trait, gene, mvp_id, rsid, trait_cut, vep_annot, grouped_vep, mapped_ancestries, max_overall_pip, 
                implication_type, max_abc_score, num_cell_types, cell_type,
                AFR.overall_pip, AFR.susie_mu, beta_AFR, pvalue_AFR, maf_AFR, 
                AMR.overall_pip, AMR.susie_mu, beta_AMR, pvalue_AMR, maf_AMR,
                EAS.overall_pip, EAS.susie_mu, beta_EAS, pvalue_EAS, maf_EAS, 
                EUR.overall_pip, EUR.susie_mu, beta_EUR, pvalue_EUR, maf_EUR)
#Write to file
write.table(supplement_df, file = paste0(gene_nom_dir, "gene_nominations_supplement.tsv"), col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

# 4) Create trait pruning functions ====

#Make functions to do the trait pruning by genetic correlations
#This is a function that identifies the best trait remaining based on max pip and then random selection in case of tie
select_best_trait <- function(trait_pip_df){ 
  temp <- trait_pip_df[trait_pip_df$max_pip == max(trait_pip_df$max_pip),]
  return(as.character(temp[sample(nrow(temp), size=1),"trait"]))
}
#This function checks correlations of current trait with remaining traits and tosses any that are correlated above 0.2
check_correlations <- function(iteration_trait, remaining_traits, coef_d, corr_thresh=0.2){ 
  #Check which remaining traits are in the colnames for coef_d
  checkable_traits <- remaining_traits[remaining_traits %in% colnames(coef_d)]
  #Check the correlations for those traits
  test_corr <- coef_d[iteration_trait, checkable_traits]
  return(checkable_traits[!(abs(test_corr) > corr_thresh)])
}
#Make a function to get number of independent traits for the gene
get_indep_traits <- function(gene, filtered_merged_nominations_df, coef_p, coef_m, coef_f){
  temp <- filtered_merged_nominations_df[filtered_merged_nominations_df$gene==gene,]
  traits <- unique(temp$trait)
  if (length(traits) == 1) { #Check if only 1 trait was identified for the gene
    trait_pip_df <- temp %>% group_by(trait) %>% summarize(mvp_id=mvp_id[which.max(max_overall_pip)], max_pip=max(max_overall_pip), max_anc_maf=max(max_anc_maf)) %>% as.data.frame()
    rownames(trait_pip_df) <- trait_pip_df$trait
    return(c(gene, trait_pip_df[traits,]))
  } else { #If multiple traits were identified, then we need to see how many were independent
    #Get maximum PIP for each trait
    trait_pip_df <- temp %>% group_by(trait) %>% summarize(mvp_id=mvp_id[which.max(max_overall_pip)], max_pip=max(max_overall_pip), max_anc_maf=max(max_anc_maf)) %>% as.data.frame()
    rownames(trait_pip_df) <- trait_pip_df$trait
    #Create a vector to hold all independent traits that have been selected for a given phenotype and those remaining 
    selected_indep_traits <- vector()
    remaining_traits <- unique(trait_pip_df$trait)
    #Loop until the remaining traits are empty
    while (length(remaining_traits) > 0) {
      iteration_trait = select_best_trait(trait_pip_df)
      selected_indep_traits = append(selected_indep_traits, c(gene, trait_pip_df[iteration_trait,]))
      remaining_traits = remaining_traits[remaining_traits != iteration_trait]
      #Check which of the matrices need to be used
      if (iteration_trait %in% rownames(coef_p)) {
        remaining_traits <- check_correlations(iteration_trait = iteration_trait, remaining_traits = remaining_traits, coef_d = coef_p)
      } else if(iteration_trait %in% rownames(coef_m)){
        remaining_traits <- check_correlations(iteration_trait = iteration_trait, remaining_traits = remaining_traits, coef_d = coef_m)
      } else {
        remaining_traits <- check_correlations(iteration_trait = iteration_trait, remaining_traits = remaining_traits, coef_d = coef_f)
      }
      #Reset the trait pip df
      trait_pip_df <- trait_pip_df[remaining_traits,]
    }
    #Return the number of independent traits
    return(selected_indep_traits)
  }
}

#Load file of genetic correlations
load(paste0(gene_nom_dir, "ALL.corr.RData"))

# 5) Analyze pleiotropy specifically in EUR ancestry ====

#Filter for nominations made in EUR
EUR_nominations_df <- nominations_w_anc %>% filter(str_detect(pattern = "EUR", mapped_ancestries)) %>%
  mutate(max_anc_maf = maf_EUR)

#Get unique genes
EUR_genes_with_noms <- unique(EUR_nominations_df$gene)

#Apply the trait-pruning functions
EUR_indep_trait_per_gene <- lapply(EUR_genes_with_noms, get_indep_traits, 
                                   filtered_merged_nominations_df=EUR_nominations_df, 
                                   coef_p=coef_p, coef_m=coef_m, coef_f=coef_f)
#Convert it to a data frame
EUR_indep_trait_per_gene_df <- as.data.frame(matrix(unlist(EUR_indep_trait_per_gene), ncol = 5, byrow = TRUE))
colnames(EUR_indep_trait_per_gene_df) <- c("gene", "trait", "mvp_id", "max_overall_pip", "max_anc_maf")
EUR_indep_trait_per_gene_df$max_overall_pip <- as.numeric(EUR_indep_trait_per_gene_df$max_overall_pip)
EUR_indep_trait_per_gene_df$max_anc_maf <- as.numeric(EUR_indep_trait_per_gene_df$max_anc_maf)
split_colon_take_x <- function(string_inp, x){return(str_split(string_inp, pattern = ":")[[1]][x])}
EUR_indep_trait_per_gene_df <- cbind.data.frame(EUR_indep_trait_per_gene_df, 
                                                chromo=as.numeric(vapply(EUR_indep_trait_per_gene_df$mvp_id, FUN = split_colon_take_x, FUN.VALUE = character(1), x=1)), 
                                                pos=as.numeric(vapply(EUR_indep_trait_per_gene_df$mvp_id, FUN = split_colon_take_x, FUN.VALUE = character(1), x=2)))

#convert the table to a gene-based table
EUR_trait_count_df <- EUR_indep_trait_per_gene_df %>% group_by(gene) %>% summarize(traits=n(), max_anc_maf=max(max_anc_maf), trait_names=paste(trait,collapse = ","))

# 6) Analyze pleiotropy specifically in AFR ancestry ====

#Filter for nominations made in AFR
AFR_nominations_df <- nominations_w_anc %>% filter(str_detect(pattern = "AFR", mapped_ancestries)) %>%
  mutate(max_anc_maf = maf_AFR)

#Calculate some statistics
length(unique(AFR_nominations_df$trait_cut_variant_gene)) #4,135 AFR trait-variant-gene combinations
length(unique(AFR_nominations_df$trait_cut_gene)) #2,452 AFR trait-gene nominations
AFR_temp <- AFR_nominations_df %>% dplyr::select(trait_cut_gene) %>% unique()
AFR_temp <- AFR_temp[,1]
AFR_temp2 <- ifelse(AFR_temp %in% eur_df$trait_cut_gene, FALSE, TRUE)
names(AFR_temp2) <- AFR_temp
sum(AFR_temp2) #392 (16.0%) AFR trait-gene nominations could be found only in non-EUR

#Get unique genes
AFR_genes_with_noms <- unique(AFR_nominations_df$gene)

#Apply the trait-pruning functions
AFR_indep_trait_per_gene <- lapply(AFR_genes_with_noms, get_indep_traits, 
                               filtered_merged_nominations_df=AFR_nominations_df, 
                               coef_p=coef_p, coef_m=coef_m, coef_f=coef_f)
#Convert it to a data frame
AFR_indep_trait_per_gene_df <- as.data.frame(matrix(unlist(AFR_indep_trait_per_gene), ncol = 5, byrow = TRUE))
colnames(AFR_indep_trait_per_gene_df) <- c("gene", "trait", "mvp_id", "max_overall_pip", "max_anc_maf")
AFR_indep_trait_per_gene_df$max_overall_pip <- as.numeric(AFR_indep_trait_per_gene_df$max_overall_pip)
AFR_indep_trait_per_gene_df$max_anc_maf <- as.numeric(AFR_indep_trait_per_gene_df$max_anc_maf)
split_colon_take_x <- function(string_inp, x){return(str_split(string_inp, pattern = ":")[[1]][x])}
AFR_indep_trait_per_gene_df <- cbind.data.frame(AFR_indep_trait_per_gene_df, 
                                            chromo=as.numeric(vapply(AFR_indep_trait_per_gene_df$mvp_id, FUN = split_colon_take_x, FUN.VALUE = character(1), x=1)), 
                                            pos=as.numeric(vapply(AFR_indep_trait_per_gene_df$mvp_id, FUN = split_colon_take_x, FUN.VALUE = character(1), x=2)))

#convert the table to a gene-based table
AFR_trait_count_df <- AFR_indep_trait_per_gene_df %>% group_by(gene) %>% summarize(traits=n(), max_anc_maf=max(max_anc_maf), trait_names=paste(trait,collapse = ","))

#Check number of genes with >=2 independent traits as a % of total
(AFR_trait_count_df %>% filter(traits > 1) %>% nrow) #240 genes nominated in AFR are pleiotropic
(AFR_trait_count_df %>% filter(traits > 1) %>% nrow)/nrow(AFR_trait_count_df)# It's 16.9% of the total

#Join data frames
AFR_EUR_trait_count_df <- merge(x = EUR_trait_count_df, y=AFR_trait_count_df, by="gene", all=TRUE)
AFR_EUR_trait_count_df$traits.x <- ifelse(is.na(AFR_EUR_trait_count_df$traits.x), 0, AFR_EUR_trait_count_df$traits.x)
AFR_EUR_trait_count_df$traits.y <- ifelse(is.na(AFR_EUR_trait_count_df$traits.y), 0, AFR_EUR_trait_count_df$traits.y)

#Run regression and calculate outliers
AFR_EUR_regress <- glm(traits.y ~ traits.x, data = AFR_EUR_trait_count_df, family = "poisson")
AFR_EUR_trait_count_df <- AFR_EUR_trait_count_df %>% mutate(residuals = AFR_EUR_regress$residuals)
AFR_EUR_trait_count_df <- AFR_EUR_trait_count_df %>% mutate(outliers=ifelse(n() - rank(residuals) < 3, gene, ""))

#Make line plotting dataframe
regression_line_df <- seq(0, max(AFR_EUR_trait_count_df$traits.x), 0.01)
regression_line_df <- cbind.data.frame(traits.x=regression_line_df, traits.y=predict(AFR_EUR_regress, newdata = data.frame(traits.x = regression_line_df), type = "response"))
#Get plot limits
axis_lim = c(0, max(AFR_EUR_trait_count_df$traits.x, AFR_EUR_trait_count_df$traits.y))
#Make plot size
plot_width = (max(AFR_EUR_trait_count_df$traits.x)/max(AFR_EUR_trait_count_df$traits.x, AFR_EUR_trait_count_df$traits.y)) * 720
plot_height = (max(AFR_EUR_trait_count_df$traits.y)/max(AFR_EUR_trait_count_df$traits.x, AFR_EUR_trait_count_df$traits.y)) * 720 * 1.1

#Create Scatterplot
jpeg(paste0(gene_nom_dir,"AFR_v_EUR.pleiotropy.scatter.jpeg"), width = plot_width*10, height = plot_height*10, res = 1000)
ggplot(data=AFR_EUR_trait_count_df, mapping = aes(x=traits.x, y=traits.y)) + 
  geom_text_repel(size = 6, mapping = aes(label=outliers, color=residuals), nudge_y = 2.5, nudge_x = 0.4) + 
  scale_color_gradient2(low="blue", mid="white", high="red", oob = scales::squish, guide = "none") +
  geom_point(shape=21, size= 4, mapping = aes(fill=residuals)) + scale_fill_gradient2(low="blue", mid="white", high="red", oob = scales::squish) +
  geom_line(data = regression_line_df, aes(x=traits.x, y=traits.y)) + 
  xlab("# Traits EUR") + ylab("# Traits AFR") + 
  xlim(axis_lim) + ylim(axis_lim) + 
  guides(fill=guide_legend(title="Residuals")) + 
  theme(axis.title.x = element_text(size = 24), legend.position = "bottom", 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA), legend.title = element_text(size = 24), 
        axis.line = element_line(colour = "black"), legend.key = element_rect(fill = "white"), legend.text = element_text(size=20),
        plot.background = element_rect(fill = "transparent",colour = NA),
        axis.title.y = element_text(size = 24), axis.text.y = element_text(size = 14), axis.text.x = element_text(size = 14))
dev.off()

#Write table to file
AFR_EUR_trait_count_df <- AFR_EUR_trait_count_df[order(AFR_EUR_trait_count_df$residuals, decreasing = TRUE),]
write.table(AFR_EUR_trait_count_df, file=paste0(gene_nom_dir, "AFR_EUR_traits_per_gene.tsv"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# 7) Analyze pleiotropy specifically in AMR ancestry ====

#Filter for nominations made in AMR
AMR_nominations_df <- nominations_w_anc %>% filter(str_detect(pattern = "AMR", mapped_ancestries)) %>%
  mutate(max_anc_maf = maf_AMR)

#Calculate some statistics
length(unique(AMR_nominations_df$trait_cut_variant_gene)) #4,135 AMR trait-variant-gene combinations
length(unique(AMR_nominations_df$trait_cut_gene)) #2,452 AMR trait-gene nominations
AMR_temp <- AMR_nominations_df %>% dplyr::select(trait_cut_gene) %>% unique()
AMR_temp <- AMR_temp[,1]
AMR_temp2 <- ifelse(AMR_temp %in% eur_df$trait_cut_gene, FALSE, TRUE)
names(AMR_temp2) <- AMR_temp
sum(AMR_temp2) #392 (16.0%) AMR trait-gene nominations could be found only in non-EUR

#Get unique genes
AMR_genes_with_noms <- unique(AMR_nominations_df$gene)

#Apply the trait-pruning functions
AMR_indep_trait_per_gene <- lapply(AMR_genes_with_noms, get_indep_traits, 
                                   filtered_merged_nominations_df=AMR_nominations_df, 
                                   coef_p=coef_p, coef_m=coef_m, coef_f=coef_f)
#Convert it to a data frame
AMR_indep_trait_per_gene_df <- as.data.frame(matrix(unlist(AMR_indep_trait_per_gene), ncol = 5, byrow = TRUE))
colnames(AMR_indep_trait_per_gene_df) <- c("gene", "trait", "mvp_id", "max_overall_pip", "max_anc_maf")
AMR_indep_trait_per_gene_df$max_overall_pip <- as.numeric(AMR_indep_trait_per_gene_df$max_overall_pip)
AMR_indep_trait_per_gene_df$max_anc_maf <- as.numeric(AMR_indep_trait_per_gene_df$max_anc_maf)
split_colon_take_x <- function(string_inp, x){return(str_split(string_inp, pattern = ":")[[1]][x])}
AMR_indep_trait_per_gene_df <- cbind.data.frame(AMR_indep_trait_per_gene_df, 
                                                chromo=as.numeric(vapply(AMR_indep_trait_per_gene_df$mvp_id, FUN = split_colon_take_x, FUN.VALUE = character(1), x=1)), 
                                                pos=as.numeric(vapply(AMR_indep_trait_per_gene_df$mvp_id, FUN = split_colon_take_x, FUN.VALUE = character(1), x=2)))

#convert the table to a gene-based table
AMR_trait_count_df <- AMR_indep_trait_per_gene_df %>% group_by(gene) %>% summarize(traits=n(), max_anc_maf=max(max_anc_maf), trait_names=paste(trait,collapse = ","))

#Check number of genes with >=2 independent traits as a % of total
(AMR_trait_count_df %>% filter(traits > 1) %>% nrow) #240 genes nominated in AMR are pleiotropic
(AMR_trait_count_df %>% filter(traits > 1) %>% nrow)/nrow(AMR_trait_count_df)# It's 16.9% of the total

#Join data frames
AMR_EUR_trait_count_df <- merge(x = EUR_trait_count_df, y=AMR_trait_count_df, by="gene", all=TRUE)
AMR_EUR_trait_count_df$traits.x <- ifelse(is.na(AMR_EUR_trait_count_df$traits.x), 0, AMR_EUR_trait_count_df$traits.x)
AMR_EUR_trait_count_df$traits.y <- ifelse(is.na(AMR_EUR_trait_count_df$traits.y), 0, AMR_EUR_trait_count_df$traits.y)

#Run regression and calculate outliers
AMR_EUR_regress <- glm(traits.y ~ traits.x, data = AMR_EUR_trait_count_df, family = "poisson")
AMR_EUR_trait_count_df <- AMR_EUR_trait_count_df %>% mutate(residuals = AMR_EUR_regress$residuals)
AMR_EUR_trait_count_df <- AMR_EUR_trait_count_df %>% mutate(outliers=ifelse(n() - rank(residuals) < 1, gene, ""))

#Make line plotting datAMRame
regression_line_df <- seq(0, max(AMR_EUR_trait_count_df$traits.x), 0.01)
regression_line_df <- cbind.data.frame(traits.x=regression_line_df, traits.y=predict(AMR_EUR_regress, newdata = data.frame(traits.x = regression_line_df), type = "response"))
#Get plot limits
axis_lim = c(0, max(AMR_EUR_trait_count_df$traits.x, AMR_EUR_trait_count_df$traits.y))
#Make plot size
#plot_width = (max(AMR_EUR_trait_count_df$traits.x)/max(AMR_EUR_trait_count_df$traits.x, AMR_EUR_trait_count_df$traits.y)) * 720
#plot_height = (max(AMR_EUR_trait_count_df$traits.y)/max(AMR_EUR_trait_count_df$traits.x, AMR_EUR_trait_count_df$traits.y)) * 720 * 1.1

#Create Scatterplot
jpeg(paste0(gene_nom_dir,"AMR_v_EUR.pleiotropy.scatter.jpeg"), width = plot_width*10, height = plot_height*10, res = 1000)
ggplot(data=AMR_EUR_trait_count_df, mapping = aes(x=traits.x, y=traits.y)) + 
  geom_text_repel(size = 6, mapping = aes(label=outliers, color=residuals), nudge_y = 3) + 
  scale_color_gradient2(low="blue", mid="white", high="red", oob = scales::squish, guide = "none") +
  geom_point(shape=21, size= 4, mapping = aes(fill=residuals)) + scale_fill_gradient2(low="blue", mid="white", high="red", oob = scales::squish) +
  geom_line(data = regression_line_df, aes(x=traits.x, y=traits.y)) + 
  xlab("# Traits EUR") + ylab("# Traits AMR") + 
  xlim(axis_lim) + ylim(axis_lim) + 
  guides(fill=guide_legend(title="Residuals")) + 
  theme(axis.title.x = element_text(size = 24), legend.position = "bottom", 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA), legend.title = element_text(size = 24), 
        axis.line = element_line(colour = "black"), legend.key = element_rect(fill = "white"), legend.text = element_text(size=20),
        plot.background = element_rect(fill = "transparent",colour = NA),
        axis.title.y = element_text(size = 24), axis.text.y = element_text(size = 14), axis.text.x = element_text(size = 14))
dev.off()

#Write table to file
AMR_EUR_trait_count_df <- AMR_EUR_trait_count_df[order(AMR_EUR_trait_count_df$residuals, decreasing = TRUE),]
write.table(AMR_EUR_trait_count_df, file=paste0(gene_nom_dir, "AMR_EUR_traits_per_gene.tsv"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# 8) Create a supplemental table of pleiotropy across ancestries ====

#Rename columns in each dataframe
colnames(AFR_EUR_trait_count_df) <- str_replace(colnames(AFR_EUR_trait_count_df), "[.]x", ".EUR") %>% str_replace(pattern = "[.]y", replacement = ".AFR") %>% 
  str_replace(pattern = "residuals", replacement = "AFR_v_EUR_residuals")
colnames(AMR_EUR_trait_count_df) <- str_replace(colnames(AMR_EUR_trait_count_df), "[.]x", ".EUR") %>% str_replace(pattern = "[.]y", replacement = ".AMR") %>% 
  str_replace(pattern = "residuals", replacement = "AMR_v_EUR_residuals")
#Outer join the two non-eur vs eur trait_count_dfs
export_table <- merge(AFR_EUR_trait_count_df, AMR_EUR_trait_count_df, by = "gene", all = TRUE)
export_table <- export_table %>% dplyr::select(gene, traits.AFR, trait_names.AFR, traits.AMR, trait_names.AMR, traits.EUR.x, trait_names.EUR.x, AFR_v_EUR_residuals, AMR_v_EUR_residuals)
colnames(export_table) <- colnames(export_table) %>% str_replace("[.]x", "")
#Write to file
write.table(export_table, file=paste0(gene_nom_dir, "cross_population_pleiotropy_comparison.supplement.tsv"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
