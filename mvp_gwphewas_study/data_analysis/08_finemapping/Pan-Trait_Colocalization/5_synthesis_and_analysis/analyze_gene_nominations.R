################################################################################

#analyze_gene_nominations.R

#The purpose of this script is to study the genes nominated by the ABC/coding
#analysis and identify relationships that will cause a gene to be pleiotropic.

#This code runs locally

################################################################################

# 0) Call libraries and set variables and directories ====

#Call libraries
library(tidyverse)
library(org.Hs.eg.db)
library(GO.db)
library(biomaRt)
library(pbapply)
library(GOSemSim)

#Set directory
gene_nom_dir <- "C:/Users/mitch/Documents/UPenn/Voight_Lab/MVP_Finemapping/gene_nominations/"

#Set trait file
trait_file <- "C:/Users/mitch/Documents/UPenn/Voight_Lab/MVP_Finemapping/Figure_Building_and_Synthesis/MVP_R4.1000G_AGR.GIA.DataDictionary.revised.022023.txt"

#Set a random seed
random_seed=5
set.seed(random_seed)

# 1) Read in files and prune traits for independence ====

#Load file of genetic correlations
load(paste0(gene_nom_dir, "ALL.corr.RData"))

#Read in gene nominations
filtered_merged_nominations_df <- read.table(paste0(gene_nom_dir, "filtered_merged_nominations.tsv"), header = TRUE, sep = "\t")
#Get unique genes
genes_with_noms <- unique(filtered_merged_nominations_df$gene)

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
    trait_pip_df <- temp %>% group_by(trait) %>% summarize(mvp_id=mvp_id[which.max(max_overall_pip)], max_pip=max(max_overall_pip)) %>% as.data.frame()
    rownames(trait_pip_df) <- trait_pip_df$trait
    return(c(gene, trait_pip_df[traits,]))
  } else { #If multiple traits were identified, then we need to see how many were independent
    #Get maximum PIP for each trait
    trait_pip_df <- temp %>% group_by(trait) %>% summarize(mvp_id=mvp_id[which.max(max_overall_pip)], max_pip=max(max_overall_pip)) %>% as.data.frame()
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
#Apply the function
indep_trait_per_gene <- lapply(genes_with_noms, get_indep_traits, 
                                   filtered_merged_nominations_df=filtered_merged_nominations_df, 
                                   coef_p=coef_p, coef_m=coef_m, coef_f=coef_f)
#Convert it to a data frame
indep_trait_per_gene_df <- as.data.frame(matrix(unlist(indep_trait_per_gene), ncol = 4, byrow = TRUE))
colnames(indep_trait_per_gene_df) <- c("gene", "trait", "mvp_id", "max_overall_pip")
split_colon_take_x <- function(string_inp, x){return(str_split(string_inp, pattern = ":")[[1]][x])}
indep_trait_per_gene_df <- cbind.data.frame(indep_trait_per_gene_df, 
                                            chromo=as.numeric(vapply(indep_trait_per_gene_df$mvp_id, FUN = split_colon_take_x, FUN.VALUE = character(1), x=1)), 
                                            pos=as.numeric(vapply(indep_trait_per_gene_df$mvp_id, FUN = split_colon_take_x, FUN.VALUE = character(1), x=2)))
#Write to file
write.table(indep_trait_per_gene_df, file = paste0(gene_nom_dir, "indep_trait_per_gene.tsv"), sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

#convert the table to a gene-based table
trait_count_df <- indep_trait_per_gene_df %>% group_by(gene) %>% summarize(traits=n())

#Make a histogram of the number of independent traits per gene
jpeg(paste0(gene_nom_dir,"traits_per_gene.independent_trait.hist.jpeg"))
trait_count_df %>% 
  ggplot(aes(x=traits)) + geom_histogram(fill = "#E69F00") + 
  scale_y_continuous(breaks = c(1,10,100,1000), labels = c(1,10,100,"1,000"), trans = "pseudo_log") + 
  xlab("Independent Traits per Gene") + ylab("Number of Genes") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), legend.title = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(color="black", size=16, angle=90, hjust = 1), axis.title=element_text(size=24), 
        legend.position = "bottom", axis.text.y = element_text(color="black", size=16), 
        legend.text = element_text(size=13), legend.key=element_blank()) 
dev.off()

#Check number of genes with >=2 independent traits as a % of total
(trait_count_df %>% filter(traits > 1) %>% nrow)/nrow(trait_count_df)
#Check number of genes with >=3 independent traits
trait_count_df %>% filter(traits >= 3) %>% nrow #1,044 genes
trait_count_df %>% filter(traits >= 3) %>% dplyr::select(traits) %>% sum #4241

# 2) Check Pleiotropy of Gene Nominations for Individual GO Term Enrichments Across All Traits ====

# Specify the Ensembl dataset and organism
ensembl_dataset <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Define a function to retrieve GO annotations for a list of ensembl gene ids
get_go_annotations <- function(genes) {
  getBM(attributes = c("ensembl_gene_id", "external_gene_name", "go_id", "description"),
        filters = "external_gene_name",
        values = genes,
        mart = ensembl_dataset)
}

# Set gene list
my_gene_list <- trait_count_df$gene  # Gene IDs in Ensembl format

# Retrieve GO annotations for the specified genes
go_annotations <- get_go_annotations(my_gene_list)
#Filter out lines from go_annotations with missing go categories
go_annotations <- go_annotations %>% filter(go_id != "")
#get unique list of go annotations
unique_go_terms <- go_annotations %>% dplyr::select(go_id) %>% filter(go_id != "") %>% unique %>% as.vector()
unique_go_terms <- unique_go_terms$go_id

#Remove any genes that don't appear at least once from the annotations table
trait_count_df <- trait_count_df %>% filter(gene %in% unique(go_annotations$external_gene_name))

#Make a function that regresses inclusion in a go_term against number of genes
regress_go_on_gene_count <- function(go_term, trait_count_df, go_annotations){
  temp <- go_annotations[go_annotations$go_id == go_term,]
  temp2 <- cbind.data.frame(trait_count_df, in_go=ifelse(trait_count_df$gene %in% temp$external_gene_name, TRUE, FALSE))
  temp3 <- summary(glm(traits ~ in_go, data = temp2, family = "poisson"))
  return(c(temp3$coefficients[2, "Estimate"], temp3$coefficients[2, "Pr(>|z|)"], paste0(temp2[temp2$in_go == TRUE, "gene"], collapse = ",")))
}
#Call function to run regressions
regress_list <- pbvapply(unique_go_terms, FUN = regress_go_on_gene_count, FUN.VALUE = character(3),trait_count_df=trait_count_df, go_annotations=go_annotations)
#convert to data frame
regress_df <- cbind.data.frame(unique_go_terms, t(regress_list))
colnames(regress_df) <- c("GO_Term", "Coefficient", "P-Value", "Genes_in_GO")
regress_df$Coefficient <- as.numeric(regress_df$Coefficient)
regress_df$`P-Value` <- as.numeric(regress_df$`P-Value`)
#Get the names for the go terms
go_terms <- biomaRt::select(GO.db, keys = regress_df$GO_Term, columns = c("GOID", "TERM", "ONTOLOGY"))
#Add the adjusted p-values and append the GO Term Names to the regression_df
regress_df <- cbind.data.frame(regress_df, GO_Desc=go_terms$TERM, padj=p.adjust(regress_df$`P-Value`, method = "BH"))
colnames(regress_df) <- c("GO_ID" , "Coefficient", "P-Value", "Genes_in_GO", "GO_TERM", "P-Value_BH")
regress_df <- regress_df[,c("GO_ID" , "GO_TERM", "Coefficient", "P-Value", "P-Value_BH", "Genes_in_GO")]
#Write data frame to file
write.table(x = regress_df, file = paste0(gene_nom_dir, "pleiotropy_regression_dataframe.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

#Filter for significance results
sig_regress_df <- regress_df[regress_df[,"P-Value_BH"] < 0.05,]

#Write an output file useable for Go-Figure
go_figure_out <- sig_regress_df[,c("GO_ID", "P-Value_BH", "Coefficient")]
colnames(go_figure_out) <- c("GOterm", "regression_P-value", "Beta")
write.table(go_figure_out, file = paste0(gene_nom_dir, "go-figure_input.tsv"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

#Write an output suitable for ClueGO
cluego_out <- sig_regress_df[,c("GO_ID", "P-Value_BH")]
colnames(cluego_out) <- c("GOID.PathwayID", "P value")
write.table(cluego_out, file = paste0(gene_nom_dir, "cluego.pleiotropy.input.tsv"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

#Execute the GO-Figure Mapping in the command line
#python C:/Users/mitch/Documents/UPenn/Voight_Lab/MVP_Finemapping/GO-Figure/gofigure.py 
# -i C:/Users/mitch/Documents/UPenn/Voight_Lab/MVP_Finemapping/gene_nominations/go-figure_input.tsv 
# -j standard-plus
# -si 0.1
# -p Blues
# -e 60
# -o pleiotropy_go

#python C:/Users/mitch/Documents/UPenn/Voight_Lab/MVP_Finemapping/GO-Figure/gofigure.py 
# -i C:/Users/mitch/Documents/UPenn/Voight_Lab/MVP_Finemapping/gene_nominations/go-figure_input.tsv 
# -j standard-plus 
# -c user
# -s pval
# -si 0.1
# -p bwr
# -e 60
# -o pleiotropy_go_beta

# 3) Check Pleiotropy of Gene Nominations for Individual GO Term Enrichments Across Binary Traits ====

#Read in the file of trait info
trait_raw <- read.csv(trait_file, header = TRUE, sep = "\t")
rownames(trait_raw) <- trait_raw$Trait
#Split the filtered_merged nominations by trait type
binary_nominations_df <- filtered_merged_nominations_df %>% filter(trait_raw[trait, "Trait_Type"] == "binary")
quantitative_nominations_df <- filtered_merged_nominations_df %>% filter(trait_raw[trait, "Trait_Type"] == "quantitative")

#Prepare the table of independent traits for the binary traits
genes_with_noms_binary <- unique(binary_nominations_df$gene)
indep_trait_per_gene_binary <- lapply(genes_with_noms_binary, get_indep_traits, 
                               filtered_merged_nominations_df=binary_nominations_df, 
                               coef_p=coef_p, coef_m=coef_m, coef_f=coef_f)
#Convert it to a data frame
indep_trait_per_gene_df_binary <- as.data.frame(matrix(unlist(indep_trait_per_gene_binary), ncol = 4, byrow = TRUE))
colnames(indep_trait_per_gene_df_binary) <- c("gene", "trait", "mvp_id", "max_overall_pip")
indep_trait_per_gene_df_binary <- cbind.data.frame(indep_trait_per_gene_df_binary, 
                                            chromo=as.numeric(vapply(indep_trait_per_gene_df_binary$mvp_id, FUN = split_colon_take_x, FUN.VALUE = character(1), x=1)), 
                                            pos=as.numeric(vapply(indep_trait_per_gene_df_binary$mvp_id, FUN = split_colon_take_x, FUN.VALUE = character(1), x=2)))
#convert the table to a gene-based table
binary_trait_count_df <- indep_trait_per_gene_df_binary %>% group_by(gene) %>% summarize(traits=n())

# Set gene list
binary_gene_list <- binary_trait_count_df$gene  # Gene IDs in Ensembl format
# Retrieve GO annotations for the specified genes
binary_go_annotations <- get_go_annotations(binary_gene_list)
#Filter out lines from go_annotations with missing go categories
binary_go_annotations <- binary_go_annotations %>% filter(go_id != "")
#get unique list of go annotations
binary_unique_go_terms <- binary_go_annotations %>% dplyr::select(go_id) %>% filter(go_id != "") %>% unique %>% as.vector()
binary_unique_go_terms <- binary_unique_go_terms$go_id
#Remove any genes that don't appear at least once from the annotations table
binary_trait_count_df <- binary_trait_count_df %>% filter(gene %in% unique(binary_go_annotations$external_gene_name))

#Call function to run regressions
binary_regress_list <- pbvapply(binary_unique_go_terms, FUN = regress_go_on_gene_count, FUN.VALUE = character(3), trait_count_df=binary_trait_count_df, go_annotations=binary_go_annotations)
#convert to data frame
binary_regress_df <- cbind.data.frame(binary_unique_go_terms, t(binary_regress_list))
colnames(binary_regress_df) <- c("GO_Term", "Coefficient", "P-Value", "Genes_in_GO")
binary_regress_df$Coefficient <- as.numeric(binary_regress_df$Coefficient)
binary_regress_df$`P-Value` <- as.numeric(binary_regress_df$`P-Value`)
#Get the names for the go terms
binary_go_terms <- biomaRt::select(GO.db, keys = binary_regress_df$GO_Term, columns = c("GOID", "TERM", "ONTOLOGY"))
#Add the adjusted p-values and append the GO Term Names to the regression_df
binary_regress_df <- cbind.data.frame(binary_regress_df, GO_Desc=binary_go_terms$TERM, padj=p.adjust(binary_regress_df$`P-Value`, method = "BH"))
colnames(binary_regress_df) <- c("GO_ID" , "Coefficient", "Genes_in_GO", "P-Value", "GO_TERM", "P-Value_BH")
binary_regress_df <- binary_regress_df[,c("GO_ID" , "GO_TERM", "Coefficient", "P-Value", "P-Value_BH", "Genes_in_GO")]
#Write data frame to file
write.table(x = binary_regress_df, file = paste0(gene_nom_dir, "binary_pleiotropy_regression_dataframe.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
#Filter for significance results
binary_sig_regress_df <- binary_regress_df %>% filter(`P-Value_BH` < 0.05) %>% as.data.frame()

#Write an output file useable for Go-Figure
binary_go_figure_out <- binary_sig_regress_df[,c("GO_ID", "P-Value_BH", "Coefficient")]
colnames(binary_go_figure_out) <- c("GOterm", "regression_P-value", "Beta")
write.table(binary_go_figure_out, file = paste0(gene_nom_dir, "binary.go-figure_input.tsv"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

#Execute the GO-Figure Mapping in the command line
#python C:/Users/mitch/Documents/UPenn/Voight_Lab/MVP_Finemapping/GO-Figure/gofigure.py 
# -i C:/Users/mitch/Documents/UPenn/Voight_Lab/MVP_Finemapping/gene_nominations/binary.go-figure_input.tsv 
# -j standard-plus
# -si 0.1
# -p Blues
# -e 60
# -o pleiotropy_go_binary

# 4) Check Pleiotropy of Gene Nominations for Individual GO Term Enrichments Across Quantitative Traits ====

#Prepare the table of independent traits for the quantitative traits
genes_with_noms_quantitative <- unique(quantitative_nominations_df$gene)
indep_trait_per_gene_quantitative <- lapply(genes_with_noms_quantitative, get_indep_traits, 
                                      filtered_merged_nominations_df=quantitative_nominations_df, 
                                      coef_p=coef_p, coef_m=coef_m, coef_f=coef_f)
#Convert it to a data frame
indep_trait_per_gene_df_quantitative <- as.data.frame(matrix(unlist(indep_trait_per_gene_quantitative), ncol = 4, byrow = TRUE))
colnames(indep_trait_per_gene_df_quantitative) <- c("gene", "trait", "mvp_id", "max_overall_pip")
indep_trait_per_gene_df_quantitative <- cbind.data.frame(indep_trait_per_gene_df_quantitative, 
                                                   chromo=as.numeric(vapply(indep_trait_per_gene_df_quantitative$mvp_id, FUN = split_colon_take_x, FUN.VALUE = character(1), x=1)), 
                                                   pos=as.numeric(vapply(indep_trait_per_gene_df_quantitative$mvp_id, FUN = split_colon_take_x, FUN.VALUE = character(1), x=2)))
#convert the table to a gene-based table
quantitative_trait_count_df <- indep_trait_per_gene_df_quantitative %>% group_by(gene) %>% summarize(traits=n())

# Set gene list
quantitative_gene_list <- quantitative_trait_count_df$gene  # Gene IDs in Ensembl format
# Retrieve GO annotations for the specified genes
quantitative_go_annotations <- get_go_annotations(quantitative_gene_list)
#Filter out lines from go_annotations with missing go categories
quantitative_go_annotations <- quantitative_go_annotations %>% filter(go_id != "")
#get unique list of go annotations
quantitative_unique_go_terms <- quantitative_go_annotations %>% dplyr::select(go_id) %>% filter(go_id != "") %>% unique %>% as.vector()
quantitative_unique_go_terms <- quantitative_unique_go_terms$go_id
#Remove any genes that don't appear at least once from the annotations table
quantitative_trait_count_df <- quantitative_trait_count_df %>% filter(gene %in% unique(quantitative_go_annotations$external_gene_name))

#Call function to run regressions
quantitative_regress_list <- pbvapply(quantitative_unique_go_terms, FUN = regress_go_on_gene_count, FUN.VALUE = character(3), trait_count_df=quantitative_trait_count_df, go_annotations=quantitative_go_annotations)
#convert to data frame
quantitative_regress_df <- cbind.data.frame(quantitative_unique_go_terms, t(quantitative_regress_list))
colnames(quantitative_regress_df) <- c("GO_Term", "Coefficient", "P-Value", "Genes_in_GO")
quantitative_regress_df$Coefficient <- as.numeric(quantitative_regress_df$Coefficient)
quantitative_regress_df$`P-Value` <- as.numeric(quantitative_regress_df$`P-Value`)
#Get the names for the go terms
quantitative_go_terms <- biomaRt::select(GO.db, keys = quantitative_regress_df$GO_Term, columns = c("GOID", "TERM", "ONTOLOGY"))
#Add the adjusted p-values and append the GO Term Names to the regression_df
quantitative_regress_df <- cbind.data.frame(quantitative_regress_df, GO_Desc=quantitative_go_terms$TERM, padj=p.adjust(quantitative_regress_df$`P-Value`, method = "BH"))
colnames(quantitative_regress_df) <- c("GO_ID" , "Coefficient", "Genes_in_GO", "P-Value", "GO_TERM", "P-Value_BH")
quantitative_regress_df <- quantitative_regress_df[,c("GO_ID" , "GO_TERM", "Coefficient", "P-Value", "P-Value_BH", "Genes_in_GO")]
#Write data frame to file
write.table(x = quantitative_regress_df, file = paste0(gene_nom_dir, "quantitative_pleiotropy_regression_dataframe.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
#Filter for significance results
quantitative_sig_regress_df <- quantitative_regress_df %>% filter(`P-Value_BH` < 0.05) %>% as.data.frame()

#Write an output file useable for Go-Figure
quantitative_go_figure_out <- quantitative_sig_regress_df[,c("GO_ID", "P-Value_BH", "Coefficient")]
colnames(quantitative_go_figure_out) <- c("GOterm", "regression_P-value", "Beta")
write.table(quantitative_go_figure_out, file = paste0(gene_nom_dir, "quantitative.go-figure_input.tsv"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

#Execute the GO-Figure Mapping in the command line
#python C:/Users/mitch/Documents/UPenn/Voight_Lab/MVP_Finemapping/GO-Figure/gofigure.py 
# -i C:/Users/mitch/Documents/UPenn/Voight_Lab/MVP_Finemapping/gene_nominations/quantitative.go-figure_input.tsv 
# -j standard-plus
# -si 0.1
# -p Blues
# -e 60
# -o pleiotropy_go_quantitative

# 5) Rerun Pleiotropy Checks Across All Traits after Removing APOE ====

#Remove APOE from trait_count dataframe
trait_count_df_minus_apoe <- trait_count_df %>% filter(gene != "APOE")

# Set gene list
my_gene_list_minus_apoe <- trait_count_df_minus_apoe$gene  # Gene IDs in Ensembl format

# Retrieve GO annotations for the specified genes
go_annotations_minus_apoe <- get_go_annotations(my_gene_list_minus_apoe)
#Filter out lines from go_annotations with missing go categories
go_annotations_minus_apoe <- go_annotations_minus_apoe %>% filter(go_id != "")
#get unique list of go annotations
unique_go_terms_minus_apoe <- go_annotations_minus_apoe %>% dplyr::select(go_id) %>% filter(go_id != "") %>% unique %>% as.vector()
unique_go_terms_minus_apoe <- unique_go_terms_minus_apoe$go_id

#Call function to run regressions
regress_list_minus_apoe <- pbvapply(unique_go_terms_minus_apoe, FUN = regress_go_on_gene_count, FUN.VALUE = character(3),trait_count_df=trait_count_df_minus_apoe, go_annotations=go_annotations_minus_apoe)
#convert to data frame
regress_df_minus_apoe <- cbind.data.frame(unique_go_terms_minus_apoe, t(regress_list_minus_apoe))
colnames(regress_df_minus_apoe) <- c("GO_Term", "Coefficient", "P-Value", "Genes_in_GO")
regress_df_minus_apoe$Coefficient <- as.numeric(regress_df_minus_apoe$Coefficient)
regress_df_minus_apoe$`P-Value` <- as.numeric(regress_df_minus_apoe$`P-Value`)
#Get the names for the go terms
go_terms_minus_apoe <- biomaRt::select(GO.db, keys = regress_df_minus_apoe$GO_Term, columns = c("GOID", "TERM", "ONTOLOGY"))
#Add the adjusted p-values and append the GO Term Names to the regression_df
regress_df_minus_apoe <- cbind.data.frame(regress_df_minus_apoe, GO_Desc=go_terms_minus_apoe$TERM, padj=p.adjust(regress_df_minus_apoe$`P-Value`, method = "BH"))
colnames(regress_df_minus_apoe) <- c("GO_ID" , "Coefficient", "P-Value", "Genes_in_GO", "GO_TERM", "P-Value_BH")
regress_df_minus_apoe <- regress_df_minus_apoe[,c("GO_ID" , "GO_TERM", "Coefficient", "P-Value", "P-Value_BH", "Genes_in_GO")]
#Write data frame to file
write.table(x = regress_df_minus_apoe, file = paste0(gene_nom_dir, "pleiotropy_regression_dataframe.minus_apoe.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

#Filter for significance results
sig_regress_df_minus_apoe <- regress_df_minus_apoe[regress_df_minus_apoe[,"P-Value_BH"] < 0.05,]

#Write an output file useable for Go-Figure
go_figure_out_minus_apoe <- sig_regress_df_minus_apoe[,c("GO_ID", "P-Value_BH", "Coefficient")]
colnames(go_figure_out_minus_apoe) <- c("GOterm", "regression_P-value", "Beta")
write.table(go_figure_out_minus_apoe, file = paste0(gene_nom_dir, "go-figure_input.minus_apoe.tsv"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

#Execute the GO-Figure Mapping in the command line
#python C:/Users/mitch/Documents/UPenn/Voight_Lab/MVP_Finemapping/GO-Figure/gofigure.py 
# -i C:/Users/mitch/Documents/UPenn/Voight_Lab/MVP_Finemapping/gene_nominations/go-figure_input.minus_apoe.tsv 
# -j standard-plus
# -si 0.1
# -p Blues
# -e 60
# -o pleiotropy_go_minus_apoe

# 6) Check Pleiotropy Relationship with Total Number of GO Term Annotations ====

#Filter out empty go categories from the list of results pulled from BiomaRt and collapse counts per gene
count_go_df <- go_annotations %>% filter(go_id != "") %>% group_by(external_gene_name) %>% 
  summarise(num_go_terms = n())
colnames(count_go_df) <- c("gene", "num_go_terms")
count_go_df <- inner_join(count_go_df, trait_count_df, by = "gene")

# Function to convert a string to Unicode superscript
convert_to_superscript <- function(input_str) {
  digits <- strsplit(as.character(input_str), '')[[1]]
  superscript_map <- c('0' = '⁰', '1' = '¹', '2' = '²', '3' = '³', '4' = '⁴', 
                       '5' = '⁵', '6' = '⁶', '7' = '⁷', '8' = '⁸', '9' = '⁹', '-' = '\U207B')
  superscripted_digits <- sapply(digits, function(digit) ifelse(digit %in% names(superscript_map), superscript_map[digit], digit))
  return(paste(superscripted_digits, collapse = ''))
}


#Regress traits back onto go terms
trait_on_go_lm <- glm(traits ~ num_go_terms, data = count_go_df, family = "poisson")
#Make a dataframe for the regression line on a log scale
regression_line_df <- seq(min(count_go_df$num_go_terms), max(count_go_df$num_go_terms), 0.01)
regression_line_df <- cbind.data.frame(num_go_terms=regression_line_df, traits=predict(trait_on_go_lm, newdata = data.frame(num_go_terms = regression_line_df), type = "response"))
#print p-value
summary_lm <- summary(trait_on_go_lm)
regression_pval <- summary_lm$coefficients[2,"Pr(>|z|)"]
regression_pval_formatted <- format(regression_pval, digits = 3, scientific = TRUE, style = "f")
# Extract the components of scientific notation
pval_components <- strsplit(regression_pval_formatted, "e")[[1]]
mantissa <- pval_components[1]
exponent <- as.numeric(pval_components[2])
formatted_exponent <- convert_to_superscript(exponent)
# Create the formatted p-value string with superscript exponent
formatted_pval <- paste("P-Value = ", mantissa, " x 10", formatted_exponent, sep = "")
#Create annotation data frame
coefficients <- coef(trait_on_go_lm)
equation <- paste("  ln(Traits) =", round(coefficients["(Intercept)"], 2), "+", round(coefficients["num_go_terms"], 3), " * (GO Terms)\n ", formatted_pval)
#plot Regression
jpeg(paste0(gene_nom_dir,"traits_v_go_per_gene.scatter.jpeg"), width = 720, height=720)
ggplot(count_go_df, mapping = aes(y = traits, x = num_go_terms)) + geom_point() +
  geom_line(data = regression_line_df, mapping = aes(y = traits, x=num_go_terms), color = "Red",  linewidth = 3) + 
  annotate("text", x = -Inf, y = Inf, label = equation, hjust = 0, vjust = 1, size = 10, color="Red") +
  xlab("GO Terms per Gene") + ylab("Indep. Traits per Gene") +
  scale_y_continuous(breaks = c(1,10,100), limits = c(0, max(count_go_df$traits)), labels = c(1,10,100), trans = "pseudo_log") + 
  scale_x_continuous(breaks = c(2,20,200), limits = c(0, max(count_go_df$num_go_terms)), labels = c(2,20,200), trans = "pseudo_log") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), legend.title = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(color="black", size=18, angle=90, hjust = 1), axis.title=element_text(size=28), 
        legend.position = "bottom", axis.text.y = element_text(color="black", size=18), 
        legend.text = element_text(size=13), legend.key=element_blank()) 
dev.off()



