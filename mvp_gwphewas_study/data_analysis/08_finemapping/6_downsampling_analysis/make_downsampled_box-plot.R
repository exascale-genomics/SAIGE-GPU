################################################################################

#make_downsampled_box-plot.R

#The purpose of this script is to analyze the results of our down-sampled 
#analysis in which we compared the EUR and AFR credible set sizes for signals 
#mapped in both ancestries with matched sample-sizes of ~120K:

################################################################################

# 0) Call libraries and set directories ====

#Call libraries
library(ggplot2)
library(tidyverse)
library(scales)
library(ggpmisc)
library(ggpubr)

#Set directory
work_dir <- "C:/Users/mitch/Documents/UPenn/Voight_Lab/MVP_Finemapping/downsampling/"

# 1) Read in merged results file ====

#Read if file of merged signals
merged_raw <- read.csv(paste0(work_dir, "master.signals.downsampled.txt"), sep = "\t", header = TRUE)

#Create columns indicating whether the signal was mapped in each ancestry
merged_append <- merged_raw %>% mutate(AFR=ifelse(is.na(AFR.num_variants) == FALSE, 1, 0),
                                       AMR=ifelse(is.na(AMR.num_variants) == FALSE, 1, 0),
                                       EAS=ifelse(is.na(EAS.num_variants) == FALSE, 1, 0),
                                       EUR=ifelse(is.na(EUR.num_variants) == FALSE, 1, 0))

# 2) Compare Number of Variants in EUR/AFR Shared Signals After EUR downsampling ====

#Make a function that takes 2 ancestries and outputs a comparison box plot
compare_variant_count <- function(merged_append, ANC1, ANC2, work_dir){
  #Filter down for signals found in both ancestries
  merged_ANC1_ANC2 <- merged_append[which(merged_append[,ANC1] == 1 & merged_append[,ANC2] == 1),]
  merged_ANC1_ANC2 <- cbind.data.frame('Number.Variants'=c(merged_ANC1_ANC2[,paste0(ANC1,".num_variants")], merged_ANC1_ANC2[,paste0(ANC2,".num_variants")]),
                                       Ancestry=c(rep(ANC1, nrow(merged_ANC1_ANC2)), rep(ANC2, nrow(merged_ANC1_ANC2))), 
                                       Unique.ID=c(paste0(merged_ANC1_ANC2[,"trait"],".",merged_ANC1_ANC2[,"locus"],".",merged_ANC1_ANC2[,"signal"]), paste0(merged_ANC1_ANC2[,"trait"],".",merged_ANC1_ANC2[,"locus"])))
  merged_ANC1_ANC2 <- merged_ANC1_ANC2[order(merged_ANC1_ANC2$Unique.ID),]
  #Calculate Differences
  ANC1_subset <- subset(merged_ANC1_ANC2,  Ancestry == ANC1, Number.Variants, drop = TRUE)
  ANC2_subset <- subset(merged_ANC1_ANC2,  Ancestry == ANC2, Number.Variants, drop = TRUE)
  differences <- cbind.data.frame(Unique.ID=merged_ANC1_ANC2[which(merged_ANC1_ANC2$Ancestry == ANC1),"Unique.ID"], ANC1_min_ANC2=ANC1_subset - ANC2_subset)
  res <- wilcox.test(differences$ANC1_min_ANC2)
  #Create annotation df
  annotations <- data.frame(xpos = c(-Inf), ypos =  c(Inf), annotateText = c(paste0("Wilcoxon*','~~italic(p)==", str_replace(format(res$p.value, scientific = TRUE, digits = 3), pattern = "e", replacement = "~~x~~10^"))), hjustvar = c(0), vjustvar = c(1))
  #Make difference plot
  temp <- ggplot(data = differences, mapping = aes(y = ANC1_min_ANC2)) + geom_boxplot() + scale_y_continuous(trans=pseudo_log_trans(base=exp(10))) + 
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


