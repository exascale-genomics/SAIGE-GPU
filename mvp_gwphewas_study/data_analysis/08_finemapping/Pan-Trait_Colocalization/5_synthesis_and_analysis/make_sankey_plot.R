################################################################################

#make_sankey_plot.R

#The purpose of this script is to make the pseudo Sankey plot that summarizes 
#the fine-mapping pipeline


################################################################################

# 0) Call libraries and directories ====

#Call libraries
library(ggplot2)
library(stringr)
library(pbapply)
library(ComplexUpset)
library(plyr)
library(ggrepel)

#Set working directory
work_dir <- "C:/Users/mitch/Documents/UPenn/Voight_Lab/MVP_Finemapping/Figure_Building_and_Synthesis/"

# 1) Read in files for the first column of the plot ====

#set file locations
locus_file <- paste0(work_dir, "final_output.grch37.high_thresh.broad.csv")
trait_file <- paste0(work_dir, "MVP_R4.1000G_AGR.GIA.DataDictionary.revised.022023.txt")
neg_control_file <- paste0(work_dir, "negative_control_list.txt")

#Read in the files
locus_raw <- read.table(locus_file, header = FALSE, sep = ",")
colnames(locus_raw) <- c("chr", "start", "end", "length", "traits")
trait_raw <- read.csv(trait_file, header = TRUE, sep = "\t")
neg_control_raw <- read.csv(neg_control_file, header = TRUE, sep = "\t")

#Check for intersection between mapped traits and expected negative controls
mapped_neg_controls <- intersect(unique(locus_raw$traits), neg_control_raw$Trait) #There are 56 of these. This is a problem!

#Remove the HIV traits we dropped from the trait_raw and extract the complete list of traits
all_traits <- trait_raw[which(!(trait_raw$Trait %in% c("Phe_071", "Phe_071_1"))), "Trait"]
#trait_category <- ifelse(all_traits %in% unique(unlist(str_split(locus_raw$traits, "/"))), "1 or More Loci Defined", ifelse(all_traits %in% neg_control_raw$Trait, "Expected Negative Control", "No Loci Detected"))
trait_category <- ifelse(all_traits %in% unique(unlist(str_split(locus_raw$traits, "/"))), "1 or More Loci Defined", "No Loci\nDetected")

#Get unique categories for first column of plot
loci_categories <- unique(trait_category)

#Set colors
shades_o_red <- c("firebrick1", "firebrick2", "firebrick3", "firebrick4")
shades_o_green <- c("#03AC13", "#3DED97")
shades_o_blue <- c("steelblue1", "steelblue2", "steelblue3", "steelblue4")
#Create data frame for loci column
loci_df=data.frame(x1=c(1), x2=c(2), y1=c(0), y2=c(sum(trait_category=="1 or More Loci Defined")), t=c(shades_o_green[1]), r=c("1 or More\nLoci Defined"))
for (i in seq(1,length(loci_categories[loci_categories != "1 or More Loci Defined"]),1)) {
  category = loci_categories[loci_categories != "1 or More Loci Defined"][i]
  loci_df <- rbind.data.frame(loci_df,
                              data.frame(x1=c(1), x2=c(2), y1=c(max(loci_df$y2)), y2=c(max(loci_df$y2) + sum(trait_category==category)), t=c(shades_o_red[i]), r=c(category)))
}

# 2) Get Data for Second Column of Figure ====

#Read in merged fine-mapping results file to get the counts of mapped loci-trait pairs
merged_raw <- read.csv(paste0(work_dir, "master.signals.txt"), sep = "\t", header = TRUE)

#Input completion data from Alex (Make sure Mapped Signals are up front!!!!)
trait_loci_raw <- c(length(unique(paste0(merged_raw$trait, ".", merged_raw$locus))), 
                    3016, 
                    71)
names(trait_loci_raw) <- c("Mapped /\n1 or More\nSignals", 
                           "Mapped / No\nSets Found",  
                           "Unable to Map")
#Create dataframe
loci_trait_df=data.frame(x1=c(3), x2=c(4), y1=c(0), y2=c(trait_loci_raw[1]), t=c(shades_o_green[1]), r=c(names(trait_loci_raw)[1]))
loci_trait_df=rbind.data.frame(loci_trait_df,
                               data.frame(x1=c(3), x2=c(4), y1=c(max(loci_trait_df$y2)), y2=c(max(loci_trait_df$y2) + trait_loci_raw[2]), t=c(shades_o_green[2]), r=c(names(trait_loci_raw)[2])))
for (i in seq(1,length(trait_loci_raw[substr(names(trait_loci_raw), 1, 6) != "Mapped"]),1)) {
  j = i + 2
  loci_trait_df <- rbind.data.frame(loci_trait_df,
                              data.frame(x1=c(3), x2=c(4), y1=c(max(loci_trait_df$y2)), y2=c(max(loci_trait_df$y2) + trait_loci_raw[j]), t=c(shades_o_red[i]), r=c(names(trait_loci_raw)[j])))
}

#Calculate number of loci_trait pairs
num_loci_trait <- sum(trait_loci_raw)

#Modify labels
loci_trait_df[3,"r"] <- "Unable\nto Map"

# 3) Calculate the sharing data ====

#Clean-up
remove(trait_raw, neg_control_raw, locus_raw)

#Calculate number of signals first
num_signals <- nrow(merged_raw)
#Count ancestries
merged_raw <- cbind.data.frame(merged_raw, num_anc=4-rowSums(is.na(merged_raw[,c("AFR.num_variants", "AMR.num_variants", "EAS.num_variants", "EUR.num_variants")])))
anc_raw <- table(merged_raw$num_anc)
remove(merged_raw)

#Make plotting df
anc_df=data.frame(x1=c(5), x2=c(6), y1=c(0), y2=c(anc_raw[1]), t=c(shades_o_blue[1]), r=c(paste0(1, " Population")))
for (i in seq(2, length(anc_raw), 1)) {
  anc_df <- rbind.data.frame(anc_df,
                             data.frame(x1=c(5), x2=c(6), y1=c(max(anc_df$y2)), y2=c(max(anc_df$y2) + anc_raw[i]), t=c(shades_o_blue[i]), r=c(paste0(i, " Populations"))))
}
#Modify labels
anc_df[3,"r"] <- "3 or 4\nPopulations"

# 4) Make the Pseudo Sankey Plot ====

#Rescale the y-axis in the data frames
loci_df$y1 = loci_df$y1/max(loci_df$y2)
loci_df$y2 = loci_df$y2/max(loci_df$y2)
loci_trait_df$y1 = loci_trait_df$y1/max(loci_trait_df$y2)
loci_trait_df$y2 = loci_trait_df$y2/max(loci_trait_df$y2)
anc_df$y1 = anc_df$y1/max(anc_df$y2)
anc_df$y2 = anc_df$y2/max(anc_df$y2)

#Extract info for polygons
trapezoid_df <- data.frame(id=rep("loci_to_loci_trait", 4), x=c(2,2,3,3), y=c(loci_df[1,"y2"],0,0,1), fill=shades_o_green[1])
trapezoid_df <- rbind.data.frame(trapezoid_df, data.frame(id=rep("loci_trait_to_anc", 4), x=c(4,4,5,5), y=c(loci_trait_df[1,"y2"],0,0,1), fill=shades_o_green[1]))

#Make figure
jpeg(paste0(work_dir, "finemapping_overview.sankey.jpeg"), width = 8, height = 8, units = 'in', res = 500)
ggplot() + 
  scale_x_continuous(name="x") + 
  scale_y_continuous(name="y") +
  geom_rect(data=loci_df, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), fill=loci_df$t, color="black") +
  geom_text(data=loci_df, aes(x=(x1+(x2-x1)/2), y=y1+(y2-y1)/2, label=r), size=5) + 
  geom_rect(data=loci_trait_df, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), fill=loci_trait_df$t, color="black") +
  geom_text(data=loci_trait_df[1:2,], aes(x=(x1+(x2-x1)/2), y=y1+(y2-y1)/2, label=r), size=5) + 
  geom_text(data=loci_trait_df[3,], aes(x=(x1+(x2-x1)/2)-1, y=y1+(y2-y1)/2, label=r), size=5, color=loci_trait_df[3,"t"]) +
  geom_segment(data=loci_trait_df, x=2.9, xend=3, y=(loci_trait_df[3,"y1"]+(loci_trait_df[3,"y2"]-loci_trait_df[3,"y1"])/2), yend=(loci_trait_df[3,"y1"]+(loci_trait_df[3,"y2"]-loci_trait_df[3,"y1"])/2), color=shades_o_red[2]) + 
  geom_rect(data=anc_df, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), fill=anc_df$t, color="black") +
  geom_text(data=anc_df[1:2,], aes(x=(x1+(x2-x1)/2), y=y1+(y2-y1)/2, label=r), size=5) + 
  geom_text(data=anc_df[3,], aes(x=(x1+(x2-x1)/2)-0.95, y=y1+(y2-y1)/2, label=r), size=5, color=anc_df[3,"t"]) +
  geom_segment(data=anc_df, x=4.9, xend=5, y=(anc_df[3,"y1"]+(anc_df[3,"y2"]-anc_df[3,"y1"])/2), yend=(anc_df[3,"y1"]+(anc_df[3,"y2"]-anc_df[3,"y1"])/2), color=shades_o_blue[3]) + 
  geom_polygon(data=trapezoid_df, aes(x=x, y=y, group = id), fill=trapezoid_df$fill, alpha=0.2) + 
  annotate("text", x = 1.5, y = 1.07, label = paste0(format(length(all_traits), big.mark=","), " traits"), size=6) + 
  annotate("text", x = 3.5, y = 1.07, label = paste0(format(num_loci_trait, big.mark=","), " locus-\ntrait pairs"), size=6, lineheight=0.9) + 
  annotate("text", x = 5.5, y = 1.07, label = paste0(format(num_signals, big.mark=","), " signals"), size=6) + 
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y = element_blank(),
        axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.x = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())
dev.off()

