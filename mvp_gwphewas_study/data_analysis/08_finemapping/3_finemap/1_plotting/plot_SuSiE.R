################################################################################

#plot_SuSiE.R

#This script is designed to plot the results of fine mapping a particular locus 
#with SuSiE. It will take a .rds file output from SuSiE and make a p-value plot.
#If the option is selected, it will instead make absolute or preceding residual
#association plots.

################################################################################

# 0) Load Needed Libraries ====

library(coloc)
library(ggplot2)
library(stringr)

# 1) Read in Needed Arguments from Command ====

#Need to read in the following: 
  # 1) File of SuSiE fine-mapping results
  # 2) output file directory
  # 3) plot type (non-residual, preceding_residual, absolute_residual)
  # 4) association type for non-residual plot threshold (gwas, absolute, preceding)
  # 5) purity threshold for plotting non-residual plots
  # 6) Necessary max association p-value for plotting credible sets on non-residual plots (defaults to 1, i.e. no maximum value)


args <- commandArgs(trailingOnly = TRUE);
if (length(args) == 6) {
  #Print confirmation output
  print("Necessary Parameters Input")
  #Set variables from the arguments in this case
  susie_results <- args[1]
  out_dir <- args[2]
  plot_type <- args[3]
  assoc_type <- args[4]
  purity <- as.numeric(args[5])
  max_assoc <- as.numeric(args[6])
} else if(length(args) == 3) {
  susie_results <- args[1]
  out_dir <- args[2]
  plot_type <- args[3]
  assoc_type <- "gwas"
  purity <- 0.5
  max_assoc <- 1
} else {
  print("Usage: %> Rscript call_susie.R susie_results output_folder plot_type [assoc_type_for_non-residual] [purity_for_non-residual] [max_assoc_non-residual]");
  print("WARNING: If any optional parameters are input, then all must be.")
  quit(save="no");
}

# 2) Read in file ====

#Define out_loc
if(substr(out_dir, nchar(out_dir), nchar(out_dir)) != "/"){
  out_dir <- paste0(out_dir, "/")
} else {
  out_dir <- out_dir
}
if (file.exists(out_dir) == FALSE) {
  print("ERROR: Invalid output folder given.")
  quit(save="no");
}

#Get file prefix
file_prefix = as.vector(strsplit(susie_results, "/"))
file_prefix = file_prefix[[1]][length(file_prefix[[1]])]
file_prefix = gsub(".rds", "", file_prefix)

#Read in rds file
fitted_rss1 <- readRDS(susie_results)

#Check whether RDS file was a dummy RDS file for an umappable region
if (is.null(fitted_rss1)) {
  print(paste0("WARNING: attempted to plot an unmappable region, ", file_prefix))
  quit(save="no");
}

# 3) Make Plots ====

#Set colors for plot
#Set colors
coloor <- c(
  "dodgerblue2",
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown"
)
#Get credible sets
cred_sets=names(fitted_rss1$sets$cs)
cred_sets=vapply(cred_sets, str_replace, character(1), pattern="L", replacement="")
cred_sets=vapply(cred_sets, as.numeric, 0)
cred_sets=sort(cred_sets)
#Verify that at least some credible sets were identified
no_cred_sets <- FALSE
if (length(cred_sets) == 0){
  no_cred_sets <- TRUE
}

## Make plot(s) ##
if (plot_type == "non-residual") {
  print("NOTE: Non-residual plot type selected.")
  plot_loc <- paste0(out_dir, file_prefix, ".purity-", purity, ".", assoc_type, "-", max_assoc, ".pvalue.jpg")
  #Check for no credible sets possibility
  if (no_cred_sets == TRUE) {
    print(paste0("WARNING: attempted to plot a region with no identified credible sets, ", file_prefix))
    #Make dataframe for plotting
    plot_df <- cbind.data.frame(SNP=names(fitted_rss1$pip), CHR=rep(fitted_rss1$chr, length(fitted_rss1$bp)), BP=fitted_rss1$bp, neglogp_gwas=fitted_rss1$neglogp_gwas)
    #Make P-Value Plot
    temp_plot <- ggplot(plot_df) + xlab(paste0("Chromosome ", as.character(median(plot_df$CHR)))) + 
      ylab("-log(P-Value)") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                    panel.background = element_blank(), axis.line = element_line(colour = "black"),
                                    legend.position = "none") + 
      geom_point(aes(x = BP, y = neglogp_gwas), color = "gray87")
    jpeg(file=plot_loc)
    print(temp_plot)
    dev.off()
    quit(save="no");
  }
  #Make dataframe for plotting
  plot_df <- cbind.data.frame(SNP=names(fitted_rss1$pip), CHR=rep(fitted_rss1$chr, length(fitted_rss1$bp)), BP=fitted_rss1$bp, log10p=fitted_rss1$neglogp_gwas, 
                              fitted_rss1$neglogp_absolute, fitted_rss1$neglogp_precede)
  #Change colnames
  colnames(plot_df) <- c("SNP", "CHR", "BP", "neglogp_gwas", paste0("neglogp_absolute_", cred_sets), paste0("neglogp_precede_", cred_sets))
  #Identify pure credible sets
  pure_cred_sets <- cred_sets[fitted_rss1$sets$purity[paste0("L",cred_sets),"min.abs.corr"] >= purity]
  #Identify the significant credible sets
  temp_sets <- fitted_rss1$sets$cs[names(pure_cred_sets)]
  if(assoc_type == 'gwas'){
    check_sig_gwas <- function(temp_set, neglogp_gwas, max_assoc){
      return(max(neglogp_gwas[temp_set] >= -log10(max_assoc)))
    }
    sig_cred_sets = pure_cred_sets[unlist(lapply(temp_sets, check_sig_gwas, neglogp_gwas=fitted_rss1$neglogp_gwas, max_assoc=max_assoc))[names(pure_cred_sets)] == 1]
  } else if (assoc_type == 'absolute'){
    check_sig_residual <- function(set_name, neglogp_residual, max_assoc, temp_sets){
      temp_set <- temp_sets[[set_name]]
      set_name <- gsub(pattern = "L", replacement = "", x = set_name)
      return(max(neglogp_residual[temp_set, set_name] >= -log10(max_assoc)))
    }
    sig_cred_sets = pure_cred_sets[unlist(lapply(names(temp_sets), check_sig_residual, neglogp_residual=fitted_rss1$neglogp_absolute, max_assoc=max_assoc, temp_sets=temp_sets)) == 1]
  } else {
    check_sig_residual <- function(set_name, neglogp_residual, max_assoc, temp_sets){
      temp_set <- temp_sets[[set_name]]
      set_name <- gsub(pattern = "L", replacement = "", x = set_name)
      return(max(neglogp_residual[temp_set, set_name] >= -log10(max_assoc)))
    }
    sig_cred_sets = pure_cred_sets[unlist(lapply(names(temp_sets), check_sig_residual, neglogp_residual=fitted_rss1$neglogp_precede, max_assoc=max_assoc, temp_sets=temp_sets)) == 1]
  }
  #Bind on credible set info to plot_df
  sig_cred_set_list <- fitted_rss1$sets$cs[names(sig_cred_sets)]
  #Check first if any credible sets were identified that pass the thresholds
  if(length(sig_cred_set_list) > 0){
    temp <- strsplit(names(unlist(sig_cred_set_list)), split = "[.]")
    temp_list <- list()
    for (i in 1:length(temp)) {
      temp_list[temp[[i]][2]] <- gsub(x = temp[[i]][1], pattern = "L", replacement = "")
    } 
  } else { #Set a dummy temp_list variable
    temp_list <- list()
    print(paste0("WARNING: No credible sets identified passing both purity and significance thresholds for ", file_prefix))
  }
  plot_df <- cbind.data.frame(plot_df, cs=factor(ifelse(plot_df$SNP %in% names(temp_list), temp_list[plot_df$SNP], 0), levels = as.character(seq(0, max(cred_sets), 1))))
  #Make P-Value Plot
  temp_plot <- ggplot(plot_df) + xlab(paste0("Chromosome ", as.character(median(plot_df$CHR)))) + 
    ylab("-log(P-Value)") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                 panel.background = element_blank(), axis.line = element_line(colour = "black"),
                                 legend.position = "none") + 
    scale_color_manual(values=coloor[1:length(levels(plot_df$cs))]) + 
    geom_point(aes(x = BP, y = neglogp_gwas), color = "gray87") +
    geom_point(data = subset(plot_df, cs != 0), aes(x = BP, y = neglogp_gwas, color = cs), shape = 21, fill = "gray87", size = 1, stroke = 5) + 
  jpeg(file=plot_loc)
  print(temp_plot)
  dev.off()
} else if (no_cred_sets == FALSE) {
  #Get the correct set of residuals:
  if (plot_type == 'absolute_residual') {
    print("NOTE: Absolute residual plot type selected.")
    type_abbr <- 'absolute'
  } else if (plot_type == 'preceding_residual'){
    print("NOTE: Preceding residual plot type selected.")
    type_abbr <- 'precede'
  } else{
    print("ERROR: Invalid plot type selected. Please select non-residual, preceding_residual, or absolute_residual")
    quit(save="no");
  }
  #Loop over credible sets and make residual plots
  for (i in cred_sets) {
    #Set plot name
    plot_loc <- paste0(out_dir, file_prefix, ".cs-", i, ".", plot_type, ".jpg")
    #Make dataframe for plotting
    plot_df <- cbind.data.frame(SNP=names(fitted_rss1$pip), CHR=rep(fitted_rss1$chr, length(fitted_rss1$bp)), BP=fitted_rss1$bp,  
                                neglogp_residual=fitted_rss1[[paste0("neglogp_", type_abbr)]][as.character(i)])
    colnames(plot_df)[ncol(plot_df)] = "neglogp_residual"
    plot_df <- cbind.data.frame(plot_df, cs=factor(ifelse(plot_df$SNP %in% names(fitted_rss1$sets$cs[[paste0("L",i)]]), as.numeric(i), 0), levels = c(0, i)))
    #Make plot
    temp_plot <- ggplot(plot_df) + xlab(paste0("Chromosome ", as.character(median(plot_df$CHR)))) + 
      ylab("-log(P-Value)") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                   panel.background = element_blank(), axis.line = element_line(colour = "black"),
                                   legend.position = "none") + 
      scale_color_manual(values=coloor[i]) + 
      geom_point(aes(x = BP, y = neglogp_residual), color = "gray87") +
      geom_point(data = subset(plot_df, cs != 0), aes(x = BP, y = neglogp_residual, color = cs), shape = 21, fill = "gray87", size = 1, stroke = 5)
    jpeg(file=plot_loc)
    print(temp_plot)
    dev.off()
  }
} else {
  #This is the case of a residual plot where there are no credible sets
  print(paste0("ERROR: attempted to plot residuals for a region with no identified credible sets, ", file_prefix))
}





