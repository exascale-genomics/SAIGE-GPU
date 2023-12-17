################################################################################

#run_deming_high_pip.R

#The purpose of this script is to run the Deming regression analysis of variant 
#effect sizes for the high-PIP variants found in multiple ancestries.

#This code runs locally.

################################################################################

# 0) Call libraries and set variables and directories ====

# Load necessary libraries
library(ggplot2)
library(plotly)
library(deming)
library(dplyr)

#Set Output (and input) directory
out_dir="C:/Users/mitch/Documents/UPenn/Voight_Lab/MVP_Finemapping/Figure_Building_and_Synthesis/"
downsample_dir="C:/Users/mitch/Documents/UPenn/Voight_Lab/MVP_Finemapping/downsampling/"

# 1) Create a function to create the high-PIP Deming regressions for a given population combination ====


make_deming <- function(ANC1, ANC2, out_dir, file_prefix=NULL, type="mu"){
  #Set input file location
  if (is.null(file_prefix)) {
    ANC1_ANC2_raw <- read.csv(file = paste0(out_dir, ANC1, ".", ANC2, ".demming_input_table.tsv"), header=TRUE, sep = "\t")
  } else {
    ANC1_ANC2_raw <- read.csv(file = paste0(out_dir, ANC1, ".", ANC2, ".", file_prefix, ".demming_input_table.tsv"), header=TRUE, sep = "\t")
  }
  
  # Subset the data frame using the input columns
  data <- ANC1_ANC2_raw %>% mutate(mu.se.ANC1=sqrt(mu2.ANC1 - (mu.ANC1)^2), mu.se.ANC2=sqrt(mu2.ANC2 - (mu.ANC2)^2)) 
  if (type == "mu") {
    data <- data %>% select(mu.ANC1, mu.ANC2, mu.se.ANC1, mu.se.ANC2, trait)
  } else if(type == "beta") {
    data <- data %>% select(beta.ANC.ANC1, beta.ANC.ANC2, se.ANC.ANC1, se.ANC.ANC2, trait)
  } else {
    print("ERROR: Invalid effect estimate selected")
    return(NULL)
  }
  
  # Rename the columns to match the formula
  colnames(data) <- c("x", "y", "xstd", "ystd","trait")
  
  # Call the deming function with the updated formula and data
  deming_result <- deming(formula = y ~ x, data = data, xstd = xstd, ystd = ystd)
  least_squares_result <- lm(y ~x, data = data)
  
  # Extract regression parameters
  intercept <- deming_result$coefficients[1]
  slope <- deming_result$coefficients[2]
  lower_intercept <- deming_result$ci[1]
  lower_slope <- deming_result$ci[2]
  upper_intercept <- deming_result$ci[3]
  upper_slope <- deming_result$ci[4]
  confidence_int = paste0("(", round(lower_slope, 2), ", ", round(upper_slope, 2), ")")
  
  # Calculate confidence intervals for the range of x values
  x_range <- seq(min(data$x), max(data$x), length.out = nrow(data))
  y_fit <- intercept + slope * x_range
  y_lower <- lower_intercept + lower_slope * x_range
  y_upper <- upper_intercept + upper_slope * x_range
  
  # Calculate midpoint for label placement
  label_x <- mean(data$x) 
  label_y <- intercept + slope * label_x
  
  # Generate the ggplot
  p <- ggplot(data, aes(x=x, y=y)) +
    geom_point(aes(color=trait)) +
    geom_path(aes(x=x_range, y=y_fit), color="blue", linetype="dotted") +
    geom_ribbon(aes(x=x_range, ymin=y_lower, ymax=y_upper), fill="blue", alpha=0.3) +
    geom_text(aes(x=label_x - 0.1*(max(x)-min(x)), y=label_y, 
                  label=paste0("Slope: ", round(slope, 3))), 
              vjust=-1.5, hjust=0.5, size=4, color="black") +
    annotate("text", x = -Inf, y = Inf, label = confidence_int, hjust = 0, vjust = 1, size = 10, color="black") +
    xlab(paste0("mu ", ANC1)) +
    ylab(paste0("mu ", ANC2)) +
    theme_bw() +
    theme(legend.position = "none")
  
  #Output the plot to a jpeg
  if (is.null(file_prefix)) {
    jpeg(paste0(out_dir, ANC1, ".", ANC2, ".deming_regression.", type, ".scatter.jpeg"), width = 720, height=720)
  } else {
    jpeg(paste0(out_dir, ANC1, ".", ANC2, ".deming_regression.", file_prefix, ".", type, ".scatter.jpeg"), width = 720, height=720)
  }
    print(p)
  dev.off()
}

# 2) Call Deming Function ====

#Call function
make_deming("AFR", "AMR", out_dir)
make_deming("AFR", "EUR", out_dir)
make_deming("AMR", "EUR", out_dir)

# 3) Read in Downsampled Results and Re-Run for AFR/EUR ====

#Read in downsampled files
signal_downsample_raw <- read.csv(paste0(downsample_dir, "master.signals.downsampled.txt"), sep = "\t", header = TRUE)
variant_downsample_raw <- read.csv(paste0(downsample_dir, "master.variants.downsampled.txt"), sep = "\t", header = TRUE)

#Filter for case where best variants are the same from each population
signal_downsample_filt <- signal_downsample_raw %>% filter(EUR.best_variant == AFR.best_variant) %>% 
  filter(EUR.max_overall_pip >= 0.95 & AFR.max_overall_pip >= 0.95)

#Make columns to inner join on
signal_downsample_filt <- signal_downsample_filt %>% mutate(trait_variant=paste0(trait, ".", EUR.best_variant))
variant_downsample_raw <- variant_downsample_raw %>% mutate(trait_variant=paste0(phenotype, ".", mvp_id))

#Join the dataframes
signal_downsample_append <- signal_downsample_filt %>% left_join((variant_downsample_raw %>% filter(ancestry == "AFR")), by="trait_variant") %>%
  left_join((variant_downsample_raw %>% filter(ancestry == "EUR")), by="trait_variant")
colnames(signal_downsample_append) <- str_replace_all(colnames(signal_downsample_append), pattern = "[.]x", replacement = ".ANC1")
colnames(signal_downsample_append) <- str_replace_all(colnames(signal_downsample_append), pattern = "[.]y", replacement = ".ANC2")

#Write to file
write.table(signal_downsample_append, paste0(downsample_dir, "AFR", ".", "EUR", ".", "downsampled", ".demming_input_table.tsv"), sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

#Call Deming function
make_deming("AFR", "EUR", downsample_dir, file_prefix="downsampled")

# 4) Read in Full Results to match variants and traits with downsampled analysis ====

#Read in full files
signal_full_raw <- read.csv(paste0(out_dir, "master.signals.txt"), sep = "\t", header = TRUE)
variant_full_raw <- read.csv(paste0(out_dir, "master.variants.txt"), sep = "\t", header = TRUE)

#Filter for case where best variants are the same from each population
signal_full_filt <- signal_full_raw %>% filter(EUR.best_variant == AFR.best_variant) %>%
  mutate(trait_variant=paste0(trait, ".", EUR.best_variant)) %>% 
  filter(trait_variant %in% signal_downsample_append$trait_variant)

#Make column to inner join on
variant_full_raw <- variant_full_raw %>% mutate(trait_variant=paste0(phenotype, ".", mvp_id))

#Join the dataframes
signal_full_append <- signal_full_filt %>% left_join((variant_full_raw %>% filter(ancestry == "AFR")), by="trait_variant") %>%
  left_join((variant_full_raw %>% filter(ancestry == "EUR")), by="trait_variant")
colnames(signal_full_append) <- str_replace_all(colnames(signal_full_append), pattern = "[.]x", replacement = ".ANC1")
colnames(signal_full_append) <- str_replace_all(colnames(signal_full_append), pattern = "[.]y", replacement = ".ANC2")

#Write to file
write.table(signal_full_append, paste0(downsample_dir, "AFR", ".", "EUR", ".", "matched_full", ".demming_input_table.tsv"), sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

#Call Deming function
make_deming("AFR", "EUR", downsample_dir, file_prefix="matched_full")

