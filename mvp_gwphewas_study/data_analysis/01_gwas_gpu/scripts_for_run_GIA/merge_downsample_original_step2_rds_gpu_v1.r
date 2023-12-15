args = commandArgs(trailingOnly=TRUE)

path=".//output/GIA_ANC_Run"
group = args[1]
pattern = paste("step2.downsampled", '*', group, "original.rds", sep=".")
#rds_files <- list.files(path, pattern = pattern, recursive = TRUE)
rds_files <- system(sprintf('find "%s" -name "%s"', path, pattern), intern=TRUE)

complete_rds <- NULL
out_rds <- paste(".//output/GIA_ANC_Run/SUBMIT/complete_step2.downsample_original.gpu", group, "rds", sep=".")
for (rds in rds_files)
{
    var_files <- list.files(path=paste(dirname(rds), "../step1", sep="/"), pattern = "\\.gpu.varianceRatio.txt$")
    if (length(var_files) > 0 )
    {
      if (file.info(paste(dirname(rds), "../step1", var_files[1], sep="/"))$size > 0)
      {
    	x <- as.data.frame(readRDS(paste(rds, sep="/")))
        x$output_col <- gsub("assoc.txt", "assoc.gpu.txt", x$output_col)
	x$outputlog_col <- x$output_col
	x$outputlog_col <- gsub("assoc.gpu.txt", "assoc.gpu.txt.log", x$outputlog_col)
	x$step2_complete <- file.exists(x[,'outputlog_col'])
	x$var_col <- gsub("out.varianceRatio.txt", "out.gpu.varianceRatio.txt", x$var_col)
	x$GMMAT_col <- gsub("out.rda", "out.gpu.rda", x$GMMAT_col)
	complete_rds <- rbind(complete_rds, x[x$'step2_complete' == FALSE,])
      }
    }
}
complete_rds <- as.data.frame(subset(complete_rds, select = -step2_complete))
complete_rds$outputlog_col <-  trimws( complete_rds$output_col)
saveRDS(complete_rds, file=out_rds)
