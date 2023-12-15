args = commandArgs(trailingOnly=TRUE)

path="./output/HARE_ANC_Run"
group = args[1]
pattern = paste("step2", '*', group, "rds$", sep=".")
rds_files <- list.files(path, pattern = pattern, recursive = TRUE)
complete_rds <- NULL
out_rds <- paste("./output/HARE_ANC_Run/SUBMIT/complete_step2.finishedstep1tasks", group, "rds", sep=".")
for (rds in rds_files)
{
    if (length(list.files(path=paste(path, dirname(rds), "../step1", sep="/"), pattern = "\\.omp.varianceRatio.txt$")) > 0)
    {
    	x <- as.data.frame(readRDS(paste(path, rds, sep="/")))
	x$step2_complete <- file.exists(x[,'output_col'])
	x$var_col <- gsub("out.varianceRatio.txt", "out.omp.varianceRatio.txt", x$var_col)
	x$GMMAT_col <- gsub("out.rda", "out.omp.rda", x$GMMAT_col)
	complete_rds <- rbind(complete_rds, x[x$'step2_complete' == FALSE,])
    }
}
complete_rds <- as.data.frame(subset(complete_rds, select = -step2_complete))
complete_rds$outputlog_col <-  trimws( complete_rds$output_col)
saveRDS(complete_rds, file=out_rds)
