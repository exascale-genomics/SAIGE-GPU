args = commandArgs(trailingOnly=TRUE)

path="/ccs/home/arodriguez/med112/task0101113/output/GIA_ANC_Run"
group = args[1]
pattern = paste("*step2.downsampled", '*', group, "rds", sep=".")
pattern_orig = paste("*step2.downsampled", '*', group, "rds.original.rds", sep=".")
#rds_files <- list.files(path, pattern = pattern, recursive = TRUE)
rds_files <- system(sprintf('find "%s" -name "%s"', path, pattern), intern=TRUE)
rds_files_orig <- system(sprintf('find "%s" -name "%s"', path, pattern_orig), intern=TRUE)

complete_rds <- NULL
complete_orig_rds <- NULL
out_rds <- paste("/ccs/home/arodriguez/med112/task0101113/output/GIA_ANC_Run/SUBMIT/complete_step2.downsample.gpu", group, "rds", sep=".")
out_rds_orig <- paste("/ccs/home/arodriguez/med112/task0101113/output/GIA_ANC_Run/SUBMIT/complete_step2.downsample.original.gpu", group, "rds", sep=".")

for (rds in rds_files_orig)
{
    var_files <- list.files(path=paste(dirname(rds), "../../step1", sep="/"), pattern = "^phewas.ld.maf.*out.gpu.varianceRatio.txt$")
    if (length(var_files) > 0 )
    {
      if (file.info(paste(dirname(rds), "../../step1", var_files[1], sep="/"))$size > 0)
      {
        x <- as.data.frame(readRDS(paste(rds, sep="/")))
        x$orig_output_col <- gsub("assoc.txt", "assoc.gpu.txt", x$orig_output_col)
        x$orig_outputlog_col <- x$orig_output_col
        x$orig_outputlog_col <- gsub("assoc.gpu.txt", "assoc.gpu.txt.log", x$orig_outputlog_col)
        x$step2_complete <- file.exists(x[,'orig_outputlog_col'])
        x$orig_var_col <- gsub("out.varianceRatio.txt", "out.gpu.varianceRatio.txt", x$orig_var_col)
        x$orig_GMMAT_col <- gsub("out.rda", "out.gpu.rda", x$orig_GMMAT_col)
        complete_orig_rds <- rbind(complete_orig_rds, x[x$'step2_complete' == FALSE,])
      }
    }
}
colnames(complete_orig_rds) <- c("bgen_col", "bgen_index_col", "samp_col", "GMMAT_col", "var_col", "output_col", "minMAF_col", "minMAC_col", "numLinesOutput_col", "IsOutputNinCaseCtrl_col", "IsOutputAFinCaseCtrl_col", "LOCO_col", "outputlog_col", "include_col", "step2_complete" )
complete_orig_rds <- as.data.frame(subset(complete_orig_rds, select = -step2_complete))
complete_orig_rds$outputlog_col <-  trimws( complete_orig_rds$output_col)
saveRDS(complete_orig_rds, file=out_rds_orig)

for (rds in rds_files)
{
    var_files <- list.files(path=paste(dirname(rds), "../../step1", sep="/"), pattern = "^downsample.*.gpu.varianceRatio.txt$")
    if (length(var_files) > 0 )
    {
      if (file.info(paste(dirname(rds), "../../step1", var_files[1], sep="/"))$size > 0)
      {
    	x <- as.data.frame(readRDS(paste(rds, sep="/")))
        x$output_col <- gsub("assoc.txt", "assoc.gpu.txt", x$output_col)
	x$outputlog_col <- x$output_col
	x$outputlog_col <- gsub("assoc.gpu.txt", "assoc.gpu.txt.log", x$outputlog_col)
	x$step2_complete <- file.exists(x[,'outputlog_col'])
	x$var_col <- gsub("out.varianceRatio.txt", "out.gpu.varianceRatio.txt", x$var_col)
	x$GMMAT_col <- gsub("out.rda", "out.gpu.rda", x$GMMAT_col)
	#complete_rds <- rbind(complete_rds, x[x$'step2_complete' == FALSE,])
	complete_rds <- rbind(complete_rds, x)
      }
    }
}
complete_rds <- as.data.frame(subset(complete_rds, select = -step2_complete))
complete_rds$outputlog_col <-  trimws( complete_rds$output_col)
saveRDS(complete_rds, file=out_rds)
