args = commandArgs(trailingOnly=TRUE)

path="./output/HARE_ANC_Run"
pheno_file = args[1]
pheno_list <- read.table(pheno_file, header=F)
rds_files <- c()
groups <- c("EUR", "AFR", "HIS", "ASN")

# loop phenotypes list 
for (pheno in pheno_list$V1) {
  for (group in groups){
    tmp <- paste(path, pheno, group, "step2", paste("step2", pheno, group, "rds", sep="."), sep="/")
    if (file.exists(tmp)){
      rds_files = c(rds_files, tmp)
    }
  }
}

complete_rds <- NULL
out_rds <- args[2]

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
	x$group <- data.frame(do.call('rbind', strsplit(as.character(basename(x$output_col)),'.',fixed=TRUE)))$X6
	x$include_file <- paste("./data/genetic_annotations", paste("R4.include_Imp3MAC20.HARE", x$group, "txt.gz", sep="."), sep="/")
	#complete_rds <- rbind(complete_rds, x[x$'step2_complete' == FALSE,])
	complete_rds <- rbind(complete_rds, x)
      }
    }
}
complete_rds <- as.data.frame(subset(complete_rds, select = -c(step2_complete,group)))
complete_rds$outputlog_col <-  trimws( complete_rds$output_col)
saveRDS(complete_rds, file=out_rds)
