args = commandArgs(trailingOnly=TRUE)

path="/ccs/home/arodriguez/med112/task0101113/output/HARE_ANC_Run"
group = args[1]
#group="EUR"
pattern = paste("step1", '*', group, "rds", sep=".")
rds_files <- list.files(path, pattern = pattern, recursive = TRUE)
complete_rds <- NULL
out_rds <- paste("/ccs/home/arodriguez/med112/task0101113/output/HARE_ANC_Run/SUBMIT/complete_step1.tasks", group, "rds", sep=".")
for (rds in rds_files)
{
    x <- readRDS(paste(path, rds, sep="/"))
    #if (dim(x)[1] > 1){
    #	    print(rds)
    #}
    rda_file <- paste(x[1,'output_col'], 'rda', sep=".")
    if (!file.exists(rda_file)){
        complete_rds <- rbind(complete_rds, x)
    }
}
complete_rds <- as.data.frame(complete_rds)
complete_rds$outputlog_col <-  trimws( complete_rds$output_col)
saveRDS(complete_rds, file=out_rds)
