args = commandArgs(trailingOnly=TRUE)

path=".//output/HARE_ANC_Run"
group = args[1]
#group="EUR"
pattern = paste("step1", '*', group, "rds", sep=".")
print(length(args))
if (length(args) == 1){
    print("IAMHERE")
    rds_files <- paste(path, list.files(path, pattern = pattern, recursive = TRUE), sep="/")
    out_rds <- paste(path, "SUBMIT", paste("complete_step1.tasks.omp", group, "rds", sep="."), sep="/")
} else {
    print("ISHOULDNOT BEHERE")
    phecode_list <- read.table(args[2])
    rds_files <- paste(path, phecode_list$V1, group, "step1", paste("step1", phecode_list$V1, group, "rds", sep="."), sep="/")
    out_rds <- paste(path, "SUBMIT", paste("complete_step1.tasks.omp.list", length(rds_files),  group, "rds", sep="."), sep="/")
}
complete_rds <- NULL
for (rds in rds_files)
{
    if (!grepl("SUBMIT", rds, fixed=TRUE)){
        x <- readRDS(rds)
        #if (dim(x)[1] > 1){
        #	    print(rds)
        #}
        rda_file <- paste(x[1,'output_col'], 'omp.rda', sep=".")
        if (!file.exists(rda_file)){
            complete_rds <- rbind(complete_rds, x)
        }
    }
}
complete_rds <- as.data.frame(complete_rds)
complete_rds$plinkFile <- "./output/pheCodes/inputs/genotypes/re-run/20200917.GenotypeData.Release4.mvpcoreid.ld.maf.05"
complete_rds$output_col <- paste(complete_rds$output_col, "omp", sep=".")
complete_rds$outputlog_col <-  trimws( complete_rds$output_col)
saveRDS(complete_rds, file=out_rds)
