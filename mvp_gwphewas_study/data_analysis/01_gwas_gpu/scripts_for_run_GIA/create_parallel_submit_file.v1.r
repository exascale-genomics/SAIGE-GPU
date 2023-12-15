#Rscript /ccs/home/arodriguez/med112/task0101113/batch/pre-gwas/gwPheWAS-Summit/scripts_for_run_GIA/create_parallel_submit_file.v1.r EAS > ./submit.step1.parallel.ASN.gpu.txt

args = commandArgs(trailingOnly=TRUE)

path="/ccs/home/arodriguez/med112/task0101113/output/GIA_ANC_Run"
saige_path="/gpfs/alpine/proj-shared/med112/task0101113/tools/saige_20220326/SAIGE-DOE/extdata"
group = args[1]
gpu_qty = 1
if (group == "EAS") {
	gpu_qty = 1
} else if (group == "AMR") {
        gpu_qty = 1 # 2 for regular batch
} else if (group == "AFR") {
        gpu_qty = 2 # 4 for regular batch
} else if (group == "EUR") {
        gpu_qty = 8 # 16 for regular batch
}
#group="EUR"
pattern = paste("step1", '*', group, "rds", sep=".")
#cat(length(args))
if (length(args) == 1){
    #cat("IAMHERE")
    #rds_files <- paste(path, list.files(path, pattern = pattern, recursive = TRUE), sep="/")
    rds_files <- system(sprintf('find "%s" -name "%s"', path, pattern), intern=TRUE)
    out_rds <- paste(path, "SUBMIT", paste("complete_step1.tasks.gpu", group, "rds", sep="."), sep="/")
} else {
    #cat("ISHOULDNOT BEHERE")
    phecode_list <- read.table(args[2])
    rds_files <- paste(path, phecode_list$V1, group, "step1", paste("step1", phecode_list$V1, group, "rds", sep="."), sep="/")
    out_rds <- paste(path, "SUBMIT", paste("complete_step1.tasks.gpu.list", length(rds_files),  group, "rds", sep="."), sep="/")
}
complete_rds <- NULL
for (rds in rds_files)
{
    if (!grepl("SUBMIT", rds, fixed=TRUE)){
        x <- readRDS(rds)
        #if (dim(x)[1] > 1){
        #	    cat(rds)
        #}
        rda_file <- paste(x[1,'output_col'], 'gpu.rda', sep=".")
	variance_file <- paste(x[1,'output_col'], 'varianceRatio.txt', sep=".")
        if (!file.exists(rda_file) ){
            complete_rds <- rbind(complete_rds, x)
        }
    }
}
complete_rds <- as.data.frame(complete_rds)
complete_rds$plinkFile <- "/gpfs/alpine/proj-shared/med112/task0101113/output/pheCodes/inputs/genotypes/re-run/20200917.GenotypeData.Release4.mvpcoreid.ld.maf.05"
complete_rds$output_col <- paste(complete_rds$output_col, "gpu", sep=".")
complete_rds$outputlog_col <-  paste(trimws( complete_rds$output_col), "log", sep=".")

total_nodes = paste("#BSUB -nnodes", dim(complete_rds)[1], sep=" ")
wall_time = paste("#BSUB -W", "24:00", sep=" ")
queue_name = paste("#BSUB -q", "batch-hm", sep=" ")
proj_name = paste("#BSUB -P", "MED112", sep=" ")
stdout_name = paste("#BSUB -o", paste("all.s1", group, "gpu.stdout", sep="."), sep=" ")
stderr_name = paste("#BSUB -e", paste("all.s1", group, "gpu.stderr", sep="."), sep=" ")
run_name = paste("#BSUB -J", paste("all.s1", group, sep="."), sep=" ")
bb = paste("#BSUB -alloc_flags", "nvme", sep=" ")
thread_count = "1"
step1_script = paste(saige_path, "step1_fitNULLGLMM.R", sep="/")

module_list = c("module load python/2.7.15-anaconda2-5.3.0",
		"module load r/4.0.5",
		"module load gcc/9.1.0",
		"module load cmake",
		"module load openblas",
		"module load cuda/11.0.2",
		"module load spectrum-mpi/10.4.0.3-20210112")


complete_rds$command <- paste("jsrun -n", gpu_qty, "-a1 -c1 -g1", 
			      paste("Rscript", step1_script, sep=" "),
			      paste("--plinkFile=", complete_rds$plinkFile, sep=""),
			      paste("--phenoFile=", complete_rds$phenoFile, sep=""),
			      paste("--outputPrefix=", complete_rds$output_col, sep=""),
			      paste("--phenoCol=", "Phenotype", sep=""),
			      paste("--covarColList=", "Age,Gender,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10", sep=""),
			      paste("--qCovarColList=Gender", sep=""),
			      paste("--sampleIDColinphenoFile=", "MVPCore_ID", sep=""),
			      paste("--traitType=", complete_rds$traitType, sep=""),
                              paste("--nThreads=", thread_count, sep=""),
                              paste("--minMAFforGRM", "0.05", sep=" "),
                              paste("--maxiterPCG", "500", sep=" "),
                              paste("--maxiter", "25", sep=" "),
                              paste("--LOCO", "TRUE", sep=" "),
                              paste(complete_rds$gender_col, sep=" "),
			      #paste("--is_Firth_beta", "TRUE", sep="="),
			      #paste("--pCutoffforFirth", "0.05", sep="="),
                              paste(">", complete_rds$outputlog_col, "&", sep=" "),
			      sep = " ")

### write the submit file
cat("#!/bin/bash\n")
cat(total_nodes, "\n")
cat(wall_time, "\n")
cat(queue_name, "\n")
cat(proj_name, "\n")
cat(stdout_name, "\n")
cat(stderr_name, "\n")
cat(run_name, "\n")
cat(bb, "\n\n")
writeLines(paste(module_list, sep="\n"))
writeLines(complete_rds$command);cat("\nwait\n")

