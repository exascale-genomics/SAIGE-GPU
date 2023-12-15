args = commandArgs(trailingOnly=TRUE)

path="/ccs/home/arodriguez/med112/task0101113/output/HARE_ANC_Run"
group = args[1]
pattern = paste("step2", '*', group, "rds", sep=".")
#rds_files <- list.files(path, pattern = pattern, recursive = TRUE)
rds_files <- system(sprintf('find "%s" -name "%s"', path, pattern), intern=TRUE)

total_nodes = paste("#BSUB -nnodes", length(rds_files)[1], sep=" ")
wall_time = paste("#BSUB -W", "2:00", sep=" ")
queue_name = paste("#BSUB -q", "batch", sep=" ")
proj_name = paste("#BSUB -P", "MED112", sep=" ")
stdout_name = paste("#BSUB -o", paste("all.s2", group, "gpu.stdout", sep="."), sep=" ")
stderr_name = paste("#BSUB -e", paste("all.s2", group, "gpu.stderr", sep="."), sep=" ")
run_name = paste("#BSUB -J", paste("all.s2", group, sep="."), sep=" ")
bb = paste("#BSUB -alloc_flags", "nvme", sep=" ")
thread_count = "1"

module_list = c("module load python/2.7.15-anaconda2-5.3.0",
                "module load r/4.0.5",
                "module load gcc/9.1.0",
                "module load cmake",
                "module load openblas",
                "module load cuda/11.0.2",
                "module load spectrum-mpi/10.4.0.3-20210112")

complete_rds <- NULL
out_rds <- paste("/ccs/home/arodriguez/med112/task0101113/output/HARE_ANC_Run/SUBMIT/complete_step2.finishedstep1tasks.gpu", group, "rds", sep=".")
count = 0
submit_count = 1
for (rds in rds_files)
{
    if (count == 0){
        # write the header
	submit_file <- paste("/ccs/home/arodriguez/med112/task0101113/output/HARE_ANC_Run/SUBMIT", paste("submit.step2.all", group, "nvlm.includesnps.parallel.gpu", submit_count, "txt", sep="."), sep="/")
        cat("#!/bin/bash\n", file=submit_file, append=FALSE)
        cat(paste(total_nodes, "\n", sep=""), file=submit_file, append=TRUE)
        cat(paste(wall_time, "\n", sep=""), file=submit_file, append=TRUE)
        cat(paste(queue_name, "\n", sep=""), file=submit_file, append=TRUE)
        cat(paste(proj_name, "\n", sep=""), file=submit_file, append=TRUE)
        cat(paste(paste(stdout_name, submit_count, sep="."), "\n", sep=""), file=submit_file, append=TRUE)
        cat(paste(paste(stderr_name, submit_count, sep="."), "\n", sep=""), file=submit_file, append=TRUE)
        cat(paste(paste(run_name, submit_count, sep="."), "\n", sep=""), file=submit_file, append=TRUE)
        cat(paste(bb, "\n\n", sep=""), file=submit_file, append=TRUE)
        cat(paste(module_list), sep="\n", file=submit_file, append=TRUE)
	#cat("World",file="outfile.txt",append=TRUE)
    }
    var_files <- system(sprintf('find "%s" -name "%s"', paste(dirname(rds), "../step1", sep="/"), "*.gpu.varianceRatio.txt"), intern=TRUE)
    if (length(var_files) > 0 )
    {
      if (file.info(var_files[1])$size > 0)
      {
    	x <- as.data.frame(readRDS(paste(rds, sep="/")))
        x$output_col <- gsub("assoc.txt", "assoc.gpu.txt", x$output_col)
	x$outputlog_col <- x$output_col
	x$outputlog_col <- gsub("assoc.gpu.txt", "assoc.gpu.txt.log", x$outputlog_col)
	x$step2_complete <- file.exists(x[,'outputlog_col'])
	x$var_col <- gsub("out.varianceRatio.txt", "out.gpu.varianceRatio.txt", x$var_col)
	x$GMMAT_col <- gsub("out.rda", "out.gpu.rda", x$GMMAT_col)

	# get run_name
        x$varianceRatioFile_base = basename(x$var_col)
        x$run_name = paste(unlist(strsplit(x$varianceRatioFile_base, ".", fixed=T))[4], unlist(strsplit(x$varianceRatioFile_base, ".", fixed=T))[5], basename(x$bgen_col), sep=".")
        x$nvlm_path = paste("/mnt/bb/arodriguez", x$run_name, sep="/")
        ####dir.create(nvlm_path)
	x$cmd <- paste("jsrun -n1 -c1 csh -c \'mkdir", "-p", x$nvlm_path)

	# get the include list file
        x$local_group = unlist(strsplit(x$varianceRatioFile_base, ".", fixed=T))[5]
        x$include_list_file = NULL
        x$include_path = "/gpfs/alpine/med112/proj-shared/data/genetic_annotations/"
        #include_path = "/ccs/home/arodriguez/med112/task0101113/batch/pre-gwas/scripts_for_run/"
	x$include_list_file = paste(x$include_path, paste("R4.include_Imp3MAC20.HARE", x$local_group, "txt.gz", sep="."), sep="/")
        x$nvlm_include_path = paste(x$nvlm_path, paste("R4.include_Imp3MAC20.HARE", x$local_group, "txt.gz", sep="."), sep="/")

	####file.copy(include_list_file, nvlm_include_path)
	x$cmd <- paste(x$cmd, paste("cp", x$include_list_file, x$nvlm_include_path), sep="; ")
        x$outFile_base = basename(x$output_col)
	out <- strsplit(as.character(x$outFile_base), ".", fixed=TRUE)
	tmp <- do.call(rbind, out)[,1]
	tmp <- gsub("chr", "", tmp)
	x$chromosome <- tmp
        x$node_outfile_path = paste(x$nvlm_path, x$outFile_base, sep="/")
        x$nvlm_logfile <- paste(x$nvlm_path, paste(x$outFile_base, "log", sep="."), sep="/")
        ####print(paste("START", node_outfile_path, Sys.time(), sep=" "))
	x$cmd <- paste(x$cmd, paste("echo", "START", x$node_outfile_path, "`date`"), sep="; ")
        ####system(paste("Rscript ", "/ccs/home/arodriguez/med112/task0101113/tools/saige_20220326/SAIGE-DOE/extdata/step2_SPAtests.R", " --bgenFile ",  bgenFile, " --bgenFileIndex ", bgenFileIndex, " --minMAF=0.0001 --minMAC=20 --markers_per_chunk=10000 --LOCO=TRUE --chrom=", chromosome, " --GMMATmodelFile ", GMMATmodelFile, " --varianceRatioFile ", varianceRatioFile, " --sampleFile ", sampleFile, " --SAIGEOutputFile ", node_outfile_path, " --idstoIncludeFile ", nvlm_include_path, " --is_Firth_beta=TRUE --pCutoffforFirth=0.05 > ", nvlm_logfile, sep=" "))
	x$cmd <- paste(x$cmd, paste("Rscript", "/ccs/home/arodriguez/med112/task0101113/tools/saige_20220326/SAIGE-DOE/extdata/step2_SPAtests.R", " --bgenFile ",  x$bgen_col, " --bgenFileIndex ", x$bgen_index_col, " --minMAF=0.01 --minMAC=20 --markers_per_chunk=10000 --LOCO=TRUE", paste("--chrom=", x$chromosome, sep=""), " --GMMATmodelFile ", x$GMMAT_col, " --varianceRatioFile ", x$var_col, " --sampleFile ", x$samp_col, " --SAIGEOutputFile ", x$node_outfile_path, " --idstoIncludeFile ", x$nvlm_include_path, " --is_Firth_beta=TRUE --pCutoffforFirth=0.05 > ", x$nvlm_logfile), sep="; ")
        ####file.copy(node_outfile_path, SAIGEOutputFile, overwrite=TRUE)
	x$cmd <- paste(x$cmd, paste("cp", x$node_outfile_path, x$output_col), sep="; ")
        ####file.copy(nvlm_logfile, paste(logfile, "log", sep="."), overwrite=TRUE)
	x$cmd <- paste(x$cmd, paste("cp", x$nvlm_logfile, x$outputlog_col), sep="; ")
        ####unlink(nvlm_path, recursive = TRUE)
	x$cmd <- paste(x$cmd, paste("rm -rf", x$nvlm_path), sep="; ")
        ####print(paste("END", node_outfile_path, Sys.time(), sep=" "))
	x$cmd <- paste(x$cmd, paste("echo", "END", x$node_outfile_path, "`date`"), sep="; ")

	x$cmd <- paste(x$cmd, "\'&", sep=" ")
        x <- x[x$'step2_complete' == FALSE,]
	cat(as.character(x$cmd), sep="\n", file=submit_file, append=TRUE)
        count <- count + dim(x)[1]
	if (count >= 100000){
	    cat("\nwait\n", file=submit_file, append=TRUE)
	    count = 0
	    submit_count = submit_count + 1
	}
      }
    }
}
cat("\nwait\n", file=submit_file, append=TRUE)
#complete_rds <- as.data.frame(subset(complete_rds, select = -step2_complete))
#complete_rds$outputlog_col <-  trimws( complete_rds$output_col)
#saveRDS(complete_rds, file=out_rds)



