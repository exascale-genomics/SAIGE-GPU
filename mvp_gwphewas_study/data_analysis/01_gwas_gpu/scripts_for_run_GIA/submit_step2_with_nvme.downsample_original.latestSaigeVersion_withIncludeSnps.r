mylibpath <- ".//tools/saige_20220326/SAIGE-DOE/extdata/R_lib"
tasktools_lib <- ".//tools/R-2/R-4.0.3/library/"
#.libPaths( c( .libPaths(), mylibpath) )
#suppressMessages(library(SAIGE, lib.loc=".//tools/saige_20220326/SAIGE-DOE/extdata/R_lib"))
suppressMessages(library(SAIGE, lib.loc=".//tools/R-2/R-4.0.3/library/"))
suppressMessages(library(pbdMPI, lib.loc=tasktools_lib))
suppressMessages(library(tasktools, lib.loc=tasktools_lib))
library(parallel)

path = function(...) paste0(list(...), collapse="/")
root = "./output/GIA_ANC_Run"
#params = path(root, "output/phe454_1.submit_df.rda")
args = commandArgs(trailingOnly=TRUE)
params = args[1]
phecode = args[2]
group = args[3]
core_count = args[4]
subsample_file = "./results/DOWNSAMPLING/Samples.EUR.AFR_downsample.txt"

#run_name = paste(phecode, group, sep=".")
#nvlm_path = paste("/mnt/bb/", run_name, sep="/")
#dir.create(nvlm_path)

checkpoint_path = path(root, "checkpoints", phecode, group)
##checkpoint_path = NULL
#if (!dir.exists(checkpoint_path)) dir.create(checkpoint_path, recursive=TRUE)
# -------------------------------------------------------------------------
p = readRDS(params)
#p$bin <- (seq(nrow(p))-1)%%as.numeric(core_count)
#p <- p[order(as.numeric(as.character(p$bin))), ]
wrapper = function(i)
{
  bgenFile = p[i, "bgen_col"]
  bgenFileIndex = p[i, "bgen_index_col"]
  sampleFile = p[i, "samp_col"]
  GMMATmodelFile = p[i, "GMMAT_col"]
  varianceRatioFile = p[i, "var_col"]
  SAIGEOutputFile = p[i, "output_col"]
  minMAF = p[i, "minMAF_col"]
  minMAC = p[i, "minMAC_col"]
  numLinesOutput = p[i, "numLinesOutput_col"]
  LOCO = p[i, "LOCO_col"]
  logfile = p[i, "outputlog_col"]
  idstoIncludeFile = p[i, "include_col"]

  # get run_name
  varianceRatioFile_base = basename(varianceRatioFile)
  run_name = paste(unlist(strsplit(varianceRatioFile_base, ".", fixed=T))[4], unlist(strsplit(varianceRatioFile_base, ".", fixed=T))[5], basename(bgenFile), sep=".")
  nvlm_path = paste("/mnt/bb/", run_name, sep="/")
  dir.create(nvlm_path)

  ##nvlm_bgenFile = paste(nvlm_path, basename(bgenFile), sep="/")
  ##file.copy(bgenFile, nvlm_bgenFile)
  ##nvlm_bgenFileIndex = paste(nvlm_path, basename(bgenFileIndex), sep="/")
  ##file.copy(bgenFileIndex, nvlm_bgenFileIndex)
  ##nvlm_sampleFile = paste(nvlm_path, basename(sampleFile), sep="/")
  ##file.copy(sampleFile, nvlm_sampleFile)
  ##nvlm_GMMATmodelFile = paste(nvlm_path, basename(GMMATmodelFile), sep="/")
  ##file.copy(GMMATmodelFile, nvlm_GMMATmodelFile)
  ##nvlm_varianceRatioFile = paste(nvlm_path, basename(varianceRatioFile), sep="/")
  ##file.copy(varianceRatioFile, nvlm_varianceRatioFile)

  outFile_base = basename(SAIGEOutputFile)
  chromosome = gsub("chr", "", unlist(strsplit(outFile_base, ".", fixed=TRUE))[1])
  node_outfile_path = paste(nvlm_path, outFile_base, sep="/")
  nvlm_logfile <- paste(nvlm_path, paste(outFile_base, "log", sep="."), sep="/")
  start_time <- Sys.time()
  print(paste("START", node_outfile_path, start_time, sep=" "))
  print(paste("COPY:", node_outfile_path, SAIGEOutputFile, sep=" "))
  #capture.output(file=logfile, append=FALSE,
    ##SPAGMMATtest_new(
    #SPAGMMATtest( 
    #  bgenFile=bgenFile,
    #  bgenFileIndex=bgenFileIndex,
    #  sampleFile=sampleFile,
    #  GMMATmodelFile=GMMATmodelFile,
    #  varianceRatioFile=varianceRatioFile,
    #  SAIGEOutputFile=node_outfile_path,
    #  min_MAF=0.0001,
    #  min_MAC=20,
    #  markers_per_chunk=10000,
    #  LOCO=FALSE
    #)
  #)

  print (paste("Rscript ", ".//tools/saige_step1GPU_step2wei/SAIGE/extdata/step2_SPAtests.R", " --bgenFile ",  bgenFile, " --bgenFileIndex ", bgenFileIndex, " --minMAF=0.0001 --minMAC=20 --markers_per_chunk=10000 --LOCO=TRUE --chrom", chromosome, " --GMMATmodelFile ", GMMATmodelFile, " --varianceRatioFile ", varianceRatioFile, " --sampleFile ", sampleFile, " --SAIGEOutputFile ", node_outfile_path, " --idstoIncludeFile ", idstoIncludeFile, " --is_Firth_beta=TRUE --pCutoffforFirth=0.05 --subSampleFile ",  subsample_file, " > ", nvlm_logfile, sep=" "))
  system(paste("Rscript ", ".//tools/saige_step1GPU_step2wei/SAIGE/extdata/step2_SPAtests.R", " --bgenFile ",  bgenFile, " --bgenFileIndex ", bgenFileIndex, " --minMAF=0.0001 --minMAC=20 --markers_per_chunk=10000 --LOCO=TRUE --chrom ", chromosome, " --GMMATmodelFile ", GMMATmodelFile, " --varianceRatioFile ", varianceRatioFile, " --sampleFile ", sampleFile, " --SAIGEOutputFile ", node_outfile_path, " --idstoIncludeFile ", idstoIncludeFile, " --is_Firth_beta=TRUE --pCutoffforFirth=0.05 --subSampleFile ", subsample_file, " > ", nvlm_logfile, sep=" "))
#  print(paste("Rscript ", "./latest_2022030318/SAIGE/extdata/step2_SPAtests.R --bgenFile ",  bgenFile, " --bgenFileIndex ", bgenFileIndex, " --minMAF=0.0001 --minMAC=20 --markers_per_chunk=10000 --LOCO=FALSE --GMMATmodelFile ", GMMATmodelFile, " --varianceRatioFile ", varianceRatioFile, " --sampleFile ", sampleFile, " --SAIGEOutputFile ", node_outfile_path, " > ", nvlm_logfile, sep=" "))
     #idstoIncludeFile=nvlm_include_path,

  # copy back to filesystem the 3 output files
  file.copy(node_outfile_path, SAIGEOutputFile, overwrite=TRUE)
  file.copy(nvlm_logfile, paste(logfile, "log", sep="."), overwrite=TRUE)
  #Sys.sleep(60)
  unlink(nvlm_path, recursive = TRUE)  
  end_time <- Sys.time()
  print(paste("END", node_outfile_path, start_time, end_time, sep=" "))
}

n = NROW(p)
#ret <- mpi_napply(n, checkpoint_path=checkpoint_path,
#  function(i) tryCatch(error=identity, wrapper(i))
#)
#ret <- mpi_napply(n, preschedule=FALSE, 
#  function(i) tryCatch(error=identity, wrapper(i))
#)
ret <- mpi_napply(n,
  function(i) tryCatch(error=identity, wrapper(i))
)


err = sapply(ret, function(x) inherits(x, "simpleError"))
if (any(err)){
  cat(paste("An error occurred in jobs: ", paste(which(err), collapse=",")))
  cat("\n")
}
finalize()
#memuse::Sys.procmem()
