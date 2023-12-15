mylibpath <- "/ccs/home/arodriguez/med112/task0101113/tools/R-2/R-4.0.3/library/"
tasktools_lib <- "/ccs/home/arodriguez/med112/task0101113/tools/R-2/R-4.0.3/library/"
print("hello")
.libPaths( c( .libPaths(), mylibpath) )
suppressMessages(library(pbdMPI, lib.loc=tasktools_lib))
suppressMessages(library(tasktools, lib.loc=tasktools_lib))
library(parallel)
suppressMessages(library(SAIGE, lib.loc="/ccs/home/arodriguez/med112/task0101113/tools/R-2/R-4.0.3/library/"))

path = function(...) paste0(list(...), collapse="/")
root = "/gpfs/alpine/proj-shared/med112/task0101113/output/GIA_ANC_Run"
#params = path(root, "output/phe454_1.submit_df.rda")
args = commandArgs(trailingOnly=TRUE)
params = args[1]
phecode = args[2]
group = args[3]
core_count = args[4]

#run_name = paste(phecode, group, sep=".")
#nvlm_path = paste("/mnt/bb/arodriguez", run_name, sep="/")
#dir.create(nvlm_path)

checkpoint_path = path(root, "checkpoints", phecode, group)
##checkpoint_path = NULL
#if (!dir.exists(checkpoint_path)) dir.create(checkpoint_path, recursive=TRUE)
# -------------------------------------------------------------------------
p = readRDS(params)
p$bin <- (seq(nrow(p))-1)%%as.numeric(core_count)
p <- p[order(as.numeric(as.character(p$bin))), ]
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

  # get run_name
  varianceRatioFile_base = basename(varianceRatioFile)
  run_name = paste(unlist(strsplit(varianceRatioFile_base, ".", fixed=T))[4], unlist(strsplit(varianceRatioFile_base, ".", fixed=T))[5], basename(bgenFile), sep=".")
  nvlm_path = paste("/mnt/bb/arodriguez", run_name, sep="/")
  dir.create(nvlm_path)

  # get the include list file
  local_group = unlist(strsplit(varianceRatioFile_base, ".", fixed=T))[5]
  include_list_file = NULL
  include_path = "/gpfs/alpine/med112/proj-shared/GIA/output/gia_impute-r2"
  if (local_group == "EAS") {
    include_list_file = paste(include_path, "R4.include_Imp3MAC20.GIA.EAS.txt.gz", sep="/")
    nvlm_include_path = paste(nvlm_path, "R4.include_Imp3MAC20.GIA.EAS.txt.gz", sep="/")
  } else if (local_group == "AFR") {
    include_list_file = paste(include_path, "R4.include_Imp3MAC20.GIA.AFR.txt.gz", sep="/")
    nvlm_include_path = paste(nvlm_path, "R4.include_Imp3MAC20.GIA.AFR.txt.gz", sep="/")
  } else if (local_group == "AMR") {
    include_list_file = paste(include_path, "R4.include_Imp3MAC20.GIA.AMR.txt.gz", sep="/")
    nvlm_include_path = paste(nvlm_path, "R4.include_Imp3MAC20.GIA.AMR.txt.gz", sep="/")
  } else if (local_group == "EUR") {
    include_list_file = paste(include_path, "R4.include_Imp3MAC20.GIA.EUR.txt.gz", sep="/")
    nvlm_include_path = paste(nvlm_path, "R4.include_Imp3MAC20.GIA.EUR.txt.gz", sep="/")
  }

  file.copy(include_list_file, nvlm_include_path)

  outFile_base = basename(SAIGEOutputFile)
  chromosome = gsub("chr", "", unlist(strsplit(outFile_base, ".", fixed=TRUE))[1])
  node_outfile_path = paste(nvlm_path, outFile_base, sep="/")
  nvlm_logfile <- paste(nvlm_path, paste(outFile_base, "log", sep="."), sep="/")
  start_time <- Sys.time()
  print(paste("START", node_outfile_path, start_time, sep=" "))
  print(paste("COPY:", node_outfile_path, SAIGEOutputFile, sep=" "))

  print (paste("Rscript ", "/ccs/home/arodriguez/med112/task0101113/tools/saige_yougdae_20221010/SAIGE/extdata/step2_SPAtests.R", " --bgenFile ",  bgenFile, " --bgenFileIndex ", bgenFileIndex, " --minMAF=0.0001 --minMAC=20 --markers_per_chunk=10000 --LOCO=TRUE --chrom", chromosome, " --GMMATmodelFile ", GMMATmodelFile, " --varianceRatioFile ", varianceRatioFile, " --sampleFile ", sampleFile, " --SAIGEOutputFile ", node_outfile_path, " --idstoIncludeFile ", nvlm_include_path, " --is_Firth_beta=TRUE --pCutoffforFirth=0.05 > ", nvlm_logfile, sep=" "))
  system(paste("Rscript ", "/ccs/home/arodriguez/med112/task0101113/tools/saige_yougdae_20221010/SAIGE/extdata/step2_SPAtests.R", " --bgenFile ",  bgenFile, " --bgenFileIndex ", bgenFileIndex, " --minMAF=0.0001 --minMAC=20 --markers_per_chunk=10000 --LOCO=TRUE --chrom ", chromosome, " --GMMATmodelFile ", GMMATmodelFile, " --varianceRatioFile ", varianceRatioFile, " --sampleFile ", sampleFile, " --SAIGEOutputFile ", node_outfile_path, " --idstoIncludeFile ", nvlm_include_path, " --is_Firth_beta=TRUE --pCutoffforFirth=0.05 > ", nvlm_logfile, sep=" "))

  # copy back to filesystem the 3 output files
  file.copy(node_outfile_path, paste(SAIGEOutputFile, "youngdae.txt", sep="."), overwrite=TRUE)
  file.copy(nvlm_logfile, paste(logfile, "youngdae.log", sep="."), overwrite=TRUE)
  unlink(nvlm_path, recursive = TRUE)  
  end_time <- Sys.time()
  print(paste("END", node_outfile_path, start_time, end_time, sep=" "))
}

n = NROW(p)
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
