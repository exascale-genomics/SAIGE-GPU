mylibpath <- "/ccs/home/arodriguez/med112/task0101113/tools/saige_20220326/SAIGE-DOE/extdata/R_lib"
tasktools_lib <- "/ccs/home/arodriguez/med112/task0101113/tools/R-2/R-4.0.3/library/"
.libPaths( c( .libPaths(), mylibpath) )
suppressMessages(library(pbdMPI, lib.loc=tasktools_lib))
suppressMessages(library(tasktools, lib.loc=tasktools_lib))
suppressMessages(library(SAIGE, lib.loc=mylibpath))
path = function(...) paste0(list(...), collapse="/")
root = "/gpfs/alpine/proj-shared/med112/task0101113/output/HARE_ANC_Run"
#params = path(root, "output/phe454_1.submit_df.rda")
args = commandArgs(trailingOnly=TRUE)
params = args[1]
phecode = args[2]
group = args[3]
#run_name = paste(phecode, group, sep=".")
#nvlm_path = paste("/mnt/bb/arodriguez", run_name, sep="/")
#if (comm.localrank() == 0) dir.create(nvlm_path)

checkpoint_path = path(root, "checkpoints", "step1", phecode, group)
#checkpoint_path = NULL
if (!dir.exists(checkpoint_path)) dir.create(checkpoint_path, recursive=TRUE)
# -------------------------------------------------------------------------
p = readRDS(params)
#print(p)
wrapper = function(i)
{
  plinkFile = p[i, "plinkFile"]
  phenoFile = p[i, "phenoFile"]
  outputPrefix = p[i, "output_col"]
  logfile = p[i, "outputlog_col"]
  traitType = p[i, "traitType"]

  # get run_name
  phenoFile_base = basename(phenoFile)
  run_name = paste(unlist(strsplit(phenoFile_base, ".", fixed=T))[2], unlist(strsplit(phenoFile_base, ".", fixed=T))[3], sep=".")
  nvlm_path = paste("/mnt/bb/arodriguez", run_name, sep="/")
  dir.create(nvlm_path)
  print(phenoFile_base)

  # copy genotype files
  plinkFile_base = basename(plinkFile)
  node_plink_path = paste(nvlm_path, plinkFile_base, sep="/") 
  file.copy(paste(plinkFile, "fam", sep="."), paste(node_plink_path, "fam", sep="."))
  file.copy(paste(plinkFile, "bim", sep="."), paste(node_plink_path, "bim", sep="."))
  file.copy(paste(plinkFile, "bed", sep="."), paste(node_plink_path, "bed", sep="."))

  # copy phenotype file
  node_phenoFile_path = paste(nvlm_path, phenoFile_base, sep="/")
  file.copy(phenoFile, node_phenoFile_path)

  # set names for outputs
  outputPrefix_base = paste(basename(outputPrefix), sep=".")
  node_outputPrefix = paste(nvlm_path, outputPrefix_base, sep="/")

  outputLog_base = paste(basename(logfile), sep=".")
  node_outputLog = paste(nvlm_path, outputLog_base, sep="/")

  covars <- strsplit("Age,Gender,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10", ",")[[1]]
  qcovarList <- c('Gender')

  convertoNumeric = function(x,stringOutput){
          y= tryCatch(expr = as.numeric(x),warning = function(w) {return(NULL)})
          if(is.null(y)){
                  stop(stringOutput, " is not numeric\n")
          }else{
                  cat(stringOutput, " is ", y, "\n")
          }
          return(y)
  }

  tauInit <- convertoNumeric(strsplit("0,0", ",")[[1]], "tauInit")
  cateVarRatioMinMACVecExclude <- convertoNumeric(x=strsplit("0.5,1.5,2.5,3.5,4.5,5.5,10.5,20.5",",")[[1]], "cateVarRatioMinMACVecExclude")
  cateVarRatioMaxMACVecInclude <- convertoNumeric(x=strsplit("1.5,2.5,3.5,4.5,5.5,10.5,20.5",",")[[1]], "cateVarRatioMaxMACVecInclude")
  print(paste("Phenotype for this run:", node_phenoFile_path, sep=" "))
  print(paste("START", outputPrefix_base, Sys.time(), sep=" "))

  #capture.output(file=logfile, append=FALSE,
    fitNULLGLMM(
      plinkFile=node_plink_path,
      phenoFile=node_phenoFile_path,
      phenoCol = "Phenotype",
      outputPrefix=node_outputPrefix,
      traitType = traitType,
      covarColList = covars,
      qCovarCol = qcovarList,
      sampleIDColinphenoFile = "MVPCore_ID",
      maxiter = 25,
      maxiterPCG = 500,
      LOCO = FALSE,
      minMAFforGRM = 0.01,
      nThreads=1
    )
  #)

  #capture.output(file=logfile, append=FALSE,
  #)

  # copy back to filesystem the 3 output files
  #print("PATHS:")
  #print(paste(node_outputPrefix, "rda", sep="."))
  #print(paste(outputPrefix, "rda", sep="."))
  if (file.exists(paste(node_outputPrefix, "rda", sep="."))) {
      file.copy(paste(node_outputPrefix, "rda", sep="."), paste(outputPrefix, "rda", sep="."))
      file.copy(paste(node_outputPrefix, "_30markers.SAIGE.results.txt", sep=""), paste(outputPrefix, "_30markers.SAIGE.results.txt", sep=""))
      file.copy(paste(node_outputPrefix, "varianceRatio.txt", sep="."), paste(outputPrefix, "varianceRatio.txt", sep="."))
      file.copy(node_outputLog, paste(logfile, "log", sep="."))
  }
  unlink(nvlm_path, recursive = TRUE)
  rm(list = ls())
  gc()
  print(paste("END", outputPrefix_base, Sys.time(), sep=" "))
}
#test_case_index = 1 # just an example, replace with the appropriate number
#wrapper(test_case_index)
n = NROW(p)
#ret <- mpi_napply(n, checkpoint_path=checkpoint_path,
#  function(i) tryCatch(error=identity, wrapper(i))
#)
ret <- mpi_napply(n,
  function(i) tryCatch(error=identity, wrapper(i))
)
comm.print(unlist(ret))

##err = sapply(ret, function(x) inherits(x, "simpleError"))
###print(ret)
##if (any(err)){
##  cat(paste("An error occurred in jobs: ", paste(which(err), ret, collapse=",")))
##  cat("\n")
##}
finalize()
#memuse::Sys.procmem()
