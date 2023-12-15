.libPaths( c( .libPaths(), "/ccs/home/arodriguez/med112/task0101113/YoungDae_work/R_libs/") )
suppressMessages(library(pbdMPI, lib.loc="/ccs/home/arodriguez/med112/task0101113/tools/R-2/R-4.0.3/library/"))
suppressMessages(library(tasktools, lib.loc="/ccs/home/arodriguez/med112/task0101113/tools/R-2/R-4.0.3/library/"))
suppressMessages(library(SAIGE, lib.loc="/ccs/home/arodriguez/med112/task0101113/YoungDae_work/R_libs/"))
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
  outputPrefix_base = basename(outputPrefix)
  node_outputPrefix = paste(nvlm_path, outputPrefix_base, sep="/")

  outputLog_base = basename(logfile)
  node_outputLog = paste(nvlm_path, outputLog_base, sep="/")

  covars <- strsplit("Age,Gender,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10", ",")[[1]]

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

  #capture.output(file=logfile, append=FALSE,
    fitNULLGLMM(
      plinkFile=node_plink_path,
      phenoFile=node_phenoFile_path,
      phenoCol = "Phenotype",
      outputPrefix=node_outputPrefix,
      traitType = traitType,
      invNormalize = FALSE,
      covarColList = covars,
      qCovarCol = NULL,
      sampleIDColinphenoFile = "MVPCore_ID",
      tol = 0.02,
      maxiter = 25,
      tolPCG = 1e-5,
      maxiterPCG = 500,
      nThreads = 42,
      SPAcutoff = 2,
      numMarkers = 30,
      skipModelFitting = FALSE,
      memoryChunk = 2,
      tauInit = tauInit,
      LOCO = FALSE,
      traceCVcutoff = 0.0025,
      ratioCVcutoff = 0.001,
      outputPrefix_varRatio = NULL,
      IsOverwriteVarianceRatioFile = FALSE,
      IsSparseKin = FALSE,
      sparseGRMFile = NULL,
      sparseGRMSampleIDFile=NULL,
      numRandomMarkerforSparseKin = 2000,
      relatednessCutoff = 0.125,
      isCateVarianceRatio = FALSE,
      cateVarRatioMinMACVecExclude = cateVarRatioMinMACVecExclude,
      cateVarRatioMaxMACVecInclude = cateVarRatioMaxMACVecInclude,
      isCovariateTransform = TRUE,
      isDiagofKinSetAsOne = FALSE,
      useSparseSigmaConditionerforPCG = FALSE,
      useSparseSigmaforInitTau = FALSE,
      minMAFforGRM = 0.0001,
      minCovariateCount=-1,
      includeNonautoMarkersforVarRatio=FALSE,
      sexCol="",
      FemaleCode=1,
      FemaleOnly=FALSE,
      MaleCode=0,
      MaleOnly=FALSE,
      noEstFixedEff=FALSE
    )
  #)

  #capture.output(file=logfile, append=FALSE,
##  system(paste("export OMP_NUM_THREADS=42; Rscript /gpfs/alpine/proj-shared/med112/task0101113/YoungDae_work/src/SAIGE/extdata/step1_fitNULLGLMM.R --plinkFile=", node_plink_path, " --phenoFile=", node_phenoFile_path, " --outputPrefix=", node_outputPrefix, " --phenoCol=Phenotype --covarColList=Age,Gender,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 --sampleIDColinphenoFile=MVPCore_ID --traitType=binary --nThreads=42 --minMAFforGRM 0.0001 --maxiterPCG 500 --maxiter 25 --LOCO FALSE > ", node_outputLog, sep=""))
  #system(paste("jsrun -n1 -c42 --smpiargs='-disable_gpu_hooks' -bpacked:42 -EOMP_NUM_THREADS=42 Rscript /gpfs/alpine/proj-shared/med112/task0101113/YoungDae_work/src/SAIGE/extdata/step1_fitNULLGLMM.R --plinkFile=", node_plink_path, " --phenoFile=", node_phenoFile_path, " --outputPrefix=", node_outputPrefix, " --phenoCol=Phenotype --covarColList=Age,Gender,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 --sampleIDColinphenoFile=MVPCore_ID --traitType=binary --nThreads=42 --minMAFforGRM 0.0001 --maxiterPCG 500 --maxiter 25 --LOCO FALSE > ", node_outputLog, sep=""))    
  #)

  # copy back to filesystem the 3 output files
  #print("PATHS:")
  #print(paste(node_outputPrefix, "rda", sep="."))
  #print(paste(outputPrefix, "rda", sep="."))
  if (file.exists(paste(node_outputPrefix, "rda", sep="."))) {
      file.copy(paste(node_outputPrefix, "rda", sep="."), paste(outputPrefix, "rda", sep="."))
      file.copy(paste(node_outputPrefix, "_30markers.SAIGE.results.txt", sep=""), paste(outputPrefix, "_30markers.SAIGE.results.txt", sep=""))
      file.copy(paste(node_outputPrefix, "varianceRatio.txt", sep="."), paste(outputPrefix, "varianceRatio.txt", sep="."))
      file.copy(node_outputLog, logfile)
  }
  unlink(nvlm_path, recursive = TRUE)
  rm(list = ls())
  gc()
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
memuse::Sys.procmem()
