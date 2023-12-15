#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

# script will create the df for a particular phenotype and group

colClasses = c("character", "character", "character", "character", "character")
col.names = c("--plinkFile=", "--phenoFile=", "--outputPrefix=", "--outputLog=", "--traitType=")
df <- read.table(text = "", colClasses = colClasses, col.names = col.names)

phenoID = args[1]
group = args[2]
outrds = args[3]
traitType = args[4]
#step1_output = args[4]

output_path= paste("./output/HARE_ANC_Run", phenoID, group, "step1", sep="/")
plinkFile="./output/pheCodes/inputs/genotypes/re-run/20200917.GenotypeData.Release4.mvpcoreid.ld.maf"
phenoFile=paste("./output/HARE_ANC_Run", phenoID, group, "inputs", paste("PhenoFile", phenoID, group, "txt", sep="."), sep="/")
output_col = paste(output_path, paste("phewas.ld.maf", phenoID, group, "out", sep="."), sep="/")
outputlog_col = paste(output_path, paste("phewas.ld.maf", phenoID, group, "log", sep="."), sep="/")

df <- rbind(df, as.matrix(data.frame(plinkFile, phenoFile, output_col, outputlog_col, traitType)))
df <- as.matrix(df)
saveRDS(object=df, file=outrds)

