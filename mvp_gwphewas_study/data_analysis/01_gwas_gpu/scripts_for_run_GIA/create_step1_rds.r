#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

# script will create the df for a particular phenotype and group

colClasses = c("character", "character", "character", "character", "character", "character")
col.names = c("--plinkFile=", "--phenoFile=", "--outputPrefix=", "--outputLog=", "--traitType=", "--genderCode")
df <- read.table(text = "", colClasses = colClasses, col.names = col.names)

phenoID = args[1]
group = args[2]
outrds = args[3]
traitType = args[4]
#step1_output = args[4]
gender_code = args[5]
base_path = "/gpfs/alpine/med112/proj-shared/task0101113/output/GIA_ANC_Run"

output_path= paste(base_path, phenoID, group, "step1", sep="/")
plinkFile="/gpfs/alpine/proj-shared/med112/task0101113/output/pheCodes/inputs/genotypes/re-run/20200917.GenotypeData.Release4.mvpcoreid.ld.maf"
phenoFile=paste(base_path, phenoID, group, "inputs", paste("PhenoFile", phenoID, group, "txt", sep="."), sep="/")
output_col = paste(output_path, paste("phewas.ld.maf", phenoID, group, "out", sep="."), sep="/")
outputlog_col = paste(output_path, paste("phewas.ld.maf", phenoID, group, "log", sep="."), sep="/")
gender_col = ""

if (gender_code == "F"){
    gender_col = " --FemaleOnly TRUE --FemaleCode 2 --sexCol Gender "
} else if (gender_code == "M") {
    gender_col = " --MaleOnly TRUE --MaleCode 1 --sexCol Gender "
}

df <- rbind(df, as.matrix(data.frame(plinkFile, phenoFile, output_col, outputlog_col, traitType, gender_col)))
df <- as.matrix(df)
saveRDS(object=df, file=outrds)

