################################################################################

#save_empty_rds.R

#This script is designed to take in a file location for an rds file and save a
#blank rds file in that location

################################################################################

# 0) Load Needed Libraries ====

# 1) Read in Needed Arguments from Command ====

#Need to read in the following: 
  # 1) Location to save blank rds file REQUIRED

args <- commandArgs(trailingOnly = TRUE);
if (length(args) == 1) {
  #Print confirmation output
  print("Necessary Parameters Input")
  #Set variables from the arguments in this case
  rds_file_loc <- args[1]
} else {
  print("Usage: %> Rscript save_empty_rds.R Blank_RDS_Location");
  quit(save="no");
}

# 2) Save empty file ====

if(file.exists(dirname(rds_file_loc))){
  temp = strsplit(basename(rds_file_loc), "[.]")
  #Check that extension is rds
  if (temp[[1]][length(temp[[1]])] != "rds") {
    print("WARNING: Saving rds object to location with non-rds file extension: ", rds_file_loc)
  }
  #Save empty file
  empty_object = NULL
  saveRDS(empty_object, file = rds_file_loc)
} else{
  print(paste0("ERROR: Attempted to save blank rds file to invalid location: ", rds_file_loc));
  quit(save="no");
}

