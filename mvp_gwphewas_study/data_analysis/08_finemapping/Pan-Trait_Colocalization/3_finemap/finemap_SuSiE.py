'''
finemap_SuSiE.py

This script that will fine map a given locus using SuSiE. It will call an 
identically named Rscript and will return as it's output a .rds file that 
contains the output of fine mapping plus the negative log P-values, the base 
pair positions of the SNPs, and the rsids. The output will be in a standard 
format so that it can be used for plotting results, testing the effects of 
different purity/association thresholds, and any other downstream analyses 
that arise.

***This needs to be run on Python/3.9 and R/4.0 as there an incompatability ***
***with other versions of the coloc package.                                ***
'''

###############################################################################

#Load libraries
import getopt
import numpy as np
import pandas as pd
import sys
import os
from os.path import exists
from scipy.stats import norm

#Define help function to be called if there is a problem
def help(exit_num=1):
    print("""-----------------------------------------------------------------
ARGUMENTS
    -p => <txt> file prefix for files should consist of trait.chr.start.stop.ancestry REQUIRED
    -o => <txt> output folder location REQUIRED
    -m => <txt> munged summary statistics directory (default is /gpfs/alpine/med112/proj-shared/results/FOR_SUSIE) OPTIONAL
    -l => <txt> Storage location for map and matrix files (default is /gpfs/alpine/med112/proj-shared/results/PostGWAS_Analysis/finemaping/GIA/matrix_dir) OPTIONAL
    -c => <num> desired confidence level for credible sets (default is 0.95) OPTIONAL
    -s => <num> random seed (default is 5) OPTIONAL
    -f => <num> minimum minor allele frequency (default is 0 / no-minimum) OPTIONAL
    -x => <txt> prefix for map/matrix files if different from file prefix / Use chr.start.stop.ancestry format (Default is file prefix) OPTIONAL
    -a => <str> flip orientation of variants with allele frequencies > 0.5 (default is False) OPTIONAL
    -n => <str> summary stats file naming convention showing where ANCestry and TRAIT are found relative to periods 
                (default is TRAIT.ANC.GIA.KDI.txt.gz) OPTIONAL
    -i => <str> file of case/control counts (default is /gpfs/alpine/med112/proj-shared/results/MVP_R4.1000G_AGR.GIA.DataDictionary.txt) OPTIONAL
    -r => <txt> location of finemap_SuSiE.R (default is finemap_SuSiE.R in this script's directory) OPTIONAL
""")
    sys.exit(exit_num)

###############################################################################
###########################  COLLECT AND CHECK ARGUMENTS  #####################
###############################################################################

## GLOBAL VARS

## MAIN ##

def main(argv): 
    try: 
        opts, args = getopt.getopt(sys.argv[1:], "p:o:m:l:c:s:f:x:a:n:i:r:nh")
    except getopt.GetoptError:
        print("ERROR: Incorrect usage of getopts flags!")
        help()

    options_dict = dict(opts)
    
    ## Check for help flag and return empty string to break function
    help_ind = options_dict.get('-h', False)
    if help_ind != False:
        help()
    
    ## Required arguments
    try:
        file_prefix = options_dict['-p']
        out_dir = options_dict['-o']

    except KeyError:
        print("ERROR: One of your required arguments does not exist.")
        help()
    
    munge_dir = options_dict.get('-m', "/gpfs/alpine/med112/proj-shared/results/FOR_SUSIE")
    ld_loc = options_dict.get('-l', "/gpfs/alpine/med112/proj-shared/results/PostGWAS_Analysis/finemaping/GIA/matrix_dir")
    confidence = options_dict.get('-c', 0.95)
    random_seed = options_dict.get('-s', 5)
    min_maf = options_dict.get('-f', 0)
    ld_prefix = options_dict.get('-x', file_prefix)
    swap_af_ind = options_dict.get('-a', 'False')
    file_format = options_dict.get('-n', 'TRAIT.ANC.GIA.KDI.txt.gz')
    count_loc = options_dict.get('-i', '/gpfs/alpine/med112/proj-shared/results/MVP_R4.1000G_AGR.GIA.DataDictionary.txt')
    temp = os.path.realpath(__file__)
    rscript_loc = options_dict.get('-r', temp[:temp.find('finemap_SuSiE.py')] + 'finemap_SuSiE.R')
    #Confirm that the file/folder locs exist
    if exists(munge_dir) == False:
        print("ERROR: Directory for munged summary statistics does not exist.")
        sys.exit(1)
    if exists(ld_loc) == False:
        print("ERROR: Storage location for LD matrices and map files does not exist.")
        sys.exit(1)
    if exists(out_dir) == False:
        print("ERROR: Output directory does not exist.")
        sys.exit(1)
    if exists(count_loc) == False:
        print("ERROR: File of case/control counts does not exist at specified location.")
        sys.exit(1)
    if exists(rscript_loc) == False:
        print("ERROR: finemap_SuSiE.R does not exist at specified location.")
        sys.exit(1)
    #Recast confidence value 
    try:
        confidence = float(confidence)
        if confidence >= 1 or confidence <= 0: 
            print("ERROR: Confidence level should be a proportion in (0,1).")
            sys.exit(1)
    except ValueError:
        print("ERROR: Confidence level is not coercible to a float. It should be a proportion in (0,1).")
        sys.exit(1)
    #Recast random_seed 
    try:
        random_seed = int(random_seed)
    except ValueError:
        print("ERROR: Random seed must be an integer.")
        sys.exit(1)
    #Recast minimum allele frequency 
    try:
        min_maf = float(min_maf)
        if min_maf >= 1 or min_maf < 0: 
            print("ERROR: Minimum allele frequency for variants should be a proportion in [0,1).")
            sys.exit(1)
    except ValueError:
        print("ERROR: Minimum allele frequency is not coercible to a float. It should be a proportion in [0,1).")
        sys.exit(1)
    #Verify that ANC and TRAIT are present in file format
    if 'ANC' not in file_format or 'TRAIT' not in file_format:
        print("ERROR: Unable to identify position of ANCestry or TRAIT in file name format.")
        sys.exit(1)
    #Recast swap_af_ind
    if(type(swap_af_ind)) == str:
        if(swap_af_ind.lower() == 'false' or swap_af_ind.lower() == 'f' or swap_af_ind.lower() == 'no' or swap_af_ind.lower() == 'n'):
            swap_af_ind = False
        elif(swap_af_ind.lower() == 'true' or swap_af_ind.lower() == 't' or swap_af_ind.lower() == 'yes' or swap_af_ind.lower() == 'y'):
            swap_af_ind = True
        else:
            print("ERROR: Unrecognized option selected for whether to swap alleles for variants with frequencies > 0.5. Select True/False.")
            sys.exit(1)
    elif(type(swap_af_ind) == bool):
        pass #In this case the indicator is already a boolean, so there's no need to recast.
    else:
        print("ERROR: Unrecognized option selected for whether to swap alleles for variants with frequencies > 0.5. Select True/False.")
        sys.exit(1)
    print("Acceptable Inputs Given for " + file_prefix)
    
    
    #Call driver function
    try:
        driver(file_prefix, out_dir, munge_dir, ld_loc, confidence, random_seed, min_maf, ld_prefix, swap_af_ind, file_format, count_loc, rscript_loc)
    except MemoryError:
        ### Output a blank rds file ###
        out_dir = add_slash(out_dir) #Add slash to out directory if it needs it
        temp = os.path.realpath(__file__)
        empty_script_loc = options_dict.get('-r', temp[:temp.find('finemap_SuSiE.py')] + 'save_empty_rds.R') #Get location of Rscript that saves empty rds files
        filename = out_dir + file_prefix + '.rds' #Set the empty file name
        os.system('Rscript ' + empty_script_loc + ' ' + filename) #Call the rscript

###############################################################################
#############################  HELPFUL FUNCTIONS ##############################
###############################################################################

#Function that adds slashes to directory names
def add_slash(directory):
    if directory[-1] != '/':
        return directory + '/'
    else:
        return directory

###############################################################################
#############################  DRIVER  ########################################
###############################################################################

def driver(file_prefix, out_dir, munge_dir, ld_loc, confidence, random_seed, min_maf, ld_prefix, swap_af_ind, file_format, count_loc, rscript_loc): 
    
    #Add final slash to out_dir and ld_loc if it's missing
    out_dir = add_slash(out_dir)
    ld_loc = add_slash(ld_loc)
    munge_dir = add_slash(munge_dir)
    
    #Extract info from the file_prefix
    temp = file_prefix.split(".")
    trait = temp[0]
    chromo = int(temp[1][3:])
    bp_start =  int(temp[2])
    bp_end = int(temp[3])
    ancestry = temp[4]
    
    #Identify ancestry and trait position from file_format
    name_conv = file_format.split(".")
    anc_pos = name_conv.index("ANC")
    trait_pos = name_conv.index("TRAIT")
    remain_name_indices=list(range(len(name_conv)))
    remain_name_indices.remove(anc_pos)
    remain_name_indices.remove(trait_pos)
    
    #Get map and matrix file names from other variables
    int_map_name = ld_loc + ld_prefix + ".map"
    int_matrix_name = ld_loc + ld_prefix + ".ld"
    #Double check that ld files exist
    if not exists(int_map_name) or not exists(int_matrix_name):
        print("ERROR: Specified LD files not found.")
        sys.exit(1)
    
    #Get munged summary stats file loc
    temp = os.listdir(munge_dir)
    temp = [test.split(".") for test in temp]
    temp = [test for test in temp if len(test) == len(name_conv)]
    temp = [test for test in temp if test[anc_pos] == ancestry and test[trait_pos] == trait]
    for i in remain_name_indices:
        temp = [test for test in temp if test[i] == name_conv[i]]
    #Check if we have an error
    if len(temp) > 1: 
        print("ERROR: Multiple files with correct ancestry, trait, and file name format found in directory.")
        sys.exit(1)
    else:
        munge_loc = munge_dir + '.'.join(temp[0])
       
    #Read munged summary stats into a pandas df 
    munge_raw = pd.read_table(munge_loc, sep="\t")
    #Filter summary stats for those in the locus
    temp = munge_raw['ID'].astype(str).str.split(":", expand=True)
    munge_raw.rename(columns = {'pos':'pos38'}, inplace = True)
    munge_raw['pos37'] = temp[1]
    munge_raw['pos37'] = munge_raw['pos37'].astype(int)
    munge_filt = munge_raw[(munge_raw['chrom'] == chromo) & (munge_raw['pos37'] >= bp_start) & (munge_raw['pos37'] <= bp_end) & (munge_raw['maf'] >= min_maf)]
    #Flip the data for any rows where the allele frequency is > 0.5 if the setting is on
    if swap_af_ind == True:
        munge_temp = munge_filt.copy() #Create a version of munge_filt that we can modify so we can save it again later
        #Handle fields in both binary and quantitative traits
        munge_temp.loc[munge_temp["af"] > 0.5, 'ac'] = 2 * munge_temp.loc[munge_temp["af"] > 0.5, 'num_samples'] - munge_temp.loc[munge_temp["af"] > 0.5, 'ac']
        munge_temp.loc[munge_temp["af"] > 0.5, 'beta'] = -1 * munge_temp.loc[munge_temp["af"] > 0.5, 'beta']
        temp = munge_temp.loc[munge_temp["af"] > 0.5, 'ref']
        munge_temp.loc[munge_temp["af"] > 0.5, 'ref'] = munge_temp.loc[munge_temp["af"] > 0.5, 'alt']
        munge_temp.loc[munge_temp["af"] > 0.5, 'alt'] = temp
        #Check for binary trait fields for the rest
        if 'case_af' in munge_temp.columns and 'control_af' in munge_temp.columns and 'or' in munge_temp.columns:
            munge_temp.loc[munge_temp["af"] > 0.5, 'case_af'] = 1 - munge_temp.loc[munge_temp["af"] > 0.5, 'case_af']
            munge_temp.loc[munge_temp["af"] > 0.5, 'control_af'] = 1 - munge_temp.loc[munge_temp["af"] > 0.5, 'control_af']
            munge_temp.loc[munge_temp["af"] > 0.5, 'or'] = 1 / munge_temp.loc[munge_temp["af"] > 0.5, 'or']            
        #Handle af last or the rest won't work
        munge_temp.loc[munge_temp["af"] > 0.5, 'af'] = 1 - munge_temp.loc[munge_temp["af"] > 0.5, 'af'] 
        #Set munge_filt equal to the new munge_temp
        munge_filt = munge_temp.copy()
        del munge_temp
    #Save a copy of the munge_filt to later track which variants were removed due to missing data
    munge_copy = munge_filt.copy()
    del munge_raw
    
    #Read in the ld variant index information
    map_raw = pd.read_table(int_map_name, sep="\t", header=None)
    
    #Make a dictionary that relates the direction of a variants alleles to a +/-1 that can be used to check the munged allele direction
    temp = map_raw[1].tolist()
    #temp2 = map_raw[1].astype(str).str.split(":", expand=True) #Commenting out allele reversal logic
    #temp2 = temp2[0] + ":" + temp2[1] + ":" + temp2[3] + ":" + temp2[2] #Commenting out allele reversal logic
    #temp2 = temp2.tolist() #Commenting out allele reversal logic
    temp3 = [1 for x in range(len(temp))]
    #temp3.extend([-1 for x in range(len(temp2))]) #Commenting out allele reversal logic
    #temp.extend(temp2) #Commenting out allele reversal logic
    map_var_dict = {temp[i] : temp3[i] for i in range(len(temp3))}
    #Append a column onto the munge dataframe reflecting whether the variant is in the map file
    temp = munge_filt['ID']
    temp = temp.tolist()
    temp2 = [map_var_dict[i] if i in map_var_dict.keys() else 0 for i in temp]
    munge_append = munge_filt.copy()
    munge_append['in_map'] = temp2
    #Filter for the variants in the map file
    munge_filt = munge_append[munge_append['in_map'] != 0]
    
    #Next make a dictionary from the variants in the munge_filt df
    temp = munge_filt['ID'].astype(str)
    temp = temp.tolist()
    #temp3 = munge_filt['ID'].astype(str).str.split(":", expand=True) #Commenting out allele reversal logic
    #temp3 = temp3[0] + ":" + temp3[1] + ":" + temp3[3] + ":" + temp3[2] #Commenting out allele reversal logic
    #temp.extend(temp3.tolist()) #Commenting out allele reversal logic
    temp2 = munge_filt['in_map'].tolist()
    #temp2.extend(temp2) #We want to add back the same orientation here for the reversed alleles since we already reveresed orientation when making the other dict #Commenting out allele reversal logic
    munge_var_dict = {temp[i] : temp2[i] for i in range(len(temp))}
    #Append a column back onto the map file that gives allele alignment/whether the variant is present
    temp = map_raw[1]
    temp = temp.tolist()
    temp2 = [munge_var_dict[i] if i in munge_var_dict.keys() else 0 for i in temp]
    map_raw['in_munge'] = temp2
    #Clean up
    del munge_append
    
    #Read in the ld_table
    ld_raw = np.loadtxt(int_matrix_name, dtype=float)
    
    #Filter down the ld matrix and map file
    good_map_idx = map_raw[map_raw['in_munge'] != 0].index
    ld_filt = ld_raw[list(good_map_idx)]
    del ld_raw #Reduce pull on memory
    ld_filt = ld_filt[:,list(good_map_idx)]
    map_filt = map_raw.iloc[good_map_idx, ]
    #Reset the pandas indices
    map_filt.reset_index(inplace = True, drop = True)
    munge_filt.reset_index(inplace = True, drop = True)
    
    #Flip the signs in the LD matrix where needed
    flip_sign_idx = map_filt[map_filt['in_munge'] == -1].index
    ld_filt[flip_sign_idx] = -1 * ld_filt[flip_sign_idx] #Flips the rows first
    ld_filt[:,list(flip_sign_idx)] = -1 * ld_filt[:,list(flip_sign_idx)] #Then flip the column
    
    #Overwrite map file unique ids with IDs from the munge file
    map_filt = map_filt.copy()
    map_filt.iloc[:,1] = munge_filt['ID']
    
    ##### Filter out completely empty columns/rows from LD matrix #####
    #Convert numpy array to pandas dataframe
    ld_filt_df = pd.DataFrame(ld_filt, columns = map_filt.iloc[:,1])
    ld_filt_df.index = ld_filt_df.columns
    #Drop completely blank columns and rows
    ld_filt_df.dropna(how = 'all', axis = 1, inplace=True)
    ld_filt_df.dropna(how = 'all', axis = 0, inplace=True)
    #Filter down the map file and summary stats too
    map_filt.index = map_filt.iloc[:,1]
    map_filt = map_filt.filter(items = ld_filt_df.columns, axis=0)
    munge_filt.index = munge_filt['ID']
    munge_filt = munge_filt.filter(items = ld_filt_df.columns, axis=0)
    #Send back to numpy array
    ld_filt = ld_filt_df.to_numpy(dtype=float)
    #Return to pandas df
    ld_filt_df = pd.DataFrame(ld_filt, columns = map_filt.iloc[:,1])
    ld_filt_df.index = ld_filt_df.columns
    
    ##### Filter out partially empty columns/rows from LD matrix in order #####
    #Get positions of all nans
    temp = np.argwhere(np.isnan(ld_filt))
    temp = np.asarray(np.unique(temp, return_counts=True))
    temp = pd.DataFrame(temp.transpose(), columns = ['Index', 'Count'])
    temp['Index'] = ld_filt_df.columns[temp['Index']]
    temp.sort_values(by=['Count'], inplace=True, ascending=False) #At this point we have a sorted list of the most problematic row/column indices to loop over
    #Define a copy of ld_filt_temp outside of loop to avoid later deletion error
    ld_filt_temp = ld_filt_df.to_numpy(dtype=float)
    #Sequentially remove the most problematic variant until matrix is complete
    while len(temp) != 0:
        var_index=temp.iloc[0,0]
        #Drop entries off tof ld_filt_df
        ld_filt_df.drop(index=var_index, inplace=True)
        ld_filt_df.drop(columns=var_index, inplace=True)
        #Send back to a new numpy array
        ld_filt_temp = ld_filt_df.to_numpy(dtype=float)
        #Check if all missing values are gone
        if np.sum(np.isnan(ld_filt_temp)) == 0:
            break #Don't need to remove any more variants
        else: #Reset the prioritized list of problematic variants
            temp = np.argwhere(np.isnan(ld_filt_temp))
            temp = np.asarray(np.unique(temp, return_counts=True))
            temp = pd.DataFrame(temp.transpose(), columns = ['Index', 'Count'])
            temp['Index'] = ld_filt_df.columns[temp['Index']]
            temp.sort_values(by=['Count'], inplace=True, ascending=False)
    #Down subset the summary stats and map file
    map_filt = map_filt.filter(items = ld_filt_df.columns, axis=0)
    munge_filt = munge_filt.filter(items = ld_filt_df.columns, axis=0)
    #Identify the SNPs that are being tossed by the filtering
    munge_copy.index = munge_copy['ID']
    temp = [x for x in munge_copy['ID'].tolist() if x not in munge_filt['ID'].tolist()]
    munge_copy = munge_copy.loc[temp,['ID', 'rsid', 'pos38', 'pval']]
    
    ##### Prep variables/files for Fine-Mapping #####
    #Check what column names we have and process z-scores accordingly
    munge_columns = munge_filt.columns.to_list()
    if 'beta' in munge_columns and 'sebeta' in munge_columns:
        #Calculate z-scores
        munge_filt['z'] = munge_filt['beta']/munge_filt['sebeta']      
    else: 
        print("ERROR: Missing betas and odds ratios/confidence intervals in summary statistics. Exiting.")
        sys.exit(1)
    #Use z-scores to calculate -logP values
    munge_filt['-logP'] = (-np.log(2) - norm.logcdf(-np.abs(munge_filt['z'])))/np.log(10)
        
    #Set file paths for temp files
    munge_temp = out_dir + file_prefix + '.temp.txt'
    map_temp = out_dir + file_prefix + '.temp.map'
    ld_temp = out_dir + file_prefix + '.temp.ld'
    #Write temp files for map_filt, ld_filt, and munge_filt to be used in Rscript
    munge_filt.to_csv(munge_temp, sep = "\t", header = True, index = False, na_rep='NA')
    map_filt.to_csv(map_temp, sep = "\t", header = False, index = False)
    ld_filt_df.to_csv(ld_temp, sep = "\t", header = False, index = False, na_rep='nan')
    #Set file path for filtered out SNPs and export
    filt_temp = out_dir + file_prefix + '.filtered_out.txt'
    munge_copy.to_csv(filt_temp, sep = "\t", header = True, index = False, na_rep='NA')
    
    #Clean up
    del map_filt, munge_filt, ld_filt, ld_filt_temp, ld_filt_df, map_raw, flip_sign_idx, good_map_idx, munge_copy, temp
        
    #Call Rscript
    os.system('Rscript ' + rscript_loc + ' ' + file_prefix + ' ' + str(random_seed) + ' ' + str(confidence) + ' ' + munge_temp + ' ' + map_temp + ' ' + ld_temp + ' ' + count_loc + ' ' + out_dir)
    
    #Delete temp files
    os.remove(map_temp)
    os.remove(ld_temp)
    os.remove(munge_temp)
    
###############################################################################


#cProfile.run('main(sys.argv[1:])')
         
if __name__ == "__main__":
    main(sys.argv[1:])
