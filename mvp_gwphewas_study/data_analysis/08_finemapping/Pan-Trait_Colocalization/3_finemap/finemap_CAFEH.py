'''
finemap_CAFEH.py

This script that will fine map a given locus using CAFEH. It will return as its
output a .pkl pickle file contains the default output of fine mapping plus the
negative log P-values, the base pair positions of the SNPs, and the rsids. The 
output will be in a standard format so that it can be used for plotting 
results, testing the effects of different purity/association thresholds, and 
any other downstream analyses that arise.

***This needs to be run on Python/3.9 as there an incompatability with other ***
***versions of the .to_list() function and the CAFEH package for python.    ***
'''

###############################################################################

#Load libraries
import getopt
import numpy as np
import pandas as pd
import sys
from os.path import exists
import os
from cafeh.cafeh_summary import fit_cafeh_summary, fit_cafeh_z
from cafeh.model_queries import *
from cafeh.cafeh_summary import CAFEHSummary
from scipy.stats import t, norm
import copy
import pickle

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
""")
    sys.exit(exit_num)

###############################################################################
###########################  COLLECT AND CHECK ARGUMENTS  #####################
###############################################################################

## GLOBAL VARS

## MAIN ##

def main(argv): 
    try: 
        opts, args = getopt.getopt(sys.argv[1:], "p:o:m:l:c:s:f:x:a:n:nh")
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
        driver(file_prefix, out_dir, munge_dir, ld_loc, confidence, random_seed, min_maf, ld_prefix, swap_af_ind, file_format)
    except MemoryError:
        #Output a blank pickle file
        out_dir = add_slash(out_dir) #Add slash to out directory if it needs it
        #Save empty pickle file
        filename = out_dir + file_prefix + '.pkl'
        outfile = open(filename,'wb')
        pickle.dump(None,outfile)
        outfile.close()
        
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

def driver(file_prefix, out_dir, munge_dir, ld_loc, confidence, random_seed, min_maf, ld_prefix, swap_af_ind, file_format): 
    
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
    temp = [test for test in temp if test[anc_pos] == ancestry and test[trait_pos] == trait and len(test) == len(name_conv)]
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
    #Also reset the ld_filt
    ld_filt = ld_filt_temp
    del ld_filt_temp, ld_filt_df
    #Identify the SNPs that are being tossed by the filtering
    munge_copy.index = munge_copy['ID']
    temp = [x for x in munge_copy['ID'].tolist() if x not in munge_filt['ID'].tolist()]
    munge_copy = munge_copy.loc[temp,['ID', 'rsid', 'pos38', 'pval']]
    #Set file path for filtered out SNPs and export
    filt_temp = out_dir + file_prefix + '.filtered_out.txt'
    munge_copy.to_csv(filt_temp, sep = "\t", header = True, index = False, na_rep='NA')
    
    ##### Prep variables/files for CAFEH #####
    study = ['study0']
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
    
    ## Finish out z-score data
    z = munge_filt[['ID','z']]
    z = z.T
    z.columns = z.iloc[0]
    z = z[1:] 
    z.columns.name=''
    z.index= study
    z = z.astype(np.float32)
    
    ## sample size - snps as columns, studies (1) as rows
    sample = munge_filt[['ID','num_samples']]
    sample = sample.T
    sample.columns = sample.iloc[0]
    sample = sample[1:] 
    sample.columns.name=''
    sample.index= study
    
    ##fit with betas and standard errors
    np.random.seed(random_seed) # Set random seed to ensure reproducible results
    n = np.max(sample.max())
    cafehs = fit_cafeh_z(ld_filt, z, n=n)
    #Get number of components
    K = len(cafehs.purity)
    
    #Store extra data in the cafehs object
    cafehs.neglogP = np.array(munge_filt['-logP'])
    cafehs.bp = np.array(munge_filt['pos37'])
    cafehs.bp38 = np.array(munge_filt['pos38'])
    cafehs.chr = np.array(munge_filt['chrom'])
    cafehs.rsid = np.array(munge_filt['rsid'])
    cafehs.ID = np.array(munge_filt['ID'])
    
    #Hard print confidence level and purity in the object as well
    cafehs.realpure = copy.deepcopy(cafehs.purity) #Need to prevent recalculation every time!!!
    cafehs.conf = confidence
    
    ##### Calculate residuals and add them to the cafehs object #####
    for i in range(K):
        #Calculate test statistics for preceding residuals
        if i > 0:
            precede_prediction = np.sum([cafehs.compute_first_moment(l) for l in range(cafehs.dims['K']) if l < i], axis=0)
            precede_test = cafehs.B - precede_prediction #Residuals when predicting betas given the preceding components
            precede_test2 = precede_test/cafehs.S #Estimated Z-scores when predicting given the preceding components
        else:
            precede_test2 = cafehs.B/cafehs.S #Actual Z-scores when looking at the first component
        #Calculate test statistics for absolute residuals
        abs_prediction = np.sum([cafehs.compute_first_moment(l) for l in range(cafehs.dims['K']) if l != i], axis=0)
        abs_test = cafehs.B - abs_prediction #Residuals when predicting betas given all other components
        abs_test2 = abs_test/cafehs.S #Estimated Z-scores when predicting given all other components
        #Calculate pvals for the residuals
        abs_neglogp_residual = (-np.log(2) - t.logcdf(-np.abs(abs_test2), df=n-1))/np.log(10)
        precede_neglogp_residual = (-np.log(2) - t.logcdf(-np.abs(precede_test2), df=n-1))/np.log(10)
        if i == 0:
            cafehs.abs_neglogp_resid = abs_neglogp_residual
            cafehs.precede_neglogp_resid = precede_neglogp_residual
        else:
            cafehs.abs_neglogp_resid = np.concatenate((cafehs.abs_neglogp_resid, abs_neglogp_residual), axis=0)
            cafehs.precede_neglogp_resid = np.concatenate((cafehs.precede_neglogp_resid, precede_neglogp_residual), axis=0)
    
    #Save pickle file
    cafehs.save(f'{out_dir}/{file_prefix}.pkl', save_data=True)
    
###############################################################################


#cProfile.run('main(sys.argv[1:])')
         
if __name__ == "__main__":
    main(sys.argv[1:])
