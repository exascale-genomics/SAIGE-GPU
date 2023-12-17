'''
compare_matrices.py

This script is designed to check a matrix that is a subset of another matrix 
against the matrix it was filtered from. It's designed to find any errors 
created during the processing

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
from math import isnan

#Define help function to be called if there is a problem
def help(exit_num=1):
    print("""-----------------------------------------------------------------
ARGUMENTS
    -a => <txt> location of original matrix (ld) file REQUIRED
    -b => <txt> location of original map file REQUIRED
    -x => <txt> location of subset matrix (ld) file REQUIRED
    -y => <txt> location of subset map file REQUIRED
""")
    sys.exit(exit_num)

###############################################################################
###########################  COLLECT AND CHECK ARGUMENTS  #####################
###############################################################################

## GLOBAL VARS

## MAIN ##

def main(argv): 
    try: 
        opts, args = getopt.getopt(sys.argv[1:], "a:b:x:y:nh")
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
        orig_ld_loc = options_dict['-a']
        orig_map_loc = options_dict['-b']
        temp_ld_loc = options_dict['-x']
        temp_map_loc = options_dict['-y']

    except KeyError:
        print("ERROR: One of your required arguments does not exist.")
        help()
    
    #Confirm that the file/folder locs exist
    if exists(orig_ld_loc) == False:
        print("ERROR: Original matrix file location does not exist.")
        sys.exit(1)
    if exists(orig_map_loc) == False:
        print("ERROR: Original map file location does not exist.")
        sys.exit(1)
    if exists(temp_ld_loc) == False:
        print("ERROR: Temporary matrix file location does not exist.")
        sys.exit(1)
    if exists(temp_map_loc) == False:
        print("ERROR: Temporary map file location does not exist.")
        sys.exit(1)
    
    print("Acceptable Inputs Given")
    #Call driver function
    driver(orig_ld_loc, orig_map_loc, temp_ld_loc, temp_map_loc)

###############################################################################
#############################  Testing Locations ##############################
###############################################################################

'''

orig_ld_loc = "/home/mconery/Matrix_Building/plink_testing_dir/chr1.500001.1250000.EUR.ld"
orig_map_loc = "/home/mconery/Matrix_Building/plink_testing_dir/chr1.500001.1250000.EUR.map"
temp_ld_loc = "/home/mconery/Matrix_Building/plink_testing_dir/chr1.500001.1250000.EUR.matrix"
temp_map_loc = "/home/mconery/Matrix_Building/plink_testing_dir/chr1.500001.1250000.EUR.map"

'''

###############################################################################
#############################  DRIVER  ########################################
###############################################################################

def driver(orig_ld_loc, orig_map_loc, temp_ld_loc, temp_map_loc): 
    
    #Set an indicator variable to check for errors
    error_ind = False
    
    #Read in the ld variant index information
    orig_map_raw = pd.read_table(orig_map_loc, sep="\t", header=None)
    temp_map_raw = pd.read_table(temp_map_loc, sep="\t", header=None)
    
    #Set the original map file index to a column and reset it with the variant ids
    orig_map_raw['orig_index'] = orig_map_raw.index
    orig_map_raw.index = orig_map_raw[1]
    #Filter the original map file for the variants in the temp file and verify the order is preserved
    orig_map_filt = orig_map_raw.filter(items=temp_map_raw[1].tolist(), axis=0)
    if sorted(orig_map_filt['orig_index'].tolist()) != orig_map_filt['orig_index'].tolist():
        print("WARNING: Order of map file has gotten mixed.")
        error_ind = True
    else:
        print("NOTE: Map file order check passes successfully.")
    
    #Set the temp file index equal to the ids
    temp_map_raw['orig_index'] = temp_map_raw.index
    temp_map_raw.index = temp_map_raw[1]
    
    #Next read in the ld matrices
    orig_ld_raw = np.loadtxt(orig_ld_loc, dtype=float)
    temp_ld_raw = np.loadtxt(temp_ld_loc, dtype=float)
    #Convert them to pandas dfs
    orig_ld_df = pd.DataFrame(orig_ld_raw, columns = orig_map_raw[1])
    orig_ld_df.index = orig_ld_df.columns
    temp_ld_df = pd.DataFrame(temp_ld_raw, columns = temp_map_raw[1])
    temp_ld_df.index = temp_ld_df.columns
    
    #Iterate over indices in temp matrix and make sure that values match
    print("NOTE: Beginning individual value checks.")
    for i in temp_ld_df.index:
        for j in temp_ld_df.index:
            if isnan(temp_ld_df.loc[i,j]) and isnan(orig_ld_df.loc[i,j]):
                pass
            elif temp_ld_df.loc[i,j] != orig_ld_df.loc[i,j]:
                print("WARNING: " + i + ', ' + j + ' does not match between the original and temp matrices.')
                print(('Orig value is ' + str(orig_ld_df.loc[i,j])))
                print(('Subset value is ' + str(temp_ld_df.loc[i,j])))
                error_ind = True
        print("NOTE: Finished row " + i + " value checks")
    
    #Print out message if any errors have occured
    if error_ind == True:
        print("WARNING: One or more inconsistencies detected between the original and subset matrices.")
    else:
        print("SUCCESS: No errors detected.")
    
###############################################################################


#cProfile.run('main(sys.argv[1:])')
         
if __name__ == "__main__":
    main(sys.argv[1:])
