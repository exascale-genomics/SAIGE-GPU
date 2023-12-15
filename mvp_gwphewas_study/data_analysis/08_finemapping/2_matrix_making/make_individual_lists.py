'''
make_individual_lists.py

The purpose of this script is to read over the phenotype registry files
and make lists of individuals for each ancestry who were phenotyped for the 
given trait/disease.

***This needs to be run on Python/3.9 as there's an incompatability with ***
***other versions of python.                                             ***
'''

#import modules that will be needed
import sys #Taking command line parameters as inputs
import getopt
from os.path import exists
import os
import pandas as pd

#Define help function to be called if there is a problem
def help(exit_num=1):
    print("""-----------------------------------------------------------------
ARGUMENTS
    -o => <txt> Output directory for lists of individuals REQUIRED
    -b => <txt> Input phenotype file of binary phenotypes (default is /gpfs/alpine/med112/proj-shared/data/phenotype_data/processed/GIA/gwphewas_binary_traits.pheno) OPTIONAL
    -q => <txt> Input phenotype file of quantitative phenotypes (default is /gpfs/alpine/med112/proj-shared/data/phenotype_data/processed/GIA/gwphewas_quantitative_traits.pheno) OPTIONAL
    -d => <txt> Input directory of phenotype-agnostic ancestry lists (default is /gpfs/alpine/med112/proj-shared/GIA/output/projections/smartpca_loadings_projection) OPTIONAL
    -f => <txt> Ancestry file naming convention showing where ANCestry is found relative to periods 
                (default release4.mvp.gia.ANC.txt) OPTIONAL
""")
    sys.exit(exit_num)

###############################################################################
###########################  COLLECT AND CHECK ARGUMENTS  #####################
###############################################################################

## GLOBAL VARS

## MAIN ##

def main(argv): 
    try: 
        opts, args = getopt.getopt(sys.argv[1:], "o:b:q:d:f:nh")
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
        out_dir = options_dict['-o']
    except KeyError:
        print("ERROR: One of your required arguments does not exist.")
        help()
    
    #Optional arguments
    binary_file = options_dict.get('-b', "/gpfs/alpine/med112/proj-shared/data/phenotype_data/processed/GIA/gwphewas_binary_traits.pheno")
    quant_file = options_dict.get('-q', "/gpfs/alpine/med112/proj-shared/data/phenotype_data/processed/GIA/gwphewas_quantitative_traits.pheno")
    anc_dir = options_dict.get('-d', "/gpfs/alpine/med112/proj-shared/GIA/output/projections/smartpca_loadings_projection")
    name_conv = options_dict.get('-f', 'release4.mvp.gia.ANC.txt')
    
    #Confirm that the file locs exist
    if exists(out_dir) == False:
        print("ERROR: Input file of binary phenotypes does not exist.")
        sys.exit(1)
    if exists(binary_file) == False:
        print("ERROR: Input file of binary phenotypes does not exist.")
        sys.exit(1)
    if exists(quant_file) == False:
        print("ERROR: Input file of quantitative phenotypes does not exist.")
        sys.exit(1)
    if exists(anc_dir) == False:
        print("ERROR: Specified directory holding files of individuals by ancestry does not exist.")
        sys.exit(1)
    #Verify that ANC is present in naming convention
    if 'ANC' not in name_conv:
        print("ERROR: Unable to identify position of ANCestry in file name convention.")
        sys.exit(1)
    print("Acceptable Inputs Given")
    
    #Call driver function
    driver(out_dir, binary_file, quant_file, anc_dir, name_conv)

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

def driver(out_dir, binary_file, quant_file, anc_dir, name_conv): 
    
    #Add final slash to out_dir and anc_dir if it's missing
    out_dir = add_slash(out_dir)
    anc_dir = add_slash(anc_dir)
    
    #Set list of ancestries
    ancestries = ['AFR', 'AMR', 'EAS', 'EUR']    
    
    #Identify ancestry position from name_conv
    name_conv_list = name_conv.split(".")
    anc_pos = name_conv_list.index("ANC")
    remain_name_indices=list(range(len(name_conv_list)))
    remain_name_indices.remove(anc_pos)
        
    #Get munged summary stats file loc
    anc_files = os.listdir(anc_dir)
    anc_files = [test.split(".") for test in anc_files]
    anc_files = [test for test in anc_files if len(test) == len(name_conv_list)]
    anc_files = [test for test in anc_files if test[anc_pos] in ancestries]
    for i in remain_name_indices:
        anc_files = [test for test in anc_files if test[i] == name_conv_list[i]]
    #Check if we have an error
    if len(anc_files) != 4: 
        print("ERROR: Incorrect number of files matching name format found in ancestry directory.")
        sys.exit(1)
    
    #Read in ancestry file
    anc_dict = {}
    for anc in ancestries:
        anc_dict[anc] = pd.read_table(anc_dir + name_conv.replace("ANC", anc), header=None, sep=" ").iloc[:,0]
    
    #Create a function to read over one of the phenotype files and make the output lists
    def make_type_lists(pheno_file, anc_dict, out_dir):
        #Read in pheno_file
        pheno_raw = pd.read_table(pheno_file, sep="\t")
        #Loop over phenotypes in the pheno_raw (skip the ID column)
        for pheno in pheno_raw.columns[1:]:
            #Extract list of individuals for phenotype
            pheno_ind = pheno_raw.iloc[pheno_raw[pheno].dropna(how = 'all').index, 0] #Remove non-phenotyped individuals
            #Loop over ancestries
            for anc in anc_dict.keys():
                temp = pd.DataFrame(list(set(anc_dict[anc]) & set(pheno_ind))) #Find the intersection between the phenotyped individuals and those of the ancestry
                #check length of intersection
                if len(temp) > 0:
                    temp[1] = temp[0] #Duplicate iid for famid
                    #Create file name and output file
                    out_file = out_dir + pheno + "." + anc + ".individuals.txt"
                    temp.to_csv(out_file, sep = " ", header = False, index = False, na_rep='NA')
                else: 
                    print("WARNING: No phenotyped " + anc + " individuals detected for " + pheno + ".")
    
    #Call function
    make_type_lists(binary_file, anc_dict, out_dir)
    make_type_lists(quant_file, anc_dict, out_dir)
                
###############################################################################


#cProfile.run('main(sys.argv[1:])')
         
if __name__ == "__main__":
    main(sys.argv[1:])
