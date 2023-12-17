'''
define_downsampling_loci.py

The purpose of this script is to make a list of loci for the downsampling
fine-mapping experiment that is compatible with our pipeline. 

***This needs to be run on Python/3.9 as there's an incompatability with ***
***other versions of python.                                             ***
'''

#import modules that will be needed
import sys #Taking command line parameters as inputs
import getopt
from os.path import exists
import os
import pandas as pd
import numpy as np

#Define help function to be called if there is a problem
def help(exit_num=1):
    print("""-----------------------------------------------------------------
ARGUMENTS
    -i => <txt> Input file of fine-mapped signals REQUIRED
    -o => <txt> Output location for locus file REQUIRED
    -p => <txt> Output location for list of traits REQUIRED
    -n => <txt> Ancestry that will be downsampled (default is EUR) OPTIONAL
    -t => <txt> Ancestry that will remain at full strength (default is AFR) OPTIONAL
""")
    sys.exit(exit_num)

###############################################################################
###########################  COLLECT AND CHECK ARGUMENTS  #####################
###############################################################################

## GLOBAL VARS

## MAIN ##

def main(argv): 
    try: 
        opts, args = getopt.getopt(sys.argv[1:], "i:o:p:n:t:nh")
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
        signal_file = options_dict['-i']
        output_loc = options_dict['-o']
        trait_loc = options_dict['-p']
    except KeyError:
        print("ERROR: One of your required arguments does not exist.")
        help()
    
    #Optional arguments
    anc_1 = options_dict.get('-n', "EUR")
    anc_2 = options_dict.get('-t', "AFR")
    
    #Confirm that ancestry is a valid option
    if type(anc_1) == str and type(anc_2) == str:
        anc_1 = anc_1.upper()
        anc_2 = anc_2.upper()
        if anc_1 not in ['AFR', 'AMR', 'EAS', 'EUR'] or anc_2 not in ['AFR', 'AMR', 'EAS', 'EUR']:
            print("ERROR: Ancestries must be AFR, AMR, EAS, or EUR.")
            sys.exit(1)
    else:
        print("ERROR: Ancestries must be strings equivalent to AFR/AMR/EAS/EUR.")
        sys.exit(1)
    #Confirm that the file locs exist
    if exists(os.path.dirname(output_loc)) == False:
        print("ERROR: Locus file output location is unwritable.")
        sys.exit(1)
    if exists(os.path.dirname(trait_loc)) == False:
        print("ERROR: Trait output file location is unwritable.")
        sys.exit(1)
    if exists(signal_file) == False:
        print("ERROR: Input file of merged signals does not exist.")
        sys.exit(1)

    print("Acceptable Inputs Given")
    #Call driver function
    driver(signal_file, output_loc, trait_loc, anc_1, anc_2)

###############################################################################
#############################  HELPFUL FUNCTIONS ##############################
###############################################################################

#Function that adds slashes to directory names
def add_slash(directory):
    if directory[-1] != '/':
        return directory + '/'
    else:
        return directory
    
#Function that encodes commands for decoder
def encode_command(locus_cmd):
    locus_cmd = locus_cmd.replace(" ", "___")
    locus_cmd = locus_cmd.replace(";", "^^^")
    return(locus_cmd)

###############################################################################
#############################  DRIVER  ########################################
###############################################################################

## drive the script ##
## ONE-TIME CALL -- called by main
def driver(signal_file, output_loc, trait_loc, anc_1, anc_2):

    #Read in the registry file
    signal_raw = pd.read_csv(signal_file, sep="\t")
    #Filter for the loci that were mapped in both ancestries
    signal_filt = signal_raw[~(signal_raw[anc_1 + '.num_variants'].isnull()) &  ~(signal_raw[anc_2 + '.num_variants'].isnull())]
    
    signal_filt['chr'] = signal_filt['chr'].str.replace('chr', '').astype(int)
    #Scrub duplicates from the list of loci and adjust range boundaries where needed
    signal_filt.drop_duplicates(subset=['trait','locus'], keep='first', inplace=True)
    signal_filt = signal_filt.sort_values(['chr', 'start', 'end', 'trait'], ascending=[True, True, True, True])
    final_ranges = [[signal_filt['chr'].iloc[1], signal_filt['start'].iloc[1], signal_filt['end'].iloc[1], signal_filt['end'].iloc[1] - signal_filt['start'].iloc[1] + 1, signal_filt['trait'].iloc[1]]]
    for i in range(0,len(signal_filt)):
        temp = [signal_filt['chr'].iloc[i], signal_filt['start'].iloc[i], signal_filt['end'].iloc[i], signal_filt['end'].iloc[i] - signal_filt['start'].iloc[i] + 1, signal_filt['trait'].iloc[i]]
        #Check if range already exists yet
        if temp[:3] != final_ranges[-1][:3]:
            final_ranges.append(temp)
        else:
            final_ranges[-1][-1] = final_ranges[-1][-1] + '/' + temp[-1]
    
    #Final sort
    final_ranges = sorted(final_ranges, key = lambda x: (x[3]))
    
    #Write to file
    import csv
    with open(output_loc, "w") as out_path:
        writer = csv.writer(out_path, delimiter = ',')
        for i in final_ranges:
            writer.writerow(i)
        out_path.close()
    
    #Get list of unique traits and write to file
    traits = signal_filt['trait'].unique()
    traits = np.append(traits, "")
    with open(trait_loc, "w") as out_path:
        out_path.write('\n'.join(traits))
    


###############################################################################


#cProfile.run('main(sys.argv[1:])')
         
if __name__ == "__main__":
    main(sys.argv[1:])

