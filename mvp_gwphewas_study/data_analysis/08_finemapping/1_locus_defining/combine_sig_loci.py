'''
combine_sig_loci.py

The purpose of this script is to combine significant loci that are identified 
by define_sig_loci_for_trait.py. The script takes in a directory that contains 
the output files as well as a file format and an ancestry that specifies which
files in the directoryshould be used. It will then output a combined list to 
the given location.

***This needs to be run on Python/3.9 as there's an incompatability with ***
***other versions of python.                                             ***
'''

#import modules that will be needed
import sys #Taking command line parameters as inputs
import getopt
import os
import pandas as pd
import numpy as np
from os.path import exists

#Define help function to be called if there is a problem
def help(exit_num=1):
    print("""-----------------------------------------------------------------
ARGUMENTS
    -l => <txt> Loci File Directory REQUIRED
    -o => <txt> output file location  REQUIRED
    -a => <txt> ancestry of selected files if multiple found in directory OPTIONAL (default is none)
    -f => <str> file name format for loci file showing where ANCestry and TRAIT are found relative to periods 
                (default is XXX.XXXX.TRAIT.ANC.XXX.csv)
""")
    sys.exit(exit_num)

###############################################################################
###########################  COLLECT AND CHECK ARGUMENTS  #####################
###############################################################################

## GLOBAL VARS

## MAIN ##

def main(argv): 
    try: 
        opts, args = getopt.getopt(sys.argv[1:], "l:o:a:f:nh")
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
        loci_folder = options_dict['-l']
        output_loc = options_dict['-o']

    except KeyError:
        print("ERROR: One of your required arguments does not exist.")
        help()
    
    ancestry = options_dict.get('-a', False)
    file_format = options_dict.get('-f', 'xxxx.xxxx.TRAIT.ANC.xxx.csv')    
    #Verify that ANC and TRAIT are present in file format
    if 'ANC' not in file_format or 'TRAIT' not in file_format:
        print("ERROR: Unable to identify position of ANCestry or TRAIT in file name format.")
        sys.exit(1)
    print("Acceptable Inputs Given")
    #Call driver function
    driver(loci_folder, output_loc, ancestry, file_format)

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

## drive the script ##
## ONE-TIME CALL -- called by main
def driver(loci_folder, output_loc, ancestry, file_format):

    #Identify ancestry and trait position from file_format
    temp = file_format.split(".")
    anc_pos = temp.index("ANC")
    trait_pos = temp.index("TRAIT")
    
    #Get list of all traits in the directory
    loci_files = os.listdir(loci_folder)
    #Filter for those of a particular ancestry
    if ancestry:
        loci_files = [file for file in loci_files if ancestry == file.split(".")[anc_pos]]
    
    #Check if the munged_folder ends in a slash
    loci_folder = add_slash(loci_folder)
    
    #Make a list to hold the output data
    out_ranges = [[]]
    #Read in and process the files one at a time
    for file in loci_files:
        #Print trait name
        trait = file.split(".")[trait_pos]
        print(trait)
        #Read in chromosome lengths file and store rows in out_ranges
        with open(loci_folder + file) as file_read:
            for line in file_read.readlines():
                temp = line.strip().split(",") 
                temp[0] = int(temp[0])
                temp[1] = int(temp[1])
                temp[2] = int(temp[2])
                temp[3] = int(temp[3])
                out_ranges.append(temp)
    del out_ranges[0]
    
    #Scrub duplicates from the list of loci and adjust range boundaries where needed
    out_ranges = sorted(out_ranges, key = lambda x: (x[0], x[1], x[2]))
    final_ranges = [out_ranges[0][:]]
    for i in range(1,len(out_ranges)):
        #Check if range already exists yet
        if out_ranges[i][:3] != final_ranges[-1][:3]:
            final_ranges.append(out_ranges[i])
        else:
            final_ranges[-1][-1] = final_ranges[-1][-1] + '/' + out_ranges[i][-1]
    
    #Final sort
    final_ranges = sorted(final_ranges, key = lambda x: (x[3]))
    
    #Write to file
    import csv
    with open(output_loc, "w") as out_path:
        writer = csv.writer(out_path, delimiter = ',')
        for i in final_ranges:
            writer.writerow(i)
        out_path.close()
        

###############################################################################


#cProfile.run('main(sys.argv[1:])')
         
if __name__ == "__main__":
    main(sys.argv[1:])
