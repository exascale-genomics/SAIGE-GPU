'''
remove_grch38_loci.py

The purpose of this script is to remove loci that were only defined in GRCh38
and don't share borders with a locus defined in GRCh37. This script is needed
to clean-up a bug in the define_loci.py code that made loci from the position
column (GRCh38) instead of the ID column (GRCh37).

***This needs to be run on Python/3.9 as there's an incompatability with ***
***other versions of python.                                             ***
'''

#import modules that will be needed
import sys #Taking command line parameters as inputs
import getopt
from os.path import exists
import os

#Define help function to be called if there is a problem
def help(exit_num=1):
    print("""-----------------------------------------------------------------
ARGUMENTS
    -d => <txt> Input file of trait-combined significant loci REQUIRED
    -l => <txt> Folder in which to output pvar files REQUIRED
    -n => <txt> Naming Convention (default is "CHR.START.END.ANC.temp.pvar") OPTIONAL
""")
    sys.exit(exit_num)

###############################################################################
###########################  COLLECT AND CHECK ARGUMENTS  #####################
###############################################################################

## GLOBAL VARS

## MAIN ##

def main(argv): 
    try: 
        opts, args = getopt.getopt(sys.argv[1:], "d:l:n:nh")
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
        del_dir = options_dict['-d']
        locus_file = options_dict['-l']
    except KeyError:
        print("ERROR: One of your required arguments does not exist.")
        help()
    
    #Confirm that the file loc and directory exist
    if exists(locus_file) == False:
        print("ERROR: Input file of merged significant loci does not exist.")
        sys.exit(1)
    if exists(del_dir) == False:
        print("ERROR: Directory with files to be deleted does not exist.")
        sys.exit(1)
    
    #Optional arguments
    suffix = options_dict.get('-n', "CHR.START.END.ANC.temp.pvar")
    print("Acceptable Inputs Given")
    #Call driver function
    driver(del_dir, locus_file, suffix)

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
def driver(del_dir, locus_file, name_conv):
    
    #Add slash
    del_dir = add_slash(del_dir)
    #Find positions in naming convention
    temp = name_conv.split(".")
    chr_pos = temp.index("CHR")
    start_pos = temp.index("START")
    end_pos = temp.index("END")
    suffix = temp[-1]
    name_len = len(temp)
    
    #Read in list of loci
    loci_list = []  
    with open(locus_file) as locus_read:
        for line in locus_read.readlines():
            temp = line.strip().split(",") 
            loci_list.append("chr" + temp[0] + "." + temp[1] + "." + temp[2])
    
    #List all files in directory and filter for those whose names have the right length and suffix 
    all_files = os.listdir(del_dir)
    all_files = [x for x in all_files if len(x.split(".")) == name_len and x.split(".")[-1] == suffix]
    
    #loop over the remaining files and test each for inclusion in the loci_list or delete if not
    for file in all_files:
        #Check if file is in loci list
        temp = file.split(".")
        temp2 = temp[chr_pos] + "." + temp[start_pos] + "." + temp[end_pos]
        if temp2 not in loci_list:
            os.remove(del_dir + file)
    

###############################################################################


#cProfile.run('main(sys.argv[1:])')
         
if __name__ == "__main__":
    main(sys.argv[1:])

