'''
save_empty_pkl.py

This script will take a location and save an empty pickle file there.

***This needs to be run on Python 3.9 as there an incompatability with other ***
***versions of the .to_list() function  and the CAFEH package for python.    ***
'''

###############################################################################

#Load libraries
import pickle
import sys
from os.path import exists, dirname, basename, splitext
import getopt

#Define help function to be called if there is a problem
def help(exit_num=1):
    print("""-----------------------------------------------------------------
ARGUMENTS
    -p => <txt> Location to save empty pickle file REQUIRED
""")
    sys.exit(exit_num)

###############################################################################
###########################  COLLECT AND CHECK ARGUMENTS  #####################
###############################################################################

## GLOBAL VARS

## MAIN ##

def main(argv): 
    try: 
        opts, args = getopt.getopt(sys.argv[1:], "p:nh")
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
        file_loc = options_dict['-p']

    except KeyError:
        print("ERROR: One of your required arguments does not exist.")
        help()
    
    #Confirm that the file/folder locs exist
    if exists(dirname(file_loc)) == False:
        print("ERROR: Directory for specified file location does not exist.")
        sys.exit(1)
    print("Acceptable Inputs Given")
    
    #Call driver function
    driver(file_loc)    
    
###############################################################################
#############################  DRIVER  ########################################
###############################################################################

def driver(file_loc): 
    
    #Check file extension
    if basename(splitext(file_loc)) != '.pkl':
        print("WARNING: Non .pkl file extension given for pickle file. Creating " + file_loc + " anyways.")
    #Save empty pickle file
    outfile = open(file_loc,'wb')
    pickle.dump(None,outfile)
    outfile.close()
    
    
    