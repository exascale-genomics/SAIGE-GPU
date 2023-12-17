'''
convert_to_mappable_list.py

The purpose of this script is to take the list of defined signficant loci that 
result from define_sig_loci_for_trait.py and combine_sig_loci.py and separate
out the loci with multiple traits onto distinct lines, one for each trait. It 
takes as input a single results file and outputs a single results file as well.

***This needs to be run on Python/3.9 as there's an incompatability with ***
***other versions of python.                                             ***
'''

#import modules that will be needed
import sys #Taking command line parameters as inputs
import getopt
from os.path import exists
import gzip
import csv
from itertools import islice

#Define help function to be called if there is a problem
def help(exit_num=1):
    print("""-----------------------------------------------------------------
ARGUMENTS
    -i => <txt> Input file of trait-combined significant loci REQUIRED
    -o => <txt> Output file of trait-separated significant loci REQUIRED
    -m => <txt> Only output multi-trait loci True/False (default is False) OPTIONAL
""")
    sys.exit(exit_num)

###############################################################################
###########################  COLLECT AND CHECK ARGUMENTS  #####################
###############################################################################

## GLOBAL VARS

## MAIN ##

def main(argv): 
    try: 
        opts, args = getopt.getopt(sys.argv[1:], "i:o:m:nh")
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
        input_loc = options_dict['-i']
        output_loc = options_dict['-o']
    except KeyError:
        print("ERROR: One of your required arguments does not exist.")
        help()
    
    #Optional arguments
    multi_trait = options_dict.get('-m', False)
    #Check that multi_trait came through correctly
    if isinstance(multi_trait, str):
        if multi_trait[0].upper() == 'T':
            multi_trait = True
        else:
            multi_trait = False
    if exists(input_loc) == False:
        print("ERROR: Input GWAS results file does not exist.")
        sys.exit(1)
    print("Acceptable Inputs Given")
    #Call driver function
    driver(input_loc, output_loc, multi_trait)

###############################################################################
#############################  HELPFUL FUNCTIONS ##############################
###############################################################################

#Function that adds slashes to directory names
def add_slash(directory):
    if directory[-1] != '/':
        return directory + '/'
    else:
        return directory

#Make a function to open files
def opener(filename):
    f = open(filename, 'rb')
    if (f.read(2) == b'\x1f\x8b'):
        return gzip.open(filename, mode='rb')
    else:
        f.seek(0)
        return f

###############################################################################
#############################  DRIVER  ########################################
###############################################################################

## drive the script ##
## ONE-TIME CALL -- called by main
def driver(input_loc, output_loc, multi_trait):

    #Open file and write location
    input_file = opener(input_loc)
    with open(output_loc, "w") as out_file:
        writer = csv.writer(out_file, delimiter = ',')
        #Loop over lines in file
        for line in islice(input_file, 0, None, 1):
            line = str(line)
            line = line.replace('\\n\'', '')
            line = line.replace('b\'', '')
            input_row = line.split(",")
            if '/' in input_row[-1]:
                temp = input_row[-1].split('/')
                for trait in temp:
                    output_row = input_row[:-1]
                    output_row.append(trait)
                    writer.writerow(output_row)
            elif multi_trait == False:
                writer.writerow(input_row)
    #Close file locations
    out_file.close()
    input_file.close()


###############################################################################


#cProfile.run('main(sys.argv[1:])')
         
if __name__ == "__main__":
    main(sys.argv[1:])

