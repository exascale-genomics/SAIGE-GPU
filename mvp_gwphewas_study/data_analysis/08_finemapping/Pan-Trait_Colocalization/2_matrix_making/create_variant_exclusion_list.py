'''
create_variant_exclusion_list.py

The purpose of this script is to create a list of the multiallelic SNPs and 
indels to exlcude from the plink files. These need to be excluded as they 
create problems for merging plink files and could make issue for fine-mapping
due to strand flip errors during impuation.

***This needs to be run on Python/3.9 as there's an incompatability with ***
***other versions of python.                                             ***
'''

#import modules that will be needed
import sys #Taking command line parameters as inputs
import getopt
from os.path import exists
import gzip
import pandas as pd

#Define help function to be called if there is a problem
def help(exit_num=1):
    print("""-----------------------------------------------------------------
ARGUMENTS
    -i => <txt> Input file of variant annotations REQUIRED
    -o => <txt> Output directory for files of variants to be excluded REQUIRED
""")
    sys.exit(exit_num)

###############################################################################
#############################  TEST FILE LOCATIONS ############################
###############################################################################

'''

input_file = "/home/mconery/Matrix_Building/Annotation_gwPheWAS_dbGaP_b38.HARE.txt.gz"
output_dir = "/home/mconery/Matrix_Building"

'''

###############################################################################
###########################  COLLECT AND CHECK ARGUMENTS  #####################
###############################################################################

## GLOBAL VARS

## MAIN ##

def main(argv): 
    try: 
        opts, args = getopt.getopt(sys.argv[1:], "i:o:nh")
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
        input_file = options_dict['-i']
        output_dir = options_dict['-o']
    except KeyError:
        print("ERROR: One of your required arguments does not exist.")
        help()
    
    #Confirm that the file loc and directory exist
    if exists(input_file) == False:
        print("ERROR: Input file of variant annotations does not exist.")
        sys.exit(1)
    if exists(output_dir) == False:
        print("ERROR: Output directory for files of excluded variants does not exist.")
        sys.exit(1)
    print("Acceptable Inputs Given")
    #Call driver function
    driver(input_file, output_dir)

###############################################################################
#############################  HELPFUL FUNCTIONS ##############################
###############################################################################

#Function that adds slashes to directory names
def add_slash(directory):
    if directory[-1] != '/':
        return directory + '/'
    else:
        return directory

#Define a function that opens a file as a python read object
def opener(filename):
    f = open(filename, 'rb')
    if (f.read(2) == b'\x1f\x8b'):
        return gzip.open(filename, mode='rb')
    else:
        f.seek(0)
        return f

#Define a function that figures out how to split a line and returns the separator
def line_splitter(line, inp_separator=False):
    #Test different splits
    space_count = len(line.split(" "))
    none_count = len(line.split()) 
    comma_count = len(line.split(","))
    line_count = len(line.split("|"))
    extra_slash_tab_count = len(line.split("\\t"))
    #Check input_separator if it's not false
    if inp_separator != False:
        custom_count = len(line.split(inp_separator))
    else:
        custom_count = 0 
    #Get max count value, then it's index in a dict, and then return the dict key 
    #i.e. the separator
    count_dict = {" ": space_count, None:none_count, ",":comma_count, "|":line_count, "\\t":extra_slash_tab_count, inp_separator: custom_count}
    max_count_index = list(count_dict.values()).index(max(count_dict.values()))
    output = list(count_dict.keys())[max_count_index]
    return output  

###############################################################################
#############################  DRIVER  ########################################
###############################################################################

## drive the script ##
## ONE-TIME CALL -- called by main
def driver(input_file, output_dir):

    #Add slash to output_dir
    output_dir = add_slash(output_dir)
    #Get file separator
    temp = opener(input_file)
    line = str(temp.readline())
    file_sep = line_splitter(line)
    #Read in file as pandas df
    annotate_raw = pd.read_csv(input_file, sep=file_sep)
    
    #Split the id column
    temp = [x.split(':') for x in annotate_raw['MVPID'].tolist()]
    annotate_raw['CHROM_b37'] = [x[0] for x in temp]
    annotate_raw['POS_b37'] = [x[1] for x in temp]
    annotate_raw['REF'] = [len(x[2]) for x in temp]
    annotate_raw['ALT'] = [len(x[3]) for x in temp]
    
    #Create function that identifies duplicated IDs and exports them to files
    def export_junk(annotate_raw, output_dir, ANC=False):
        if ANC:
            annotate_sort= annotate_raw.sort_values(by = 'AF_' + ANC)
            out_file = output_dir + 'variant_exclude_list.' + ANC + '.txt'
        else:
            annotate_sort= annotate_raw.sort_values(by = 'AF')
            out_file = output_dir + 'variant_exclude_list.GLOBAL.txt'
        #Get multi-allelic exclusion variantss
        annotate_dups = annotate_sort.duplicated(subset = ['CHROM_b37', 'POS_b37'], keep=False)
        exclude_vars = annotate_raw.loc[annotate_dups[annotate_dups].index,'MVPID']
        #Get indel exclusion variants
        indel_1 = annotate_sort[(annotate_sort['REF'] != 1) | (annotate_sort['ALT'] != 1)]
        exclude_vars = pd.concat([exclude_vars, indel_1['MVPID']], axis = 0) 
        exclude_vars = exclude_vars.sort_values()
        exclude_vars = exclude_vars[exclude_vars.duplicated(keep = 'first') == False]
        exclude_vars.to_csv(out_file, index=False, header=False)
    #Call function for single ancestry
    export_junk(annotate_raw, output_dir)

###############################################################################


#cProfile.run('main(sys.argv[1:])')
         
if __name__ == "__main__":
    main(sys.argv[1:])

