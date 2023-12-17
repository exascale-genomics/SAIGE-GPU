'''
make_downsampling_variant_lists_wrapper.py

The purpose of this script is to be a wrapper that determines for each trait 
which loci need to be rerun for the down-sampling analysis. It prepares a file 
for each trait of awk commands which will then be submitted by the shell 
script.

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
    -i => <txt> Input file of fine-mapped signals REQUIRED
    -c => <txt> Directory of command files to be submitted to the console REQUIRED
    -o => <txt> Directory for outputting variant lists REQUIRED
    -s => <txt> Directory of summary statistics for identifying variants (default is /gpfs/alpine/med112/proj-shared/results/FOR_SUSIE)
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
        opts, args = getopt.getopt(sys.argv[1:], "i:c:o:s:n:t:nh")
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
        command_dir = options_dict['-c']
        output_dir = options_dict['-o']
    except KeyError:
        print("ERROR: One of your required arguments does not exist.")
        help()
    
    #Optional arguments
    sum_stats_dir = options_dict.get('-s', "/gpfs/alpine/med112/proj-shared/results/FOR_SUSIE")
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
    if exists(signal_file) == False:
        print("ERROR: Input file of merged signals does not exist.")
        sys.exit(1)
    #Check that directories exist
    if exists(output_dir) == False:
        print("ERROR: Output directory does not exist.")
        sys.exit(1)
    if exists(command_dir) == False:
        print("ERROR: Directory for command files does not exist.")
        sys.exit(1)
    if exists(sum_stats_dir) == False:
        print("ERROR: Summary statistics directory does not exist.")
        sys.exit(1)

    print("Acceptable Inputs Given")
    #Call driver function
    driver(signal_file, command_dir, output_dir, sum_stats_dir, anc_1, anc_2)

###############################################################################
#############################  TEST LOCATIONS ##############################
###############################################################################

signal_file = "C:/Users/mitch/Documents/UPenn/Voight_Lab/MVP_Finemapping/Figure_Building_and_Synthesis/master.signals.txt"
command_dir = "C:/Users/mitch/Documents/UPenn/Voight_Lab/MVP_Finemapping/testing/downsampling_commands"
output_dir = "C:/Users/mitch/Documents/UPenn/Voight_Lab/MVP_Finemapping/testing/downsampling_variant_lists"
sum_stats_dir = "C:/Users/mitch/Documents/UPenn/Voight_Lab/MVP_Finemapping/testing/downsampling_sumstats"
anc_1 = "EUR"
anc_2 =  "AFR"

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
def driver(signal_file, command_dir, output_dir, sum_stats_dir, anc_1, anc_2):
    
    #Add slashes to directories
    output_dir = add_slash(output_dir)
    command_dir = add_slash(command_dir)
    sum_stats_dir = add_slash(sum_stats_dir)
    
    #Read in the registry file
    signal_raw = pd.read_csv(signal_file, sep="\t")
    #Filter for the loci that were mapped in both ancestries
    signal_filt = signal_raw[~(signal_raw[anc_1 + '.num_variants'].isnull()) &  ~(signal_raw[anc_2 + '.num_variants'].isnull())]
    
    #Get list of unique traits
    traits = signal_filt['trait'].unique()
    
    #Loop over the traits and prepare command files for each
    for trait in traits:
        #Create a blank list that will hold the commands
        trait_coms = []
        #Pull out the signals for the trait
        signal_trait = signal_filt[signal_filt['trait'] == trait]
        #Create file names for command files
        trait_command_file = command_dir + trait + '.downsampling_variant_list_commands.txt'
        #Extrapolate file names for the summary statistics
        sum_stats_file_1 = sum_stats_dir + trait + '.' + anc_1 + ".GIA.KDI.txt.gz"
        #Extrapolate file names for the summary statistics
        variant_out_file_1 = output_dir + trait + '.' + anc_1 + ".downsampling_variants.txt"
        #Extract the list of loci for the trait
        signal_trait_loci = signal_trait['locus'].unique()
        #Write command for clearing the variant files
        trait_coms.append(encode_command(">" + variant_out_file_1))
        #Loop over the loci for the trait
        for locus in signal_trait_loci:
            #Split up the locus
            locus_pieces = locus.split(".")
            chromo = locus_pieces[0].replace("chr", "")
            #Write awk commands for locus to file
            trait_coms.append(encode_command("zcat " + sum_stats_file_1 + " | awk -F \"\\t\" '{print $1}' | awk -F \":\" 'NR > 1 && $1 == " + chromo + " && " + locus_pieces[1] + " <= $2 && " +  locus_pieces[2] + " >= $2 {print $0}' >> " + variant_out_file_1))
        #Add a do-nothing line at the end of the file
        trait_coms.append("echo___\"empty\"")
        #Print the list of commands to a file
        with open(trait_command_file, "w") as out_path:
            out_path.write('\n'.join(trait_coms))
    

###############################################################################


#cProfile.run('main(sys.argv[1:])')
         
if __name__ == "__main__":
    main(sys.argv[1:])

