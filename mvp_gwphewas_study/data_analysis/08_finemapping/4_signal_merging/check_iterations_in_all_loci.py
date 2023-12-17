'''
check_iterations_in_all_loci.py

The purpose of this script is to iterate over the list of loci/trait pairs that
were fine-mapped and call an rscript to output the number of iterations run for 
each locus.

***This needs to be run on Python/3.9 and R/4.0 as there's an               ***
***incompatability with other versions of python.                           ***
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
    -a => <str> Directory for AFR fine-mapping results (default is /gpfs/alpine/med112/proj-shared/results/PostGWAS_Analysis/finemaping/GIA/susie_dir_AFR) OPTIONAL
    -b => <str> Directory for AMR fine-mapping results (default is /gpfs/alpine/med112/proj-shared/results/PostGWAS_Analysis/finemaping/GIA/susie_dir_AMR) OPTIONAL
    -c => <str> Directory for EAS fine-mapping results (default is /gpfs/alpine/med112/proj-shared/results/PostGWAS_Analysis/finemaping/GIA/susie_dir_EAS) OPTIONAL
    -d => <str> Directory for EUR fine-mapping results (default is /gpfs/alpine/med112/proj-shared/results/PostGWAS_Analysis/finemaping/GIA/susie_dir_EUR) OPTIONAL
""")
    sys.exit(exit_num)

###############################################################################
###########################  COLLECT AND CHECK ARGUMENTS  #####################
###############################################################################

## GLOBAL VARS

## MAIN ##

def main(argv): 
    try: 
        opts, args = getopt.getopt(sys.argv[1:], "m:a:nh")
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
        command_file = options_dict['-m']
        ANC_dir = options_dict['-a']
    except KeyError:
        print("ERROR: One of your required arguments does not exist.")
        help()
    
    #Check that directories exist
    if exists(ANC_dir) == False:
        print("ERROR: ANC directory does not exist.")
        sys.exit(1)
    
    #Call driver function
    driver(command_file, ANC_dir)

###############################################################################
#############################  HELPFUL FUNCTIONS ##############################
###############################################################################

#Function that adds slashes to directory names
def add_slash(directory):
    if directory[-1] != '/':
        return directory + '/'
    else:
        return directory
    
#Create a function that gets the identities of all files in a directory
def make_commands_for_dir(ANC_dir, iter_script, locus_coms):
    #Add slash
    ANC_dir = add_slash(ANC_dir)
    #List files
    all_files = os.listdir(ANC_dir)
    #Get rds files
    rds_files = [x for x in all_files if x.split(".")[-1] == "rds"]
    #Loop over rds file and make commands
    for file in rds_files:
        locus_cmd = "Rscript " + iter_script + " " + ANC_dir + file
        locus_cmd = locus_cmd.replace(" ", "___")
        locus_cmd = locus_cmd.replace(";", "^^^")
        locus_coms.append(locus_cmd)
    #Return full list of locus_coms
    return locus_coms

###############################################################################
#############################  DRIVER  ########################################
###############################################################################

## drive the script ##
## ONE-TIME CALL -- called by main
def driver(command_file, ANC_dir):
    
    #Get location of iteration counting script and verify that exists
    temp = os.path.realpath(__file__)
    iter_script = temp[:temp.find('check_iterations_in_all_loci.py')] + 'check_iterations_and_log_pvals.R'
    if exists(iter_script) == False:
        print("ERROR: check_iterations_and_log_pvals.R script not found in same folder as check_iterations_in_all_loci.py. Consider repulling github.")
        sys.exit(1)

    #Create a blank list that will hold the commands
    locus_coms = []
    
    #Make list of commands
    locus_coms = make_commands_for_dir(ANC_dir, iter_script, locus_coms)
    
    #Verify that list of locus_coms is non-empty (i.e. that at least one locus/trait/anc combo needs mapping)
    if len(locus_coms) > 0:
        #Add a do-nothing line at the end of the file
        locus_coms.append("echo___\"empty\"")
        #Print the list of commands to a file
        with open(command_file, "w") as out_path:
            out_path.write('\n'.join(locus_coms))
    else:
        print("WARNING: No RDS files detected in specified directories.")
    

###############################################################################


#cProfile.run('main(sys.argv[1:])')
         
if __name__ == "__main__":
    main(sys.argv[1:])

