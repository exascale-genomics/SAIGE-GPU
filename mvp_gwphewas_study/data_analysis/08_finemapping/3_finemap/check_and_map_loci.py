'''
check_and_map_loci.py

The purpose of this script is to iterate over the list of loci/trait pairs to 
finemap, and decide which need to be fine-mapped based on whether the matrices 
have been successfully generated, whether any of the ancestries lacked data, 
and whether the mapping has already been completed. If all three conditions are 
met successfully, this script outputs an encoded command to map the pair to a 
text file that can be decoded by execute_command.py.

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
    -l => <str> List of trait-combined loci output from combine_loci.py REQUIRED
    -c => <str> Output location for command file REQUIRED
    -o => <str> Output directory for fine-mapping results (default is /gpfs/alpine/med112/proj-shared/results/PostGWAS_Analysis/finemaping/GIA/susie_dir) OPTIONAL
    -m => <str> Folder in which matrix/map files are stored (default is /gpfs/alpine/med112/proj-shared/results/PostGWAS_Analysis/finemaping/GIA/matrix_dir) OPTIONAL
    -s => <str> Directory of summary statistics (default is /gpfs/alpine/med112/proj-shared/results/FOR_SUSIE) OPTIONAL
    -a => <str> Ancestry choice of AFR, AMR, EAS, EUR, or ALL (default is ALL) OPTIONAL
    -r => <str> flip orientation of variants with allele frequencies > 0.5 (default is False) OPTIONAL
    -t => <str> Type of fine-mapping (SuSiE or CAFEH; default is SuSiE) OPTIONAL
    -f => <num> minimum minor allele frequency (default is 0 / no-minimum) OPTIONAL
    -n => <str> Summary stats file naming convention showing where ANCestry and TRAIT are found relative to periods 
                (default is TRAIT.ANC.GIA.KDI.txt.gz) OPTIONAL
""")
    sys.exit(exit_num)

###############################################################################
###########################  COLLECT AND CHECK ARGUMENTS  #####################
###############################################################################

## GLOBAL VARS

## MAIN ##

def main(argv): 
    try: 
        opts, args = getopt.getopt(sys.argv[1:], "l:c:o:m:s:a:r:t:f:n:nh")
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
        locus_file = options_dict['-l']
        command_file = options_dict['-c']        
    except KeyError:
        print("ERROR: One of your required arguments does not exist.")
        help()
    
    #Optional arguments
    output_dir = options_dict.get('-o', "/gpfs/alpine/med112/proj-shared/results/PostGWAS_Analysis/finemaping/GIA/susie_dir")
    matrix_dir = options_dict.get('-m', "/gpfs/alpine/med112/proj-shared/results/PostGWAS_Analysis/finemaping/GIA/matrix_dir")
    stats_dir = options_dict.get('-s', "/gpfs/alpine/med112/proj-shared/results/FOR_SUSIE")
    ancestries = options_dict.get('-a', "ALL")
    swap_af_ind = options_dict.get('-r', 'False')
    map_type = options_dict.get('-t', 'susie')
    min_maf = options_dict.get('-f', 0)
    file_format = options_dict.get('-n', 'TRAIT.ANC.GIA.KDI.txt.gz')
    #Confirm that output_type is a valid option
    if type(map_type) == str:
        map_type = map_type.lower()
        if map_type not in ['susie', 'cafeh']:
            print("ERROR: Mapping Type must be SuSiE or CAFEH.")
            sys.exit(1)
    else:
        print("ERROR: Mapping Type must be a string equivalent to SuSiE or CAFEH.")
        sys.exit(1)
    #Confirm that the file locs exist
    if exists(locus_file) == False:
        print("ERROR: Input file of merged significant loci does not exist.")
        sys.exit(1)
    #Check that directories exist
    if exists(output_dir) == False:
        print("ERROR: Output directory for fine-mapping results does not exist.")
        sys.exit(1)
    if exists(matrix_dir) == False:
        print("ERROR: Map/matrix directory does not exist.")
        sys.exit(1)
    if exists(stats_dir) == False:
        print("ERROR: Invalid directory given for summary statistics location.")
        sys.exit(1)
    #Verify that ANC and TRAIT are present in naming convention
    if 'ANC' not in file_format:
        print("ERROR: Unable to identify position of ANCestry in file name convention.")
        sys.exit(1)
    elif 'TRAIT' not in file_format:
        print("ERROR: Unable to identify position of TRAIT in file name convention.")
        sys.exit(1)
    #Verify that ancestry option is acceptable
    if type(ancestries) == str:
        if ancestries.upper() in ['AFR', 'AMR', 'EAS', 'EUR']:
            ancestries = [ancestries.upper()]             
        elif ancestries.upper() == 'ALL':
            ancestries = ['AFR', 'AMR', 'EAS', 'EUR']
        else:
            print("ERROR: Ancestry be equivalent to AMR, AFR, EAS, EUR, or ALL.")
            sys.exit(1)
    else:
        print("ERROR: Ancestry must be a string equivalent to AMR, AFR, EAS, EUR, or ALL.")
        sys.exit(1)
    #Recast minimum allele frequency 
    try:
        min_maf = float(min_maf)
        if min_maf >= 1 or min_maf < 0: 
            print("ERROR: Minimum allele frequency for variants should be a proportion in [0,1).")
            sys.exit(1)
    except ValueError:
        print("ERROR: Minimum allele frequency is not coercible to a float. It should be a proportion in [0,1).")
        sys.exit(1)
    #Recast swap_af_ind
    if type(swap_af_ind) == str:
        if swap_af_ind.lower() == 'false' or swap_af_ind.lower() == 'f' or swap_af_ind.lower() == 'no' or swap_af_ind.lower() == 'n':
            swap_af_ind = False
        elif swap_af_ind.lower() == 'true' or swap_af_ind.lower() == 't' or swap_af_ind.lower() == 'yes' or swap_af_ind.lower() == 'y':
            swap_af_ind = True
        else:
            print("ERROR: Unrecognized option selected for whether to swap alleles for variants with frequencies > 0.5. Select True/False.")
            sys.exit(1)
    elif type(swap_af_ind) == bool:
        pass #In this case the indicator is already a boolean, so there's no need to recast.
    else:
        print("ERROR: Unrecognized option selected for whether to swap alleles for variants with frequencies > 0.5. Select True/False.")
        sys.exit(1)
    print("Acceptable Inputs Given")
    #Call driver function
    driver(locus_file, command_file, output_dir, matrix_dir, stats_dir, ancestries, swap_af_ind, map_type, min_maf, file_format)

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
def driver(locus_file, command_file, output_dir, matrix_dir, stats_dir, ancestries, swap_af_ind, map_type, min_maf, file_format):
    
    #Create a blank list that will hold the commands
    locus_coms = []
    
    #Set fine-mapping script paths
    temp = os.path.realpath(__file__)
    if map_type == 'susie':
        map_script = temp[:temp.find('check_and_map_loci.py')] + 'finemap_SuSiE.py'
        empty_script = 'Rscript ' + temp[:temp.find('check_and_map_loci.py')] + 'save_empty_rds.R'
        out_file_ext = ".rds"
    elif map_type == 'cafeh':
        map_script = temp[:temp.find('check_and_map_loci.py')] + 'finemap_CAFEH.py'
        empty_script = 'python3 ' + temp[:temp.find('check_and_map_loci.py')] + 'save_empty_pkl.py' + ' -p'
        out_file_ext = ".pkl"
    else:
        print("ERROR: Invalid mapping-type selected.")
        sys.exit(1)
    
    #Create a second output dir that will actually receive the output (Files will be manually transferred after each job to the other output_dir)
    output_dir2 = "%s.2" % (output_dir) 
    
    #Add final slash to directories if they're missing
    output_dir = add_slash(output_dir)
    matrix_dir = add_slash(matrix_dir)
    stats_dir = add_slash(stats_dir)
    output_dir2 = add_slash(output_dir2)

    #Identify ancestry and trait position from file_format
    name_conv = file_format.split(".")
    anc_pos = name_conv.index("ANC")
    trait_pos = name_conv.index("TRAIT")
    
    #Read in the list of loci
    locus_raw = []  
    with open(locus_file) as locus_read:
        for line in locus_read.readlines():
            temp = line.strip().split(",") 
            temp[0] = int(temp[0])
            temp[1] = int(temp[1])
            temp[2] = int(temp[2])
            temp[3] = int(temp[3])
            temp[4] = temp[4].split("/")
            locus_raw.append(temp)
      
    #Create a list of all the matrix/map files in the output directory
    out_files = os.listdir(matrix_dir)
    map_files = [x.replace('.map','') for x in out_files if x.split(".")[-1] == 'map']
    matrix_files = [x.replace('.ld','') for x in out_files if x.split(".")[-1] == 'ld']   
    #Get only map_files and matrix files that are matched
    map_files = [x for x in map_files if x in matrix_files]
    matrix_files = [x for x in matrix_files if x in map_files]
    map_files = [x.split(".") for x in map_files]
    matrix_files = [x.split(".") for x in matrix_files]
    
    #Loop over the locus file
    for locus in locus_raw:
        #Loop over traits
        for trait in locus[4]:
            #Loop over the ancestries
            for ancestry in ancestries:
                #Create skip indicator
                ld_skip_ind = False
                #Set file prefix
                file_prefix = trait + '.chr' + str(locus[0]) + '.' + str(locus[1]) + '.' + str(locus[2]) + '.' + ancestry
                #Check for matching map and matrix files
                locus_map_files = [x for x in map_files if trait == x[0] and 'chr' + str(locus[0]) == x[1] and locus[1] >= int(x[2]) and locus[2] <= int(x[3]) and  ancestry == x[4]]
                locus_matrix_files = [x for x in matrix_files if trait == x[0] and 'chr' + str(locus[0]) == x[1] and locus[1] >= int(x[2]) and locus[2] <= int(x[3]) and ancestry == x[4]]
                if len(locus_map_files) > 0 and len(locus_matrix_files) > 0:
                    #Check for intersections between the two lists
                    locus_ld_files = [x for x in locus_map_files if x in locus_matrix_files]
                    if len(locus_ld_files) == 1:
                        ld_prefix = '.'.join(locus_ld_files[0])
                    elif len(locus_ld_files) > 1:                    
                        temp = [int(x[3]) - int(x[2]) + 1 for x in locus_ld_files]
                        ld_prefix = '.'.join(locus_ld_files[temp.index(min(temp))]) #Finds smallest overlapping matrix
                    else:
                        print("WARNING: No acceptable matrix/map fair found for " + file_prefix + ". Locus will be skipped for this ancestry.")
                        ld_skip_ind = True
                else: #Event where no acceptable matrix/map combo is found
                    print("WARNING: No acceptable matrix/map fair found for " + file_prefix + ". Locus will be skipped for this ancestry.")
                    ld_skip_ind = True
                if ld_skip_ind == False: #Case where we won't skip over the locus for the given ancestry
                    #Get location of corresponding summary statistics
                    stat_loc = name_conv
                    stat_loc[anc_pos] = ancestry
                    stat_loc[trait_pos] = trait
                    stat_loc = stats_dir + '.'.join(stat_loc)
                    #Check that the output file has not already been generated.
                    out_loc = output_dir + file_prefix + out_file_ext
                    if exists(out_loc):
                        print("NOTE: " + file_prefix + " has already been mapped. Skipping.")
                        continue
                    #Check that the summary statistics exist for the given trait for the given ancestry
                    if not exists(stat_loc):
                        print("WARNING: Summary statistics missing for " + file_prefix + ". Locus/trait will be skipped for this ancestry.")
                        locus_cmd = empty_script + ' ' + output_dir2 + file_prefix + out_file_ext
                        #Encode the command
                        locus_cmd = locus_cmd.replace(" ", "___")
                        locus_cmd = locus_cmd.replace(";", "^^^")
                        #Append the command to the list
                        locus_coms.append(locus_cmd)
                    else:
                        #This is the condition where we finally have to do the mapping
                        #Create a string to hold the command we're going to submit
                        locus_cmd = 'python3 ' + map_script + ' -p ' + file_prefix + ' -o ' + output_dir2 + ' -m ' + stats_dir + ' -l ' + matrix_dir + ' -x ' + ld_prefix + ' -n ' + file_format + ' -f ' + str(min_maf) + ' -a ' + str(swap_af_ind)
                        #Encode the command
                        locus_cmd = locus_cmd.replace(" ", "___")
                        locus_cmd = locus_cmd.replace(";", "^^^")
                        #Append the command to the list
                        locus_coms.append(locus_cmd)
    
    #Verify that list of locus_coms is non-empty (i.e. that at least one locus/trait/anc combo needs mapping)
    if len(locus_coms) > 0:
        #Add a do-nothing line at the end of the file
        locus_coms.append("echo___\"empty\"")
        #Print the list of commands to a file
        with open(command_file, "w") as out_path:
            out_path.write('\n'.join(locus_coms))  

###############################################################################


#cProfile.run('main(sys.argv[1:])')
         
if __name__ == "__main__":
    main(sys.argv[1:])

