'''
make_matrices.py

The purpose of this script is to be a wrapper that determines for each locus
which plink files are needed to map the locus. It submits separate srun 
commands to combine the plink files and then ultimately makes map 
files, or both map files and matrices depending upon the input parameters

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
    -l => <txt> Input file of trait-combined significant loci REQUIRED
    -o => <txt> Folder in which to output matrix/map files REQUIRED
    -a => <txt> ancestry (options are AFR/AMR/EAS/EUR) REQUIRED
    -i => <txt> Directory of files specifying individuals for each ancestry/phenotype pair REQUIRED
    -r => <txt> Registry file that lists out the ranges found in each plink segment REQUIRED
    -c => <txt> File of commands to be submitted to the console REQUIRED
    -t => <txt> Output type (options are Matrix/Map/Both; default is Both) Optional
    -d => <Num> Maximum distance in bp to make matrices for (default is Inf) Optional
    -e => <txt> Variant exclusion file (default is variant_exclude_list.GLOBAL.txt in this script's directory) OPTIONAL
    -p => <txt> Directory that houses plink files (default is /gpfs/alpine/med112/proj-shared/task0101113/output/pheCodes/inputs/bgen) OPTIONAL
    -f => <txt> Ancestry/Phenotype file naming convention showing where ANCestry and PHENO are found relative to periods 
                (default is PHENO.ANC.individuals.txt) OPTIONAL
""")
    sys.exit(exit_num)

###############################################################################
###########################  COLLECT AND CHECK ARGUMENTS  #####################
###############################################################################

## GLOBAL VARS

## MAIN ##

def main(argv): 
    try: 
        opts, args = getopt.getopt(sys.argv[1:], "l:o:a:i:r:c:t:d:e:p:f:nh")
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
        output_dir = options_dict['-o']
        ancestry = options_dict['-a']
        individual_dir = options_dict['-i']
        register_file = options_dict['-r']
        command_file = options_dict['-c']
    except KeyError:
        print("ERROR: One of your required arguments does not exist.")
        help()
    
    #Optional arguments
    output_type = options_dict.get('-t', "BOTH")
    max_dist = options_dict.get('-d', 1e300)
    plink_dir = options_dict.get('-p', "/gpfs/alpine/med112/proj-shared/task0101113/output/pheCodes/inputs/bgen")
    name_conv = options_dict.get('-f', 'PHENO.ANC.individuals.txt')
    temp = os.path.realpath(__file__)
    exclude_file = options_dict.get('-e', temp[:temp.find('make_matrices.py')] + 'variant_exclude_list.GLOBAL.txt')
    #Confirm that output_type is a valid option
    if type(output_type) == str:
        output_type = output_type.upper()
        if output_type not in ['MATRIX', 'MAP', 'BOTH']:
            print("ERROR: Output Type must be Matrix, Map, or Both.")
            sys.exit(1)
    else:
        print("ERROR: Output Type must be a string equivalent to Matrix/Map/Both.")
        sys.exit(1)
    #Confirm that ancestry is a valid option
    if type(ancestry) == str:
        ancestry = ancestry.upper()
        if ancestry not in ['AFR', 'AMR', 'EAS', 'EUR']:
            print("ERROR: Ancestry must be AFR, AMR, EAS, or EUR.")
            sys.exit(1)
    else:
        print("ERROR: Ancestry must be a string equivalent to AFR/AMR/EAS/EUR.")
        sys.exit(1)
    #Confirm that the file locs exist
    if exists(locus_file) == False:
        print("ERROR: Input file of merged significant loci does not exist.")
        sys.exit(1)
    if exists(individual_dir) == False:
        print("ERROR: Input directory of files listing individuals by ancestry/phenotype does not exist.")
        sys.exit(1)
    if exists(register_file) == False:
        print("ERROR: Input file detailing breakdown of plink file ranges does not exist.")
        sys.exit(1)
    if exists(exclude_file) == False:
        print("ERROR: Input file of variant ids to exclude does not exist.")
        sys.exit(1)
    #Check that directories exist
    if exists(output_dir) == False:
        print("ERROR: Output directory does not exist.")
        sys.exit(1)
    if exists(plink_dir) == False:
        print("ERROR: Invalid directory given for plink files' location.")
        sys.exit(1)
    #Confirm that the maximum distance is a positive value
    try:
        max_dist = int(max_dist)
    except ValueError:
        print("ERROR: Maximum bp distance must be a positive integer.")
        sys.exit(1)
    if max_dist <= 0:
        print("ERROR: Maximum bp distance must be positive.")
        sys.exit(1)
    #Verify that ANC is present in naming convention
    if 'ANC' not in name_conv:
        print("ERROR: Unable to identify position of ANCestry in file name convention.")
        sys.exit(1)
    print("Acceptable Inputs Given")
    #Call driver function
    driver(locus_file, output_dir, output_type, max_dist, ancestry, individual_dir, register_file, command_file, exclude_file, plink_dir, name_conv)

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
def driver(locus_file, output_dir, output_type, max_dist, ancestry, individual_dir, register_file, command_file, exclude_file, plink_dir, name_conv):
    
    #Set plink directories
    plink2 = "/gpfs/alpine/proj-shared/med112/GIA/bin/plink2"
    plink1 = "/gpfs/alpine/proj-shared/med112/GIA/bin/plink"
        
    #Add slashes to directories
    output_dir = add_slash(output_dir)
    plink_dir = add_slash(plink_dir)
    individual_dir = add_slash(individual_dir)
    #Get list of all plink files in the directory
    plink_files = os.listdir(plink_dir)
    plink_files = [file.replace(".pgen","") for file in plink_files if "pgen" == file.split(".")[-1]]
    plink_files = [file for file in plink_files if "dose" == file.split(".")[-1]] #Filter for the dose files
    
    #Read in the registry file
    register_raw = pd.read_csv(register_file, sep="\t")
    #Add dose to the file names in the register_raw
    register_raw["Chunk"] = register_raw["Chunk"].astype(str) + '.dose'
    #Read in the locus file
    locus_raw = []  
    with open(locus_file) as locus_read:
        for line in locus_read.readlines():
            temp = line.strip().split(",")
            temp[3] = int(temp[3])
            if temp[3] <= max_dist:
                temp[0] = int(temp[0])
                temp[1] = int(temp[1])
                temp[2] = int(temp[2])
                temp[4] = temp[4].split("/")
                locus_raw.append(temp)
    
    #Create a list of all the matrix/map files in the output directory
    out_files = os.listdir(output_dir) 
    map_files = [x.replace('.map','') for x in out_files if x.split(".")[-1] == 'map']
    matrix_files = [x.replace('.ld','') for x in out_files if x.split(".")[-1] == 'ld']
    #Check if output_type is 'BOTH' and if so, filter the lists of matrix and map files for only the intersection
    if output_type == 'BOTH':
        map_files = [x for x in map_files if x in matrix_files]
        matrix_files = [x for x in matrix_files if x in map_files]
    matrix_files = [x.split(".") for x in matrix_files]
    map_files = [x.split(".") for x in map_files]
    
    #Create a blank list that will hold the commands
    locus_coms = []
    #Loop over the locus_raw list and start making matrices
    for locus in locus_raw:
        #Loop over the traits for the locus
        for trait in locus[4]:
            #Get ancestry_pheno_file
            ancestry_pheno_file = individual_dir + name_conv.replace('ANC', ancestry).replace('PHENO', trait)
            #Reset execute indicator
            execute_ind = False
            #Create a blank string to hold the command we're going to submit
            locus_cmd = ''
            #Get locus prefix
            file_prefix = trait + '.chr' + str(locus[0]) + '.' + str(locus[1]) + '.' + str(locus[2]) + '.' + ancestry
            #Identify the plink files that are needed for the locus 
            temp1 = register_raw[register_raw['Chr'] == locus[0]]
            temp2 = temp1[temp1['Start'] <= locus[1]]
            if temp2.empty: #Handle chromosome start issues
                temp2 = temp1.iloc[0,:]
            else:
                temp2 = temp2.loc[max(temp2.index),]
            temp3 = temp1[temp1['End'] >= locus[2]]
            if temp3.empty: #Handle chromosome end issues
                temp3 = temp1.iloc[len(temp1)-1,:]
            else:
                temp3 = temp3.loc[min(temp3.index),]
            #Create indices
            if temp2['Start'] <= locus[1] and temp2['End'] >= locus[1]:
                lower_index = temp2.name
            elif temp2['Start'] > locus[1]: #This will only be possible in the chromosome start case
                lower_index = temp2.name
            else:
                lower_index = temp2.name + 1
            if temp3['Start'] <= locus[2] and temp3['End'] >= locus[2]:
                upper_index = temp3.name + 1
            elif  temp3['End'] < locus[2]: #This will only be possible in the chromosome end case
                upper_index = temp3.name + 1
            else:
                upper_index = temp3.name
            plink_count = upper_index - lower_index
            #Print update message and append to locus_cmd
            locus_cmd = locus_cmd + 'echo \"START: ' + file_prefix + ' `date`\";'
            locus_cmd = locus_cmd + 'echo \"HOST: ' + file_prefix + ' \$HOSTNAME\";'
            if plink_count > 1:
                print('NOTE: Merging ' + str(plink_count) + ' plink files ranging from ' + 
                      register_raw.loc[lower_index, 'Chunk'] + ' to ' + register_raw.loc[upper_index-1, 'Chunk'] + 
                      ' for ' + file_prefix)
                #Export files to a merge list
                temp_chunks=register_raw.loc[range(lower_index, upper_index), 'Chunk']
                temp_chunks=plink_dir + temp_chunks
                temp_chunks.to_csv(plink_dir + file_prefix + '.temp.merge_list.txt', sep = "\t", header = False, index = False, na_rep='NA')
                #Add to the command
                locus_cmd = locus_cmd + plink2 + ' --pmerge-list ' + plink_dir + file_prefix + '.temp.merge_list.txt --exclude ' + exclude_file + ' --keep ' + ancestry_pheno_file + ' --chr ' + str(locus[0]) + ' --from-bp ' + str(locus[1]) + ' --to-bp ' + str(locus[2]) + ' --make-bed --out ' + plink_dir + file_prefix + '.temp --memory 200000;'
            else: 
                print('NOTE: Merging unnecessary for ' + file_prefix)
                #Make plink1.9 compatible files
                locus_cmd  = locus_cmd + plink2 + ' --pfile ' + plink_dir + str(register_raw.loc[lower_index, 'Chunk']) + ' --exclude ' + exclude_file + ' --keep ' + ancestry_pheno_file + ' --chr ' + str(locus[0]) + ' --from-bp ' + str(locus[1]) + ' --to-bp ' + str(locus[2]) + ' --make-bed --out ' + plink_dir + file_prefix + '.temp --memory 200000;'
            #Process locus according to desired output
            if output_type == 'MAP' or output_type == 'BOTH':
                #Check map_files for matching files
                locus_map_files = [x for x in map_files if trait == x[0] and 'chr' + str(locus[0]) == x[1] and locus[1] >= int(x[2]) and locus[2] <= int(x[3]) and ancestry == x[4]]
                #Check length of matching map files
                if len(locus_map_files) == 1:
                    print("NOTE: Map file " + locus_map_files[0][0] + "." + locus_map_files[0][1] + "." + locus_map_files[0][2] + "." + locus_map_files[0][3] + " already created for " + file_prefix +'. Skipping locus/ancestry pair.')
                elif len(locus_map_files) > 1:
                    print("NOTE: Multiple useable map files exist for " + file_prefix +'. Skipping locus/ancestry pair.')
                else:
                    locus_cmd  = locus_cmd + ' ' + plink1 + ' --bfile ' + plink_dir + file_prefix + '.temp --recode --out ' + output_dir + file_prefix + ' --memory 200000;'
                    execute_ind = True
            if output_type == 'MATRIX' or output_type == 'BOTH':
                #Check matrix_files for matching files
                locus_matrix_files = [x for x in matrix_files if trait == x[0] and 'chr' + str(locus[0]) == x[1] and locus[1] >= int(x[2]) and locus[2] <= int(x[3]) and ancestry == x[4]]
                #Check length of matching matrix files
                if len(locus_matrix_files) == 1:
                    print("NOTE: Matrix file " + locus_matrix_files[0][0] + "." + locus_matrix_files[0][1] + "." + locus_matrix_files[0][2] + "." + locus_matrix_files[0][3] + " already created for " + file_prefix +'. Skipping locus/ancestry pair.')
                elif len(locus_matrix_files) > 1:
                    print("NOTE: Multiple useable matrix files exist for " + file_prefix +'. Skipping locus/ancestry pair.')
                else:
                    locus_cmd  = locus_cmd + ' ' + plink1 + ' --bfile ' + plink_dir + file_prefix + '.temp --r square --keep-allele-order --out ' + output_dir + file_prefix + ' --memory 200000;'
                    execute_ind = True
            #Delete the temporary files
            locus_cmd = locus_cmd + 'rm ' + plink_dir + file_prefix + '.temp.*;'
            locus_cmd = locus_cmd + 'rm ' + output_dir + file_prefix + '.log;'
            locus_cmd = locus_cmd + 'rm ' + output_dir + file_prefix + '.nosex;'
            locus_cmd = locus_cmd + 'rm ' + output_dir + file_prefix + '.ped;'
            locus_cmd = locus_cmd + 'echo \"END: ' + file_prefix +  ' `date`\"'
            #Replace all the spaces and semi-colons for bash interpretation
            locus_cmd = locus_cmd.replace(" ", "___")
            locus_cmd = locus_cmd.replace(";", "^^^")
            #Append the commands to the locus_coms list
            if execute_ind == True:
                locus_coms.append(locus_cmd)
    
    #Add a do-nothing line at the end of the file
    locus_coms.append("echo___\"empty\"")
    #Print the list of commands to a file
    with open(command_file, "w") as out_path:
        out_path.write('\n'.join(locus_coms))
    

###############################################################################


#cProfile.run('main(sys.argv[1:])')
         
if __name__ == "__main__":
    main(sys.argv[1:])

