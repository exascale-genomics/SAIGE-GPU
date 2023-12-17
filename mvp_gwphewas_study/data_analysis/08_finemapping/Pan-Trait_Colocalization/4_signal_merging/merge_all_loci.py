'''
merge_all_loci.py

The purpose of this script is to iterate over the list of loci/trait pairs that
were fine-mapped.

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
    -m => <str> Output location for command file REQUIRED
    -o => <str> Output directory for merged fine-mapping results (default is /gpfs/alpine/med112/proj-shared/results/PostGWAS_Analysis/finemaping/GIA/susie_dir_MERGE) OPTIONAL
    -a => <str> Directory for AFR fine-mapping results (default is /gpfs/alpine/med112/proj-shared/results/PostGWAS_Analysis/finemaping/GIA/susie_dir_AFR) OPTIONAL
    -b => <str> Directory for AMR fine-mapping results (default is /gpfs/alpine/med112/proj-shared/results/PostGWAS_Analysis/finemaping/GIA/susie_dir_AMR) OPTIONAL
    -c => <str> Directory for EAS fine-mapping results (default is /gpfs/alpine/med112/proj-shared/results/PostGWAS_Analysis/finemaping/GIA/susie_dir_EAS) OPTIONAL
    -d => <str> Directory for EUR fine-mapping results (default is /gpfs/alpine/med112/proj-shared/results/PostGWAS_Analysis/finemaping/GIA/susie_dir_EUR) OPTIONAL
    -t => <str> Significance type to use for thresholding (Options are gwas/absolute_residual/preceding_residual; default is gwas) OPTIONAL
    -s => <num> Significance threshold (default is 5e-8) OPTIONAL
    -p => <num> Purity threshold (default is 0.1) OPTIONAL
    -u => <num> Cutting height for Jaccard distance dendrogram OPTIONAL (default is 0.9; which is equivalent to merging all credible sets with Jaccard indices >= 0.1)
    -n => <str> SuSiE Results file naming convention showing where ANCestry, TRAIT, CHROmosome, START, and END are found relative to periods 
                (default is TRAIT.CHR.START.END.ANC.rds) OPTIONAL
    -r => <num> Random seed (default is 5) OPTIONAL
    -e => <boo> Force merge partially mapped locus/trait pairs (True/False; default is False) OPTIONAL
""")
    sys.exit(exit_num)

###############################################################################
###########################  COLLECT AND CHECK ARGUMENTS  #####################
###############################################################################

## GLOBAL VARS

## MAIN ##

def main(argv): 
    try: 
        opts, args = getopt.getopt(sys.argv[1:], "l:m:o:a:b:c:d:t:s:p:u:n:r:e:nh")
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
        command_file = options_dict['-m']        
    except KeyError:
        print("ERROR: One of your required arguments does not exist.")
        help()
    
    #Optional arguments
    output_dir = options_dict.get('-o', "/gpfs/alpine/med112/proj-shared/results/PostGWAS_Analysis/finemaping/GIA/susie_dir_MERGE")
    AFR_dir = options_dict.get('-a', "/gpfs/alpine/med112/proj-shared/results/PostGWAS_Analysis/finemaping/GIA/susie_dir_AFR")
    AMR_dir = options_dict.get('-b', "/gpfs/alpine/med112/proj-shared/results/PostGWAS_Analysis/finemaping/GIA/susie_dir_AMR")
    EAS_dir = options_dict.get('-c', "/gpfs/alpine/med112/proj-shared/results/PostGWAS_Analysis/finemaping/GIA/susie_dir_EAS")
    EUR_dir = options_dict.get('-d', "/gpfs/alpine/med112/proj-shared/results/PostGWAS_Analysis/finemaping/GIA/susie_dir_EUR")
    sig_type = options_dict.get('-t', 'gwas')
    sig_thresh = options_dict.get('-s', 5e-8)
    purity_thresh = options_dict.get('-p', 0.1)
    cutting_height = options_dict.get('-u', 0.9)
    name_conv = options_dict.get('-n', 'TRAIT.CHR.START.END.ANC.rds')
    random_seed = options_dict.get('-r', 5)
    empty_ind = options_dict.get('-e', False)
    
    #Recast empty_ind
    if(type(empty_ind)) == str:
        if(empty_ind.lower() == 'false' or empty_ind.lower() == 'f' or empty_ind.lower() == 'no' or empty_ind.lower() == 'n'):
            empty_ind = False
        elif(empty_ind.lower() == 'true' or empty_ind.lower() == 't' or empty_ind.lower() == 'yes' or empty_ind.lower() == 'y'):
            empty_ind = True
        else:
            print("ERROR: Unrecognized option selected for whether to force merge partially mapped loci-trait pairs. Select True/False.")
            sys.exit(1)
    elif(type(empty_ind) == bool):
        pass #In this case the indicator is already a boolean, so there's no need to recast.
    else:
        print("ERROR: Unrecognized option selected for whether to force merge partially mapped loci-trait pairs. Select True/False.")
        sys.exit(1)
    #Confirm that output_type is a valid option
    if type(sig_type) == str:
        sig_type = sig_type.lower()
        if sig_type not in ['gwas', 'absolute_residual', 'preceding_residual']:
            print("ERROR: Significance Type must be gwas, absolute_residual, preceding_residual.")
            sys.exit(1)
    else:
        print("ERROR: Mapping Type must be a string equivalent to gwas, absolue_residual, or preceding_residual.")
        sys.exit(1)
    #Confirm that the file locs exist
    if exists(locus_file) == False:
        print("ERROR: Input file of merged significant loci does not exist.")
        sys.exit(1)
    #Check that directories exist
    if exists(output_dir) == False:
        print("ERROR: Output directory for merged fine-mapping results does not exist.")
        sys.exit(1)
    if exists(AFR_dir) == False:
        print("ERROR: AFR directory does not exist.")
        sys.exit(1)
    if exists(AMR_dir) == False:
        print("ERROR: AMR directory does not exist.")
        sys.exit(1)
    if exists(EAS_dir) == False:
        print("ERROR: EAS directory does not exist.")
        sys.exit(1)
    if exists(EUR_dir) == False:
        print("ERROR: EUR directory does not exist.")
        sys.exit(1)
    if exists(os.path.dirname(command_file)) == False:
        print("ERROR: Invalid location specified for command file.")
        sys.exit(1)
    #Verify that ANC and TRAIT are present in naming convention
    if 'ANC' not in name_conv:
        print("ERROR: Unable to identify position of ANCestry in file name convention.")
        sys.exit(1)
    elif 'TRAIT' not in name_conv:
        print("ERROR: Unable to identify position of TRAIT in file name convention.")
        sys.exit(1)
    #Recast significance/purity thresholds and cutting height
    try:
        sig_thresh = float(sig_thresh)
        if sig_thresh > 1 or sig_thresh <= 0: 
            print("ERROR: Signficance threshold for credible sets should be a proportion in (0,1].")
            sys.exit(1)
    except ValueError:
        print("ERROR: Signficance threshold for credible sets should be a float in (0,1].")
        sys.exit(1)
    try:
        purity_thresh = float(purity_thresh)
        if purity_thresh >= 1 or purity_thresh < 0: 
            print("ERROR: Purity threshold for credible sets should be a proportion in [0,1).")
            sys.exit(1)
    except ValueError:
        print("ERROR: Purity threshold for credible sets should be a float in [0,1).")
        sys.exit(1)
    try:
        cutting_height = float(cutting_height)
        if cutting_height > 1 or cutting_height < 0: 
            print("ERROR: Cutting height for merging should be a float in (0,1).")
            sys.exit(1)
    except ValueError:
        print("ERROR: Cutting height for merging should be a float in (0,1).")
        sys.exit(1)
    #Recast random_seed 
    try:
        random_seed = int(random_seed)
    except ValueError:
        print("ERROR: Random seed must be an integer.")
        sys.exit(1)
    print("Acceptable Inputs Given")
    #Call driver function
    driver(locus_file, command_file, output_dir, AFR_dir, AMR_dir, EAS_dir, EUR_dir, sig_type, sig_thresh, purity_thresh, cutting_height, name_conv, random_seed, empty_ind)

###############################################################################
#############################  HELPFUL FUNCTIONS ##############################
###############################################################################

#Function that adds slashes to directory names
def add_slash(directory):
    if directory[-1] != '/':
        return directory + '/'
    else:
        return directory

#Function that lists SuSiE results for a given ancestry, checks matching format against the naming convention, and replaces the ancestry with ANC
def get_anc_files(ANC_dir, name_conv, ancestry):
    name_conv = name_conv.split(".")
    anc_pos = name_conv.index("ANC")
    anc_files = os.listdir(ANC_dir)
    anc_files = [x.split(".") for x in anc_files]
    anc_match_files = [x for x in anc_files if x[anc_pos] == ancestry and len(x) == len(name_conv) and x[-1] == name_conv[-1]]
    anc_comb_match_files = ['.'.join(x) for x in anc_match_files]
    anc_comb_match_files = [x.replace("." + ancestry + ".", ".ANC.") for x in anc_comb_match_files]
    return anc_comb_match_files

###############################################################################
#############################  DRIVER  ########################################
###############################################################################

## drive the script ##
## ONE-TIME CALL -- called by main
def driver(locus_file, command_file, output_dir, AFR_dir, AMR_dir, EAS_dir, EUR_dir, sig_type, sig_thresh, purity_thresh, cutting_height, name_conv, random_seed, empty_ind):
    
    #Get location of locus merging script and verify that it exists
    temp = os.path.realpath(__file__)
    merge_script = temp[:temp.find('merge_all_loci.py')] + 'merge_locus.R'
    if exists(merge_script) == False:
        print("ERROR: merge_locus.R merging script not found in same folder as merge_all_loci.py. Consider repulling github.")
        sys.exit(1)

    #Create a blank list that will hold the commands
    locus_coms = []
    
    #Add final slash to directories if they're missing
    output_dir = add_slash(output_dir)
    AFR_dir = add_slash(AFR_dir)
    AMR_dir = add_slash(AMR_dir)
    EAS_dir = add_slash(EAS_dir)
    EUR_dir = add_slash(EUR_dir)
    
    #Read in the list of loci
    locus_raw = []  
    with open(locus_file) as locus_read:
        for line in locus_read.readlines():
            temp = line.strip().split(",") 
            locus_name = name_conv.replace("CHR", 'chr' + temp[0]).replace("START", temp[1]).replace("END", temp[2])
            traits = temp[4].split("/")
            for trait in traits:
                locus_raw.append(locus_name.replace("TRAIT", trait))
    
    #Create a list of all the rds files in the output directories
    ancestry_dirs = {'AFR':AFR_dir, 'AMR':AMR_dir, 'EAS':EAS_dir, 'EUR':EUR_dir}
    ancestry_files = {}
    for ancestry in ancestry_dirs.keys():
        ancestry_files[ancestry] = get_anc_files(ancestry_dirs[ancestry], name_conv, ancestry)
    
    #Get intersection and union of lists
    intersected_files = [x for x in locus_raw if x in ancestry_files['AFR'] and x in ancestry_files['AMR'] and x in ancestry_files['EAS'] and x in ancestry_files['EUR']]
    union_files = list(set(locus_raw) | set(ancestry_files['AFR']) | set(ancestry_files['AMR']) | set(ancestry_files['EAS']) | set(ancestry_files['EUR']))
    missing_files = list(set([x for x in union_files if x not in intersected_files]))
    
    #Loop over the missing files and create ancestry-specific lists of missing files
    missing_ancestry_files = {}
    partially_missing_files = []
    for file in missing_files:
        if file not in locus_raw: #Case where an rds file was made even though it shouldn't have been - this theoretically shouldn't happen
            print("WARNING: RDS file(s) exist for " + file.replace(".rds", "") + " even though this locus/trait combo is not in the locus file")
            continue
        missing_ancestry_files[file] = {}
        for ancestry in ancestry_dirs.keys(): # Add ancestry specific missingness info to the dictionary
            if file not in ancestry_files[ancestry]:
                missing_ancestry_files[file][ancestry] = False
            else:
                missing_ancestry_files[file][ancestry] = True
        if sum(missing_ancestry_files[file].values()) == 0:
            print("WARNING: " + file.replace(".ANC.rds", "") + " is unmapped in all ancestries.")
        else: 
            if empty_ind == False:
                print("WARNING: " + file.replace(".ANC.rds", "") + " is mapped in some ancestries. Files will remain unmerged")
            else:
                print("WARNING: " + file.replace(".ANC.rds", "") + " is mapped in some ancestries. Force merging completed files.")
                partially_missing_files.append(file)
    
    #Loop over the intersected files 
    for trait_locus in intersected_files:
        #Check if the output file already exists for the locus
        if exists(output_dir + trait_locus.replace('.ANC.', '.MERGE.').replace('.rds', '.txt')) == False:
            locus_cmd = 'Rscript ' + merge_script + ' '
            #Add on ancestry-specific results files
            for ancestry in ancestry_dirs.keys(): #This will only work because the ancestries are in alphabetical order for inputs in the Rscript
                locus_cmd = locus_cmd + ancestry_dirs[ancestry] + trait_locus.replace('.ANC.', '.' + ancestry + '.') + ' '
            #Add on output dir and the thresholds
            locus_cmd = locus_cmd + output_dir + ' ' + str(sig_thresh) + ' ' + sig_type + ' ' + str(purity_thresh) + ' ' + str(cutting_height) + ' ' + name_conv + ' ' + str(random_seed)
            #Encode the command
            locus_cmd = locus_cmd.replace(" ", "___")
            locus_cmd = locus_cmd.replace(";", "^^^")
            #Append the command to the list
            locus_coms.append(locus_cmd)
        else:
            #print a note that the locus already exists
            print('NOTE: Merging already complete for ' + trait_locus.replace('.ANC.', '').replace('.rds', ''))
    
    #If empty_ind is set to true also force merge the partially_missing_files
    if empty_ind == True:
        for file in partially_missing_files:
            #Similar to intersected_files except don't check for file existing. Always force remerge these in case new files were created.            
            #Start command
            locus_cmd = 'Rscript ' + merge_script + ' '
            #Add on ancestry-specific results files
            for ancestry in ancestry_dirs.keys(): #This will only work because the ancestries are in alphabetical order for inputs in the Rscript
                if missing_ancestry_files[file][ancestry] == True: #Case where the ancestry has been mapped
                    locus_cmd = locus_cmd + ancestry_dirs[ancestry] + file.replace('.ANC.', '.' + ancestry + '.') + ' '
                else: #Case where the ancestry hasn't been mapped
                    locus_cmd = locus_cmd + "NA" + ' '
            #Add on output dir and the thresholds
            locus_cmd = locus_cmd + output_dir + ' ' + str(sig_thresh) + ' ' + sig_type + ' ' + str(purity_thresh) + ' ' + str(cutting_height) + ' ' + name_conv + ' ' + str(random_seed)
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
    else:
        print("NOTE: All files with complete rds files have already been merged.")
    

###############################################################################


#cProfile.run('main(sys.argv[1:])')
         
if __name__ == "__main__":
    main(sys.argv[1:])

