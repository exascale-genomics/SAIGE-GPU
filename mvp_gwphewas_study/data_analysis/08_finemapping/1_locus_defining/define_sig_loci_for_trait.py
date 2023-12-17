'''
define_sig_loci_for_trait.py

The purpose of this script is to read over the GWAS results for a trait and 
define the significant loci bounds. Significant loci are defined as ranges that 
are multiples of a fixed distance tile size where each consecutive tile of 
the stagger size in the middle of the interval contains one or more SNP 
passing a set secondardy p-value threshold and overall there is at least one 
variant in the locus passing a stricter primary threshold. Each tailing 
staggered tile contains no signficant variants. 

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
    -g => <txt> GWAS summary statistics location for selected phenotype REQUIRED
    -o => <txt> output file location  REQUIRED
    -p => <flt> primary p-value significance threshold (default is 4.6e-11) OPTIONAL
    -s => <flt> secondary p-value significance threshold (default is 5e-8) OPTIONAL
    -t => <int> tile size (default is 250kb) OPTIONAL
    -c => <str> chromosome ranges file (default is grch38_chr_bounds_mhc_only.bed in this script's directory) OPTIONAL
    -n => <str> summary stats file naming convention showing where ANCestry and TRAIT are found relative to periods 
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
        opts, args = getopt.getopt(sys.argv[1:], "g:o:p:s:t:c:n:nh")
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
        gwas_file = options_dict['-g']
        output_loc = options_dict['-o']

    except KeyError:
        print("ERROR: One of your required arguments does not exist.")
        help()
    
    prime_p_thresh = options_dict.get('-p', 4.6e-11)
    sec_p_thresh = options_dict.get('-s', 5e-8)
    stagger_bp = options_dict.get('-t', 250000)
    temp = os.path.realpath(__file__)
    chromo_lengths = options_dict.get('-c', temp[:temp.find('define_sig_loci_for_trait.py')] + 'grch37_chr_bounds_mhc_only.bed')
    file_format = options_dict.get('-n', 'TRAIT.ANC.GIA.KDI.txt.gz')
    #Confirm that integers were given for optional parameters
    try:
        stagger_bp = int(stagger_bp)
        prime_p_thresh = float(prime_p_thresh)
        sec_p_thresh = float(sec_p_thresh)
    except ValueError:
        print("ERROR: Stagger amount must be a positive integer and p-value thresholds must be floats.")
        sys.exit(1)
    #Verify that the stagger amount is a positive integer
    if stagger_bp <= 0:
        print("ERROR: Stagger amount must be positive.")
        sys.exit(1)
    #Verify that the p-value thresholds are between 0 and 1
    if prime_p_thresh <= 0 or prime_p_thresh > 1 or sec_p_thresh <= 0 or sec_p_thresh > 1 or prime_p_thresh > sec_p_thresh:
        print("ERROR: p-value thresholds must be in (0,1] and primary threshold must be <= secondary threshold.")
        sys.exit(1)
    #Confirm that the file locs exist
    if exists(chromo_lengths) == False:
        print("ERROR: Input chromosome lengths file does not exist.")
        sys.exit(1)
    if exists(gwas_file) == False:
        print("ERROR: Input GWAS results file does not exist.")
        sys.exit(1)
    #Verify that ANC and TRAIT are present in file format
    if 'ANC' not in file_format or 'TRAIT' not in file_format:
        print("ERROR: Unable to identify position of ANCestry or TRAIT in file name format.")
        sys.exit(1)
    print("Acceptable Inputs Given")
    #Call driver function
    driver(gwas_file, output_loc, prime_p_thresh, sec_p_thresh, stagger_bp, chromo_lengths, file_format)

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
def driver(gwas_file, output_loc, prime_p_thresh, sec_p_thresh, stagger_bp, chromo_lengths, file_format):

    #Identify ancestry and trait position from file_format
    temp = file_format.split(".")
    trait_pos = temp.index("TRAIT")
    
    #Make a list to hold the output data
    out_ranges = []

    #Get trait name
    trait = gwas_file.split("/")[-1].split(".")[trait_pos]
    cur_file = pd.read_csv(gwas_file, index_col=0, sep="\t")
    print(trait)
    if cur_file.shape[1] <= 1:
        cur_file = pd.read_csv(gwas_file, index_col=0, sep=" ")
    #Filter for the SNPs below the secondary threshold
    filt_file = cur_file[cur_file['pval'] <= sec_p_thresh] #May need to modify row if file changes
    #Figure out which multiple of the stagger_bp each SNP is in
    temp = filt_file.index.str.split(":")
    sig_tiles = [int(x[0]) + (int(x[1])//stagger_bp)/1000 for x in temp]
    #Sort the list of significant tiles and remove duplicates
    sig_tiles_no_dups = list(set(sig_tiles))
    sig_tiles_no_dups = np.sort(sig_tiles_no_dups)
    #Check to make sure there's at least one significant tile
    if len(sig_tiles_no_dups) == 0:
        print("NOTE: No significant loci identified for " + trait + ".")
        sys.exit(0)
    #Loop over the significant tiles and append information to the final_ranges
    cur_locus = []
    for i in range(len(sig_tiles_no_dups)):
        cur_tile = sig_tiles_no_dups[i]
        chromo = int(cur_tile//1)
        start_bp = int(round(1000*(cur_tile%1)) * stagger_bp) #This gets the start position of the current tile (ending in 0)
        if cur_locus == []: #Handle the first interval
            cur_locus = [chromo, start_bp+1, start_bp+stagger_bp]
        elif start_bp == cur_locus[2] and chromo == cur_locus[0]: #Extend the current locus
            cur_locus[2] = cur_locus[2] + stagger_bp
        else: #We need to start a new locus in this case
            #Add a length to the cur_locus
            cur_locus.append(cur_locus[2] - cur_locus[1] + 1)
            cur_locus.append(trait)
            out_ranges.append(cur_locus) #Append the current locus to the output list
            cur_locus = [chromo, start_bp+1, start_bp+stagger_bp]
    #Handle the final sig tile
    #Add a length to the cur_locus
    cur_locus.append(cur_locus[2] - cur_locus[1] + 1)
    cur_locus.append(trait)
    out_ranges.append(cur_locus)
    
    #Filter down the filt_file again for the variants significant at the primary threshold
    temp = filt_file[filt_file['pval'] <= prime_p_thresh].index.str.split(":")
    prime_sig_vars = [int(x[0]) + float(x[1])/1e9 for x in temp]
    
    #Read in chromosome lengths file
    chr_lens = [] #This file should already be ordered by chromosome 
    with open(chromo_lengths) as chromo_len_read:
        for line in chromo_len_read.readlines():
            temp = line.strip().split("\t") 
            temp[0] = temp[0][3:] #This will remove the "chr"s in the bed file
            chr_lens.append(temp)
    #Make a dictionary out of the chromosome lengths
    chr_lens_dict = {}
    for i in range(len(chr_lens)):
        if int(chr_lens[i][0]) in chr_lens_dict.keys():
            chr_lens_dict[int(chr_lens[i][0])].append([int(chr_lens[i][1]), int(chr_lens[i][2])])
        else: 
            chr_lens_dict[int(chr_lens[i][0])] = [[int(chr_lens[i][1]), int(chr_lens[i][2])]]

    #Filter for loci containing a primary significan variant and adjust range boundaries where needed
    out_ranges = sorted(out_ranges, key = lambda x: (x[0], x[1], x[2]))
    final_ranges = ["temp"]
    for i in range(len(out_ranges)):
        #Adjust the range boundaries if need be
        cur_chromo = out_ranges[i][0]
        for cur_chromo_range in chr_lens_dict[cur_chromo]:
            if (out_ranges[i][1] >= cur_chromo_range[0] and out_ranges[i][1] < cur_chromo_range[1]) or (out_ranges[i][2] > cur_chromo_range[0] and out_ranges[i][2] <= cur_chromo_range[1]):
                range_start = max(out_ranges[i][1] - stagger_bp, cur_chromo_range[0]) #Added padding equal to the stagger bp
                range_end = min(out_ranges[i][2] + stagger_bp, cur_chromo_range[1]) #Added padding equal to the stagger bp
                range_len = range_end - range_start + 1
                range_temp = [cur_chromo, range_start, range_end]
                #Verify that the locus contains at least one primary significant variant before appending
                if len([x for x in prime_sig_vars if x >= range_temp[0] + float(range_temp[1])/1e9 and x <= range_temp[0] + float(range_temp[2])/1e9]) > 0:
                    #Append new locus to output file list
                    range_temp.extend([range_len, out_ranges[i][4]])
                    final_ranges.append(range_temp)
    del final_ranges[0]
    #Output message if no rows are left in final_ranges
    if len(final_ranges) == 0:
        print("NOTE: No significant loci identified for " + trait + ".")
        sys.exit(0)
    
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
