'''
plot_CAFEH.py

This script will plot the results of fine mapping a particular locus with 
CAFEH. It will take a pickle file output from CAFEH and make a p-value plot.
If the option is selected, it will also make absolute or preceding residual
association plots.

***This needs to be run on Python 3.9 as there an incompatability with other ***
***versions of the .to_list() function  and the CAFEH package for python.    ***
'''

###############################################################################

#Load libraries
import getopt
import numpy as np
import pandas as pd
import sys
from os.path import exists
import matplotlib.pyplot as plt
import seaborn as sns
from cafeh.cafeh_summary import fit_cafeh_summary, fit_cafeh_z
from cafeh.model_queries import *
from cafeh.cafeh_summary import CAFEHSummary
import pickle

#Define help function to be called if there is a problem
def help(exit_num=1):
    print("""-----------------------------------------------------------------
ARGUMENTS
    -p => <txt> Pickled CAFEH output REQUIRED
    -o => <txt> output folder location REQUIRED
    -t => <txt> plot type: non-residual, absolute residual, or preceding residual  (default is non-residual) OPTIONAL
    -a => <txt> association type for non-residual plot threshold: gwas, absolute_residual, preceding_residual (default is gwas) OPTIONAL
    -u => <num> purity threshold for plotting non-residual plots (default is 0.5) OPTIONAL
    -m => <num> Necessary max association p-value for plotting credible sets on non-residual plots (default is 1) OPTIONAL
""")
    sys.exit(exit_num)

###############################################################################
###########################  COLLECT AND CHECK ARGUMENTS  #####################
###############################################################################

## GLOBAL VARS

## MAIN ##

def main(argv): 
    try: 
        opts, args = getopt.getopt(sys.argv[1:], "p:o:t:a:u:m:nh")
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
        pickle_file = options_dict['-p']
        out_dir = options_dict['-o']

    except KeyError:
        print("ERROR: One of your required arguments does not exist.")
        help()
    
    #Get optional inputs
    plot_type = options_dict.get('-t', 'non-residual')
    assoc_type = options_dict.get('-a', 'gwas')
    purity = options_dict.get('-u', 0.5)
    max_assoc = options_dict.get('-m', 1)
    #Confirm that the file/folder locs exist
    if exists(out_dir) == False:
        print("ERROR: Output directory does not exist.")
        sys.exit(1)
    if exists(pickle_file) == False:
        print("ERROR: Pickle file of CAFEH output doesn't exist.")
        sys.exit(1)
    #Check that the plot_type makes sense
    if type(plot_type) != str:
        print("ERROR: Plot type is invalid. Select non-residual, absolute_residual, or preceding_residual.")
        sys.exit(1)
    else:
        if plot_type.lower()[0] == 'n':
            plot_type = 'non-residual'
            print("NOTE: Non-residual plot type selected.")
        elif plot_type.lower()[0] == 'a':
            plot_type = 'absolute_residual'
            print("NOTE: Absolute residual plot type selected.")
        elif plot_type.lower()[0] == 'p':
            plot_type = 'preceding_residual'
            print("NOTE: Preceding residual plot type selected.")
        else: 
            print("ERROR: Plot type is invalid. Select non-residual, absolute_residual, or preceding_residual.")
            sys.exit(1)
    #Check that the assoc_type makes sense
    if type(plot_type) == 'non-residual' and type(assoc_type) != str:
        print("ERROR: Association type is invalid for a non-residual plot. Select gwas, absolute, or preceding.")
        sys.exit(1)
    elif type(plot_type) == 'non-residual':
        if assoc_type.lower()[0] == 'g':
            assoc_type = 'gwas'
            print("NOTE: GWAS association selected for non-residual plot.")
        elif assoc_type.lower()[0] == 'a':
            assoc_type = 'absolute'
            print("NOTE: Absolute residual association selected for non-residual plot.")
        elif assoc_type.lower()[0] == 'p':
            assoc_type = 'preceding'
            print("NOTE: Preceding residual association selected for non-residual plot.")
        else: 
            print("ERROR: Association type is invalid for a non-residual plot. Select gwas, absolute, or preceding.")
            sys.exit(1)
    #Recast purity value 
    try:
        purity = float(purity)
        if purity >= 1 or purity < 0: 
            print("ERROR: Purity level should be  in [0,1).")
            sys.exit(1)
    except ValueError:
        print("ERROR: Purity level is not coercible to a float. It should be a proportion in [0,1).")
        sys.exit(1)
    #Recast maximum association value 
    try:
        max_assoc = float(max_assoc)
        if max_assoc > 1 or purity <= 0: 
            print("ERROR: Maximum association level should be in (0,1].")
            sys.exit(1)
    except ValueError:
        print("ERROR: Maximum association level is not coercible to a float. It should be a p-value in (0,1].")
        sys.exit(1)
    print("Acceptable Inputs Given")
    #Call driver function
    driver(pickle_file, out_dir, plot_type, assoc_type, purity, max_assoc)

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

def driver(pickle_file, out_dir, plot_type, assoc_type, purity, max_assoc): 

    #Extract file prefix from pickle_file
    file_prefix = pickle_file.split("/")[-1]
    file_prefix = file_prefix[:file_prefix.rfind(".")]
    
    #Add slash to out_dir if needed
    out_dir = add_slash(out_dir)
    
    #Read in pickle file
    with open(pickle_file, 'rb') as f:
        cafehs = pickle.load(f)
    
    #Check if the pickle file is empty
    if cafehs is None:
        print("WARNING: Attempted to plot unmapped locus that encountered RAM issue, " + file_prefix + '.')
        sys.exit(0)
    
    #Check for convergence
    if cafehs.check_convergence() == False:
        print("WARNING: Plotting results for a locus that failed to converge, " + file_prefix + '.')

    #Get credible sets
    credible_sets = cafehs.credible_sets

    #Check what plot_type we have an print the corresponding plots
    #Non-residual case
    if plot_type == 'non-residual':
        active_components = cafehs.active.max(0) > 0.5
        plt.clf()
        #Set an array of SNPs that aren't in any credible sets
        remain_snps = cafehs.snp_ids
        #Iterate over credible sets and make the non-residual plot
        for k in np.arange(cafehs.dims['K'])[active_components]:
            #Get the correct assoc_type data
            if assoc_type == 'gwas':
                assoc_array = cafehs.neglogP
            elif assoc_type == 'absolute':
                assoc_array = cafehs.abs_neglogp_resid.T[:, k]
            else:
                assoc_array = cafehs.precede_neglogp_resid.T[:,k]
            #Check that the credible set passes the purity and association thresholds
            if np.max(assoc_array[credible_sets[k]]) >= -np.log10(max_assoc) and cafehs.realpure[k] >= purity: 
                #Remove credible set SNPs from remaining snps
                remain_snps = np.delete(remain_snps, np.argwhere(np.isin(remain_snps,credible_sets[k])))
                #Plot the cred set
                temp = plt.scatter(
                    cafehs.bp[credible_sets[k]],
                    cafehs.neglogP[credible_sets[k]],
                    label='k{}'.format(k))
                temp.set_zorder(cafehs.dims['K'] + 1 - k)
        rest = plt.scatter(cafehs.bp[remain_snps], cafehs.neglogP[remain_snps], color="#DEDEDE")
        rest.set_zorder(1)
        plt.title(file_prefix)
        plt.xlabel('chr' + str(np.max(cafehs.chr)))
        plt.ylabel('-log(P-value)')
        plt.ticklabel_format(axis="x", style='plain')
        plt.xticks(rotation = 45) 
        plt.legend(bbox_to_anchor=(1.04,1), loc="upper left")
        plt.savefig(f'{out_dir}{file_prefix}.purity-{purity}.{assoc_type}-{max_assoc}.pvalue.png', bbox_inches='tight')
    
    #Residual case
    else:
        #Get the correct set of residuals
        if plot_type == 'absolute_residual':
            assoc_array = cafehs.abs_neglogp_resid.T
        else:
            assoc_array = cafehs.precede_neglogp_resid.T
        #Iterate over credible sets and make the residual plots
        for i in np.arange(cafehs.dims['K']):
            plt.clf()
            logp_test = assoc_array[:,i]
            cred_set = credible_sets[i]
            rest_set = np.array([i for i in cafehs.snp_ids if i not in cred_set])
            logp_test_cred = logp_test[np.sort(cred_set)]
            logp_test_rest = logp_test[np.sort(rest_set)]
            plt.scatter(cafehs.bp[np.sort(cred_set)], logp_test_cred, color=sns.color_palette("tab10")[i])
            plt.scatter(cafehs.bp[np.sort(rest_set)], logp_test_rest, color=sns.color_palette(['#DEDEDE']))
            plt.title(f"{file_prefix}/component-{i} {plot_type}")
            plt.xlabel('SNP')
            plt.ticklabel_format(axis="x", style='plain')
            plt.xticks(rotation = 45) 
            plt.ylabel('-log(residual association)')
            plt.savefig(f'{out_dir}{file_prefix}.component-{i}.{plot_type}.png', bbox_inches='tight')
            
###############################################################################

#cProfile.run('main(sys.argv[1:])')
         
if __name__ == "__main__":
    main(sys.argv[1:])