'''
execute_command_file.py

The purpose of this script is to decode and execute the commands in an input 
command file. It takes in a single command file parameter and executes the 
steps listed in that file in sequence. It can accept multiple commands per line
separated by ";"s as well.

***This needs to be run on Python/3.9 as there's an incompatability with ***
***other versions of python.                                             ***
'''

#import modules that will be needed
import sys #Taking command line parameters as inputs
import os
from os.path import exists

###############################################################################
###########################  COLLECT AND CHECK ARGUMENTS  #####################
###############################################################################

## GLOBAL VARS

## MAIN ##

def main(argv): 

    ## Required arguments
    try:
        command_file = sys.argv[1]
    except KeyError:
        print("ERROR: One of your required arguments does not exist.")
        help()
    #Confirm that the file locs exist
    if exists(command_file) == False:
        print("ERROR: Input file of commands does not exist.")
        sys.exit(1)
    print("Acceptable Inputs Given")
    #Call driver function
    driver(command_file)

###############################################################################
#############################  TEST LOCATIONS ##############################
###############################################################################

command_file = "C:/Users/mitch/Documents/UPenn/Voight_Lab/MVP_Finemapping/testing/downsampling_commands/A1C_Max_INT.downsampling_variant_list_commands.txt"

###############################################################################
#############################  HELPFUL FUNCTIONS ##############################
###############################################################################
    
#Function that encodes commands for decoder
def decode_command(locus_cmd):
    locus_cmd = locus_cmd.replace("___", " ")
    locus_cmd = locus_cmd.replace("^^^", ";")
    return(locus_cmd)

###############################################################################
#############################  DRIVER  ########################################
###############################################################################

## drive the script ##
## ONE-TIME CALL -- called by main
def driver(command_file):
    
    #Loop over the lines in the file
    with open(command_file) as file:
        for line in file:
            line = line.strip() #preprocess line
            command_string = decode_command(line)
            #Split command
            command_list = command_string.split(";")
            #Loop over the command_list and submit the commands
            for command in command_list:
                os.system(command)

###############################################################################


#cProfile.run('main(sys.argv[1:])')
         
if __name__ == "__main__":
    main(sys.argv[1:])

