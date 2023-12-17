'''
execute_command.py

The purpose of this script is to make the matrix/map files determined by 
make_matrices.py. It takes in a single command parameter and executes the 
steps of that command in sequence. 

***This needs to be run on Python/3.9 as there's an incompatability with ***
***other versions of python.                                             ***
'''

#import modules that will be needed
import sys #Taking command line parameters as inputs
import os


###############################################################################
###########################  COLLECT AND CHECK ARGUMENTS  #####################
###############################################################################

## GLOBAL VARS

## MAIN ##

def main(argv): 

    ## Required arguments
    try:
        command_string = sys.argv[1]
    except KeyError:
        print("ERROR: One of your required arguments does not exist.")
        help()
    print("Acceptable Inputs Given")
    #Call driver function
    driver(command_string)

###############################################################################
#############################  DRIVER  ########################################
###############################################################################

## drive the script ##
## ONE-TIME CALL -- called by main
def driver(command_string):
    
    #Decode triple underscore and carrots
    command_string = command_string.replace("___", " ")
    command_string = command_string.replace("^^^", ";")
    #Split command
    command_list = command_string.split(";")
    #Loop over the command_list and submit the commands
    for command in command_list:
        os.system(command)

###############################################################################


#cProfile.run('main(sys.argv[1:])')
         
if __name__ == "__main__":
    main(sys.argv[1:])

