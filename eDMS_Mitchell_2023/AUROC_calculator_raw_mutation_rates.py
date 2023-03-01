###########################################
#           ROC Plot Generator            #
#     Uncorrected (raw) reactivities      #
#                                         #
#           David Mitchell III            #
#                (c) 2023                 #
###########################################

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plot
import sys
import RNAtools2 as RNAtools
import numpy as np
from sklearn.metrics import roc_curve, roc_auc_score


## Subprogram to create data for ROC plot from input files ##
def getNtLists(profile, pairedlist, ntype):

    react = []
    ispaired = []
    for i,v in enumerate(profile.rawprofile):               
        if v > -10 and profile.sequence[i] == ntype:
            react.append(v)
            ispaired.append( int((i+1) in pairedlist) )
        
    return react, ispaired


## Open and Read User Arguments from Text File ##

afile = open(sys.argv[1])             # Open the file with the file name given by the user
atext = afile.read()                  # Read the contents of the file
atext1 = atext.split('\n')            # Split the file into a list where each list element is a complete line from the file
atextlen = len(atext1)                # Give count of number of lines in original file


## User Confirmation ##
print('\nOpening file ' + sys.argv[1] + ' to generate ROC plots.\n')

##################################################################
# Program runs in a While loop from 0 to atextlen in order to    #
# further parse the line elements and operate on the components. #
##################################################################

j = 0                              # The variable 'j' is the loop index counter
line_check = 0                     # Variable 'line_check' is the loop index counter for the empty line check
line_error = 0                     # Variable 'line_error' is subtracted from 'textlen' to eliminate empty lines
loopend = atextlen

## Check for Empty Lines at the End of the Input File ##
print('Checking input file for blank lines...')

while line_check < atextlen:
    if atext1[line_check] == '':
        line_error += 1
        print('** Detected ' + str(line_error) + ' terminal blank lines in the input file ' + str(sys.argv[1]) + ' **')
    line_check += 1

loopend = loopend - line_error

## User Confirmation of Total Number of Entries ##
print('\nCalculating a total of ' + str(loopend) + ' AUCs.\n')

####################
# BEGIN WHILE LOOP #
####################

while j < loopend:                # Start of the while loop
    atext2 = str(atext1[j])        # Changes list created from file line into a string
    atext3 = atext2.split()        # Splits the string into a list containing 4 elements

    ################################################
    # From here, the correlations between variable #
    # names and input data are as follows:         #
    #                                              #
    # atext3[0] = CT file (.ct)                    #
    # atext3[1] = Profile text file (.txt)         #
    ################################################


    ## Creating Profile and CT Objects from User Input ##
    ctobj = RNAtools.CT(atext3[0])
    profobj = ReactivityProfile(atext3[1])

    ## Using RNAtools2 to generate paired residue list ##
    # Residue added to list if value in column 4 (which gives the pairing partner for a residue) = 0
    pairedlist = ctobj.pairedResidueList(False)

    ## Creating Graphical Compilation of ROC Plots Using MatPlotLib ##
    labels = ('A', 'C', 'G', 'U')
    for i,nt in enumerate(labels):
        r, p = getNtLists(profobj, pairedlist, nt)
        print('Name:', atext3[1])           
        print('  Nucleotide:', nt, '\n  AUC: {:.2f}'.format(roc_auc_score(p,r)), '\n')
    
    j = j + 1              # Iterate the loop counter

##################
# END WHILE LOOP #
##################

print('\nDone.\n')