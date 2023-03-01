############################################
#   Multiple ROC Plot Generator Program    #
#                                          #
#           David Mitchell III             #
#                (c) 2023                  #
############################################

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plot
import sys
import RNAtools2 as RNAtools
import numpy as np
from sklearn.metrics import roc_curve, roc_auc_score
from ReactivityProfile import ReactivityProfile


## Subprogram to create data for ROC plot from input files ##
def getNtLists(profile, pairedlist, ntype):

    react = []
    ispaired = []
    for i,v in enumerate(profile.subprofile):               
        if v > -10 and profile.sequence[i] == ntype:
            react.append(v)
            ispaired.append( int((i+1) in pairedlist) )
        
    return react, ispaired

#########################################################

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
print('\nGenerating a total of ' + str(loopend) + ' ROC plots.\n')


####################
# BEGIN WHILE LOOP #
####################

while j < loopend:                 # Start of the while loop
    atext2 = str(atext1[j])        # Changes list created from file line into a string
    atext3 = atext2.split()        # Splits the string into a list containing 4 elements

    ################################################
    # From here, the correlations between variable #
    # names and input data are as follows:         #
    #                                              #
    # atext3[0] = Output directory                 #
    # atext3[1] = Output file name                 #
    # atext3[2] = CT file (.ct)                    #
    # atext3[3] = Profile file 1                   #
    # atext3[4] = Profile file 2                   #
    # atext3[...] = Profile file ...               #
    # atext3[n] = Profile file n                   #
    ################################################

    ## Creating Profile and CT Objects from User Input ##
    # Runs ReactivityProfile.py on profile TXT file
    profobj1 = ReactivityProfile(atext3[3])
    profobj2 = ReactivityProfile(atext3[4])

    ctobj1 = RNAtools.CT(atext3[2]) 

    ## Using RNAtools2 to generate paired residue list ##
    # Residue added to list if value in column 4 (which gives the pairing partner for a residue) > 0
    pairedlist1 = ctobj1.pairedResidueList(False)


    ## Creating Graphical Compilation of ROC Plots Using MatPlotLib ##
    #labels = ('A', 'C', 'G', 'U')  # <-- * Normal labels list, do not remove *
    labels = ('G')

    fig,ax = plot.subplots(1, 4, figsize=(8,2))

    for i,nt in enumerate(labels):
        # Setting variables for ROC plot generation by SciKit-Learn (sklearn.metrics)
        # 'r' = Profile file background-subtracted mutation rates
        # 'p' = Paired residue list from CT file
        r1, p1 = getNtLists(profobj1, pairedlist1, nt)         
        r2, p2 = getNtLists(profobj2, pairedlist1, nt)                


        # Variable 'fpr' = false positive rate
        # Variable 'tpr' = true psitive rate
        # Variable 'thr' = thresholds for the decision function
        fpr1, tpr1, thr1 = roc_curve(p1, r1, drop_intermediate=False)
        fpr2, tpr2, thr2 = roc_curve(p2, r2, drop_intermediate=False)
        
        ax[i].plot(fpr1, tpr1, color='dimgray', lw=2)         
        ax[i].plot(fpr2, tpr2, color='violet', lw=2)                
        ax[i].text(1,0.01,'Orig={:.2f}'.format(roc_auc_score(p1,r1)) + '\nFilter={:.2f}'.format(roc_auc_score(p2,r2)), horizontalalignment='right', transform=ax[i].transAxes, fontsize=6)
        ax[i].plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
        ax[i].set_xlim([-0.02, 1.0])
        ax[i].set_ylim([0.0, 1.02])
        ax[i].set_title(nt)
    
    fig.suptitle('Multi-ROC Plot: ' + atext3[1])
    fig.tight_layout()


    ## Creating Output PDF File ##
    print('\nCreating ROC plot at output file: multiroc_plot_10_' + atext3[1] + '.pdf\n..................')
    outPath1 = atext3[0] + '/multiroc_plot_10_' + atext3[1] + '.pdf'
    fig.savefig(outPath1, dpi=100, bbox_inches="tight", transparent=True)


    j = j + 1              # Iterate the loop counter

##################
# END WHILE LOOP #
##################

print('\nCompleted ROC plot generation.\n')