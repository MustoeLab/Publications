###########################################
#    Pairmapper PPV and Sens Plotter      #
#                                         #
#           David Mitchell III            #
#                (c) 2023                 #
###########################################

from collections import UserString
import sys
import pmanalysis
import RNAtools2 as RNAtools
import matplotlib.pyplot as plot
import numpy as np


## Open and Read User Arguments from Text File ##

afile = open(sys.argv[1])             # Open the file with the file name given by the user
atext = afile.read()                  # Read the contents of the file
atext1 = atext.split('\n')            # Split the file into a list where each list element is a complete line from the file
atextlen = len(atext1)                # Give count of number of lines in original file


## User Confirmation ##
print('\nOpening file ' + sys.argv[1] + ' to calculate PPV and Sensitivity values from Pairmapper output and create plots.\n')

## Initiating Data Lists ##
namelist = []
ptype_lbls = ('Principal', 'Minor', 'Both')
data_ppv = []
data_sens = []

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

pt = int(sys.argv[2])
ptype_name = ptype_lbls[(pt - 1)]

outputname = sys.argv[3]

## User Confirmation of Total Number of Entries ##
print('\nCalculating Pairmapper PPV and Sensitivity values for ' + ptype_name + ' pairs for ' + str(loopend) + ' items.\n')


####################
# BEGIN WHILE LOOP #
####################

while j < loopend:                # Start of the while loop
    atext2 = str(atext1[j])        # Changes list created from file line into a string
    atext3 = atext2.split()        # Splits the string into a list containing 4 elements

    #####################################################
    # From here, the correlations between variable      #
    # names and input data are as follows:              #
    #                                                   #
    # atext3[0] = CT file (.ct)                         #
    # atext3[1] = Pairmap File (*-pairmap.txt)          #
    # atext3[2] = DMS File (.dms)                       #
    #####################################################

    ## Perform PPV and Sensitivity Calculations ##
    ct = RNAtools.CT(atext3[0], filterNC=True, filterSingle=True)
    pm = pmanalysis.PairMap(atext3[1])                           
    profile, seq = RNAtools.readSHAPE(atext3[2])                 

    # Calculate PPV and Sensitivity using Module from Tony
    ppv, sens = pm.ppvsens_duplex(ct, ptype = pt, exact = False, profile = profile)
    # Convert PPV/Sensitivity values from fraction to integer and round up to next integer
    ppv = int(round(100 * ppv, 0))
    sens = int(round(100 * sens, 0))
    print('Name: ', atext3[2], '\n   PPV: ', ppv, '\n   Sens: ', sens)
    # Append data to lists
    data_ppv.append(ppv)
    data_sens.append(sens)

    # Iterate Loop Counter
    j = j + 1

##################
# END WHILE LOOP #
##################

data_ppv_use = np.array(data_ppv)
data_sens_use = np.array(data_sens)
ecolors = ('cornflowerblue', 'limegreen', 'darkorange', 'cornflowerblue', 'limegreen', 'darkorange', 'cornflowerblue', 'limegreen', 'darkorange', 'cornflowerblue', 'limegreen', 'darkorange', 'cornflowerblue', 'limegreen', 'darkorange', 'cornflowerblue', 'limegreen', 'darkorange', 'cornflowerblue', 'limegreen', 'darkorange')

x_pos = (0,1,2,0,1,2,3,4,5,3,4,5,6,7,8,6,7,8,9,10,11)

## Create PPV and Sensitivity Plots ##
names = ('PPV', 'Sensitivity')
fig, ax1 = plot.subplots(2,1)

# Create Box plot #
for z, ztype in enumerate(names):
    plot.rc('xtick', labelsize=8)
    plot.rc('ytick', labelsize=8)
    plot.rc('legend', fontsize=6)
    plot.figure(figsize=(4,4))
    if ztype == 'PPV':
        y_use = data_ppv_use
    elif ztype == 'Sensitivity':
        y_use = data_sens_use
    for az in range(y_use.shape[0]):
        ax1[z].scatter(x_pos, y_use, s=60, marker='o', c = ecolors, edgecolors = ecolors, linewidths=2)
    ax1[z].set_ylabel(ztype)
    ax1[z].set_ylim(-5,105)
    
fig.tight_layout()

# Save box plot as PDF
outPath = 'Output/Pairmap_MSFIG_SI_FigS10A_Plot_' + outputname + '_' + ptype_name + '.pdf'
fig.savefig(outPath, dpi=100)

print('\n\nFinished with calculating Pairmapper PPV and Sensitivity values and creating plots.\n\n')
