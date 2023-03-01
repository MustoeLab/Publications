###########################################
#    Pairmapper PPV and Sens Plotter      #
#        Read Depth Undersampling         #
#                                         #
#           David Mitchell III            #
#                (c) 2023                 #
###########################################

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
lbls = ('Data', 'Enzyme', 'Us_count')
ptype_lbls = ('Principal', 'Minor', 'Both')
data_ppv_eDMS_50K = []
data_ppv_eDMS_200K = []
data_ppv_eDMS_600K = []
data_sens_eDMS_50K = []
data_sens_eDMS_200K = []
data_sens_eDMS_600K = []
data_ppv_origDMS_50K = []
data_ppv_origDMS_200K = []
data_ppv_origDMS_600K = []
data_sens_origDMS_50K = []
data_sens_origDMS_200K = []
data_sens_origDMS_600K = []

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
    # atext3[3] = Undersample read count                #
    # atext3[4] = Identifier                            #
    # atext3[5] = Output directory                      #
    #####################################################

    ct = RNAtools.CT(atext3[0], filterNC=True, filterSingle=True)   # Input CT file
    pm = pmanalysis.PairMap(atext3[1])                              # Input Pairmap file
    profile, seq = RNAtools.readSHAPE(atext3[2])                    # Input DMS file
    us_count = atext3[3]
    rna_info = atext3[4]

    ## Calculate PPV and Sensitivity using Module from Tony ##
    ppv, sens = pm.ppvsens_duplex(ct, ptype = pt, exact = False, profile = profile)
    # Convert PPV/Sensitivity values from fraction to integer and round up to next integer
    ppv = int(round(100 * ppv, 0))
    sens = int(round(100 * sens, 0))
    
    # Append data to lists
    if us_count == 'Us50000':
        if rna_info == 'origDMS':
            data_ppv_origDMS_50K.append(ppv)
            data_sens_origDMS_50K.append(sens)
        elif rna_info == 'eDMS':
            data_ppv_eDMS_50K.append(ppv)
            data_sens_eDMS_50K.append(sens)
    elif us_count == 'Us200000':
        if rna_info == 'origDMS':
            data_ppv_origDMS_200K.append(ppv)
            data_sens_origDMS_200K.append(sens)
        elif rna_info == 'eDMS':
            data_ppv_eDMS_200K.append(ppv)
            data_sens_eDMS_200K.append(sens)
    elif us_count == 'Us600000':
        if rna_info == 'origDMS':
            data_ppv_origDMS_600K.append(ppv)
            data_sens_origDMS_600K.append(sens)
        elif rna_info == 'eDMS':
            data_ppv_eDMS_600K.append(ppv)
            data_sens_eDMS_600K.append(sens)
  
    # Iterate Loop Counter
    j = j + 1

##################
# END WHILE LOOP #
##################

## Calculate Means ##

mean_ppv_origDMS_50K = np.mean(data_ppv_origDMS_50K)
mean_ppv_origDMS_200K = np.mean(data_ppv_origDMS_200K)
mean_ppv_origDMS_600K = np.mean(data_ppv_origDMS_600K)
mean_ppv_eDMS_50K = np.mean(data_ppv_eDMS_50K)
mean_ppv_eDMS_200K = np.mean(data_ppv_eDMS_200K)
mean_ppv_eDMS_600K = np.mean(data_ppv_eDMS_600K)

mean_sens_origDMS_50K = np.mean(data_sens_origDMS_50K)
mean_sens_origDMS_200K = np.mean(data_sens_origDMS_200K)
mean_sens_origDMS_600K = np.mean(data_sens_origDMS_600K)
mean_sens_eDMS_50K = np.mean(data_sens_eDMS_50K)
mean_sens_eDMS_200K = np.mean(data_sens_eDMS_200K)
mean_sens_eDMS_600K = np.mean(data_sens_eDMS_600K)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
print('PPV Mean:')
print('  OrigDMS US 50K:', mean_ppv_origDMS_50K)
print('  OrigDMS US 200K:', mean_ppv_origDMS_200K)
print('  OrigDMS US 600K:', mean_ppv_origDMS_600K)
print('  eDMS US 50K:', mean_ppv_eDMS_50K)
print('  eDMS US 200K:', mean_ppv_eDMS_200K)
print('  eDMS US 600K:', mean_ppv_eDMS_600K, '\n')
print('Sens Mean:')
print('  OrigDMS US 50K:', mean_sens_origDMS_50K)
print('  OrigDMS US 200K:', mean_sens_origDMS_200K)
print('  OrigDMS US 600K:', mean_sens_origDMS_600K)
print('  eDMS US 50K:', mean_sens_eDMS_50K)
print('  eDMS US 200K:', mean_sens_eDMS_200K)
print('  eDMS US 600K:', mean_sens_eDMS_600K, '\n')
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


## Create PPV and Sensitivity Plots ##
names = ('ppv (%)', 'sens (%)')

plot_ppv_mean = (mean_ppv_origDMS_50K, mean_ppv_eDMS_50K, mean_ppv_origDMS_200K, mean_ppv_eDMS_200K, mean_ppv_origDMS_600K, mean_ppv_eDMS_600K)
plot_sens_mean = (mean_sens_origDMS_50K, mean_sens_eDMS_50K, mean_sens_origDMS_200K, mean_sens_eDMS_200K, mean_sens_origDMS_600K, mean_sens_eDMS_600K)
x_data = (1,2,4,5,7,8)
color_use = ('red', 'blue', 'red', 'blue', 'red', 'blue')


fig, ax1 = plot.subplots(2,1, figsize=(6,4.5))

# Create Box plot #
for z, ztype in enumerate(names):
    plot.rc('xtick', labelsize=8)   # Change the font size of the x-axis text
    plot.rc('ytick', labelsize=8)   # Change the font size of the y-axis text
    plot.rc('legend', fontsize=6)   # Change the font size of the legend text
    if ztype == 'ppv (%)':
        y_use = plot_ppv_mean
    elif ztype == 'sens (%)':
        y_use = plot_sens_mean
    print('Using data:', ztype)
    ax1[z].scatter(x_data, y_use, c = color_use)        
    ax1[z].set_ylabel(ztype)
    ax1[z].set_ylim(-5,105)
    
#fig.suptitle('PPV and Sens Box Plot: ' + ptype_name)
fig.tight_layout()

# Save box plot as PDF
outPath_W = atext3[5] + '/PAIRppvsens70_US_eDMSvDMS_' + ptype_name + '.pdf'
fig.savefig(outPath_W, dpi=100)

print('\n\nFinished with calculating Pairmapper PPV and Sensitivity values and creating plots.\n\n')
