############################################
#   Plot Combined Mutation Rate and AUC    #
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
from sklearn.metrics import roc_auc_score
from ReactivityProfile import ReactivityProfile

####################################################################

## Subprogram to create data for ROC plot from input files ##
def getFullLists(profile, pairedlist, nt):

    react = []
    ispaired = []
    for i,v in enumerate(profile.subprofile):               
        if v > -10 and profile.sequence[i] == nt:
            react.append(v)
            ispaired.append( int((i+1) in pairedlist) )
        
    return react, ispaired

####################################################################

## Open and Read User Arguments from Text File ##

afile = open(sys.argv[1])             # Open the file with the file name given by the user
atext = afile.read()                  # Read the contents of the file
atext1 = atext.split('\n')            # Split the file into a list where each list element is a complete line from the file
atextlen = len(atext1)                # Give count of number of lines in original file


## User Confirmation ##
print('\nOpening file ' + sys.argv[1] + ' to generate AUC and mutation rate bar plots.\n')

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
    else:
        print('Line OK')
    line_check += 1

loopend = loopend - line_error

## User Confirmation of Total Number of Entries ##
print('\nGenerating a total of ' + str(loopend) + ' AUC bar plots.\n')


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
    # atext3[2] = CT file                          #
    # atext3[3] = Replicate A Profile file 1       #
    # atext3[4] = Replicate A Profile file 2       #
    # atext3[5] = Replicate A Profile file 3       #
    # atext3[6] = Replicate B Profile file 1       #
    # atext3[7] = Replicate B Profile file 2       #
    # atext3[8] = Replicate B Profile file 3       #
    ################################################

    ## Creating Profile and CT Objects from User Input ##
    profobj_A1 = ReactivityProfile(atext3[3]) # REP A unbuffered
    profobj_A2 = ReactivityProfile(atext3[4]) # REP A pH 7
    profobj_A3 = ReactivityProfile(atext3[5]) # REP A pH 8
    profobj_B1 = ReactivityProfile(atext3[6]) # REP B unbuffered
    profobj_B2 = ReactivityProfile(atext3[7]) # REP B pH 7
    profobj_B3 = ReactivityProfile(atext3[8]) # REP B pH 8

    ctobj1 = RNAtools.CT(atext3[2])    
    pairedlist1 = ctobj1.pairedResidueList(False)

    ntlist = ('A', 'C', 'G', 'U')

    grand_auc_list = []
    grand_auc_error_list = []
    grand_mutrate_list = []
    grand_mutrate_error_list = []

    ## Calculate Mutation Rates and AUC Values for each Nucleotide #    
    for i, nt in enumerate(ntlist):
        ra1, pa1 = getFullLists(profobj_A1, pairedlist1, nt)        
        ra2, pa2 = getFullLists(profobj_A2, pairedlist1, nt)               
        ra3, pa3 = getFullLists(profobj_A3, pairedlist1, nt)
        rb1, pb1 = getFullLists(profobj_B1, pairedlist1, nt)        
        rb2, pb2 = getFullLists(profobj_B2, pairedlist1, nt)               
        rb3, pb3 = getFullLists(profobj_B3, pairedlist1, nt)

        r_sA1 = []
        r_sA2 = []
        r_sA3 = []
        r_sB1 = []
        r_sB2 = []
        r_sB3 = []

        # For Rep A s1 #
        for w in range(len(pa1)):
            if pa1[w] == 1:
                r_sA1.append(ra1[w])
        # For Rep A s2 #
        for w in range(len(pa2)):
            if pa2[w] == 1:
                r_sA2.append(ra2[w])
        # For Rep A s3 #
        for w in range(len(pa3)):
            if pa3[w] == 1:
                r_sA3.append(ra3[w])
        # For Rep B s1 #
        for w in range(len(pb1)):
            if pb1[w] == 1:
                r_sB1.append(rb1[w])
        # For Rep B s2 #
        for w in range(len(pb2)):
            if pb2[w] == 1:
                r_sB2.append(rb2[w])
        # For Rep B s3 #
        for w in range(len(pb3)):
            if pb3[w] == 1:
                r_sB3.append(rb3[w])

        # Calculate mean and standard deviation for each nucleotide #

        rate_stdev_1 = []
        rate_stdev_2 = []
        rate_stdev_3 = []

        aucA1 = roc_auc_score(pa1, ra1)
        aucA2 = roc_auc_score(pa2, ra2)
        aucA3 = roc_auc_score(pa3, ra3)
        aucB1 = roc_auc_score(pb1, rb1)
        aucB2 = roc_auc_score(pb2, rb2)
        aucB3 = roc_auc_score(pb3, rb3)

        r_sA1_use = np.mean(r_sA1)
        r_sA2_use = np.mean(r_sA2)
        r_sA3_use = np.mean(r_sA3)
        r_sB1_use = np.mean(r_sB1)
        r_sB2_use = np.mean(r_sB2)
        r_sB3_use = np.mean(r_sB3)

        auc_1_use = np.mean([aucA1, aucB1])
        auc_2_use = np.mean([aucA2, aucB2])
        auc_3_use = np.mean([aucA3, aucB3])

        rate_1_use = (r_sA1_use + r_sB1_use) / 2
        rate_2_use = (r_sA2_use + r_sB2_use) / 2
        rate_3_use = (r_sA3_use + r_sB3_use) / 2

        rate_stdev_1 = np.std([r_sA1_use, r_sB1_use])
        rate_stdev_2 = np.std([r_sA2_use, r_sB2_use])
        rate_stdev_3 = np.std([r_sA3_use, r_sB3_use])

        auc_stdev_1 = np.std([aucA1, aucB1])
        auc_stdev_2 = np.std([aucA2, aucB2])
        auc_stdev_3 = np.std([aucA3, aucB3])

        auctotal = (auc_1_use, auc_2_use, auc_3_use)
        aucerror = (auc_stdev_1, auc_stdev_2, auc_stdev_3)
        mutratetotal = (rate_1_use, rate_2_use, rate_3_use)
        mutrateerror = (rate_stdev_1, rate_stdev_2, rate_stdev_3)

        grand_auc_list.append(auctotal)
        grand_auc_error_list.append(aucerror)
        grand_mutrate_list.append(mutratetotal)
        grand_mutrate_error_list.append(mutrateerror)

    xtiks = (1,1.1,1.2)
    labels = ('Unbuffered', 'pH 7.2', 'pH 8')

    # Define color dictionary #
    colorA = {'A' : 'blue', 'C' : 'limegreen', 'G' : 'orange', 'U' : 'violet'}
    colorB = {'A' : 'navy', 'C' : 'darkgreen', 'G' : 'darkorange', 'U' : 'purple'}

    ## Creating Plots ##
    names = ('AUC', 'Mutation Rate')
    fig, ax = plot.subplots(2,1, figsize=(3,4))
    plot.rc('xtick', labelsize=10)
    plot.rc('ytick', labelsize=10)

    for z, ztype in enumerate(names):
        for w, wtype in enumerate(ntlist):
            if ztype == 'AUC':
                data_use = grand_auc_list[w]
                error_use = grand_auc_error_list[w]
                ylimrange = (0.38, 1.02)
                yscale_use = 'linear'
                coloruseA = colorA.get(wtype)
                coloruseB = colorB.get(wtype)
            elif ztype == 'Mutation Rate':
                data_use = grand_mutrate_list[w]
                error_use = grand_mutrate_error_list[w]
                ylimrange = (0.0001, 0.1)
                yscale_use = 'log'
                coloruseA = colorA.get(wtype)
                coloruseB = colorB.get(wtype)

            ax[z].plot(xtiks, data_use, 'o', ms = 8 , mec = coloruseA, mfc = 'None', ls = '-', label = nt)
            ax[z].errorbar(xtiks, data_use, yerr=error_use, color=coloruseA, ecolor=coloruseB, elinewidth=0)
            ax[z].set_ylim(ylimrange)
            ax[z].set_ylabel(ztype)
            ax[z].set_yscale(yscale_use)

    plot.suptitle(atext3[1])
    plot.tight_layout()

    ## Creating Output PDF File ##
    print('\nCreating AUC and mutation rate bar plots at output file: MutRateAUC_21_' + atext3[1] + '.pdf\n..................')
    outPath1 = atext3[0] + '/MutRateAUC_21_' + atext3[1] + '.pdf'
    plot.savefig(outPath1, dpi=100, bbox_inches="tight", transparent=True)


    j = j + 1              # Iterate the loop counter

##################
# END WHILE LOOP #
##################

print('\nCompleted AUC bar plot generation.\n')