##########################################
#   RT Mutation Rate Histogram Plotter   #
#                                        #
#           David Mitchell III           #
#                (c) 2023                #
##########################################

import sys
import numpy as np
import RNAtools2 as RNAtools
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plot
from ReactivityProfile import ReactivityProfile


## Subprogram to create data for ROC plot from input files ##
def getNtLists(profile, pairedlist, ntype):

    react = []
    ispaired = []
    for i,v in enumerate(profile.subprofile):               
        if v > -20 and profile.sequence[i] == ntype:
            react.append(v)
            ispaired.append( int((i+1) in pairedlist) )
        
    return react, ispaired

######################################################################

## Open and Read User Arguments from Text File ##

afile = open(sys.argv[1])             # Open the file with the file name given by the user
atext = afile.read()                  # Read the contents of the file
atext1 = atext.split('\n')            # Split the file into a list where each list element is a complete line from the file
atextlen = len(atext1)                # Give count of number of lines in original file


## User Confirmation of File Open ##
print('\nOpening file ' + sys.argv[1] + ' to analyze mutation spectrum and generate mean mutation rate bar plots.\n')

##################################################################
# Program runs in a While loop from 0 to atextlen in order to    #
# further parse the line elements and operate on the components. #
##################################################################


q = 0                              # Variable 'q' is the primary loop index counter
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
print('\nGenerating a total of ' + str(loopend) + ' mutation spectrum plots.\n')


####################
# BEGIN WHILE LOOP #
####################

while q < loopend:                 # Start of the while loop
    atext2 = str(atext1[q])        # Changes list created from file line into a string
    atext3 = atext2.split()        # Splits the string into a list containing 6 elements

    ################################################## 
    # From here, the correlations between variable   #
    # names and input data are as follows:           #
    #                                                #
    # atext3[0] = REF0 CT file                       #
    # atext3[1] = SS2 Profile TXT file               #
    # atext3[2] = Marathon Profile TXT file          #
    # atext3[3] = TGIRT Profile TXT file             #
    # atext3[4] = Output directory                   #
    # atext3[5] = Output filename                    #
    ##################################################

    ## Define output name for histogram plot saving ##
    outPath = atext3[4] + '/MutRateViolin20_MSFig_plot_' + atext3[5] + '.pdf'

    ## List of Nucleotides ##
    nucl_tuple = ('A', 'C', 'U', 'G')

    ## CT File and Profile File Objects ##
    ctobj1 = RNAtools.CT(atext3[0])
    profobj_SS2 = ReactivityProfile(atext3[1])
    profobj_Mara = ReactivityProfile(atext3[2])
    profobj_TGIRT = ReactivityProfile(atext3[3])

    ss_list = ctobj1.pairedResidueList(False) # Returns list of ** SINGLE-STRANDED ** nucleotides


    ## Generate Lists of Single-Stranded and Base-Paired Nucleotides ##
    for i,nt in enumerate(nucl_tuple):
        # Setting variables for ROC plot generation by SciKit-Learn (sklearn.metrics)
        # 'r' = Profile file background-subtracted mutation rates
        # 'p' = Paired residue list from CT file
        r1, p1 = getNtLists(profobj_SS2, ss_list, nt)
        r2, p2 = getNtLists(profobj_Mara, ss_list, nt)
        r3, p3 = getNtLists(profobj_TGIRT, ss_list, nt)

        ss_SS2_use = []
        bp_SS2_use = []
        ss_Mara_use = []
        bp_Mara_use = []
        ss_TGIRT_use = []
        bp_TGIRT_use = []
        
        # For SS2 #
        for j in range(len(p1)):
            if p1[j] == 0:
                bp_SS2_use.append(r1[j])
            elif p1[j] == 1:
                ss_SS2_use.append(r1[j])

        # For Marathon #
        for j in range(len(p2)):
            if p2[j] == 0:
                bp_Mara_use.append(r2[j])
            elif p2[j] == 1:
                ss_Mara_use.append(r2[j])

        # For TGIRT #
        for j in range(len(p3)):
            if p3[j] == 0:
                bp_TGIRT_use.append(r3[j])
            elif p3[j] == 1:
                ss_TGIRT_use.append(r3[j])

        plot.subplot(1, 4, i+1)
        plot.rc('xtick', labelsize=12)
        plot.rc('ytick', labelsize=12)
        vparts1 = plot.violinplot(bp_SS2_use, np.arange(1)-25, widths=10, showmedians=True)
        vparts2 = plot.violinplot(ss_SS2_use, np.arange(1)-15, widths=10, showmedians=True)
        vparts3 = plot.violinplot(bp_Mara_use, np.arange(1)-5, widths=10, showmedians=True)
        vparts4 = plot.violinplot(ss_Mara_use, np.arange(1)+5, widths=10, showmedians=True)
        vparts5 = plot.violinplot(bp_TGIRT_use, np.arange(1)+15, widths=10, showmedians=True)
        vparts6 = plot.violinplot(ss_TGIRT_use, np.arange(1)+25, widths=10, showmedians=True)

        # Set properties for violin plot #
        for partname in ('cbars', 'cmins', 'cmaxes', 'cmedians'):
            pc1 = vparts1[partname]
            pc2 = vparts2[partname]
            pc3 = vparts3[partname]
            pc4 = vparts4[partname]
            pc5 = vparts5[partname]
            pc6 = vparts6[partname]
            pc1.set_edgecolor('limegreen')
            pc2.set_edgecolor('limegreen')
            pc3.set_edgecolor('cornflowerblue')
            pc4.set_edgecolor('cornflowerblue')
            pc5.set_edgecolor('darkorange')
            pc6.set_edgecolor('darkorange')
            pc1.set_linewidth(1)
            pc2.set_linewidth(1)
            pc3.set_linewidth(1)
            pc4.set_linewidth(1)
            pc5.set_linewidth(1)
            pc6.set_linewidth(1)

        # Set face color for violin plots #
        for vp in vparts1['bodies']:
            vp.set_facecolor('limegreen')
        for vp in vparts2['bodies']:
            vp.set_facecolor('limegreen')
        for vp in vparts3['bodies']:
            vp.set_facecolor('cornflowerblue')
        for vp in vparts4['bodies']:
            vp.set_facecolor('cornflowerblue')
        for vp in vparts5['bodies']:
            vp.set_facecolor('darkorange')
        for vp in vparts6['bodies']:
            vp.set_facecolor('darkorange')

        plot.title(nt)
        #plot.ylabel('Mutation Rate', fontsize=12)
        plot.ylim(-0.016, 1)
        plot.yscale('symlog', linthresh=.0099, linscale=.25, subs=[1,2,3,4,5,6,7,8,9])

    #plot.tight_layout()
    plot.savefig(outPath, dpi=100, bbox_inches="tight", transparent=True)
    plot.clf()

    ## Iterate the Loop Counter ##
    q = q + 1              

##################
# END WHILE LOOP #
##################


print('\n*** DONE ***\n')