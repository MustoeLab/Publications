###########################################
#    RT Mutation Rate Histogram Plotter   #
#                                         #
#            David Mitchell III           #
#                 (c) 2023                #
###########################################

import sys
import numpy as np
import RNAtools2 as RNAtools
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plot
from ReactivityProfile import ReactivityProfile
import xml.etree.ElementTree as ET

#############################################################################################

## Subprogram to create data for ROC plot from input files ##
def getNtLists(profile, pairedlist, ntype):
    react = []
    ispaired = []
    for i,v in enumerate(profile.subprofile):               
        if v > -20 and profile.sequence[i] == ntype:
            react.append(v)
            ispaired.append( int((i+1) in pairedlist) )
        
    return react, ispaired

## Opening XML Reactivity Files from Incarnato Lab ##
def readXML(xml):

    tree = ET.parse(xml)
    root = tree.getroot()
    s = root[0][0].text.split()
    seq = ''
    for i in s:
        seq+=i
    seq = seq.replace('T','U')

    data = [float(x) for x in root[0][1].text.split(',')]
    
    prof = ReactivityProfile()
    prof.sequence = np.array(list(seq))
    prof.subprofile = np.array(data)

    return prof


##################################################################################################

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
    else:
        print('Line OK')
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
    # atext3[0] = RMRP CT file                       #
    # atext3[1] = RnaseP CT file                     #
    # atext3[2] = HEK293 ic 18S rRNA CT file         #
    # atext3[3] = eDMS RMRP Profile TXT file         #
    # atext3[4] = 2A3 RMRP Profile TXT file          #
    # atext3[x] = eDMS RnaseP Profile TXT file       #
    # atext3[x] = 2A3 RnaseP Profile TXT file        #
    # atext3[7] = 2A3 HEK293 ic 18S rRNA Profile     #
    # atext3[8] = Output directory                   #
    # atext3[9] = Output filename                    #
    ##################################################

    ## Define output name for histogram plot saving ##
    outPath = atext3[8] + '/MutRateViolin25b_MSFig_S6_plot_' + atext3[9] + '.pdf'


    ## List of Nucleotides ##
    nucl_tuple = ('A', 'C', 'U', 'G')

    ## CT File and Profile File Objects ##
    ctobj_RMRP = RNAtools.CT(atext3[0])
    ctobj_RnaseP = RNAtools.CT(atext3[1])
    ctobj_18S_rRNA = RNAtools.CT(atext3[2])

    profobj_eDMS_RMRP = ReactivityProfile(atext3[3])
    profobj_2A3_RMRP = ReactivityProfile(atext3[4])
    profobj_eDMS_RnaseP = ReactivityProfile(atext3[5])
    profobj_2A3_RnaseP = ReactivityProfile(atext3[6])
    profobj_2A3_18S_rRNA = readXML(atext3[7])
    
    ss_list_RMRP = ctobj_RMRP.pairedResidueList(False) # Returns list of ** SINGLE-STRANDED ** nucleotides
    ss_list_RnaseP = ctobj_RnaseP.pairedResidueList(False) # Returns list of ** SINGLE-STRANDED ** nucleotides
    ss_list_18S_rRNA = ctobj_18S_rRNA.pairedResidueList(False) # Returns list of ** SINGLE-STRANDED ** nucleotides


    ## Save Output to Text File ##
    outPathData = atext3[8] + '/MutRateViolin25b_MSFig_S6_Data_' + atext3[9] + '.txt'
    VP_datafile = open(outPathData, 'w')


    ## Generate Lists of Single-Stranded and Base-Paired Nucleotides ##
    for i,nt in enumerate(nucl_tuple):
        r1, p1 = getNtLists(profobj_eDMS_RMRP, ss_list_RMRP, nt)
        r2, p2 = getNtLists(profobj_2A3_RMRP, ss_list_RMRP, nt)
        r3, p3 = getNtLists(profobj_eDMS_RnaseP, ss_list_RnaseP, nt)
        r4, p4 = getNtLists(profobj_2A3_RnaseP, ss_list_RnaseP, nt)
        r5, p5 = getNtLists(profobj_2A3_18S_rRNA, ss_list_18S_rRNA, nt)

        ss_eDMS_RMRP_use = []
        bp_eDMS_RMRP_use = []
        ss_2A3_RMRP_use = []
        bp_2A3_RMRP_use = []
        ss_eDMS_RnaseP_use = []
        bp_eDMS_RnaseP_use = []
        ss_2A3_RnaseP_use = []
        bp_2A3_RnaseP_use = []
        ss_2A3_RnaseP_use = []
        bp_2A3_RnaseP_use = []
        ss_2A3_18S_rRNA_use = []
        bp_2A3_18S_rRNA_use = []

        # For eDMS RMRP #
        for j in range(len(p1)):
            if p1[j] == 0:
                bp_eDMS_RMRP_use.append(r1[j])
            elif p1[j] == 1:
                ss_eDMS_RMRP_use.append(r1[j])

        # For 2A3 RMRP #
        for j in range(len(p2)):
            if p2[j] == 0:
                bp_2A3_RMRP_use.append(r2[j])
            elif p2[j] == 1:
                ss_2A3_RMRP_use.append(r2[j])

        # For eDMS RnaseP #
        for j in range(len(p3)):
            if p3[j] == 0:
                bp_eDMS_RnaseP_use.append(r3[j])
            elif p3[j] == 1:
                ss_eDMS_RnaseP_use.append(r3[j])

        # For 2A3 RnaseP #
        for j in range(len(p4)):
            if p4[j] == 0:
                bp_2A3_RnaseP_use.append(r4[j])
            elif p4[j] == 1:
                ss_2A3_RnaseP_use.append(r4[j])


        # For 2A3 18S rRNA #
        for j in range(len(p5)):
            if p5[j] == 0:
                bp_2A3_18S_rRNA_use.append(r5[j])
            elif p5[j] == 1:
                ss_2A3_18S_rRNA_use.append(r5[j])


        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
        VP_datafile.write('Name: ' + atext3[9] + '\n')
        VP_datafile.write('------------------------------\n')
        VP_datafile.write('Nucleotide: ' + nt + '\n')
        VP_datafile.write('  Minimum bp_eDMS_RMRP: ' + str(np.min(bp_eDMS_RMRP_use)) + '\n')
        VP_datafile.write('  Minimum ss_eDMS_RMRP: ' + str(np.min(ss_eDMS_RMRP_use)) + '\n')
        VP_datafile.write('  Minimum bp_2A3_RMRP: ' + str(np.min(bp_2A3_RMRP_use)) + '\n')
        VP_datafile.write('  Minimum ss_2A3_RMRP: ' + str(np.min(ss_2A3_RMRP_use)) + '\n')
        VP_datafile.write('  Minimum bp_eDMS_RnaseP: ' + str(np.min(bp_eDMS_RnaseP_use)) + '\n')
        VP_datafile.write('  Minimum ss_eDMS_RnaseP: ' + str(np.min(ss_eDMS_RnaseP_use)) + '\n')
        VP_datafile.write('  Minimum bp_2A3_RnaseP: ' + str(np.min(bp_2A3_RnaseP_use)) + '\n')
        VP_datafile.write('  Minimum ss_2A3_RnaseP: ' + str(np.min(ss_2A3_RnaseP_use)) + '\n')
        VP_datafile.write('  Minimum bp_2A3_18S_rRNA: ' + str(np.min(bp_2A3_18S_rRNA_use)) + '\n')
        VP_datafile.write('  Minimum ss_2A3_18S_rRNA: ' + str(np.min(ss_2A3_18S_rRNA_use)) + '\n')

        VP_datafile.write('  Median bp_eDMS_RMRP: ' + str(np.median(bp_eDMS_RMRP_use)) + '\n')
        VP_datafile.write('  Median ss_eDMS_RMRP: ' + str(np.median(ss_eDMS_RMRP_use)) + '\n')
        VP_datafile.write('  Median bp_2A3_RMRP: ' + str(np.median(bp_2A3_RMRP_use)) + '\n')
        VP_datafile.write('  Median ss_2A3_RMRP: ' + str(np.median(ss_2A3_RMRP_use)) + '\n')
        VP_datafile.write('  Median bp_eDMS_RnaseP: ' + str(np.median(bp_eDMS_RnaseP_use)) + '\n')
        VP_datafile.write('  Median ss_eDMS_RnaseP: ' + str(np.median(ss_eDMS_RnaseP_use)) + '\n')
        VP_datafile.write('  Median bp_2A3_RnaseP: ' + str(np.median(bp_2A3_RnaseP_use)) + '\n')
        VP_datafile.write('  Median ss_2A3_RnaseP: ' + str(np.median(ss_2A3_RnaseP_use)) + '\n')
        VP_datafile.write('  Median bp_2A3_18S_rRNA: ' + str(np.median(bp_2A3_18S_rRNA_use)) + '\n')
        VP_datafile.write('  Median ss_2A3_18S_rRNA: ' + str(np.median(ss_2A3_18S_rRNA_use)) + '\n')

        VP_datafile.write('  Mean bp_eDMS_RMRP: ' + str(np.mean(bp_eDMS_RMRP_use)) + '\n')
        VP_datafile.write('  Mean ss_eDMS_RMRP: ' + str(np.mean(ss_eDMS_RMRP_use)) + '\n')
        VP_datafile.write('  Mean bp_2A3_RMRP: ' + str(np.mean(bp_2A3_RMRP_use)) + '\n')
        VP_datafile.write('  Mean ss_2A3_RMRP: ' + str(np.mean(ss_2A3_RMRP_use)) + '\n')
        VP_datafile.write('  Mean bp_eDMS_RnaseP: ' + str(np.mean(bp_eDMS_RnaseP_use)) + '\n')
        VP_datafile.write('  Mean ss_eDMS_RnaseP: ' + str(np.mean(ss_eDMS_RnaseP_use)) + '\n')
        VP_datafile.write('  Mean bp_2A3_RnaseP: ' + str(np.mean(bp_2A3_RnaseP_use)) + '\n')
        VP_datafile.write('  Mean ss_2A3_RnaseP: ' + str(np.mean(ss_2A3_RnaseP_use)) + '\n')
        VP_datafile.write('  Mean bp_2A3_18S_rRNA: ' + str(np.mean(bp_2A3_18S_rRNA_use)) + '\n')
        VP_datafile.write('  Mean ss_2A3_18S_rRNA: ' + str(np.mean(ss_2A3_18S_rRNA_use)) + '\n')

        VP_datafile.write('  Maximum bp_eDMS_RMRP: ' + str(np.max(bp_eDMS_RMRP_use)) + '\n')
        VP_datafile.write('  Maximum ss_eDMS_RMRP: ' + str(np.max(ss_eDMS_RMRP_use)) + '\n')
        VP_datafile.write('  Maximum bp_2A3_RMRP: ' + str(np.max(bp_2A3_RMRP_use)) + '\n')
        VP_datafile.write('  Maximum ss_2A3_RMRP: ' + str(np.max(ss_2A3_RMRP_use)) + '\n')
        VP_datafile.write('  Maximum bp_eDMS_RnaseP: ' + str(np.max(bp_eDMS_RnaseP_use)) + '\n')
        VP_datafile.write('  Maximum ss_eDMS_RnaseP: ' + str(np.max(ss_eDMS_RnaseP_use)) + '\n')
        VP_datafile.write('  Maximum bp_2A3_RnaseP: ' + str(np.max(bp_2A3_RnaseP_use)) + '\n')
        VP_datafile.write('  Maximum ss_2A3_RnaseP: ' + str(np.max(ss_2A3_RnaseP_use)) + '\n')
        VP_datafile.write('  Maximum bp_2A3_18S_rRNA: ' + str(np.max(bp_2A3_18S_rRNA_use)) + '\n')
        VP_datafile.write('  Maximum ss_2A3_18S_rRNA: ' + str(np.max(ss_2A3_18S_rRNA_use)) + '\n\n')

        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


        plot.subplot(1, 4, i+1)
        plot.rc('xtick', labelsize=12)
        plot.rc('ytick', labelsize=12)
        vparts1 = plot.violinplot(bp_eDMS_RMRP_use, np.arange(1)-45, widths=10, showmedians=True)
        vparts2 = plot.violinplot(ss_eDMS_RMRP_use, np.arange(1)-35, widths=10, showmedians=True)
        vparts3 = plot.violinplot(bp_2A3_RMRP_use, np.arange(1)-25, widths=10, showmedians=True)
        vparts4 = plot.violinplot(ss_2A3_RMRP_use, np.arange(1)-15, widths=10, showmedians=True)
        vparts5 = plot.violinplot(bp_eDMS_RnaseP_use, np.arange(1)-5, widths=10, showmedians=True)
        vparts6 = plot.violinplot(ss_eDMS_RnaseP_use, np.arange(1)+5, widths=10, showmedians=True)
        vparts7 = plot.violinplot(bp_2A3_RnaseP_use, np.arange(1)+15, widths=10, showmedians=True)
        vparts8 = plot.violinplot(ss_2A3_RnaseP_use, np.arange(1)+25, widths=10, showmedians=True)
        vparts9 = plot.violinplot(bp_2A3_18S_rRNA_use, np.arange(1)+35, widths=10, showmedians=True)
        vparts10 = plot.violinplot(ss_2A3_18S_rRNA_use, np.arange(1)+45, widths=10, showmedians=True)


        # Set properties for violin plot #
        for partname in ('cbars', 'cmins', 'cmaxes', 'cmedians'):
            pc1 = vparts1[partname]
            pc2 = vparts2[partname]
            pc3 = vparts3[partname]
            pc4 = vparts4[partname]
            pc5 = vparts5[partname]
            pc6 = vparts6[partname]
            pc7 = vparts7[partname]
            pc8 = vparts8[partname]
            pc9 = vparts9[partname]
            pc10 = vparts10[partname]
            pc1.set_edgecolor('mediumblue')
            pc2.set_edgecolor('mediumblue')
            pc3.set_edgecolor('purple')
            pc4.set_edgecolor('purple')
            pc5.set_edgecolor('mediumblue')
            pc6.set_edgecolor('mediumblue')
            pc7.set_edgecolor('purple')
            pc8.set_edgecolor('purple')
            pc9.set_edgecolor('black')
            pc10.set_edgecolor('black')
            pc1.set_linewidth(1)
            pc2.set_linewidth(1)
            pc3.set_linewidth(1)
            pc4.set_linewidth(1)
            pc5.set_linewidth(1)
            pc6.set_linewidth(1)
            pc7.set_linewidth(1)
            pc8.set_linewidth(1)
            pc9.set_linewidth(1)
            pc10.set_linewidth(1)

        # Set face color for violin plots #
        for vp in vparts1['bodies']:
            vp.set_facecolor('cornflowerblue')
        for vp in vparts2['bodies']:
            vp.set_facecolor('cornflowerblue')
        for vp in vparts3['bodies']:
            vp.set_facecolor('violet')
        for vp in vparts4['bodies']:
            vp.set_facecolor('violet')
        for vp in vparts5['bodies']:
            vp.set_facecolor('cornflowerblue')
        for vp in vparts6['bodies']:
            vp.set_facecolor('cornflowerblue')
        for vp in vparts7['bodies']:
            vp.set_facecolor('violet')
        for vp in vparts8['bodies']:
            vp.set_facecolor('violet')
        for vp in vparts9['bodies']:
            vp.set_facecolor('gray')
        for vp in vparts10['bodies']:
            vp.set_facecolor('gray')



        plot.title(nt)
        #plot.ylabel('Mutation Rate', fontsize=12)
        plot.ylim(-0.016, 0.6)
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