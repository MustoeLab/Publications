###########################################
#    RT Mutation Rate Histogram Plotter   #
#                                         #
#           David Mitchell III            #
#                (c) 2023                 #
###########################################

import sys
import numpy as np
import RNAtools2 as RNAtools
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plot
from ReactivityProfile import ReactivityProfile

######################################################################

def getDMSNtLists(profile, pairedlist):
    react = []
    ispaired = []
    for i,v in enumerate(profile.rawprofile):               
        if v >= .0001 and (profile.sequence[i] == 'A' or profile.sequence[i] == 'C'):
            react.append(v)
            ispaired.append( int((i+1) in pairedlist) )        
    return react, ispaired

def getEtOHNtLists(profile, pairedlist):
    react = []
    ispaired = []
    for i,v in enumerate(profile.backprofile):               
        if v >= 0.0001 and (profile.sequence[i] == 'A' or profile.sequence[i] == 'C'):
            react.append(v)
            ispaired.append( int((i+1) in pairedlist) )
    return react, ispaired

def mergeprofile(plist):
    m = ReactivityProfile()
    m.sequence = np.array([])
    m.normprofile = np.array([])
    m.backprofile = np.array([])
    m.rawprofile = np.array([])
    m.subprofile = np.array([])
    for f in plist:
        p = ReactivityProfile(f)
        m.sequence = np.append(m.sequence, p.sequence)
        m.normprofile = np.append(m.normprofile, p.normprofile)
        m.backprofile = np.append(m.backprofile, p.backprofile)
        m.rawprofile = np.append(m.rawprofile, p.rawprofile)
        m.subprofile = np.append(m.subprofile, p.subprofile)
    return m

########################################################################

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
print('\nGenerating a total of ' + str(loopend) + ' plots.\n')


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
    # atext3[0] = Output directory                   #
    # atext3[1] = Output filename                    #
    # atext3[2] = SS2 Profile TXT RMRP Rep 1         #
    # atext3[3] = SS2 Profile TXT RMRP Rep 2         #
    # atext3[4] = SS2 Profile TXT RnaseP Rep 1       #
    # atext3[5] = SS2 Profile TXT tmRNA-ic Rep 1     #
    # atext3[6] = SS2 Profile TXT tmRNA-ic Rep 2     #
    # atext3[7] = SS2 Profile TXT tmRNA-cf           #
    # atext3[8] = Mara Profile TXT RMRP Rep 1        #
    # atext3[9] = Mara Profile TXT RMRP Rep 2        #
    # atext3[10] = Mara Profile TXT RnaseP Rep 1     #
    # atext3[11] = Mara Profile TXT tmRNA-ic Rep 1   #
    # atext3[12] = Mara Profile TXT tmRNA-ic Rep 2   #
    # atext3[13] = Mara Profile TXT tmRNA-cf         #
    # atext3[14] = TGIRT Profile TXT RMRP Rep 1      #
    # atext3[15] = TGIRT Profile TXT RMRP Rep 2      #
    # atext3[16] = TGIRT Profile TXT RnaseP Rep 1    # 
    # atext3[17] = TGIRT Profile TXT tmRNA-ic Rep 1  #
    # atext3[18] = TGIRT Profile TXT tmRNA-ic Rep 2  #
    # atext3[19] = TGIRT Profile TXT tmRNA-cf        #
    # atext3[20] = HIV Profile TXT RMRP Rep 1        #
    # atext3[21] = HIV Profile TXT RMRP Rep 2        #
    ##################################################

    ## Define output name for histogram plot saving ##
    outPath = atext3[0] + '/MutRateViolin21c_MSFig_plot_' + atext3[1] + '.pdf'

    ## Merge Profiles ##
    profobj_SS2 = mergeprofile([atext3[2], atext3[3], atext3[4], atext3[5], atext3[6], atext3[7]])
    profobj_Mara = mergeprofile([atext3[8], atext3[9], atext3[10], atext3[11], atext3[12], atext3[13]])
    profobj_TGIRT = mergeprofile([atext3[14], atext3[15], atext3[16], atext3[17], atext3[18], atext3[19]])
    profobj_HIV = mergeprofile([atext3[20], atext3[21]])

    seq_list = tuple(np.arange(1,50000))

    ## Generate Lists of Single-Stranded and Base-Paired Nucleotides ##
    # DMS Modified
    r_SS2_DMS, p_SS2_DMS = getDMSNtLists(profobj_SS2, seq_list)
    r_Mara_DMS, p_Mara_DMS = getDMSNtLists(profobj_Mara, seq_list)
    r_TGIRT_DMS, p_TGIRT_DMS = getDMSNtLists(profobj_TGIRT, seq_list)
    r_HIV_DMS, p_HIV_DMS = getDMSNtLists(profobj_HIV, seq_list)
    # Ethanol control
    r_SS2_EtOH, p_SS2_EtOH = getEtOHNtLists(profobj_SS2, seq_list)
    r_Mara_EtOH, p_Mara_EtOH = getEtOHNtLists(profobj_Mara, seq_list)
    r_TGIRT_EtOH, p_TGIRT_EtOH = getEtOHNtLists(profobj_TGIRT, seq_list)
    r_HIV_EtOH, p_HIV_EtOH = getEtOHNtLists(profobj_HIV, seq_list)


    DMS_SS2_use = []
    EtOH_SS2_use = []
    DMS_Mara_use = []
    EtOH_Mara_use = []
    DMS_TGIRT_use = []
    EtOH_TGIRT_use = []
    DMS_HIV_use = []
    EtOH_HIV_use = []
    
    # For SS2 #
    for j in range(len(p_SS2_DMS)):
        if r_SS2_DMS[j] > 0:
            DMS_SS2_use.append(r_SS2_DMS[j])
    for j in range(len(p_SS2_EtOH)):
        if r_SS2_EtOH[j] > 0:
            EtOH_SS2_use.append(r_SS2_EtOH[j])
    # For Marathon #
    for j in range(len(p_Mara_DMS)):
        if r_Mara_DMS[j] > 0:
            DMS_Mara_use.append(r_Mara_DMS[j])
    for j in range(len(p_Mara_EtOH)):
        if r_Mara_EtOH[j] > 0:
            EtOH_Mara_use.append(r_Mara_EtOH[j])
    # For TGIRT #
    for j in range(len(p_TGIRT_DMS)):
        if r_TGIRT_DMS[j] > 0:
            DMS_TGIRT_use.append(r_TGIRT_DMS[j])
    for j in range(len(p_TGIRT_EtOH)):
        if r_TGIRT_EtOH[j] > 0:
            EtOH_TGIRT_use.append(r_TGIRT_EtOH[j])
    # For HIV #
    for j in range(len(p_HIV_DMS)):
        if r_HIV_DMS[j] > 0:
            DMS_HIV_use.append(r_HIV_DMS[j])
    for j in range(len(p_HIV_EtOH)):
        if r_HIV_EtOH[j] > 0:
            EtOH_HIV_use.append(r_HIV_EtOH[j])
        
        
    ## Create Plot ##
    data_plot = [DMS_SS2_use, EtOH_SS2_use, DMS_Mara_use, EtOH_Mara_use, DMS_TGIRT_use, EtOH_TGIRT_use, DMS_HIV_use, EtOH_HIV_use]

    plot.rc('xtick', labelsize=12)
    plot.rc('ytick', labelsize=12)

    bp = plot.boxplot(data_plot, patch_artist=True)

    colors = ['limegreen', 'lightgreen', 'cornflowerblue', 'lightskyblue', 'darkorange', 'navajowhite', 'violet', 'pink']
    colors2 = ['limegreen', 'limegreen', 'lightgreen', 'lightgreen', 'cornflowerblue', 'cornflowerblue', 'lightskyblue', 'lightskyblue', 'darkorange', 'darkorange', 'navajowhite', 'navajowhite', 'violet', 'violet', 'pink', 'pink']
    colors3 = ['darkgreen', 'darkgreen', 'mediumblue', 'mediumblue', 'sienna', 'sienna', 'purple', 'purple']

    for patch, colorx in zip(bp['boxes'], colors):
        patch.set(facecolor = colorx)
    
    for whisker, colorx in zip(bp['whiskers'], colors2):
        whisker.set(color = colorx, linewidth = 1.5)

    for caps, colorx in zip(bp['caps'], colors2):
        caps.set(color = colorx, linewidth = 1.5)

    for median, colorx in zip(bp['medians'], colors3):
        median.set(color = colorx, linewidth = 1.5)

    for flier, colorx in zip(bp['fliers'], colors):
        flier.set(markeredgecolor = colorx, alpha = 0.5)

    plot.title(atext3[1])
    plot.ylabel('Mutation Rate', fontsize=12)
    plot.ylim(0.00008, 1)
    plot.yscale('log')

    #plot.tight_layout()
    plot.savefig(outPath, dpi=100, bbox_inches="tight", transparent=True)
    plot.clf()

    ## Iterate the Loop Counter ##
    q = q + 1              

##################
# END WHILE LOOP #
##################


print('\n*** DONE ***\n')