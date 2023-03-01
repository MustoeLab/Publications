###########################################
#      RT Mutation Spectrum Analysis      #
#                                         #
#            David Mitchell III           #
#                 (c) 2023                #
###########################################

import numpy as np
import sys
import RNAtools2 as RNAtools
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plot
from ReactivityProfile import ReactivityProfile

##################################################################

## CountedMutationFile Class to Read the Mutation Counts File and Calculate Mutation Rates ##

class CountedMutationFile:

    def __init__(self, fname, reactprofile=None):
        with open(fname) as inp:
            self.header = inp.readline().split()    
        self.data = np.genfromtxt(fname, delimiter='\t', skip_header=1, dtype=int)
        self.compute_rates(reactprofile)
   
    def compute_rates(self, reactprofile):
        RP = None
        if reactprofile is not None:                 
            RP = ReactivityProfile(reactprofile)
        with np.errstate(invalid='ignore'):
            self.rates = np.array(self.data[:,:-5], dtype=float)
            for i in range(self.data.shape[0]):
                denom = self.data[i, -4]
                if denom < 10000 or np.sum(self.rates[i]) < 10:
                    denom = 0
                if RP is not None and np.isnan(RP.normprofile[i]) or RP.normprofile[i]<-5:
                    denom = 0
                if denom == 0:
                    self.rates[i] = 0
                else:
                    self.rates[i] /= denom

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
    # atext3[1] = Modified Mutation Counts TXT file  #
    # atext3[2] = Untreated Mutation Counts TXT file #
    # atext3[3] = Profile TXT file                   #
    # atext3[4] = Nucleotide                         #
    # atext3[5] = Plot colors (color1_color2)        #
    # atext3[6] = Output directory                   #
    # atext3[7] = Output filename                    #
    ##################################################


    ## Importing CT file and performing calculations on mutation counts file using a given profile file ##
    # CountedMutationFile variables: atext3[1/2] -> 'fname' and atext3[3] -> 'reactprofile'
    ct = RNAtools.CT(atext3[0])                         
    cmfile = CountedMutationFile(atext3[1], atext3[3])  
    cmfile2 = CountedMutationFile(atext3[2], atext3[3]) 

    ## Creating lists for elements of the violin plots - .shape[1] ensures reading by column (i.e. mutation type) ##
    bp = [[] for i in range(cmfile.rates.shape[1])]     
    ss = [[] for i in range(cmfile.rates.shape[1])]     
    ubp = [[] for i in range(cmfile.rates.shape[1])]    
    uss = [[] for i in range(cmfile.rates.shape[1])]    
    ss_corr = [[] for i in range(cmfile.rates.shape[1])]
    bp_corr = [[] for i in range(cmfile.rates.shape[1])]
    
    ## Input argument for specific nucleotide ##
    nt = atext3[4]

    ## Checking for accidental 'T' instead of 'U' for nucleotide - Note: RNAtools automatically replaces 'T' with 'U' ##
    if atext3[4] == 'T':
        print('Nucleotide input is "T". Changing to "U" for analysis.')
        nt = 'U'

    ## Error checking variables ##
    goodvar = 0
    badvar = 0

    ## Reads mutation rates from self.rates (where 'self' is the user-defined 'cmfile'), reading by rows (.shape[0])
    for i in range(cmfile.rates.shape[0]):
        if np.isfinite(cmfile.rates[i,0]) and ct.seq[i]==nt:
            goodvar += 1
            if ct.ct[i] == 0:
                for j in range(cmfile.rates.shape[1]):
                    ss_corr[j].append(cmfile.rates[i,j] - cmfile2.rates[i,j])
            else:
                for j in range(cmfile.rates.shape[1]):
                    bp_corr[j].append(cmfile.rates[i,j] - cmfile2.rates[i,j])
        else:
            badvar += 1
    
    
    ## Converting lists to arrays ##
    ary_ss = np.array(ss_corr)
    ary_bp = np.array(bp_corr)

    ## Ensure correct number of X ticks ##
    xlen = len(cmfile.header)-5
    cmfileheader1 = cmfile.header[0:xlen]

    print('Header Length:', xlen)

    if xlen == 26:
        # Old ShapeMapper 2.1.5
        filter_aaa = [1,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,1,1,1,1,1]
        filter_ccc = [0,0,0,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1]
        filter_ggg = [0,0,1,0,1,1,1,1,1,0,0,0,0,0,0,1,1,1,0,0,0,1,1,1,1,1]
        filter_uuu = [0,1,0,0,1,1,1,1,1,0,0,0,1,1,1,0,0,0,0,0,0,1,1,1,1,1]

        smaller_header = []
        smaller_ss_corr = []
        smaller_bp_corr = []

    elif xlen == 30:
        # ShapeMapper 2.2
        filter_aaa = [1,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,1,1,1]
        filter_ccc = [0,0,0,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,1,1,1,1,1,0,1,0,0,1,1,1]
        filter_ggg = [0,0,1,0,1,1,1,1,1,0,0,0,0,0,0,1,1,1,0,0,0,1,1,0,0,1,0,1,1,1]
        filter_uuu = [0,1,0,0,1,1,1,1,1,0,0,0,1,1,1,0,0,0,0,0,0,1,1,0,0,0,1,1,1,1]

        smaller_header = []
        smaller_ss_corr = []
        smaller_bp_corr = []


    ## Filter impossible mutation types from selected nucleotide and save into new variables ##
    # Setting the filter based on the selected nucleotide
    if nt == 'A':
        filter_use = filter_aaa
    elif nt == 'C':
        filter_use = filter_ccc
    elif nt == 'G':
        filter_use = filter_ggg
    elif nt == 'U':
        filter_use = filter_uuu

    # Filtering processes
    for x in range(ary_bp.shape[0]):
        if filter_use[x] == 1:
            smaller_header.append(cmfileheader1[x])
            smaller_ss_corr.append(ss_corr[x])
            smaller_bp_corr.append(bp_corr[x])

    # Convert new smaller SS and BP lists into arrays #
    smaller_ss_corr_arry = np.array(smaller_ss_corr)
    smaller_bp_corr_arry = np.array(smaller_bp_corr)

    print('smaller_ss_corr_arry:\n rows:', smaller_ss_corr_arry.shape[0], '\n cols:', smaller_ss_corr_arry.shape[1])
    print('smaller_bp_corr_arry:\n rows:', smaller_bp_corr_arry.shape[0], '\n cols:', smaller_bp_corr_arry.shape[1])
    print('\n')

    ## Save Output to Text File ##
    outPath_Data = atext3[6] + '/MutSpec_Data_MSFig_3_' + atext3[7] + '_' + nt + '.txt'
    MutSpec_datafile = open(outPath_Data, 'w')
    MutSpec_datafile.write('Mutation Spectrum Rates for Each Mutation Type')
    MutSpec_datafile.write('File Name: ' + atext3[7] + '\n\n')


    ## Print values for median mutation rate on screen ##
    # Single-stranded #
    for i in range(smaller_ss_corr_arry.shape[0]):
        ss_temp = []
        for j in range(smaller_ss_corr_arry.shape[1]):
            ss_temp.append(smaller_ss_corr_arry[i,j])
        ss_med = np.median(ss_temp)
        ss_mean = np.mean(ss_temp)
        ss_95th = np.percentile(ss_temp, 95)
        print('Single-stranded')
        print('Mutation type:', smaller_header[i])
        print(' Median rate:', ss_med)
        print(' Mean rate:', ss_mean)
        print(' 95th percentile:', ss_95th, '\n')
        MutSpec_datafile.write('Single-Stranded\n')
        MutSpec_datafile.write('Mutation type: ' + smaller_header[i] + '\n')
        MutSpec_datafile.write(' Median rate: ' + str(ss_med) + '\n')
        MutSpec_datafile.write(' Mean rate: ' + str(ss_mean) + '\n')
        MutSpec_datafile.write(' 95th percentile: ' + str(ss_95th) + '\n\n')

    # Base-paired #
    for i in range(smaller_bp_corr_arry.shape[0]):
        bp_temp = []
        for j in range(smaller_bp_corr_arry.shape[1]):
            bp_temp.append(smaller_bp_corr_arry[i,j])
        bp_med = np.median(bp_temp)
        bp_mean = np.mean(bp_temp)
        bp_95th = np.percentile(bp_temp, 95)
        print('Base-paired')
        print('Mutation type:', smaller_header[i])
        print(' Median rate:', bp_med)
        print(' Mean rate:', bp_mean)
        print(' 95th percentile:', bp_95th, '\n')
        MutSpec_datafile.write('Base-Paired\n')
        MutSpec_datafile.write('Mutation type: ' + smaller_header[i] + '\n')
        MutSpec_datafile.write(' Median rate: ' + str(bp_med) + '\n')
        MutSpec_datafile.write(' Mean rate: ' + str(bp_mean) + '\n')
        MutSpec_datafile.write(' 95th percentile: ' + str(bp_95th) + '\n\n')

    ## Printing arguments for confirmation ##
    print('\nCT File: ' + atext3[0] + '\nModified Mutation Counts file: ' + atext3[1] + '\nUntreated Mutation Counts file: ' + atext3[2] + '\nProfile file: ' + atext3[3] + '\nNucleotide Chosen: ' + atext3[4] + '\nNucleotide Used: ' + nt + '\n')

    ## Determination of Integrated or Separated Violin and Bar Plots ##
    plot.rc('xtick', labelsize=8)
    plot.rc('ytick', labelsize=8)

    ## Creating violin plots for all nucleotide position modifications (mutation/deletion/insertion) ##
    fig2,ax2 = plot.subplots()

    ## Set the arrangment of the violin plot elements (ax2) ##
    vparts1 = ax2.violinplot(smaller_bp_corr, 30*np.arange(len(smaller_bp_corr))-6, widths=12, showmedians=True)
    vparts2 = ax2.violinplot(smaller_ss_corr, 30*np.arange(len(smaller_ss_corr))+6, widths=12, showmedians=True)
    ax2.set_title('BG-Corr Mutation Spectrum Output: ' + atext3[7] + '_' + nt)
    ax2.set_xticks(30*np.arange(len(smaller_ss_corr)))
    ax2.set_xticklabels(smaller_header, rotation=90)
    ax2.set_ylabel('Mutation Rate', fontsize=8)
    ax2.set_ylim(-0.02,1)
    ax2.set_yscale('symlog', linthresh=.001, linscale=0.1, subs=[1,2,3,4,5,6,7,8,9])

    # Set properties for violin plot #
    colors = atext3[5].split('_')
    coloruse = colors[0]
    coloruse2 = colors[1]

    for partname in ('cbars', 'cmins', 'cmaxes', 'cmedians'):
        pc1 = vparts1[partname]
        pc2 = vparts2[partname]
        pc1.set_edgecolor(coloruse2)
        pc2.set_edgecolor(coloruse2)
        pc1.set_linewidth(1)
        pc2.set_linewidth(1)
    # Set face color for violin plots #
    for vp in vparts1['bodies']:
        vp.set_facecolor(coloruse)
        vp.set_alpha(1)
    for vp in vparts2['bodies']:
        vp.set_facecolor(coloruse)
        vp.set_alpha(1)

    ## Creating output PDF file ##
    print('\nCreating bar plot files: MutSpec_VPlt_MSFig_3_' + atext3[7] + '_' + nt + '.pdf \n..................................\n')
    outPath2 = atext3[6] + '/MutSpec_VPlt_MSFig_3_' + atext3[7] + '_' + nt + '.pdf'
    fig2.savefig(outPath2, dpi=100, bbox_inches="tight", transparent=True)

    ## Iterate the Loop Counter ##
    q = q + 1              

##################
# END WHILE LOOP #
##################

print('\n*** DONE: Completed mutation spectrum analysis and plot generation. ***\n')