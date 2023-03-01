############################################
#  Expanded RT Mutation Spectrum Analysis  #
#                                          #
#            David Mitchell III            #
#                 (c) 2023                 #
############################################

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
print('\nOpening file ' + sys.argv[1] + ' to analyze mutation spectrum and generate violin plots.\n')

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
    # atext3[0] = REF0 CT file                       #
    # atext3[1] = Modified Mutation Counts TXT file  #
    # atext3[2] = Untreated Mutation Counts TXT file #
    # atext3[3] = Profile TXT file                   #
    # atext3[4] = Nucleotide                         #
    # atext3[5] = Plot integration (Yes/No)          #
    # atext3[6] = Output directory                   #
    # atext3[7] = Output filename                    #
    ##################################################


    ## Importing CT file and performing calculations on mutation counts file using a given profile file ##
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

    ## Error checking variablers ##
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
    
    
    ## Calculating averages for plotting next to violin plots ##
    ary_ss = np.array(ss_corr)
    ary_bp = np.array(bp_corr)

    #Variables for calculating averages: single-strand nucleotide
    ss_avg1 = [0 for x in range(ary_ss.shape[0])]
    ss_average = [0 for x in range(ary_ss.shape[0])]
    ss_denom = ary_ss.shape[1]
    
    #Variables for calculating averages: base-paired nucleotide
    bp_avg1 = [0 for x in range(ary_bp.shape[0])]
    bp_average = [0 for x in range(ary_bp.shape[0])]
    bp_denom = ary_bp.shape[1]

    # Looking through each row of the single-stranded mutation array 'ary_ss' (where row = mutation type)
    for w in range(ary_ss.shape[0]):
        for z in range(ary_ss.shape[1]):
            ss_avg1[w] = ss_avg1[w] + ary_ss[w,z]
        ss_average[w] = ss_avg1[w] / ss_denom

    # Looking through each row of the base-paired mutation array 'ary_bp' (where row = mutation type)
    for w in range(ary_bp.shape[0]):
        for z in range(ary_bp.shape[1]):
            bp_avg1[w] = bp_avg1[w] + ary_bp[w,z]
        bp_average[w] = bp_avg1[w] / bp_denom

    # Variable for calculating difference in averages
    avg_diff = [0 for y in range(len(ss_average))]

    ## Calculating the difference in single-stranded and base-paired mutation rate averages by mutation type ##
    for y in range(len(ss_average)):
        avg_diff[y] = np.subtract(ss_average[y], bp_average[y])

    ## Compilating single-stranded nucleotide average rates into single variable ##
    ss_avg_all = []
    ss_avg_all.append(ss_average)
    ary_ss_avg_all = np.array(ss_avg_all)

    diff_avg_all = []
    diff_avg_all.append(avg_diff)
    ary_diff_avg_all = np.array(diff_avg_all)


    ## Ensure correct number of X ticks ##
    xlen = len(cmfile.header)-5
    cmfileheader1 = cmfile.header[0:xlen]

    ## Create variables for filtering out impossible mutation types from a selected nucleotide ##
    filter_aaa = [1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1]
    filter_ccc = [1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1]
    filter_ggg = [1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,1,1,1,0,0,0,1,1,1,1,1,1,1,1,1]
    filter_uuu = [1,1,1,1,1,1,1,1,1,0,0,0,1,1,1,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1]

    smaller_header = []
    smaller_ss_corr = []
    smaller_bp_corr = []
    smaller_ss_average = []
    smaller_bp_average = []
    smaller_avg_diff = []
    smaller_ss_avg_all = []

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
            #print('Adjusting')
            smaller_header.append(cmfileheader1[x])
            smaller_ss_corr.append(ss_corr[x])
            smaller_bp_corr.append(bp_corr[x])
            smaller_ss_average.append(ss_average[x])
            smaller_bp_average.append(bp_average[x])
            smaller_avg_diff.append(avg_diff[x])
            

    ## Printing arguments for confirmation ##
    print('\nCT File: ' + atext3[0] + '\nModified Mutation Counts file: ' + atext3[1] + '\nUntreated Mutation Counts file: ' + atext3[2] + '\nProfile file: ' + atext3[3] + '\nNucleotide Chosen: ' + atext3[4] + '\nNucleotide Used: ' + nt + '\n')

    ## Determination of Integrated or Separated Violin and Bar Plots ##
    plotint = atext3[5]

    if plotint[0] == 'Y' or plotint[0] == 'y':
        ## Creating combined violin and bar plots for all nu
        fig,ax = plot.subplots()

        ## Set the arrangement of the violin and bar plot elements ##
        ax.violinplot(smaller_bp_corr, 25*np.arange(len(smaller_bp_corr))-5, widths=10, showmedians=True)
        ax.bar(25*np.arange(len(smaller_bp_corr)), smaller_avg_diff, width=5, color='darkorchid')        
        ax.violinplot(smaller_ss_corr, 25*np.arange(len(smaller_bp_corr))+5, widths=10, showmedians=True)
        ax.set_title('BG-Corr Mutation Spectrum Output: ' + atext3[7] + '_' + nt)
        ax.set_xticks(25*np.arange(len(smaller_bp_corr)))
        ax.set_xticklabels(smaller_header, rotation=90)
        ax.set_ylabel('Mutation rate')

        ## Creating output PDF file ##
        print('\nCreating combined violin/bar plot file: MutSpec_Int_' + atext3[7] + '_' + nt + '.pdf\n..................................\n')
        outPath = atext3[6] + '/MutSpec_Int_' + atext3[7] + '_' + nt + '.pdf'
        fig.savefig(outPath, dpi=100, bbox_inches="tight", transparent=True)

        ## Save single-stranded nucleotide background-corrected mutation rates to TXT file ##
        print ('Saving single-stranded mutation rates to file: MutSpec_data_' + atext3[7] + '_' + nt + '.txt')
        outPathData0 = atext3[6] + '/MutSpec_data_' + atext3[7] + '_' + nt + '.txt'
        np.savetxt(outPathData0, ary_ss_avg_all, header=atext3[7])
    
    elif plotint[0] == 'N' or plotint[0] == 'n':
        ## Creating violin plots for all nucleotide position modifications (mutation/deletion/insertion) ##
        fig1,ax1 = plot.subplots()
        fig2,ax2 = plot.subplots()

        ## Set the arrangement of the violin plot elements (ax1) ##
        ax1.violinplot(smaller_bp_corr, 25*np.arange(len(smaller_bp_corr))-5, widths=10, showmedians=True)
        ax1.violinplot(smaller_ss_corr, 25*np.arange(len(smaller_bp_corr))+5, widths=10, showmedians=True)
        ax1.set_title('BG-Corr Mutation Spectrum Output: ' + atext3[7] + '_' + nt)
        ax1.set_xticks(25*np.arange(len(smaller_bp_corr)))
        ax1.set_xticklabels(smaller_header, rotation=90)
        ax1.set_ylabel('Mutation rate')


        ## Set the arrangment of the bar plot elements (ax2)
        ax2.bar(25*np.arange(len(smaller_bp_corr))-5, smaller_bp_average, width=4, color='royalblue', label='Base-paired')
        ax2.bar(25*np.arange(len(smaller_bp_corr)), smaller_avg_diff, width=4, color='darkorchid', label='Difference')
        ax2.bar(25*np.arange(len(smaller_bp_corr))+5, smaller_ss_average, width=4, color='lightcoral', label='Single-stranded')
        ax2.set_title('BG-Corr Mutation Spectrum Output: ' + atext3[7] + '_' + nt)
        ax2.set_xticks(25*np.arange(len(smaller_bp_corr)))
        ax2.set_xticklabels(smaller_header, rotation=90)
        ax2.set_ylabel('Mutation rate')
        ax2.legend()

        ## Creating output PDF file ##
        print('\nCreating separate violin plot and bar plot files: MutSpec_VPlt_' + atext3[7] + '_' + nt + '.pdf and MutSpec_Bar_' + atext3[7] + '_' + nt + '.pdf \n..................................\n')
        outPath1 = atext3[6] + '/MutSpec_VPlt_' + atext3[7] + '_' + nt + '.pdf'
        outPath2 = atext3[6] + '/MutSpec_Bar_' + atext3[7] + '_' + nt + '.pdf'
        fig1.savefig(outPath1, dpi=100, bbox_inches="tight", transparent=True)
        fig2.savefig(outPath2, dpi=100, bbox_inches="tight", transparent=True)

        ## Save single-stranded nucleotide background-corrected mutation rates to TXT file ##
        print ('Saving single-stranded mutation rates to file: MutSpec_dataSS_' + atext3[7] + '_' + nt + '.txt')
        outPathData = atext3[6] + '/MutSpec_dataSS_' + atext3[7] + '_' + nt + '.txt'
        np.savetxt(outPathData, ary_ss_avg_all, header="BG-Subtracted Single-Stranded Nucleotide Mutation Rates: " + atext3[7])

        ## Save differential (SS - BP) background-corrected mutation rates to TXT file ##
        print ('Saving differential (SS - BP) mutation rates to file: MutSpec_dataDF_' + atext3[7] + '_' + nt + '.txt')
        outPathData2 = atext3[6] + '/MutSpec_dataDF_' + atext3[7] + '_' + nt + '.txt'
        np.savetxt(outPathData2, ary_diff_avg_all, header="BG-Subtracted Mutation Rate Difference (SS - BP): " + atext3[7])


    else:
        print('Error: No plot generated for ' + atext3[7])

    ## Iterate the Loop Counter ##
    q = q + 1              

##################
# END WHILE LOOP #
##################


print('\n*** DONE: Completed mutation spectrum analysis and plot generation. ***\n')