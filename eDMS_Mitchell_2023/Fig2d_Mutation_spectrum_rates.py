###########################################
#      RT Mutation Spectrum Analysis      #
#        Mutation Rate Calculator         #
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
    uni_corr = [[] for i in range(cmfile.rates.shape[1])]       
    
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
            for j in range(cmfile.rates.shape[1]):
                uni_corr[j].append(cmfile.rates[i,j] - cmfile2.rates[i,j])
        else:
            badvar += 1
    
    
    ## Converting list to arrays ##
    ary_uni = np.array(uni_corr)

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
        smaller_uni_corr = []
        smaller_bp_corr = []

    elif xlen == 30:
        # ShapeMapper 2.2
        filter_aaa = [1,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,1,1,1]
        filter_ccc = [0,0,0,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,1,1,1,1,1,0,1,0,0,1,1,1]
        filter_ggg = [0,0,1,0,1,1,1,1,1,0,0,0,0,0,0,1,1,1,0,0,0,1,1,0,0,1,0,1,1,1]
        filter_uuu = [0,1,0,0,1,1,1,1,1,0,0,0,1,1,1,0,0,0,0,0,0,1,1,0,0,0,1,1,1,1]

        smaller_header = []
        smaller_uni_corr = []


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
    for x in range(ary_uni.shape[0]):
        if filter_use[x] == 1:
            smaller_header.append(cmfileheader1[x])
            smaller_uni_corr.append(uni_corr[x])

    # Convert new smaller UNI list into array #
    smaller_uni_corr_arry = np.array(smaller_uni_corr)

    print('smaller_uni_corr_arry:\n rows:', smaller_uni_corr_arry.shape[0], '\n cols:', smaller_uni_corr_arry.shape[1])
    print('\n')

    ## Save Output to Text File ##
    outPath_Data = atext3[6] + '/MutSpec_Data_MSFig_3bb_' + atext3[7] + '_' + nt + '.txt'
    MutSpec_datafile = open(outPath_Data, 'w')
    MutSpec_datafile.write('Mutation Spectrum Rates for Each Mutation Type')
    MutSpec_datafile.write('File Name: ' + atext3[7] + '\n\n')


    ## Print values for median mutation rate on screen ##
    for i in range(smaller_uni_corr_arry.shape[0]):
        uni_temp = []
        for j in range(smaller_uni_corr_arry.shape[1]):
            uni_temp.append(smaller_uni_corr_arry[i,j])
        uni_med = np.median(uni_temp)
        uni_mean = np.mean(uni_temp)
        uni_95th = np.percentile(uni_temp, 95)
        print('Mutation type:', smaller_header[i])
        print(' Median rate:', uni_med)
        print(' Mean rate:', uni_mean, '\n')
        print(' 95th Percentile:', uni_95th)
        MutSpec_datafile.write('Mutation type: ' + smaller_header[i] + '\n')
        MutSpec_datafile.write(' Median rate: ' + str(uni_med) + '\n')
        MutSpec_datafile.write(' Mean rate: ' + str(uni_mean) + '\n')
        MutSpec_datafile.write(' 95th percentile: ' + str(uni_95th) + '\n\n')


    ## Printing arguments for confirmation ##
    print('\nCT File: ' + atext3[0] + '\nModified Mutation Counts file: ' + atext3[1] + '\nUntreated Mutation Counts file: ' + atext3[2] + '\nProfile file: ' + atext3[3] + '\nNucleotide Chosen: ' + atext3[4] + '\nNucleotide Used: ' + nt + '\n')

    ## Iterate the Loop Counter ##
    q = q + 1              

##################
# END WHILE LOOP #
##################


print('\n*** DONE: Completed mutation spectrum analysis and plot generation. ***\n')