###########################################
#      RT Mutation Spectrum Analysis      #
#      Indel vs Mismatch Calculator       #
#                                         #
#            David Mitchell III           #
#                 (c) 2023                #
###########################################

import numpy as np
import sys
import RNAtools2 as RNAtools
import matplotlib
matplotlib.use('Agg')
import math

##################################################################

class CountedMutationFile:
    def __init__(self, fname):
        with open(fname) as inp:
            self.header = inp.readline().split()
        self.data = np.genfromtxt(fname, delimiter='\t', skip_header=1, dtype=int)

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
print('\nWorking with ' + str(loopend) + ' total inputs.\n')


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
    # atext3[3] = Identifier                         #
    ##################################################


    ## Importing CT file and performing calculations on mutation counts file using a given profile file ##
    ct = RNAtools.CT(atext3[0])              
    cmfile = CountedMutationFile(atext3[1])  
    cmfile2 = CountedMutationFile(atext3[2]) 

    ## Creating lists for elements of the violin plots - .shape[1] ensures reading by column (i.e. mutation type) ##
    ss_corr = [[] for i in range(cmfile.data.shape[1])]
    
    ## Error checking variables ##
    goodvar = 0
    badvar = 0

    ## Reads mutation rates from self.rates (where 'self' is the user-defined 'cmfile'), reading by rows (.shape[0])
    for i in range(cmfile.data.shape[0]):
        if np.isfinite(cmfile.data[i,0]):
            goodvar += 1
            for j in range(cmfile.data.shape[1]):
                ss_corr[j].append(cmfile.data[i,j] - cmfile2.data[i,j])
        else:
            badvar += 1
    
    ## Converting lists to arrays ##
    ary_ss = np.array(ss_corr)

    ## Ensure correct number of X ticks ##
    xlen = len(cmfile.header)-5
    cmfileheader1 = cmfile.header[0:xlen]

    print('Length:', xlen)

    if xlen == 30:
        ## FOR NEWER SHAPEMAPPER - HEADER LENGTH 30 ##
        # Create variables for filtering out impossible mutation types from a selected nucleotide #
        # Del-A	 Del-T	Del-G	Del-C	Ins-A	Ins-T	Ins-G	Ins-C	Ins-N	AT	AG	AC	TA	TG	TC	GA	GT	GC	CA	CT	CG	multinuc_deletion	multinuc_insertion	A_multinuc_mismatch	C_multinuc_mismatch	G_multinuc_mismatch	T_multinuc_mismatch	N_multinuc_mismatch	complex_deletion	complex_insertion
        # Indel List: 0, 1, 2, 3, 4, 5, 6, 7, 8, 21, 22, 28, 29
        # Mismatch List: 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 23, 24, 25, 26, 27

        # Define variables #
        ss_sub_list = []
        ss_indel_list = []

        # Combine all indel mutation types together #
        for i in range(ary_ss.shape[1]):
            ss_indel_list.append(math.fsum([ary_ss[0,i], 
                                        ary_ss[1,i], 
                                        ary_ss[2,i], 
                                        ary_ss[3,i], 
                                        ary_ss[4,i], 
                                        ary_ss[5,i], 
                                        ary_ss[6,i], 
                                        ary_ss[7,i], 
                                        ary_ss[8,i], 
                                        ary_ss[21,i], 
                                        ary_ss[22,i], 
                                        ary_ss[28,i], 
                                        ary_ss[29,i]]))

        for i in range(ary_ss.shape[1]):
            ss_sub_list.append(math.fsum([ary_ss[9,i], 
                                        ary_ss[10,i], 
                                        ary_ss[11,i], 
                                        ary_ss[12,i], 
                                        ary_ss[13,i], 
                                        ary_ss[14,i], 
                                        ary_ss[15,i], 
                                        ary_ss[16,i], 
                                        ary_ss[17,i], 
                                        ary_ss[18,i], 
                                        ary_ss[19,i], 
                                        ary_ss[20,i],
                                        ary_ss[23,i], 
                                        ary_ss[24,i], 
                                        ary_ss[25,i], 
                                        ary_ss[26,i], 
                                        ary_ss[27,i]]))

    elif xlen == 26:
        ## FOR OLD SHAPEMAPPER - HEADER LENGTH 26 ##
        # Create variables for filtering out impossible mutation types from a selected nucleotide #
        # Del-A	 Del-T	Del-G	Del-C	Ins-A	Ins-T	Ins-G	Ins-C	Ins-N	AT	AG	AC	TA	TG	TC	GA	GT	GC	CA	CT	CG	multinuc_deletion	multinuc_insertion	multinuc_mismatch	complex_deletion	complex_insertion
        # Indel List: 0, 1, 2, 3, 4, 5, 6, 7, 8, 21, 22, 24, 25
        # Mismatch List: 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 23

        # Define variables #
        ss_sub_list = []
        ss_indel_list = []

        # Combine all indel mutation types together #
        for i in range(ary_ss.shape[1]):
            ss_indel_list.append(math.fsum([ary_ss[0,i], 
                                        ary_ss[1,i], 
                                        ary_ss[2,i], 
                                        ary_ss[3,i], 
                                        ary_ss[4,i], 
                                        ary_ss[5,i], 
                                        ary_ss[6,i], 
                                        ary_ss[7,i], 
                                        ary_ss[8,i], 
                                        ary_ss[21,i], 
                                        ary_ss[22,i], 
                                        ary_ss[24,i], 
                                        ary_ss[25,i]]))

        for i in range(ary_ss.shape[1]):
            ss_sub_list.append(math.fsum([ary_ss[9,i], 
                                        ary_ss[10,i], 
                                        ary_ss[11,i], 
                                        ary_ss[12,i], 
                                        ary_ss[13,i], 
                                        ary_ss[14,i], 
                                        ary_ss[15,i], 
                                        ary_ss[16,i], 
                                        ary_ss[17,i], 
                                        ary_ss[18,i], 
                                        ary_ss[19,i], 
                                        ary_ss[20,i],
                                        ary_ss[23,i]]))

    ## Calculate Sums ##
    indel_len_use = len(ss_indel_list)
    mm_len_use = len(ss_sub_list)
    total_indels = 0
    total_mm = 0
    grand_total = 0

    for i in range(indel_len_use):
        total_indels = total_indels + ss_indel_list[i]

    for i in range(mm_len_use):
        total_mm = total_mm + ss_sub_list[i]

    grand_total = total_indels + total_mm


    # Print Output #
    print('Name:', atext3[3])
    print('  Indels:', int(total_indels))
    print('  % Indels:', total_indels / grand_total)
    print('  Mismatches:', int(total_mm))
    print('  % Mismatches:', total_mm / grand_total)
    print('\n')


    ## Iterate the Loop Counter ##
    q = q + 1              

##################
# END WHILE LOOP #
##################


print('\n*** DONE ***\n')