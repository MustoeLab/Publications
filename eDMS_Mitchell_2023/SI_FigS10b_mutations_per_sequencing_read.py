############################################
#      Plot Histogram from Per-Read        #
#          Mutation Counts Data            #
#                                          #
#           David Mitchell III             #
#                (c) 2023                  #
############################################


import sys
import numpy as np
import ReactivityProfile
import matplotlib.pyplot as plot


## Open and Read User Arguments from Text File ##

afile = open(sys.argv[1])             # Open the file with the file name given by the user
atext = afile.read()                  # Read the contents of the file
atext1 = atext.split('\n')            # Split the file into a list where each list element is a complete line from the file
atextlen = len(atext1)                # Give count of number of lines in original file


## User Confirmation ##
print('\nOpening file ' + sys.argv[1] + ' to generate mutation count plots.\n')

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

## User Confirmation of Total Number of Entries ##
print('\nGenerating a total of ' + str(loopend) + ' data points.\n')

# Define variables #
names = []
total_nmuts = []
tot_median_nmuts = []
tot_95per_nmuts = []

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
    # atext3[0] = Parsed Mutation File (*.mut)     #
    # atext3[1] = Profile File (*.txt)             #
    ################################################


    fname = atext3[0]
    profile = ReactivityProfile.ReactivityProfile(atext3[1])
    coverage = 0.9
    outname = 'Output/Per-Read_Mut_Counts_Histogram_Orig-Untreated_Med-95th_TEST.pdf'

    with np.errstate(invalid='ignore'):
        dlen = np.sum(profile.rawprofile > 0.0001)

    cov_cut = dlen * coverage

    names.append(atext3[0])

    nmuts = []
    with open(fname) as inp:
        for line in inp:
            spl = line.split()
            if spl[0] != 'MERGED' or spl[4]!='INCLUDED':
                continue
            if spl[7].count('1') < cov_cut:
                continue
            m = spl[8].count('1')
            nmuts.append(m)


    total_nmuts.append(nmuts)

    med_nmuts = np.median(nmuts)
    mean_nmuts = np.mean(nmuts)
    p95_nmuts = np.percentile(nmuts,95)
    max_nmuts = np.max(nmuts)
    len_nmuts = len(nmuts)

    tot_median_nmuts.append(med_nmuts)
    tot_95per_nmuts.append(p95_nmuts)

    print('Name:', atext3[0])
    print('Median mutations per read: {}'.format(med_nmuts))
    print('Mean mutations per read: {}'.format(mean_nmuts))
    print('95th percentile mutations per read: {}'.format(p95_nmuts))
    print('Maximum mutations per read: {}'.format(max_nmuts))
    print('Number of reads with mutations: {}'.format(len_nmuts))
    print('\n\n')

    j = j + 1           # Iterate the loop counter

##################
# END WHILE LOOP #
##################

#print('names:', names, '\n')
#print('total_nmuts:', len(total_nmuts))

rna_names = ('RMRP_rep1_Mara', 'RMRP_rep1_SS2', 'RMRP_rep1_TGIRT', 'RnaseP_rep1_Mara', 'RnaseP_rep1_SS2', 'RnaseP_rep1_TGIRT', 'tmRNA_R1_Mara', 'tmRNA_rep1_SS2', 'tmRNA_rep1_TGIRT',
            'RMRP_rep2_Mara', 'RMRP_rep2_SS2', 'RMRP_rep2_TGIRT', 'RnaseP_rep2_Mara', 'RnaseP_rep2_SS2', 'RnaseP_rep2_TGIRT', 'tmRNA_rep2_Mara', 'tmRNA_rep2_SS2', 'tmRNA_rep2_TGIRT')

for i, iname in enumerate(rna_names):
    if names[i].find('RMRP') >= 0 and names[i].find('Mara') >= 0 and names[i].find('rep1') >=0:
        ent_RMRP_rep1_Mara = total_nmuts[i]
        ent_RMRP_rep1_Mara_med = tot_median_nmuts[i]
        ent_RMRP_rep1_Mara_95 = tot_95per_nmuts[i]
    elif names[i].find('RMRP') >= 0 and names[i].find('SS2') >= 0 and names[i].find('rep1') >=0:
        ent_RMRP_rep1_SS2_med = tot_median_nmuts[i]
        ent_RMRP_rep1_SS2_95 = tot_95per_nmuts[i]
        ent_RMRP_rep1_SS2 = total_nmuts[i]
    elif names[i].find('RMRP') >= 0 and names[i].find('TGIRT') >= 0 and names[i].find('rep1') >=0:
        ent_RMRP_rep1_TGIRT_med = tot_median_nmuts[i]
        ent_RMRP_rep1_TGIRT_95 = tot_95per_nmuts[i]
        ent_RMRP_rep1_TGIRT = total_nmuts[i]
    elif names[i].find('RnaseP') >= 0 and names[i].find('Mara') >= 0 and names[i].find('rep1') >=0:
        ent_RnaseP_rep1_Mara_med = tot_median_nmuts[i]
        ent_RnaseP_rep1_Mara_95 = tot_95per_nmuts[i]
        ent_RnaseP_rep1_Mara = total_nmuts[i]
    elif names[i].find('RnaseP') >= 0 and names[i].find('SS2') >= 0 and names[i].find('rep1') >=0:
        ent_RnaseP_rep1_SS2_med = tot_median_nmuts[i]
        ent_RnaseP_rep1_SS2_95 = tot_95per_nmuts[i]
        ent_RnaseP_rep1_SS2 = total_nmuts[i]
    elif names[i].find('RnaseP') >= 0 and names[i].find('TGIRT') >= 0 and names[i].find('rep1') >=0:
        ent_RnaseP_rep1_TGIRT_med = tot_median_nmuts[i]
        ent_RnaseP_rep1_TGIRT_95 = tot_95per_nmuts[i]
        ent_RnaseP_rep1_TGIRT = total_nmuts[i]
    elif names[i].find('tmRNA') >= 0 and names[i].find('Mara') >= 0 and names[i].find('R1') >=0:
        ent_tmRNA_rep1_Mara_med = tot_median_nmuts[i]
        ent_tmRNA_rep1_Mara_95 = tot_95per_nmuts[i]
        ent_tmRNA_rep1_Mara = total_nmuts[i]
    elif names[i].find('tmRNA') >= 0 and names[i].find('SS2') >= 0 and names[i].find('R1') >=0:
        ent_tmRNA_rep1_SS2_med = tot_median_nmuts[i]
        ent_tmRNA_rep1_SS2_95 = tot_95per_nmuts[i]
        ent_tmRNA_rep1_SS2 = total_nmuts[i]
    elif names[i].find('tmRNA') >= 0 and names[i].find('TGIRT') >= 0 and names[i].find('R1') >=0:
        ent_tmRNA_rep1_TGIRT_med = tot_median_nmuts[i]
        ent_tmRNA_rep1_TGIRT_95 = tot_95per_nmuts[i]
        ent_tmRNA_rep1_TGIRT = total_nmuts[i]
    elif names[i].find('RMRP') >= 0 and names[i].find('Mara') >= 0 and names[i].find('rep2') >=0:
        ent_RMRP_rep2_Mara = total_nmuts[i]
        ent_RMRP_rep2_Mara_med = tot_median_nmuts[i]
        ent_RMRP_rep2_Mara_95 = tot_95per_nmuts[i]
    elif names[i].find('RMRP') >= 0 and names[i].find('SS2') >= 0 and names[i].find('rep2') >=0:
        ent_RMRP_rep2_SS2_med = tot_median_nmuts[i]
        ent_RMRP_rep2_SS2_95 = tot_95per_nmuts[i]
        ent_RMRP_rep2_SS2 = total_nmuts[i]
    elif names[i].find('RMRP') >= 0 and names[i].find('TGIRT') >= 0 and names[i].find('rep2') >=0:
        ent_RMRP_rep2_TGIRT_med = tot_median_nmuts[i]
        ent_RMRP_rep2_TGIRT_95 = tot_95per_nmuts[i]
        ent_RMRP_rep2_TGIRT = total_nmuts[i]
    elif names[i].find('RnaseP') >= 0 and names[i].find('Mara') >= 0 and names[i].find('rep2') >=0:
        ent_RnaseP_rep2_Mara_med = tot_median_nmuts[i]
        ent_RnaseP_rep2_Mara_95 = tot_95per_nmuts[i]
        ent_RnaseP_rep2_Mara = total_nmuts[i]
    elif names[i].find('RnaseP') >= 0 and names[i].find('SS2') >= 0 and names[i].find('rep2') >=0:
        ent_RnaseP_rep2_SS2_med = tot_median_nmuts[i]
        ent_RnaseP_rep2_SS2_95 = tot_95per_nmuts[i]
        ent_RnaseP_rep2_SS2 = total_nmuts[i]
    elif names[i].find('RnaseP') >= 0 and names[i].find('TGIRT') >= 0 and names[i].find('rep2') >=0:
        ent_RnaseP_rep2_TGIRT_med = tot_median_nmuts[i]
        ent_RnaseP_rep2_TGIRT_95 = tot_95per_nmuts[i]
        ent_RnaseP_rep2_TGIRT = total_nmuts[i]
    elif names[i].find('tmRNA') >= 0 and names[i].find('Mara') >= 0 and names[i].find('R2') >=0:
        ent_tmRNA_rep2_Mara_med = tot_median_nmuts[i]
        ent_tmRNA_rep2_Mara_95 = tot_95per_nmuts[i]
        ent_tmRNA_rep2_Mara = total_nmuts[i]
    elif names[i].find('tmRNA') >= 0 and names[i].find('SS2') >= 0 and names[i].find('R2') >=0:
        ent_tmRNA_rep2_SS2_med = tot_median_nmuts[i]
        ent_tmRNA_rep2_SS2_95 = tot_95per_nmuts[i]
        ent_tmRNA_rep2_SS2 = total_nmuts[i]
    elif names[i].find('tmRNA') >= 0 and names[i].find('TGIRT') >= 0 and names[i].find('R2') >=0:
        ent_tmRNA_rep2_TGIRT_med = tot_median_nmuts[i]
        ent_tmRNA_rep2_TGIRT_95 = tot_95per_nmuts[i]
        ent_tmRNA_rep2_TGIRT = total_nmuts[i]


plot.rc('xtick', labelsize=4)   # Change the font size of the x-axis text
plot.rc('ytick', labelsize=4)   # Change the font size of the y-axis text

# Plot 95th percentiles #
plot.subplot(2, 1, 2)
plot.scatter(0.45, ent_RMRP_rep1_SS2_95, s=35, color='limegreen')
plot.scatter(0.50, ent_RMRP_rep1_Mara_95, s=35, color='cornflowerblue')
plot.scatter(0.55, ent_RMRP_rep1_TGIRT_95, s=35, color='darkorange')
plot.scatter(0.95, ent_RnaseP_rep1_SS2_95, s=35, color='limegreen')
plot.scatter(1.00, ent_RnaseP_rep1_Mara_95, s=35, color='cornflowerblue')
plot.scatter(1.05, ent_RnaseP_rep1_TGIRT_95, s=35, color='darkorange')
plot.scatter(1.45, ent_tmRNA_rep1_SS2_95, s=35, color='limegreen')
plot.scatter(1.50, ent_tmRNA_rep1_Mara_95, s=35, color='cornflowerblue')
plot.scatter(1.55, ent_tmRNA_rep1_TGIRT_95, s=35, color='darkorange')
plot.scatter(0.45, ent_RMRP_rep2_SS2_95, s=35, color='limegreen')
plot.scatter(0.50, ent_RMRP_rep2_Mara_95, s=35, color='cornflowerblue')
plot.scatter(0.55, ent_RMRP_rep2_TGIRT_95, s=35, color='darkorange')
plot.scatter(0.95, ent_RnaseP_rep2_SS2_95, s=35, color='limegreen')
plot.scatter(1.00, ent_RnaseP_rep2_Mara_95, s=35, color='cornflowerblue')
plot.scatter(1.05, ent_RnaseP_rep2_TGIRT_95, s=35, color='darkorange')
plot.scatter(1.45, ent_tmRNA_rep2_SS2_95, s=35, color='limegreen')
plot.scatter(1.50, ent_tmRNA_rep2_Mara_95, s=35, color='cornflowerblue')
plot.scatter(1.55, ent_tmRNA_rep2_TGIRT_95, s=35, color='darkorange')
plot.title('95th Percentile Values', fontsize=4)
plot.ylabel('Mutation Rate', fontsize=4)
plot.ylim(-0.1, 12.5)

# Plot medians #
plot.subplot(2, 1, 1)
plot.scatter(0.45, ent_RMRP_rep1_SS2_med, s=35, color='limegreen')
plot.scatter(0.50, ent_RMRP_rep1_Mara_med, s=35, color='cornflowerblue')
plot.scatter(0.55, ent_RMRP_rep1_TGIRT_med, s=35, color='darkorange')
plot.scatter(0.95, ent_RnaseP_rep1_SS2_med, s=35, color='limegreen')
plot.scatter(1.00, ent_RnaseP_rep1_Mara_med, s=35, color='cornflowerblue')
plot.scatter(1.05, ent_RnaseP_rep1_TGIRT_med, s=35, color='darkorange')
plot.scatter(1.45, ent_tmRNA_rep1_SS2_med, s=35, color='limegreen')
plot.scatter(1.50, ent_tmRNA_rep1_Mara_med, s=35, color='cornflowerblue')
plot.scatter(1.55, ent_tmRNA_rep1_TGIRT_med, s=35, color='darkorange')
plot.scatter(0.45, ent_RMRP_rep2_SS2_med, s=35, color='limegreen')
plot.scatter(0.50, ent_RMRP_rep2_Mara_med, s=35, color='cornflowerblue')
plot.scatter(0.55, ent_RMRP_rep2_TGIRT_med, s=35, color='darkorange')
plot.scatter(0.95, ent_RnaseP_rep2_SS2_med, s=35, color='limegreen')
plot.scatter(1.00, ent_RnaseP_rep2_Mara_med, s=35, color='cornflowerblue')
plot.scatter(1.05, ent_RnaseP_rep2_TGIRT_med, s=35, color='darkorange')
plot.scatter(1.45, ent_tmRNA_rep2_SS2_med, s=35, color='limegreen')
plot.scatter(1.50, ent_tmRNA_rep2_Mara_med, s=35, color='cornflowerblue')
plot.scatter(1.55, ent_tmRNA_rep2_TGIRT_med, s=35, color='darkorange')

plot.title('Median Values', fontsize=4)
plot.ylabel('Mutation Rate', fontsize=4)
plot.ylim(-0.1, 8.5)

plot.tight_layout()

plot.savefig(outname, dpi=100, bbox_inches="tight", transparent=True)


