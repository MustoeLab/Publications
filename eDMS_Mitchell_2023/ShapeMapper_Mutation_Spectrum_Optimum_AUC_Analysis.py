###########################################
#        AUC Matrix Python Program        #
#                                         #
#           David Mitchell III            #
#                (c) 2022                 #
###########################################

import sys
import numpy as np
import RNAtools2 as RNAtools
from sklearn.metrics import roc_curve, roc_auc_score
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plot
from matplotlib.colors import LinearSegmentedColormap
from ReactivityProfile import ReactivityProfile


## Function to calculate the background-subtracted mutation rate for each mutation type and for each nucleotide position ##
def getMutationDifference(dms_bb, unt_bb, len_bb):

    with np.errstate(invalid='ignore'):

        # Creates an array removing the five rightmost columns (keeping only mutation types)
        dms_cc = np.array(dms_bb[:,:-5])
        unt_cc = np.array(unt_bb[:,:-5])
    
        # Creates an array containing only the effective reads column
        # Numbers changed to account for new mutation types from Shapemapper.
        #dms_dd = np.array(dms_bb[:,28:29])
        #unt_dd = np.array(unt_bb[:,28:29])
        # Try a more Shapemapper version agnostic mode of catching the effective reads column
        dms_dd = np.array(dms_bb[:,len_bb:len_bb+1])
        unt_dd = np.array(unt_bb[:,len_bb:len_bb+1])
    
        # Define variables for mutation count sum calculation
        dms_cc_sum = [[] for q in range(dms_cc.shape[0])]
        unt_cc_sum = [[] for q in range(unt_cc.shape[0])]
    
        # Calculate the mutation count sum for each nucleotide position
        for q in range(dms_cc.shape[0]):
            dms_cc_sum[q].append(np.sum(dms_cc[q]))
            unt_cc_sum[q].append(np.sum(unt_cc[q]))
    
        # Change *_cc_sum type from a list into an array
        dms_cc_sum = np.array(dms_cc_sum)
        unt_cc_sum = np.array(unt_cc_sum)
    
        # Define variables for calculating corrected effective reads and rates per nucleotide position and mutation type
        dms_ee = [[] for y in range(dms_cc.shape[1])]       # Array for DMS-treated sample corrected effective reads
        dms_ff = [[] for n in range(dms_cc.shape[1])]       # Array for DMS-treated sample mutation rate
        unt_ee = [[] for y in range(unt_cc.shape[1])]       # Array for Untreated sample corrected effective reads
        unt_ff = [[] for n in range(unt_cc.shape[1])]       # Array for Untreated sample mutation rate
    
        # Perform the following calculation per nucleotide position (e.g. per row):
        # New effective reads [for given mutation type] = Old effective reads 
        #                                               - sum of mutation counts for nucleotide position 
        #                                               + mutation countss at position for desired mutation type
        # This subtracts from the effective read count all the mutation counts except those for the given mutation type
        # In this example, it will be dd - cc_sum + cc
    
        for x in range(dms_cc.shape[0]):
            for y in range(dms_cc.shape[1]):
                dms_ee[y].append(dms_dd[x] - dms_cc_sum[x] + dms_cc[x,y])
                unt_ee[y].append(unt_dd[x] - unt_cc_sum[x] + unt_cc[x,y])
    
        # Change *_ee type from a list into an array and transpose to restore correct array dimensions
        dms_gg = np.transpose(np.array(dms_ee))
        dms_gg = np.reshape(dms_gg, (dms_cc.shape[0],dms_cc.shape[1]))
        unt_gg = np.transpose(np.array(unt_ee))
        unt_gg = np.reshape(unt_gg, (unt_cc.shape[0],unt_cc.shape[1]))
    
        # Calculate the mutation rate for each nucleotide position and each mutation type
        for m in range(dms_cc.shape[0]):
            for n in range(dms_cc.shape[1]):
                dms_ff[n].append(dms_cc[m,n] / dms_gg[m,n])
                unt_ff[n].append(unt_cc[m,n] / unt_gg[m,n])
    
        # Change *_ff type from a list into an array and transpose to restore correct array dimensions
        dms_hh = np.transpose(np.array(dms_ff))
        unt_hh = np.transpose(np.array(unt_ff))
    
        # Check *_hh to find NaN from divide-by-zero and replace with 0 in the mutation rate
        for w in range(dms_cc.shape[0]):
            for z in range(dms_cc.shape[1]):
                if np.isnan(dms_hh[w,z]) == True:
                    dms_hh[w,z] = 0
                if np.isnan(unt_hh[w,z]) == True:
                    unt_hh[w,z] = 0
                else:
                    dms_hh[z] = dms_hh[z]
                    unt_hh[z] = unt_hh[z]
        
        # Defiine variables for calculating the difference in mutation rates between DMS-modified and untreated samples
        mutation_diff = [[] for x in range(dms_cc.shape[1])]
    
        # Calculate difference in mutation rates between DMS-modified and untreated samples for each mutation type and nucleotide position
        for x in range(dms_cc.shape[0]):
            for y in range(dms_cc.shape[1]):
                mutation_diff[y].append(dms_hh[x,y] - unt_hh[x,y])
    
    return mutation_diff


## Function to create data for AUC score calculation from original profile TXT file ##
def getNtLists1D(profile, pairedlist, ntype):

    react = []
    ispaired = []

    # .subprofle = Background-subtracted profile as defined by ReactivityProfile.py
    # enumerate gives a list of tuples -- e.g. [(0, AAA), (1, BBB), (2, CCC)]
    # Here, the second value in the tuple is the background-subtracted mutation rate
    for i,v in enumerate(profile.subprofile):               
        # Variable 'i' is the index and 'v' is the value within the index
        if v > -10 and profile.sequence[i] == ntype:
            react.append(v)
            # Looks through the CT file and appends value of '0' if given nucleotide is single-stranded
            # Otherwise, it appends '1' if the nucleotide is base-paired
            # this will append 0/1 if nt is in pairedlist or not
            ispaired.append( int((i+1) in pairedlist) )
        
    return react, ispaired


## Function to create data for AUC score calculations from 2D arrays ##
def getNtLists2D(mutations_diff, pairedlist, sequence, ntype):

    react = [[] for x in range(len(mutations_diff))]
    ispaired = [[] for x in range(len(mutations_diff))]

    for x in range(len(mutations_diff)):
        
        # enumerate gives a list of tuples -- e.g. [(0, AAA), (1, BBB), (2, CCC)]
        # Here, the second value in the tuple is the background-subtracted mutation rate
        for y,v in enumerate(mutations_diff[x]):               
            # Variable 'y' is the index and 'v' is the value within the index
            if v > -10 and sequence[y] == ntype:
                react[x].append(v)
                # Looks through the CT file and appends value of '0' if given nucleotide is single-stranded
                # Otherwise, it appends '1' if the nucleotide is base-paired
                ispaired[x].append( int((y+1) in pairedlist) )
        
    return react, ispaired


## Function to create data for AUC score calculations from grouped mutation types list ##
def getNtListsGroupMutType(mutations_diff, pairedlist, sequence, ntype):

    react = []
    ispaired = []

    # enumerate gives a list of tuples -- e.g. [(0, AAA), (1, BBB), (2, CCC)]
    # Here, the second value in the tuple is the background-subtracted mutation rate
    for x,v in enumerate(mutations_diff):               
        # Variable 'y' is the index and 'v' is the value within the index
        if v > -10 and sequence[x] == ntype:
            react.append(v)
            # Looks through the CT file and appends value of '0' if given nucleotide is single-stranded
            # Otherwise, it appends '1' if the nucleotide is base-paired
            ispaired.append( int((x+1) in pairedlist) )
        
    return react, ispaired


## Function to calculate the AUC scores ##
def getAUCscores2D(mutations_diff, true_pos_2d, score_2d):
    
    auc_scores = []
    
    for x in range(len(mutations_diff)):
        auc_scores.append(roc_auc_score(true_pos_2d[x], score_2d[x]))
        
    return auc_scores


## Function to check the AUC scores for each mutation type against the aggregate AUC score ##
def getAUCbool(auc_array, data_len):

    auc_bool_scores = [[] for x in range(len(auc_array))]
    
    for x in range(len(auc_array)):
        for y in range(len(auc_array[x])):
            if auc_array[x,y] >= auc_array[x,data_len]:
                auc_bool_scores[x].append(1)
            else:
                auc_bool_scores[x].append(0)

    return auc_bool_scores


## Function to deterime which mutation types to bundle and calculate AUC score as a group ##
def getAUCcompare(auc_array, header, mutation_difference, min_auc, diff_auc, len_dd):

    auc_compare_indices = [[] for x in range(len(auc_array))]
    header_indices = [[] for x in range(len(auc_array))]
    mutation_indices = [[] for x in range(len(auc_array))]

    # Variable 'min_auc' gives the minimum AUC score for an individual mutation type
    # Variable 'diff_auc' gives the minimum delta between the individual mutation type AUC score and aggregate AUC score

    for x in range(len(auc_array)):
        for y in range(len(auc_array[x])):
            if auc_array[x,y] > float(min_auc) and np.subtract(auc_array[x,y], auc_array[x,len_dd]) > float(diff_auc):
                auc_compare_indices[x].append(y)
                header_indices[x].append(header[y])
                mutation_indices[x].append(mutation_difference[y])

    return auc_compare_indices, header_indices, mutation_indices


## Function to create new data array for getMutationDifference() containing only filtered mutation types from getAUCcompare() ##
def getMakeNewArray(dms_data, unt_data, filter_indices, len_cc):

    # Define variables
    dms_effective_reads = []
    unt_effective_reads = []
    dms_data_filtered = []
    unt_data_filtered = []
    index_reject = []
    dms_mutations_sum = []
    unt_mutations_sum = []
    dms_full_sum = []
    unt_full_sum = []
    dms_filter_sum = []
    unt_filter_sum = []

    # Trim data from original files
    dms_data_trimmed = np.array(dms_data[:,:-5])
    unt_data_trimmed = np.array(unt_data[:,:-5])

    # Calculate overall sum of mutation counts for each nucleotide position
    for x in range(dms_data.shape[0]):
        dms_mutations_sum.append(np.sum(dms_data[x]))
        unt_mutations_sum.append(np.sum(unt_data[x]))

    # Create lists for DMS-modified and untreated effective read depths
    for x in range(dms_data.shape[0]):
        dms_effective_reads.append(dms_data[x,len_cc])
        unt_effective_reads.append(unt_data[x,len_cc])

    # Variable 'list_length' to give lengths in 'for' loops
    list_length = len(dms_effective_reads)

    dms_data_use = np.transpose(dms_data)
    unt_data_use = np.transpose(unt_data)

    # Identify indices for mutation types which passed selection criteria
    for q in range(len(filter_indices)):
        if filter_indices[q] == []:
            index_reject.append(q)
    
    # Escape function if no mutation types pass selection criteria
    if len(index_reject) == 4:
        bgsub_mutation_rate = None
        index_reject = None
        filter_indices = []
    else:
        bgsub_mutation_rate = []

    
    ## Continue only if there are mutation types which pass selection criteria 
    if bgsub_mutation_rate == None:
        ## Stopping below
        print('ERROR')

    else:
        ## Continuation below
        # Filter the original data based on which mutation types passed selection
        for x in range(len(filter_indices)):
            if x not in index_reject:
                for y in range(dms_data.shape[1]):
                    if y in filter_indices[x]:
                        dms_data_filtered.append(dms_data_use[y])
                        unt_data_filtered.append(unt_data_use[y])

        # Convert from list type to array and transpose to restore mutation types to columns 
        dms_data_filtered_array = np.transpose(np.array(dms_data_filtered))
        unt_data_filtered_array = np.transpose(np.array(unt_data_filtered))

        # Calculate sum of mutation counts for all mutation types (and change type to array)
        for x in range(dms_data_trimmed.shape[0]):
            dms_full_sum.append(np.sum(dms_data_trimmed[x]))
            unt_full_sum.append(np.sum(unt_data_trimmed[x]))

        # Calculate sum of mutation counts for filtered mutation types (and change type to array)
        for x in range(dms_data_filtered_array.shape[0]):
            dms_filter_sum.append(np.sum(dms_data_filtered_array[x]))
            unt_filter_sum.append(np.sum(unt_data_filtered_array[x]))

        # Calculate corrected effective reads
        dms_corr_effreads = []
        unt_corr_effreads = []

        for x in range(dms_data_filtered_array.shape[0]):
            dms_corr_effreads.append(dms_effective_reads[x] - dms_full_sum[x] + dms_filter_sum[x])
            unt_corr_effreads.append(unt_effective_reads[x] - unt_full_sum[x] + unt_filter_sum[x])

        # Calculate the mutation rate and background subtract #
        with np.errstate(invalid='ignore'):

            # Define variables
            dms_mutation_rate = []
            unt_mutation_rate = []

            # Perform calculation
            for x in range(list_length):
                dms_mutation_rate.append(dms_filter_sum[x] / dms_corr_effreads[x])
                unt_mutation_rate.append(unt_filter_sum[x] / unt_corr_effreads[x])

            # Change NaN values from divide-by-zero error to '0' in array
            for x in range(list_length):
                if np.isnan(dms_mutation_rate[x]) == True:
                    dms_mutation_rate[x] = 0
                if np.isnan(unt_mutation_rate[x]) == True:
                    unt_mutation_rate[x] = 0
                else:
                    dms_mutation_rate[x] = dms_mutation_rate[x]
                    unt_mutation_rate[x] = unt_mutation_rate[x]

            # Calculate background-subtracted mutation rate
            for x in range(list_length):
                bgsub_mutation_rate.append(dms_mutation_rate[x] - unt_mutation_rate[x])        

    return bgsub_mutation_rate, index_reject, filter_indices

#################################################################################################################

## Open and Read User Arguments from Text File ##

afile = open(sys.argv[1])             # Open the file with the file name given by the user
atext = afile.read()                  # Read the contents of the file
atext1 = atext.split('\n')            # Split the file into a list where each list element is a complete line from the file
atextlen = len(atext1)                # Give count of number of lines in original file


## User Confirmation ##
print('\nOpening file ' + sys.argv[1] + ' to generate ROC plots.\n')

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
print('\nCalculating AUC scores for ' + str(loopend) + ' items.\n')


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
    # atext3[1] = Modified mutation counts file (.txt)  #
    # atext3[2] = Untreated mutation counts file (.txt) #
    # atext3[3] = Profile file (.txt)                   #
    # atext3[4] = Minimum mutation type AUC score       #
    # atext3[5] = Minimum AUC score delta from aggregate#
    # atext3[6] = Output directory                      #
    # atext3[7] = Output file name                      #
    #####################################################

    ## Open mutation counts files and read data into variables ##
    cmfile1 = open(atext3[1])                                                               # Open DMS-modified mutation counts file
    cmhead = cmfile1.readline().split()                                                     # Pulls the first line from the mutation counts file as a header
    cmheader = cmhead[0:(len(cmhead)-5)]                                                    # Only all but last 5 entries of header (mutation types)
    dms_file_data = np.genfromtxt(cmfile1, delimiter='\t', skip_header=0, dtype=int)        # Obtain data from remainder of the file
    cmfile2 = open(atext3[2])                                                               # Open untreated mutation counts file
    unt_file_data = np.genfromtxt(cmfile2, delimiter='\t', skip_header=1, dtype=int)        # Obtain data from remainder of the file
    
    # Appends to the header 'Origninal AUC score' to allow title for mutation-type agnostic AUC score in spectrum plot
    cmheader.append('Original AUC score')
    data_len_a = len(cmheader)                                                              # Allows use of arbitrary number of mutation types (for future Shapemapper variants)
    data_len_b = data_len_a - 1


    ## Operations to calculate the mutation rate and mutation rate difference from complete data ##
    mut_diff = getMutationDifference(dms_file_data, unt_file_data, data_len_b)    


    ## Open profile file and obtain nucleotide sequence from first column ##
    prof_file = open(atext3[3])
    prof_seq = np.genfromtxt(prof_file, delimiter='\t', skip_header=1, dtype=str)
    prof_seq_data = np.array(prof_seq[:,1:2])
    
    ## Open the CT file and run RNAtools2 on it to generate the paired residue list ##
    ctobj = RNAtools.CT(atext3[0])
    # Residue added to list if value in column 4 (which gives the pairing partner for a residue) > 0
    pairedlist = ctobj.pairedResidueList(False)

    # Runs ReactivityProfile.py on profile TXT file
    profobj = ReactivityProfile(atext3[3])

    ## Creating Graphical Compilation of ROC Plots Using MatPlotLib ##
    labels = ('A', 'C', 'G', 'U')  
    
    ## Operations to calculate the AUC scores for each nucleotide ##
    
    auc_value_full = []
    auc_value_full_1d = []
    auc_value_groupmut_full = []

    for i,nt in enumerate(labels):
        # Setting variables for ROC plot generation by SciKit-Learn (sklearn.metrics)
        # 'r' = Profile file background-subtracted mutation rates
        # 'p' = Paired residue list from CT file
        r1, p1 = getNtLists1D(profobj, pairedlist, nt)
        r2, p2 = getNtLists2D(mut_diff, pairedlist, prof_seq_data, nt)

        # Calculating the AUC score per mutation type
        auc_value_per_mt = getAUCscores2D(mut_diff, p2, r2)
        auc_value_full.append(auc_value_per_mt)

        auc_value_1d = roc_auc_score(p1, r1)
        auc_value_full_1d.append(auc_value_1d)

    # Appending the mutation-type-agnostic AUC scores to the overall AUC score list
    for x in range(len(auc_value_full)):
        auc_value_full[x].append(auc_value_full_1d[x])

    # Changing compilated AUC values from list to array
    auc_value_full_array = np.array(auc_value_full)

    ## Calculating relative AUC scores (AUC for specific mutation type / AUC for mutation agnostic) ##
    auc_relative = [[] for x in range(auc_value_full_array.shape[0])]
    
    for x in range(auc_value_full_array.shape[0]):
        for y in range(auc_value_full_array.shape[1]):
            auc_relative[x].append(auc_value_full_array[x, y] / auc_value_full_array[x, data_len_b])
    
    # Creating array for relative AUC scores
    auc_relative_array = np.array(auc_relative)

    # Creating array of header
    cmheader_array = str(np.array(cmheader))


    ## Boolean check of per-mutation-type AUC score against aggregate AUC score ##
    auc_value_bool = getAUCbool(auc_value_full_array, data_len_b)
    auc_value_bool_array = np.array(auc_value_bool)


    ## Comparison check of per-mutation-type AUC score against aggregate AUC score ##
    auc_index_compare, cmheader_new, mut_diff_new = getAUCcompare(auc_value_full_array, cmheader, mut_diff, atext3[4], atext3[5], data_len_b)


    ## ROC data saved as text file ##
    outPath_ROCdata = atext3[6] + '/AUC_grROCdata_' + atext3[7] + '.txt'
    roc_datafile = open(outPath_ROCdata, 'w')

    roc_datafile.write('Calculating AUC scores for mutation types which pass the selection criteria for ' + atext3[7] + '.....\n')
    roc_datafile.write('Selection criteria are as follows:\n')
    roc_datafile.write('   Minimum individual mutation type AUC score: ' + str(atext3[4]) + '\n')
    roc_datafile.write('   Minimum delta between individual mutation type AUC score and aggregate AUC score:' + str(atext3[5]) + '\n\n')
    roc_datafile.write('The following mutation types were selected for A nucleotides: ' + str(cmheader_new[0]) + '\n')
    roc_datafile.write('The following mutation types were selected for C nucleotides: ' + str(cmheader_new[1]) + '\n')
    roc_datafile.write('The following mutation types were selected for G nucleotides: ' + str(cmheader_new[2]) + '\n')
    roc_datafile.write('The following mutation types were selected for U nucleotides: ' + str(cmheader_new[3]) + '\n\n')

    # Print output on screen for user
    print('Calculating AUC scores for mutation types which pass the selection criteria for ' + atext3[7] + '.....')
    print('Selection criteria are as follows:')
    print('   Minimum individual mutation type AUC score:', atext3[4])
    print('   Minimum delta between individual mutation type AUC score and aggregate AUC score:', atext3[5])
    print('')
    print('The following mutation types were selected for A nucleotides:', cmheader_new[0])
    print('The following mutation types were selected for C nucleotides:', cmheader_new[1])
    print('The following mutation types were selected for G nucleotides:', cmheader_new[2])
    print('The following mutation types were selected for U nucleotides:', cmheader_new[3], '\n')


    ## Get AUC scores for the filtered lists ##
    sub_mutrate, reject_list, number_muttypes = getMakeNewArray(dms_file_data, unt_file_data, auc_index_compare, data_len_b)

    # End if none of the mutation types pass muster
    if sub_mutrate == None:
        print('None of the mutation types passed the selection criteria for any nucleotides.')
        roc_datafile.write('None of the mutation types passed the selection criteria for any nucleotides.')


    ## Obtain nucleotide list data from grouped mutation types for calculating ROC curve ##
    fig_roc,ax_roc = plot.subplots(1, 4, figsize=(8,2))

    for i, nt in enumerate(labels):
        if sub_mutrate != None:
            r3, p3 = getNtListsGroupMutType(sub_mutrate, pairedlist, prof_seq_data, nt)
            r1a, p1a = getNtLists1D(profobj, pairedlist, nt)

        else:
            r3 = None
            p3 = None
    
        if reject_list != None and i not in reject_list and r3 != None:
            auc_value_groupmut = roc_auc_score(p3, r3)
            auc_value_groupmut_full.append(auc_value_groupmut)
            print('The AUC score for the combined grouping of mutation types for nucleotide ' + nt + ' is: ', auc_value_groupmut)
            roc_datafile.write('The AUC score for the combined grouping of mutation types for nucleotide ' + nt + ' is: ' + str(auc_value_groupmut) + '\n')

            # Create ROC plot for groupings of mutation types
            fpr, tpr, thr = roc_curve(p3, r3, drop_intermediate=False)
            fpr_a, tpr_a, thr_a = roc_curve(p1a, r1a, drop_intermediate=False)
        
            ax_roc[i].plot(fpr, tpr, color='green', lw=3)
            ax_roc[i].plot(fpr_a, tpr_a, color='orange', lw=3)
            ax_roc[i].plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
            ax_roc[i].set_xlim([-0.02, 1.0])
            ax_roc[i].set_ylim([0.0, 1.02])
            ax_roc[i].set_title(nt + ' (Mut = ' + str(len(number_muttypes[i])) + ')')
            ax_roc[i].text(1,0.01,'New AUC={:.2f}'.format(roc_auc_score(p3,r3)) + '\nOld AUC={:.2f}'.format(roc_auc_score(p1a,r1a)), horizontalalignment='right', transform=ax_roc[i].transAxes)

    print('\n. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .\n')
    roc_datafile.close()
    
    fig_roc.suptitle('ROC Plot: ' + atext3[7] + '; Min AUC=' + atext3[4] + '; delta AUC=' + atext3[5])
    fig_roc.tight_layout()


    ## Creating pcolormesh plots for raw AUC scores and relative scores ##
    # Raw AUC scores plots
    fig1,ax1 = plot.subplots()
    fig1.colorbar(ax1.pcolormesh(auc_value_full_array, cmap='RdYlBu_r', vmin=0, vmax=1))
    ax1.pcolormesh(auc_value_full_array, cmap='RdYlBu_r', vmin=0, vmax=1, edgecolors='k', linewidths=1)
    ax1.set_xticks(np.arange(0, data_len_a))
    ax1.set_yticks(np.arange(0, 4))
    ax1.set_yticklabels(('A', 'C', 'G', 'U'))
    ax1.set_xticklabels(cmheader, rotation=90)
    ax1.set_title('AUC matrix for ' + atext3[7] + '\nMin AUC=' + atext3[4] + '; delta AUC=' + atext3[5])

    # Relative AUC scores plots
    fig2,ax2 = plot.subplots()
    fig2.colorbar(ax2.pcolormesh(auc_relative_array, cmap='RdYlBu_r', vmin=0, vmax=2))    
    ax2.pcolormesh(auc_relative_array, cmap='RdYlBu_r', vmin=0, vmax=2, edgecolors='k', linewidths=1)
    ax2.set_xticks(np.arange(0, data_len_a))
    ax2.set_yticks(np.arange(0, 4))
    ax2.set_yticklabels(('A', 'C', 'G', 'U'))
    ax2.set_xticklabels(cmheader, rotation=90)
    ax2.set_title('Relative AUC matrix for ' + atext3[7] + '\nMin AUC=' + atext3[4] + '; delta AUC=' + atext3[5])

    # Boolean AUC scores plots
    fig3,ax3 = plot.subplots()
    fig3.colorbar(ax3.pcolormesh(auc_value_bool_array, cmap='Greens', vmin=0, vmax=1))    
    ax3.pcolormesh(auc_value_bool_array, cmap='Greens', vmin=0, vmax=1, edgecolors='k', linewidths=1)
    ax3.set_xticks(np.arange(0, data_len_a))
    ax3.set_yticks(np.arange(0, 4))
    ax3.set_yticklabels(('A', 'C', 'G', 'U'))
    ax3.set_xticklabels(cmheader, rotation=90)
    ax3.set_title('Boolean AUC matrix for ' + atext3[7] + '\nMin AUC=' + atext3[4] + '; delta AUC=' + atext3[5])


    ## Creating Output PDF Files ##
    # Raw AUC scores matrix plot PDF
    print('\nCreating pcolormesh matrix plot for raw AUC scores at output file: AUC_rawplot_' + atext3[7] + '.pdf\n\n................................................................\n')
    outPath1 = atext3[6] + '/AUC_rawplot_' + atext3[7] + '.pdf'
    fig1.savefig(outPath1, dpi=100, bbox_inches="tight", transparent=True)
    
    # Relative AUC scores matrix plot PDF
    print('\nCreating pcolormesh matrix plot for relative AUC scores at output file: AUC_Relplot_' + atext3[7] + '.pdf\n\n................................................................\n')
    outPath2 = atext3[6] + '/AUC_Relplot_' + atext3[7] + '.pdf'
    fig2.savefig(outPath2, dpi=100, bbox_inches="tight", transparent=True)

    # Relative AUC scores matrix plot PDF
    print('\nCreating pcolormesh matrix plot for relative AUC scores at output file: AUC_Boolplot_' + atext3[7] + '.pdf\n\n................................................................\n')
    outPath3 = atext3[6] + '/AUC_Boolplot_' + atext3[7] + '.pdf'
    fig3.savefig(outPath3, dpi=100, bbox_inches="tight", transparent=True)

    # ROC plot from groupings of mutation types PDF
    print('\nCreating ROC plot for groupings of selected mutation types at output file: AUC_grROCplot_' + atext3[7] + '.pdf\n\n................................................................\n')
    outPath_ROC = atext3[6] + '/AUC_grROCplot_' + atext3[7] + '.pdf'
    fig_roc.savefig(outPath_ROC, dpi=100, bbox_inches="tight", transparent=True)


    # AUC values saved as text file
    outPath4 = atext3[6] + '/AUC_rawdata_' + atext3[7] + '.txt'
    np.savetxt(outPath4, auc_value_full_array, header=atext3[7])
    outPath5 = atext3[6] + '/AUC_Reldata_' + atext3[7] + '.txt'
    np.savetxt(outPath5, auc_relative_array, header=atext3[7])
    outPath6 = atext3[6] + '/AUC_Booldata_' + atext3[7] + '.txt'
    np.savetxt(outPath6, auc_value_bool_array, header=atext3[7])


    ## Iterate Loop Counter ##
    j = j + 1

##################
# END WHILE LOOP #
##################

print('\n\nFinished with calculating AUC scores and creating plots.\n\n')
