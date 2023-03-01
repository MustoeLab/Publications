###########################################
#       Profile Correlation Plotter       #
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
import scipy.stats


def filterProfilesByNt(profile1, profile2, nts=None, name=None, exclhighbg=None):
    """Return matched reactivities from ReactivityProfile objects prof1 and prof2,
    filtering out any NaNs
    If nts argument is None, return all nts. Otherwise, return only these nts   
    name argument is passed to ReactivityProfile to get desired profile    
    """
    
    # lets check to make sure the sequences are the same!
    if not np.array_equal(profile1.sequence, profile2.sequence):
        raise ValueError('Sequences are not the same!')
        
    
    # initiliaze my mask (all false)
    mask = np.zeros(len(profile1.sequence), dtype=bool)
    
    if nts is None:
        mask[:] = True #reassign all values to True
    
    else:
        for n in nts:
            mask = mask | (profile1.sequence == n)
    
    if exclhighbg is not None:
        with np.errstate(invalid='ignore'):
            mask = mask & (profile1.backprofile < exclhighbg) & (profile2.backprofile < exclhighbg)
    
    
    # get the desired profiles
    r1 = profile1.profile(name)     
    r2 = profile2.profile(name)
    
    mask = mask & np.isfinite(r1) & np.isfinite(r2)
    
    return r1[mask], r2[mask]


def plotCorrelation(var1, var2, ax, colr, colr2, title='', xlabel='', ylabel=''):
    
    regress = scipy.stats.linregress(var1, var2)
    
    # get min and max values for plotting
    xmin, ymin = min(var1), min(var2)
    gmin = min(xmin, ymin)
    xmax, ymax = max(var1), max(var2)
    gmax = max(xmax, ymax)
    
    ax.scatter(var1,var2, color=colr, s=10)
    
    # plot the diagonal
    ax.plot([gmin, gmax], [gmin, gmax], 'k--', label='diagonal')

    # plot the fit
    x = np.linspace(xmin, xmax, num=100)
    y = x*regress.slope + regress.intercept
    #ax.plot(x, y, 'r', label='fit', color=colr2)
    
    # add labels and stuff
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_ylim(0,1)
    ax.set_xlim(0,1)
    ax.set_xscale('symlog', linthresh=.0001, subs=[1,2,3,4,5,6,7,8,9])
    ax.set_yscale('symlog', linthresh=.0001, subs=[1,2,3,4,5,6,7,8,9])



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
print('\nGenerating a total of ' + str(loopend) + ' DMS reactivity profile correlation plots.\n')


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
    # atext3[2] = Exclude High Background Value      #
    # atext3[3] = Profile TXT file A Rep 1           #
    # atext3[4] = Profile TXT file A Rep 2           #
    # atext3[5] = Profile TXT file B Rep 1           #
    # atext3[6] = Profile TXT file B Rep 2           #
    ##################################################

    ## Conformation of Action ##
    print('Creating Correlation plot for ' + atext3[1] + '...')

    ## Define Variables ##
    profileA1 = atext3[3]
    profileA2 = atext3[4]
    profileB1 = atext3[5]
    profileB2 = atext3[6]

    highbgfilter = float(atext3[2])
    outPath = atext3[0] + '/CorrelationPlot_' + atext3[1] + '.pdf'

    ## Profile File Objects ##
    pA1 = ReactivityProfile(profileA1)
    pA1.normalize(DMS=True)
    pA2 = ReactivityProfile(profileA2)
    pA2.normalize(DMS=True)

    pB1 = ReactivityProfile(profileA1)
    pB1.normalize(DMS=True)
    pB2 = ReactivityProfile(profileA2)
    pB2.normalize(DMS=True)
	
    fig, ax = plot.subplots(figsize=(5.1,5))

    # This loop will compute correlations for all nts combined, and then individually
    for j, nb in enumerate(('A', 'C', 'U', 'G')):
        dataA1, dataA2 = filterProfilesByNt(pA1, pA2, nts=nb, name='sub', exclhighbg=highbgfilter)
        dataB1, dataB2 = filterProfilesByNt(pA1, pA2, nts=nb, name='sub', exclhighbg=highbgfilter)
        
        # Calculate average values for replicates #


        color_dict = {'A' : 'royalblue', 'C' : 'limegreen', 'G' : 'darkgoldenrod', 'U' : 'violet'}
        coloruse = color_dict.get(nb)
        color2_dict = {'A' : 'mediumblue', 'C' : 'darkgreen', 'G' : 'darkgoldenrod', 'U' : 'purple'}
        coloruse2 = color2_dict.get(nb)
        # Call function to create correlation plot
        plotCorrelation(datax, datay, ax, coloruse, coloruse2, title='sub')
	# label the axes
    ax.set_xlabel(profile1.split('/')[-1])
    ax.set_ylabel(profile2.split('/')[-1])

    fig.tight_layout()
    fig.savefig(outPath, dpi=100, bbox_inches="tight", transparent=True)


    ## Iterate the Loop Counter ##
    q = q + 1              

##################
# END WHILE LOOP #
##################

print('\n*** DONE ***\n')