##############################
#   AUROC and  ESI Plotter   #
#        E. coli rRNA        #
#                            #
#     David Mitchell III     #
#          (c) 2023          #
##############################

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plot
import sys
import RNAtools2 as RNAtools
import numpy as np
from sklearn.metrics import roc_curve, roc_auc_score
import xml.etree.ElementTree as ET

sys.path.append('/Users/anthonymustoe/Dropbox/Code/RingMapper')
from ReactivityProfile import ReactivityProfile

from sklearn.neighbors import KernelDensity
from scipy.stats import gamma
import gammamix


#matplotlib.rcParams['xtick.major.size'] = 8
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['font.sans-serif'] = 'Arial'




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
    prof.normprofile = np.array(data)

    return prof


def background_subtract_xml(xml1, xml2, normalize=True):

    p1 = readXML(xml1)
    p2 = readXML(xml2)
    
    p1.rawprofile = p1.normprofile
    p1.backprofile = p2.normprofile

    mask = (p1.rawprofile == 0) | (p1.backprofile==0)
    p1.rawprofile[mask] = np.nan
    p1.backprofile[mask] = np.nan

    p1.backgroundSubtract(normalize=False)
    

    if normalize:
        print(p1.normMethod)
        print(p1.normalize())
    else:
        p1.normprofile = p1.subprofile


    return p1


def _si_(r, p, u):

    si = 0
    for x in range(len(r)):
        if p[x] > 0:
            xp = p[x]/(p[x]+u[x])
            si += r[x]*xp*np.log2(xp)

        if u[x] > 0:
            xu = u[x]/(p[x]+u[x])
            si += r[x]*xu*np.log2(xu)

    return 1+si


def compute_histograms(profile, ct, nts, bins, minval=0, maxval=2):
    
    # make sure we are comparing the right things!
    if len(profile.sequence) != len(ct.ct):
        raise IndexError

    allr = []
    pair = []
    unpaired =[] 

    for i,v in enumerate(profile.normprofile):
        if v>-10 and profile.sequence[i] in nts:
            
            if v <= minval:
                v = minval+0.0001
            elif v >=maxval:
                v = maxval-0.0001

            allr.append(v)
            if ct.ct[i] != 0:
                pair.append(v)
            else:
                unpaired.append(v)
            
    # compute histograms
    hr, hr_edges = np.histogram(allr, range=(minval,maxval), bins=bins)
    hr = np.array(hr, dtype=float) / np.sum(hr)
    
    pr, pr_edges = np.histogram(pair, range=(minval,maxval), bins=bins)
    pr = np.array(pr, dtype=float) / np.sum(pr)

    ur, ur_edges = np.histogram(unpaired, range=(minval,maxval), bins=bins)
    ur = np.array(ur, dtype=float) / np.sum(ur)

    return hr,pr,ur,hr_edges



def compute_si_hist(profile, ct, nts, bins=50, **kwargs):
    
    si = 0
    for n in nts:
        r,p,u,edge = compute_histograms(profile, ct, n, bins)
        
        f = float(sum(profile.sequence==n)) / len(profile.sequence)
        si += f * _si_(r,p,u)

    return si



def _getnts_(profile, ct, n, noneg=False):
    """get list of paired/unpaired nucleotides of given nt type"""
    pair = []
    unpaired =[]
    for i,v in enumerate(profile.normprofile):
        if v>-10 and profile.sequence[i] in n:
            
            if v<=0 and noneg:
                # For gamma dist fitting, values need to be >0. If all the same value, leads to 
                # fitting instability, so inject some randomness
                v = 0.0001+np.random.random()*0.01

            if ct.ct[i] != 0:
                pair.append(v)
            else:
                unpaired.append(v)

    return np.array(pair), np.array(unpaired)



def compute_si_kde(profile, ct, nts):
    """Estimated reactivity PDFs using a Gaussian Kernel Density estimator and use these for computing SI"""

    si = []
    for n in 'ACGU':
        
        p,u = _getnts_(profile, ct, n)
   
        # if we don't have data, assume SI=0 at those nts. This is for when we compute
        # SI of just DMS A/C
        if n not in nts:
            si.extend([0]*(len(p)+len(u)))
            continue

        p,u = _getnts_(profile, ct, n)
        
        p_kde = KernelDensity(kernel='gaussian',bandwidth=0.1).fit(p[:,np.newaxis])
        u_kde = KernelDensity(kernel='gaussian',bandwidth=0.1).fit(u[:,np.newaxis])

        alldata = np.concatenate((p,u))
        
        # get the probabilities
        p_prob = np.exp(p_kde.score_samples( alldata[:,np.newaxis]))
        u_prob = np.exp(u_kde.score_samples( alldata[:,np.newaxis]))
        
        tmp = []
        for i in range(len(alldata)):
            l = p_prob[i]/(p_prob[i]+u_prob[i])
            if 0<l<1:
                v = 1+l*np.log2(l)+(1-l)*np.log2(1-l)
            else:
                v = 1
            
            tmp.append(v)

        print(n, np.mean(tmp))
        si.extend(tmp)

    return np.mean(si)


def compute_si_gamma(profile, ct, nts):
    """Estimate reactivity PDF using gamma mixture model."""

    def gmixpdf(x, params):
        return params.mix_prop[0]*gamma.pdf(x,a=params.alpha[0], scale=params.invbeta[0]) \
                + params.mix_prop[1]*gamma.pdf(x,a=params.alpha[1], scale=params.invbeta[1])


    si = []
    for n in 'ACUG':
        pair,unpair = _getnts_(profile, ct, n, noneg=True)
        
        # if we don't have data, assume SI=0 at those nts. This is for when we compute
        # SI of just DMS A/C
        if n not in nts:
            si.extend([0]*(len(pair)+len(unpair)))
            continue


        pair_gamma = gammamix.gammamix_em(pair, epsilon=1e-3, mix_prop=np.array([0.5, 0.5]))
        unpair_gamma = gammamix.gammamix_em(unpair, epsilon=1e-3, mix_prop=np.array([0.5, 0.5]))
        
        tmp = []
        for r in np.concatenate((pair, unpair)):

            l = gmixpdf(r, pair_gamma.params) / (gmixpdf(r, pair_gamma.params) + gmixpdf(r, unpair_gamma.params))
            if 0<l<1:
                v = 1+l*np.log2(l)+(1-l)*np.log2(1-l)
            else:
                v = 1
        
            tmp.append(v)
        
        print(n, np.mean(tmp))
        si.extend(tmp)

    return np.mean(si)





def compute_auc(profile, ct, nts, **kwargs):
    
    if len(profile.sequence) != len(ct.ct):
        raise IndexError

    react = []
    binary_pair = []

    for i,v in enumerate(profile.normprofile):               
        # only consider nts of desired type and ensure that reactivity is defined
        if v > -10 and profile.sequence[i] in nts:
            
            react.append(v)

            if ct.ct[i] != 0:
                binary_pair.append(0)
            else:
                binary_pair.append(1)
            
    return roc_auc_score(binary_pair, react)


def mergeprofile(plist, xml=False):

    m = ReactivityProfile()
    m.sequence = np.array([])
    m.normprofile = np.array([])

    for f in plist:
        
        if xml:
            p = readXML(f)
        else:
            p = ReactivityProfile(f)
        m.sequence = np.append(m.sequence, p.sequence)
        m.normprofile = np.append(m.normprofile, p.normprofile)
    
    return m



def get_stat(func, prof, ct, nts, **kwargs):
    
    if isinstance(prof, list):
        val = []
        for prof_i in prof:
            val.append(func(prof_i, ct, nts, **kwargs))
        
        return np.mean(val), np.std(val)

    else:
        return func(prof, ct, nts, **kwargs), None



# read in the data!
rprofile = {}
rprofile['ecrrna_ic'] = {}
rprofile['ecrrna_ic']['DMS'] = mergeprofile(['Profiles_for_ESI/Ecoli_rRNA-ic/Parsed_DM_rRNA_ic_Mara_Rep1_22dev-DMS_ec16S_profile.txt', 
                                             'Profiles_for_ESI/Ecoli_rRNA-ic/Parsed_DM_rRNA_ic_Mara_Rep1_22dev-DMS_ec23S_profile.txt'])
rprofile['ecrrna_ic']['2A3x'] = mergeprofile(['Profiles_for_ESI/Ecoli_rRNA-ic/DInc_Ecoli_rRNA-ic_2A3_16S.xml', 
                                              'Profiles_for_ESI/Ecoli_rRNA-ic/DInc_Ecoli_rRNA-ic_2A3_23S.xml'], True)
rprofile['ecrrna_ic']['NAIx'] = mergeprofile(['Profiles_for_ESI/Ecoli_rRNA-ic/DInc_Ecoli_rRNA-ic_NAI_16S.xml', 
                                              'Profiles_for_ESI/Ecoli_rRNA-ic/DInc_Ecoli_rRNA-ic_NAI_23S.xml'], True)



rprofile['ecrrna_cf'] = {}
rprofile['ecrrna_cf']['DMS'] = [mergeprofile(['Profiles_for_ESI/Ecoli_rRNA-cf/Parsed_DM_rRNA_cf_Mara_Rep1_22dev-DMS_ec16S_profile.txt', 
                                              'Profiles_for_ESI/Ecoli_rRNA-cf/Parsed_DM_rRNA_cf_Mara_Rep1_22dev-DMS_ec23S_profile.txt']),
                                mergeprofile(['Profiles_for_ESI/Ecoli_rRNA-cf/Parsed_JC_rRNA_cf_Mara_Rep1_22dev-DMS_ec16S_profile.txt', 
                                              'Profiles_for_ESI/Ecoli_rRNA-cf/Parsed_JC_rRNA_cf_Mara_Rep1_22dev-DMS_ec23S_profile.txt'])]

rprofile['ecrrna_cf']['1M7'] = mergeprofile(['Profiles_for_ESI/Ecoli_rRNA-cf/Ecoli_cellfree_16S_1M7_profile.txt', 
                                             'Profiles_for_ESI/Ecoli_rRNA-cf/Ecoli_cellfree_23S_1M7_profile.txt'])
rprofile['ecrrna_cf']['2A3x'] = mergeprofile(['Profiles_for_ESI/Ecoli_rRNA-cf/DInc_Ecoli_rRNA-cf_2A3_16S.xml', 
                                              'Profiles_for_ESI/Ecoli_rRNA-cf/DInc_Ecoli_rRNA-cf_2A3_23S.xml'], True)
rprofile['ecrrna_cf']['NAIx'] = mergeprofile(['Profiles_for_ESI/Ecoli_rRNA-cf/DInc_Ecoli_rRNA-cf_NAI_16S.xml', 
                                              'Profiles_for_ESI/Ecoli_rRNA-cf/DInc_Ecoli_rRNA-cf_NAI_23S.xml'], True)



# adjust the files from Danny's lab to align with our files
# 2A3
rprofile['ecrrna_cf']['2A3x'].sequence = np.insert(rprofile['ecrrna_cf']['2A3x'].sequence, 2286, 'G') 
rprofile['ecrrna_cf']['2A3x'].normprofile = np.insert(rprofile['ecrrna_cf']['2A3x'].normprofile, 2286, np.nan) 
rprofile['ecrrna_ic']['2A3x'].sequence = np.insert(rprofile['ecrrna_ic']['2A3x'].sequence, 2286, 'G') 
rprofile['ecrrna_ic']['2A3x'].normprofile = np.insert(rprofile['ecrrna_ic']['2A3x'].normprofile, 2286, np.nan) 
# NAI
rprofile['ecrrna_cf']['NAIx'].sequence = np.insert(rprofile['ecrrna_cf']['NAIx'].sequence, 2286, 'G') 
rprofile['ecrrna_cf']['NAIx'].normprofile = np.insert(rprofile['ecrrna_cf']['NAIx'].normprofile, 2286, np.nan) 
rprofile['ecrrna_ic']['NAIx'].sequence = np.insert(rprofile['ecrrna_ic']['NAIx'].sequence, 2286, 'G') 
rprofile['ecrrna_ic']['NAIx'].normprofile = np.insert(rprofile['ecrrna_ic']['NAIx'].normprofile, 2286, np.nan) 




# read in the CT files
CTs = {}
CTs['ecrrna'] = RNAtools.CT('Files/d.16.b.E.coli.nop.mask.ct', filterNC=True, filterSingle=True)
tmp = RNAtools.CT('Files/d.23.b.E.coli.nop.mask.ct', filterNC=True, filterSingle=True)
CTs['ecrrna'].ct.extend(tmp.ct)


# make the AUC, SI plots

fig,(ax1,ax2) = plot.subplots(2)

spacer = 10

#colors = ['tan', 'salmon','firebrick']
colors = ['tan', 'firebrick']
#colors = ['cornflowerblue', 'firebrick']

colors_2a3 = ['0.2', '0.5', 'k']
colors_nai = ['0.4', '0.7', 'k']
colors_1m7 = ['0.6', '0.9', 'k']

for r_index, rna in enumerate(('ecrrna_cf', 'ecrrna_ic')):
   
    prof = rprofile[rna]
    
    for n_index, nts in enumerate(('AC', 'ACGU')):

        x = spacer*r_index + n_index
        
        label=None
        if r_index == 0:
            label = 'DMS '+nts
        
        m,e = get_stat(compute_auc, prof['DMS'], CTs[rna[:-3]], nts)
        ax1.bar(x, m, yerr=e, color=colors[n_index], label=label)
        
        m,e = get_stat(compute_si_gamma, prof['DMS'], CTs[rna[:-3]], nts)
        ax2.bar(x, m, yerr=e, color=colors[n_index])
    
    
    # increment x from what it was in above loop
    x += 1

    #2A3    
    reagent = None
    if '2A3x' in prof:
        reagent = '2A3x'
    elif '2A3' in prof:
        reagent = '2A3'
    else:
        continue

    label=None
    if r_index == 0:
        label='2A3 '+nts

    m1,e = get_stat(compute_auc, prof[reagent], CTs[rna[:-3]], 'ACGU')
    ax1.bar(x, m1, yerr=e, color=colors_2a3[n_index], label=label)

    m2,e = get_stat(compute_si_gamma, prof[reagent], CTs[rna[:-3]], 'ACGU')
    ax2.bar(x, m2, yerr=e, color=colors_2a3[n_index])

    x += 1

    # NAI
    reagent = None
    if 'NAIx' in prof:
        reagent = 'NAIx'
    elif 'NAI' in prof:
        reagent = 'NAI'
    else:
        continue

    label=None
    if r_index == 1:
        label='NAI '+nts

    m1,e = get_stat(compute_auc, prof[reagent], CTs[rna[:-3]], 'ACGU')
    ax1.bar(x, m1, yerr=e, color=colors_nai[n_index], label=label)

    m2,e = get_stat(compute_si_gamma, prof[reagent], CTs[rna[:-3]], 'ACGU')
    ax2.bar(x, m2, yerr=e, color=colors_nai[n_index])

    x += 1

    # 1M7
    reagent = None
    if '1M7x' in prof:
        reagent = '1M7x'
    elif '1M7' in prof:
        reagent = '1M7'
    else:
        continue

    label=None
    if r_index == 2:
        label='1M7 '+nts


    m1,e = get_stat(compute_auc, prof[reagent], CTs[rna[:-3]], 'ACGU')
    ax1.bar(x, m1, yerr=e, color=colors_1m7[n_index], label=label)

    m2,e = get_stat(compute_si_gamma, prof[reagent], CTs[rna[:-3]], 'ACGU')
    ax2.bar(x, m2, yerr=e, color=colors_1m7[n_index])

    x += 1
    
    print('{} {} {}'.format(rna, m1, m2))


# format the plot
ax1.set_ylim(0.5,1.0)
ax1.set_xticks([])
ax1.set_ylabel('AUROC')
ax1.legend()

ax2.set_xticks(np.arange(0,r_index*spacer+1, spacer)+1.5)
ax2.set_xticklabels(['rRNA\nCell-free','rRNA\nIn-cell',])
ax2.set_ylabel('ESI (bits)')

fig.savefig('Output/dm_FigS4_ESI_AUC_Ecoli_rRNA_GammaMix.pdf')


#################################################################################################
# Make the LL plot

fig,(ax1,ax2) = plot.subplots(2, figsize=(4,8))

r,p,u,edge = compute_histograms(rprofile['ecrrna_cf']['2A3x'], CTs['ecrrna'], 'ACGU', 10, maxval=1)

# extend the last bin for plotting purposes
hist = np.zeros(edge.shape)
hist[:-1] = p/u
hist[-1] = hist[-2]
rhist = np.zeros(edge.shape)
rhist[:-1] = r
rhist[-1] = rhist[-2]

#ax.step(edge, hist, label='2A3',where='post')
ax1.fill_between(edge, hist, np.ones(edge.shape), label='2A3',step='post', color='k', alpha=0.5)

#ax.step(edge, hist, label='2A3',where='post')
ax2.fill_between(edge, rhist, np.zeros(edge.shape),  label='2A3',step='post', color='k', alpha=0.5)


colors={'U':'orangered', 'G':'maroon', 'C':'darkorange', 'A':'darkgoldenrod'}

for n in 'ACUG':
    r,p,u,edge = compute_histograms(rprofile['ecrrna_cf']['DMS'][0], CTs['ecrrna'], n, 10, maxval=1)

    hist = np.zeros(edge.shape)
    hist[:-1] = p/u
    hist[-1] = hist[-2]
    #hist[np.isinf(hist)] = 50 #np.nan
    hist[hist==0] = np.nan

    rhist = np.zeros(edge.shape)
    rhist[:-1] = r
    rhist[-1] = rhist[-2]

    ax1.step(edge, hist, label='DMS '+n, where='post', color=colors[n], lw=2)
    ax2.step(edge, rhist, label='DMS '+n, where='post', color=colors[n], lw=2)



ax1.plot((-0.05,1.05),(1,1), color='k', linestyle='--')
ax1.set_xlim(-0.05, 1.05)
ax1.set_yscale('log')
ax1.set_ylim(0.01,10)
ax1.set_ylabel('Likelihood ratio, unpaired / paired')
ax1.set_xlabel('Normalized reactivity')

ax1.legend()


fig.savefig('Output/dm_FigS4_ESI_Histogram_Ecoli_rRNA_GammaMix.pdf')





