#!/usr/bin/env python

import numpy as np
import pickle, os, sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

modules_dir = '/afs/ifh.de/user/b/blotsd/scratch/Analysis/OscFit/releases/oscFitND_v1.0/modules'
# modules_dir = os.path.dirname(os.path.realpath(__file__)).split('/resources')[0] + '/modules'
sys.path.append(modules_dir)
import jp_mpl as jplot


####
## User defined variables
###

# Andrii wrote down a PublicDataReader class that converts the eff. area pickle files into histograms
# The code is currently not linked to oscFit. You can get it here:
# http://code.icecube.wisc.edu/svn/sandbox/terliuk/PublicDataReader/PublicDataReader.py
# You're going to need an icetray with a new-ish version of NuFlux as well.

sys.path.append('/afs/ifh.de/user/y/yanezjua/i3scripts/public_data/PublicDataReader')
import PublicDataReader as PubDataReader

oscfit_user_name = 'PRD_extended_lowdt_data'
indir = '/lustre/fs19/group/icecube/deepcore_analyses/effective_areas'

outdir = os.path.join(indir, oscfit_user_name, 'validation')
indir  = os.path.join(indir, oscfit_user_name, 'pckl')

if not os.path.exists(outdir):
    os.makedirs(outdir)


# Run the thing
from icecube import NuFlux, dataclasses
myflux =  NuFlux.makeFlux('IPhonda2014_spl_solmin')

import dataLoader
loader =  dataLoader.dataLoader(bin_edges   = [10**np.linspace(0.8,1.75,9),
                                               np.arccos(np.linspace(-1,0.,9))[::-1]],
                                observables = ['reco_energy', 'reco_zenith'],
                                user = 'PRD_extended_lowdt_data',
                                LEaxis = [],
                                weight_keys = ['weight_e', 'weight_mu'],
                                verbose = False,
                                legacy_detsys    = False,
                                detailed_detsys = False,
                                extra_cuts = {'vtxz': [-500., -250.]},
                                break_energy = 195.
                            )

# Add some oscillations
dm31 = 0.0025
dm21 = 7.53E-5
theta23 = 0.78
theta12 = np.arcsin(np.sqrt(0.846))/2.
theta13 = 0.15

#dm31 = dm21 = theta23 = theta13 = theta12 = 0.

settings = {'dm31':                      dm31, 
             'theta23':                  theta23, 
             'theta13':                  theta13,
             'mix_angle':                0., 
             'norm_e':                   1.,
             'norm':                     2.3,
             'gamma':                    0.0, 
             'domeff':                   1., 
             'hole_ice':                 0.02, 
             'pi2k':                     1., 
             'rel_qe':                   0., 
             'atmmu_template':           'data',
             'atmmu_fraction':           0.0, 
             'simulation':               'baseline', 
             'oscMode':                  'Vacuum',
             'add_detector_systematics': False}

# Converting everything to 950 days of livetime
mchist = loader.getNeutrinoHistograms(settings)
for one_key in mchist.keys():
    mchist[one_key]*= 3600.*24*950.

tot_events = 0
for one_key in mchist.keys():
    tot_events += np.sum(mchist[one_key])
print tot_events

import oscCalc
reload(oscCalc)
probCalc = oscCalc.OscCalc(oscMode = 'Vacuum',
                           doTables = False)
probCalc.setParameters(dm31=dm31, theta23=theta23, theta13=theta13)
pdreader = PubDataReader.PublicDataReader(indir, 
                                          skipTaus=False)

print 'Bins in true variables ', len(pdreader.true_cosZ_edges), len(pdreader.true_logE_edges)

settings = {'energy_oversampling':1,
            'coszen_oversampling':1,
            'lifetime': 950.,
            'flux_calc': myflux,
            'oscCalc': probCalc}
aehist = pdreader.GetHistogramFromAreas(settings)

norm = tot_events / np.sum(aehist['all'])
print 'Difference in total number of events', norm

# Do the plot per flavor
figs = []
colors = ['b','r','c','g']
fig = plt.figure(figsize=(9,18))
recoz_edges = aehist['reco_cosZ_edges']
recoe_edges = aehist['reco_logE_edges']

plotlims = [0.9, 1.1]

counter = 1
for eband in range(aehist['nue_histo'].shape[0]-1):
    
    fig.add_subplot(aehist['nue_histo'].shape[0], 2, counter)
    for i, hkey in enumerate(mchist.keys()):
        hratio = norm*aehist[hkey]/mchist[hkey][:,::-1]
        jplot.unfilledBar(recoz_edges,hratio[eband,:], color = colors[i], label = hkey)
    counter += 1
    plt.ylim(plotlims)
    plt.xlabel('cos(zenith)')
    
    fig.add_subplot(aehist['nue_histo'].shape[1], 2, counter)
    for i, hkey in enumerate(mchist.keys()):
        hratio = norm*aehist[hkey]/mchist[hkey][:,::-1]
        jplot.unfilledBar(recoe_edges,hratio[:,eband], color = colors[i], label = hkey)
    counter += 1
    #plt.legend(loc=0)
    plt.ylim(plotlims)
    plt.xlabel('log(E)')

plt.subplots_adjust(wspace = 0.3, hspace =0.7)
fig.savefig(os.path.join(outdir, 'Ratio_perFlavor.png'))

mcfull = np.zeros_like(aehist['all'])
for i, hkey in enumerate(mchist.keys()):
    mcfull += mchist[hkey]
hratio = norm*aehist['all']/mcfull[:,::-1]


# Do the plot now over the full sample
colors = ['b','r','c','g']
fig = plt.figure(figsize=(9,18))
recoz_edges = aehist['reco_cosZ_edges']
recoe_edges = aehist['reco_logE_edges']

plotlims = [0.9, 1.1]

counter = 1
for eband in range(aehist['nue_histo'].shape[0]-1):
    
    fig.add_subplot(aehist['nue_histo'].shape[0], 2, counter)
    jplot.unfilledBar(recoz_edges,hratio[eband,:], color = colors[i], label = hkey)
    counter += 1
    plt.ylim(plotlims)
    plt.xlabel('cos(zenith)')
    
    fig.add_subplot(aehist['nue_histo'].shape[1], 2, counter)
    jplot.unfilledBar(recoe_edges,hratio[:,eband], color = colors[i], label = hkey)
    counter += 1
    #plt.legend(loc=0)
    plt.ylim(plotlims)
    plt.xlabel('log(E)')

plt.subplots_adjust(wspace = 0.3, hspace =0.7)
fig.savefig(os.path.join(outdir, 'Ratio_neutrinoSum.png'))


# Plotting the deviations
x =np.linspace(0.97, 1.03, 21)
print 'Mean deviation ', hratio.mean()
print 'Std of deviations ', hratio.std()
errors2, x = np.histogram(hratio.flatten(), x)
fig = plt.figure(figsize=(8,4))
fig.add_subplot(121)
plt.plot(hratio.flatten(), 'xr')
plt.ylabel('Deviation')
plt.xlabel('Analysis bin')
fig.add_subplot(122)
jplot.unfilledBar(x, errors2, color = 'red')
plt.ylabel('Number of bins')
plt.xlabel('Deviation')
fig.savefig(os.path.join(outdir, 'Deviations_summary.png'))
