#!/usr/bin/env python

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import jp_mpl as jplot
import numpy as np
import pickle, sys, os
from itertools import cycle
# user_dir = os.path.dirname(os.path.realpath(__file__)).split('/resources')[0] + '/user_settings'
user_dir = '/afs/ifh.de/user/b/blotsd/scratch/Analysis/OscFit/releases/oscFitND_v1.0/user_settings'
sys.path.append(user_dir)

###
# User defined variables
###

'''
This script expects your pickle files to have a key with the value:
(Earth_survival_probability * MC_one_weight) / total_mcevents_generated

Where the denominator is given by:
neutrinos:      Nevents_per_file  * total_files * fraction_of_neutrinos_generated
antineutrinos:  Nevents_per_file  * total_files * fraction_of_antineutrinos_generated
'''

# Key in pickle files with the value defined above
weight_key = 'oneweight'

# oscFit user name
oscfit_user_name = 'PRD_extended_lowdt_data'

# Additional cuts used when loading the events (copy verbatim from dataLoader)
extra_cuts = {'vtxz': [-500., -250.]}

# Pick which eff. flavor you want redone (-1:all, 1:nue, 2:numu, 3:nutau, 4:nc)
do_flavor = -1

# Directory to store the output (it will be outdir/oscfit_user/)
outdir = '/lustre/fs19/group/icecube/deepcore_analyses/effective_areas'

# How many bins do you want in true logE & cosZ? Need "enough" to reproduce results.
trueE_bins =  51
trueZ_bins =  31

# Do you want to produce plots as you go?
doPlots = True

# Include all NC in calculation? This produces better eff. areas, but you need extra MC details
nc_use_all_files = True
if nc_use_all_files:
    # If you don't know these numbers, set all_nc to False
    # This is an array of the number of events produced (nFiles*nEvents) for [NuE, NuMu, NuTau]
    nfiles = np.array([2700.*300, 4000.*750, 1400.*60])


###
# Configuration
###

outdir  = os.path.join(outdir, oscfit_user_name)
pckldir = os.path.join(outdir, 'pckl')
if not os.path.exists(outdir):
    os.makedirs(outdir)
if not os.path.exists(pckldir):
    os.makedirs(pckldir)
if doPlots:
    outdir_figs = os.path.join(outdir, 'figs')
    if not os.path.exists(outdir_figs):
        os.makedirs(outdir_figs)

oscfit_user  = __import__(oscfit_user_name)

# Histogram of observables
logE_edges = np.linspace(0.8, 1.75, 9)
cosZ_edges = np.linspace(-1, 0, 9)
bin_edges  = np.meshgrid(logE_edges, cosZ_edges)

trueZ_edges   = np.linspace(-1, 1. , trueZ_bins)
trueZ_centers = (trueZ_edges[:-1] + trueZ_edges[1:])/2.
trueE_edges   = np.linspace(np.log10(3.), 3., trueE_bins)

solid_angle   = 2*np.pi*(trueZ_edges[1:] - trueZ_edges[:-1])
energy_diff   = 10**trueE_edges[1:] - 10**trueE_edges[:-1]
effarea_factor= energy_diff.reshape(len(energy_diff),1) * solid_angle


eff_area_data = {'true_logE_edges':trueE_edges,
                 'true_cosZ_edges':trueZ_edges,
                 'reco_logE_edges':logE_edges,
                 'reco_cosZ_edges':cosZ_edges}

# This 10^9 is changing ns to s - it is required
factor = (10**-9)


# Load the 3 pickle files required - produce the fake file for NC interactions
print '\nLoading all the files'
flavors  = ['nue','numu','nutau']
infiles = {}
copy_keys = ['reco_energy','reco_zenith','energy','zenith', 'interaction', 'ptype'] + extra_cuts.keys()
for i, one_flavor in enumerate(flavors):
    infiles[one_flavor] = pickle.load( open(oscfit_user.genie_p1[one_flavor] +
                                            oscfit_user.mc_sets['baseline']['baseline'] + 
                                            oscfit_user.genie_p2))
    if nc_use_all_files:
        # Start the dictionary
        if not infiles.has_key('nc'):
            infiles['nc'] = {}

        nc_factor = nfiles[i]/nfiles.sum()
        for one_key in infiles[one_flavor].keys():

            # Initialize the keys
            if not infiles['nc'].has_key(one_key):
                infiles['nc'][one_key] = np.zeros(1)

            # Need to account for having more files produced (only for weight_key)
            if one_key == weight_key:
                infiles['nc'][one_key] = np.concatenate(( infiles['nc'][one_key], 
                                                          infiles[one_flavor][one_key]*nc_factor ))
            elif one_key in copy_keys:
                infiles['nc'][one_key] = np.concatenate(( infiles['nc'][one_key], 
                                                          infiles[one_flavor][one_key] ))
# Alternatively, use NuMu for NC eff. areas
if not nc_use_all_files:
    infiles['nc'] = infiles['numu']

print '\nFinished loading the files - Calculating effective areas'


# Keys for eff areas creation
interactions  = flavors + ['nc']
flabels       = ['CC_NuE','CC_NuMu','CC_NuTau','NC_NuX']
ccnc          = [1,1,1,2]
nu_nubar_list = [1, -1]
nu_nubar_str  = ['','Bar']

if do_flavor > 0:
    interactions = [interactions[do_flavor-1]]
    flabels      = [flabels[do_flavor-1]]
    ccnc         = [ccnc[do_flavor-1]]

# Loop over the interactions
for index in range(len(interactions)):
    data       = infiles[interactions[index]]
    inter_bool = data['interaction'] == ccnc[index]

    # Get a bool array of the extra cuts already here
    cuts_bool = np.array([True]*len(data['interaction']))
    for cut_var in extra_cuts.keys():
        cuts_bool *= ((data[cut_var] > extra_cuts[cut_var][0]) * 
                      (data[cut_var] < extra_cuts[cut_var][1]))

    # Multiplying times the interaction already here
    inter_bool *= cuts_bool

    # Loop over neutrino/antineutrion
    for k, nu_nubar in enumerate(nu_nubar_list):
        nubar_bool = np.sign(data['ptype']) == nu_nubar

        flabel = flabels[index] + nu_nubar_str[k]
        outfile_name = os.path.join(pckldir, flabel)
        eff_area_data['particle']  = flabel
        eff_area_data['eff_areas'] = []

        # Loop over the reconstructed cosZ bins
        for cosZ in range(len(cosZ_edges)-1):
            cosz_bool = (np.cos(data['reco_zenith']) >= cosZ_edges[cosZ])*\
                        (np.cos(data['reco_zenith']) < cosZ_edges[cosZ+1])

            # Loop over the reconstructed logE bins
            for logE in range(len(logE_edges)-1):
                loge_bool = (np.log10(data['reco_energy']) >= logE_edges[logE])*\
                            (np.log10(data['reco_energy']) < logE_edges[logE+1])

                # All the events for this interaction (CC/NC), this (anti)neutrino, and this reco bin (E,z)
                sbool   = inter_bool*nubar_bool*cosz_bool*loge_bool
                weights = data[weight_key][sbool]
                logw    = np.log10(weights)

                # Removing extreme outliers
                wbool = logw < (logw.mean() + 3*logw.std())

                bins, x, y = np.histogram2d(np.log10(data['energy'][sbool][wbool]), 
                                            np.cos(data['zenith'][sbool][wbool]),
                                            [trueE_edges, trueZ_edges],
                                            weights = weights[wbool])
                
                effarea_onebin = (bins*factor/effarea_factor).T

                # Showing the progress
                print '\n\n******************\nDone with: ', flabel
                print 'cos(reco_zenith)', cosZ_edges[cosZ], cosZ_edges[cosZ+1]
                print 'logE(reco_E)', logE_edges[logE], logE_edges[logE+1]

                # This is a histogram of the bin weights
                if doPlots and len(logw < 10) > 0:

                    fig = plt.figure(figsize=(12,6))
                    fig_title = flabel + '_' + \
                        "%0.2f" % cosZ_edges[cosZ] + '_' +  "%0.2f" % cosZ_edges[cosZ+1] + '___' + \
                        "%0.2f" % logE_edges[logE] + '_' +  "%0.2f" % logE_edges[logE+1]

                    # Distribution of event weights
                    fig.add_subplot(1,2,1)
                    xhist = np.linspace(0, 10, 101)
                    xhist[-1] = np.inf
                    bhist, x = np.histogram(logw, xhist)
                    jplot.unfilledBar(xhist, bhist)
                    plt.axvline(x=logw.mean()+2*logw.std(), lw=2, color='r', linestyle='-')
                    plt.axvline(x=logw.mean()+3*logw.std(), lw=2, color='r', linestyle='--')
                    plt.ylabel('N events')
                    plt.xlabel('weight (a.u.)')
                    plt.yscale('log')

                    # Effective area for one bin
                    fig.add_subplot(1,2,2)
                    plt.pcolor(trueE_edges, trueZ_edges, effarea_onebin, cmap='Blues')
                    plt.xlim([trueE_edges.min(), trueE_edges.max()])
                    plt.axvspan(logE_edges[logE], logE_edges[logE+1], lw=2,
                                facecolor='0.3', alpha=0.2)
                    plt.axhspan(cosZ_edges[cosZ], cosZ_edges[cosZ+1], lw=2,
                                facecolor='0.3', alpha=0.2)

                    plt.ylabel('cos(true zenith)')
                    plt.xlabel('log10(true E/GeV)')
                    plt.colorbar()
                    plt.suptitle(fig_title)

                    plt.show()

                    fig.subplots_adjust(left=0.1, wspace=0.3, right=0.95)
                    fig.savefig(os.path.join(outdir_figs, fig_title +'.png'))
                    plt.close(fig)

                bin_data = {'logEreco':[logE_edges[logE], logE_edges[logE+1]],
                            'cosZreco':[cosZ_edges[cosZ], cosZ_edges[cosZ+1]],
                            'effArea': effarea_onebin}
                eff_area_data['eff_areas'].append(bin_data)



        pickle.dump(eff_area_data, open(outfile_name+'.pckl', 'w'))






