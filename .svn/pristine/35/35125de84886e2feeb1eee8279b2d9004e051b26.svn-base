#!/usr/bin/env python

import os, sys

# modules_dir = os.path.dirname(os.path.realpath(__file__)).split('/resources')[0] + '/modules'
modules_dir = '/afs/ifh.de/user/b/blotsd/scratch/Analysis/OscFit/releases/oscFitND_v1.0/modules'
sys.path.append(modules_dir)

###
# User defined variables
###

# oscFit user name
oscfit_user_name = 'PRD_extended_lowdt_data'

# Directory to store the output (it will be outdir/oscfit_user/)
outdir = '/lustre/fs19/group/icecube/deepcore_analyses/effective_areas'
outdir  = os.path.join(outdir, oscfit_user_name, 'txt')



import numpy as np
import dataLoader

loader =  dataLoader.dataLoader(bin_edges   = [10**np.linspace(0.8,1.75,9),
                                               np.arccos(np.linspace(-1,0.,9))[::-1]],
                                observables = ['reco_energy', 'reco_zenith'],
                                user = oscfit_user_name, 
                                LEaxis = [],
                                weight_keys = ['weight_e', 'weight_mu'], 
                                verbose = False,
                                legacy_detsys    = False,
                                detailed_detsys  = False,
                                extra_cuts = {'vtxz': [-500., -250.]},
                                break_energy = 195.
                            )

zenith_bands = loader.bin_edges[1]
energy_axis  = loader.bin_edges[0]

# Doing the DATA
data_histogram = loader.loadData(year = 'all')
txt_outfile = os.path.join(outdir,'DataCounts.txt')
txtfile = open(txt_outfile, 'w')
txtfile.write('DATA\n')
txtfile.write('cos(zenith)\tlog10(E/GeV)\tCounts\n')

for k in range(len(zenith_bands)-1)[::-1]:
    for i in range(len(energy_axis)-1):
        txtfile.write('['+"%.2f" % np.cos(zenith_bands[k+1]) + ','+ "%.2f" % np.cos(zenith_bands[k]) + ']\t' +
                      '['+"%.2f" % np.log10(energy_axis[i]) + ','+ "%.2f" % np.log10(energy_axis[i+1]) + ']\t' +
                      "%i" % data_histogram[i,k] + '\n')

txtfile.close()

# Doing the ATMOSPHERIC MUONS (from data)
data_histogram = loader.atmmu_histo['data']
print 'Total number of muons: ', np.sum(data_histogram)

txt_outfile = os.path.join(outdir,'AtmMuons_fromData.txt')
txtfile = open(txt_outfile, 'w')
txtfile.write('Atm. muons - This template is left to float freely in the final fit.\n')
txtfile.write('cos(zenith)\tlog10(E/GeV)\tCounts\n')

for k in range(len(zenith_bands)-1)[::-1]:
    for i in range(len(energy_axis)-1):
        txtfile.write('['+"%.2f" % np.cos(zenith_bands[k+1]) + ','+ "%.2f" % np.cos(zenith_bands[k]) + ']\t' +
                      '['+"%.2f" % np.log10(energy_axis[i]) + ','+ "%.2f" % np.log10(energy_axis[i+1]) + ']\t' +
                      "%i" % data_histogram[i,k] + '\n')

txtfile.close()

