#!/usr/bin/env python

import numpy as np
import pickle, os, sys

oscfit_user_name = 'PRD_extended_lowdt_data'
indir = '/lustre/fs19/group/icecube/deepcore_analyses/effective_areas'

outdir = os.path.join(indir, oscfit_user_name, 'txt')
indir  = os.path.join(indir, oscfit_user_name, 'pckl')

if not os.path.exists(outdir):
    os.makedirs(outdir)

filenames = os.listdir(indir)

# Load all of the maps
data = {}
for one_file in filenames:
    this_path = os.path.join(indir, one_file)
    if os.path.isdir(this_path):
        continue
    flavor = one_file.split('.')[0]
    data[flavor] = pickle.load(open(this_path))

for one_key in data.keys():
    txtfile = open(os.path.join(outdir, one_key + '.txt'), 'w')
    txtfile.write(one_key +'\n')
    for rbin in range(len(data[one_key]['eff_areas'])):
        bin_name = 'logEreco = [' + \
            "%.2f" % data[one_key]['eff_areas'][rbin]['logEreco'][0] + ' - '+ \
            "%.2f" % data[one_key]['eff_areas'][rbin]['logEreco'][1] + '], '+ \
            'cosZreco = [' + \
            "%.2f" % data[one_key]['eff_areas'][rbin]['cosZreco'][0] + '_'+ \
            "%.2f" % data[one_key]['eff_areas'][rbin]['cosZreco'][1] + ']\n'
        txtfile.write(bin_name)
        txtfile.write('cos(zenith)\tlog10(E/GeV)\tAeff (m^2)\n')


        for one_cosz in range(len(data[one_key]['true_cosZ_edges'])-1):
            for energy_band in range(len(data[one_key]['true_logE_edges'])-1):
                txtfile.write("[" + "%.3f" % data[one_key]['true_cosZ_edges'][one_cosz] + ',' +\
                              "%.3f" % data[one_key]['true_cosZ_edges'][one_cosz+1] + ']\t')
                txtfile.write("[" + "%.3f" % data[one_key]['true_logE_edges'][energy_band] + ',' +\
                          "%.3f" % data[one_key]['true_logE_edges'][energy_band+1] + ']\t')
                txtfile.write("%.3e" % data[one_key]['eff_areas'][rbin]['effArea'][one_cosz, energy_band] + '\n')

            
        txtfile.write('\n')
    txtfile.close()
