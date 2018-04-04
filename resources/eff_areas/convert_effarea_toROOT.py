#!/usr/bin/env python

import numpy as np
import pickle, os, sys
import ROOT

oscfit_user_name = 'PRD_extended_lowdt_data'
indir = '/lustre/fs19/group/icecube/deepcore_analyses/effective_areas'

outdir = os.path.join(indir, oscfit_user_name, 'root')
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

outfile = ROOT.TFile(os.path.join(outdir, 'IceCube_EffAreas.root'), 'new')
histos = []

for one_key in data.keys():
    for rbin in range(len(data[one_key]['eff_areas'])):

        bin_name = one_key + '_logEreco_' + \
            "%.2f" % data[one_key]['eff_areas'][rbin]['logEreco'][0] + '_'+ \
            "%.2f" % data[one_key]['eff_areas'][rbin]['logEreco'][1] + '_'+ \
            'cosZreco_' + \
            "%.2f" % data[one_key]['eff_areas'][rbin]['cosZreco'][0] + '_'+ \
            "%.2f" % data[one_key]['eff_areas'][rbin]['cosZreco'][1]
        histos.append(ROOT.TH2F(bin_name, bin_name, 
                                len(data[one_key]['true_cosZ_edges'])-1,
                                data[one_key]['true_cosZ_edges'][0],
                                data[one_key]['true_cosZ_edges'][-1],
                                len(data[one_key]['true_logE_edges'])-1,
                                data[one_key]['true_logE_edges'][0],
                                data[one_key]['true_logE_edges'][-1]))
        histos[-1].SetTitle(bin_name+';cos(zenith);log10(E/GeV)')


        for ebin in range(histos[-1].GetNbinsY()):
            for zbin in range(histos[-1].GetNbinsX()):
                histos[-1].SetBinContent(zbin+1, ebin+1, 
                                         data[one_key]['eff_areas'][rbin]['effArea'][zbin, ebin])
                
outfile.Write()
outfile.Close()
#h1.Draw('colz')
#raw_input()
