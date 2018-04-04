#!/usr/bin/env python
import os, sys
modules_dir = os.path.dirname( os.path.realpath(__file__) ).split('/scripts')[0] + '/modules'
sys.path.append(modules_dir)
print "Modules directory: \n\t ", modules_dir

import numpy as np
from copy import deepcopy
from itertools import cycle
import jp_mpl as jplot
import scipy.stats
import matplotlib.pyplot as plt
import dataLoader

def PlotSystematics(x, y, y_w2, fmatrix , x_nomind,
                    pid_bin ):
    #print "Function to plot systematic sets "
    total_chi2 = 0.0
    fig = plt.figure(figsize=(32,24))
    for i in xrange(len(fmatrix)):
        for j in xrange(len(fmatrix[0])):
            #print i, j, i*len(fmatrix[0]) + j
            curax = fig.add_subplot(len(fmatrix), len(fmatrix[0]),i*len(fmatrix[0]) + j+1 )
            
            y_values =  y[i,j,pid_bin]/y[i,j,pid_bin,x_nomind]
            y_errors = np.sqrt(y_w2[i,j,pid_bin] )/y[i,j,pid_bin,x_nomind]
            curax.errorbar(x, y_values, y_errors, 
                           fmt = "o",linestyle = '', 
                           ecolor = 'r', mec = 'r', mfc = 'r', capsize = 10)
            x_fit = np.linspace(x[0],x[-1], 101)
            fit_func = fmatrix[i,j,pid_bin]
            y_values_fit = fmatrix[i,j,pid_bin](x)
            chi2 = np.sum((y_values_fit - y_values)**2/ y_errors**2)
            total_chi2+=chi2
            curax.plot(x_fit, fit_func(x_fit), "b")
            curax.set_xlim([x[0], x[-1]] ) 
            curax.set_ylim([0.7 , 1.3] ) 
            zed = [tick.label.set_fontsize(12) for tick in curax.yaxis.get_major_ticks()]
            zed = [tick.label.set_fontsize(12) for tick in curax.xaxis.get_major_ticks()]
            x_pos_text = x.min()+0.01*(x.max() - x.min())
            curax.text(x_pos_text, 1.05, 
                       r'$\mathrm{\chi^2 = %0.2f}$ ' % chi2 +"\n"\
                     + r"$\mathrm{E \in\ [\ %0.2f,\ %0.2f\ ]\ GeV}$ \n"%(energy_bin_edges[i],energy_bin_edges[i+1]) +"\n"\
                     + r"$\mathrm{\cos(\theta)\ \in\ [\ %0.2f,\ %0.2f\ ]}$"%(np.cos(zenith_bin_edges[j]),np.cos(zenith_bin_edges[j+1])) ,
                       fontsize = 'large', )


    fig.subplots_adjust(hspace = 0.2, wspace = 0.3, right = 0.99, bottom = 0.1, left = 0.1)
    fig.text(0.08, 0.5, r'$\mathrm{high}\leftarrow\log_{10}(E_\mathrm{reco}\mathrm{GeV})\,\mathrm{bin}\rightarrow\mathrm{low}$',
                       horizontalalignment='center',
                       verticalalignment='center',
                       rotation = 'vertical',
                       fontsize = 25)
    fig.text(0.5, 0.07, r'$\mathrm{horizon}\leftarrow\cos\,(\theta_\mathrm{reco})\,\mathrm{bin}\rightarrow\mathrm{vertical-up}$',
                       horizontalalignment = 'center',
                       verticalalignment='center',
                       rotation = 'horizontal',
                       fontsize = 25)     
    fig.text(0.5, 0.93, r'$\mathrm{User\ : }$ %s$\mathrm{,\ Systematics\ :\ %s\ ,\ PID [%0.3f, %0.3f  ]}$'%(options.USER,cursys,pid_bin_edges[pid_bin],pid_bin_edges[pid_bin+1] ),
                       horizontalalignment = 'center',
                       verticalalignment='center',
                       rotation = 'horizontal',
                       fontsize = 25)    
    fig.text(0.25, 0.93, r'$\mathrm{Total\ chi2\ (in\ figure)\ : \ %0.2f}$'%total_chi2 ,
                       horizontalalignment = 'center',
                       verticalalignment='center',
                       rotation = 'horizontal',
                       fontsize = 20)       
    return fig, total_chi2

from optparse import OptionParser
usage = "usage: %prog [options]"
parser = OptionParser(usage)
parser.add_option("-u", "--user",default='jpall_lowdt',
                  dest="USER", help="Name of os")
parser.add_option("-e", "--energies",default="0.8,1.75,8",
                  dest="ENERGIES", help="Energy range and number of bins (log scale): low,high,nbins")
parser.add_option("-z", "--zeniths",default="-1.0,0.0,8",
                  dest="ZENITHS", help="Cos zenith range and number of bins: low,high,nbins")
parser.add_option("-p", "--pids",default="-inf,0.7,inf",
                  dest="PIDS", help="PID ranges. WARINING! Values of binning has to be specified, etc 0.0,1.0 means 1 bin [0.0,1.0]")
parser.add_option("-w", "--weights",default='weight_e,weight_mu',
                  dest="WEIGHTS", help="Names of weights: nue,numu")
parser.add_option("-s", "--systematics",default='domeff,holeice',
                  dest="SYSTEMATICS", help="Names of systematics to plot: domeff,holeice")
(options,args) = parser.parse_args()
if len(args) != 0:
        crap = "Got undefined options:"
        for a in args:
                crap += a
                crap += " "
        parser.error(crap)

options.ENERGIES = options.ENERGIES.split(",")
options.ZENITHS = options.ZENITHS.split(",")
options.PIDS = options.PIDS.split(",")
options.WEIGHTS = options.WEIGHTS.split(",")
options.SYSTEMATICS = options.SYSTEMATICS.split(",")
print "Will plot systematics: ", options.SYSTEMATICS
weight_keys = [options.WEIGHTS[0], options.WEIGHTS[1]]
energy_range = (float(options.ENERGIES[0]),
                float(options.ENERGIES[1]),
                  int(options.ENERGIES[2]) )

zenith_range = (float(options.ZENITHS[0]),
                float(options.ZENITHS[1]),
                  int(options.ZENITHS[2]) )
pid_bin_edges = np.zeros(len(options.PIDS), dtype=float)
for i in xrange(0, len(options.PIDS)):
    pid_bin_edges[i] = float(options.PIDS[i])
    
print "Reading user: ", options.USER
print "Log E range: \t[%0.2f, %0.2f] \t in %i bins"%energy_range
print "Cos range: \t[%0.2f, %0.2f] \t in %i bins"%zenith_range
print "PID bins: \t",  pid_bin_edges, "\t total", len(pid_bin_edges)-1,  " bins"
print "Weight keys: \t NuE : ", weight_keys[0], "\t ", "NuMu : ",weight_keys[1]
energy_bin_edges = 10.0**np.linspace(energy_range[0],energy_range[1],energy_range[2]+1)
zenith_bin_edges = np.arccos( np.linspace(zenith_range[0],zenith_range[1],zenith_range[2]+1))[::-1]

print "E: ", energy_bin_edges
print "Z: ", zenith_bin_edges
print "P: ", pid_bin_edges
#energy_bin_edges = np.linspace(energy_range[0],energy_range[1],energy_range[3]+1)


loader =  dataLoader.dataLoader(
            observables =  ['reco_energy', 
                            'reco_zenith', 
                            'pid'],
            bin_edges   =  [energy_bin_edges,
                            zenith_bin_edges,
                            pid_bin_edges],
             user = options.USER, 
             LEaxis = [],
             weight_keys = weight_keys, 
                          extra_cuts = {   
                              'pid':        [   pid_bin_edges[0] -1e-5,
                                                pid_bin_edges[-1]+1e-5], 
                              'reco_energy':[energy_bin_edges[0] -1e-5,
                                             energy_bin_edges[-1]+1e-5], 
                              'reco_zenith':[zenith_bin_edges[0] -1e-5,
                                             zenith_bin_edges[-1]+1e-5],
                              },   
             break_energy = 1000.,
             detailed_detsys =False,
             verbose = False, 
 )

settings = deepcopy(loader.default_data_settings)
histos = loader.getSumNuHistogram(settings)

for cursys in options.SYSTEMATICS:
    if not cursys in loader.user.mc_sets.keys():
        print "Error! No systematic set for ", cursys, " in data loader! Skipping"
        continue
    print "Plotting systematic parameters: ", cursys   
    x_vals = getattr(loader, "x_"+ cursys)
    if cursys == 'domeff': x_vals = np.log10(x_vals)
    x_nomind = getattr(loader, 'x_'+cursys+'_nominal')
    y_vals = getattr(loader, "y_"+ cursys)
    y_w2_vals = getattr(loader, "y_"+ cursys + "_w2")
    chi2_tot = 0.0
    for pid_i in xrange(0, len(pid_bin_edges)-1):
        #print pid_i
        figure, chi2 = PlotSystematics(x = x_vals, 
                              y = y_vals,
                              y_w2 = y_w2_vals,
                              fmatrix = getattr(loader, "pmatrix_"+ cursys),
                              x_nomind = x_nomind,
                              pid_bin = pid_i )
        figure.savefig(options.USER + "_" + cursys +"_PID_"+str(pid_i+1) + ".png", bbox_inches= "tight")
        chi2_tot +=chi2
    print "Total chi2 in all the points : ", chi2_tot
