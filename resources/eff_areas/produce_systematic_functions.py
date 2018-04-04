#!/usr/bin/env python

import os, sys, pickle
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
modules_dir = '/afs/ifh.de/user/b/blotsd/scratch/Analysis/OscFit/releases/oscFitND_v1.0/modules'
# modules_dir = os.path.dirname(os.path.realpath(__file__)).split('/resources')[0] + '/modules'
sys.path.append(modules_dir)
import dataLoader


###
# User defined variables
###

# oscFit user name
oscfit_user_name = 'PRD_extended_lowdt_data'

# Directory to store the output (it will be outdir/oscfit_user/)
outdir = '/lustre/fs19/group/icecube/deepcore_analyses/effective_areas'

# Produce figures?
doFigs = True

figsdir = os.path.join(outdir, oscfit_user_name, 'sysfigs')
outdir  = os.path.join(outdir, oscfit_user_name, 'txt')

if not os.path.exists(figsdir):
    os.makedirs(figsdir)

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

doFigs = True


##########
# HOLE ICE
##########

systematic = 'HoleIce'

fmatrix = loader.pmatrix_holeice
y       = loader.y_holeice
yw2     = loader.y_holeice_w2
x       = loader.x_holeice
basevalue = 0.02
x_nominal = np.where(x==basevalue)[0]

txtfile = open(os.path.join(outdir, systematic + '.txt'), 'w') 
    
if doFigs:
    myfig2 = plt.figure(figsize = (24,18))
    figaxes = []
    
line1 = 'All flavors\n'
line2 = 'Hole ice scattering -- Error range = [0.01, 0.03] 1/cm -- Nominal = 0.02 1/cm\n'
line3 = 'Energy bin\t\t Zenith angle bin\t\t Relative change\n'
   
txtfile.write(line1+line2+line3)
counter = 0
for i in range(0, len(fmatrix)):
    for j in range(0, len(fmatrix[0])):
        yvalues =  y[i,j]/y[i,j,x_nominal]
        yerrors = np.sqrt(yw2[i,j] )/y[i,j,x_nominal]

        newx = np.linspace(x.min(), x.max(), 100)

        bin_txt = 'logEreco=['+"%0.2f" % np.log10(loader.bin_edges[0][i])+','+\
            "%0.2f" % np.log10(loader.bin_edges[0][i+1]) +']\t'+\
            'cosZreco=['"%0.2f" % np.cos(loader.bin_edges[1][j]) +"," + \
            "%0.2f" % np.cos(loader.bin_edges[1][j+1]) + ']\t\t'
        if doFigs:
            figaxes.append(myfig2.add_subplot(len(fmatrix),len(fmatrix[0]),counter+1))
            figaxes[-1].errorbar(x, yvalues, yerrors, 
                                 linestyle = '', fmt = 'o', ecolor = 'r', 
                                 mec = 'r', mfc = 'r', capsize = 10)
            figaxes[-1].text(0.012, 1.3, r'$E=['+"%i" % loader.bin_edges[0][i]+'-'+\
                                 "%i" % loader.bin_edges[0][i+1] +']$\n$' 
                         + "%.2f" % np.cos(loader.bin_edges[1][j]) +'< \cos\,\eta <'+\
                                 "%.2f" % np.cos(loader.bin_edges[1][j+1]) + '$', fontsize = 'large')
            
            figaxes[-1].plot(newx, fmatrix[i,j](newx), 'b', linewidth = 3)

        this_line = bin_txt + "%0.3f" % (fmatrix[i][j][1]/fmatrix[i][j](basevalue)) +\
            'x + ' + "%0.3f" % (fmatrix[i][j][0]/fmatrix[i][j](basevalue)) + '\n'
        txtfile.write(this_line)
        counter += 1
        if doFigs:
            if i < (len(fmatrix)-1):
                plt.setp( figaxes[-1].get_xticklabels(), visible=False)
            else:
                figaxes[-1].axes.get_xaxis().set_ticks(np.floor(x*100.)/100.)
            plt.ylim([0., 2.])
            plt.grid(True)    
            myfig2.savefig(os.path.join(figsdir, systematic + '.png'))
txtfile.close()

##########
# HOLE ICE
##########
systematic = 'DOMeff'

fmatrix = loader.pmatrix_domeff
y       = loader.y_domeff
yw2     = loader.y_domeff_w2
x       = loader.x_domeff
basevalue = 1.
x_nominal = np.where(x==basevalue)[0]

txtfile = open(os.path.join(outdir, systematic + '.txt'), 'w') 
    
if doFigs:
    myfig2 = plt.figure(figsize = (24,18))
    figaxes = []
    
    
line1 = 'All flavors\n'
line2 = 'DOM efficiency -- Error range = [0.9, 1.1] -- Nominal = 1.\n'
line3 = 'Energy bin\t\t Zenith angle bin\t\t Relative change\n'

txtfile.write(line1+line2+line3)
counter = 0
for i in range(0, len(fmatrix)):
    for j in range(0, len(fmatrix[0])):
        yvalues =  y[i,j]/y[i,j,x_nominal]
        yerrors = np.sqrt(yw2[i,j] )/y[i,j,x_nominal]
        

        newx = np.linspace(x.min(), x.max(), 100)


        bin_txt = 'logEreco=['+"%0.2f" % np.log10(loader.bin_edges[0][i])+','+\
            "%0.2f" % np.log10(loader.bin_edges[0][i+1]) +']\t'+\
            'cosZreco=['"%0.2f" % np.cos(loader.bin_edges[1][j]) +"," +\
            "%0.2f" % np.cos(loader.bin_edges[1][j+1]) + ']\t\t'
        if doFigs:
            figaxes.append(myfig2.add_subplot(len(fmatrix),len(fmatrix[0]),counter+1))
            figaxes[-1].errorbar(x, yvalues, yerrors, 
                                 linestyle = '', fmt = 'o', ecolor = 'r', 
                                 mec = 'r', mfc = 'r', capsize = 10)
            figaxes[-1].text(0.87, 1.4, r'$E=['+"%i" % loader.bin_edges[0][i]+'-'+\
                                 "%i" % loader.bin_edges[0][i+1] +']$\n$' 
                         + "%.2f" % np.cos(loader.bin_edges[1][j]) +'< \cos\,\eta <'+\
                                 "%.2f" % np.cos(loader.bin_edges[1][j+1]) + '$', fontsize = 'large')
            
            figaxes[-1].plot(newx, fmatrix[i,j](newx), 'b', linewidth = 3)

        this_line = bin_txt + "%0.3f" % (fmatrix[i][j][1]/fmatrix[i][j](basevalue)) +\
            'x + ' + "%0.3f" % (fmatrix[i][j][0]/fmatrix[i][j](basevalue)) + '\n'
            

        txtfile.write(this_line)
        counter += 1
        if doFigs:
            if i < (len(fmatrix)-1):
                plt.setp( figaxes[-1].get_xticklabels(), visible=False)
            else:
                figaxes[-1].axes.get_xaxis().set_ticks(np.floor(x*100.)/100.)
            plt.ylim([0., 2.])
            plt.grid(True)    
            myfig2.savefig(os.path.join(figsdir, systematic + '.png'))
            
txtfile.close()



