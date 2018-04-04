import os, sys
modules_dir = os.path.dirname(os.path.realpath(__file__)).split('/plot-scripts')[0] + '/modules'
sys.path.append(modules_dir)
import matplotlib
matplotlib.use('Agg')

from numpy import *
from pylab import *
from jp_mpl import *
from numpy.polynomial.polynomial import polyfit
import numpy as np

import dataLoader as dt
loader =  dt.dataLoader(observables =  
                        ['reco_energy', 'reco_zenith', 'pid'],
                        bin_edges   =  
                        [10**np.linspace(0.8,2.,11),
                         np.arccos(np.linspace(-1,0.,9))[::-1],
                         np.array([0, 0.8, np.inf])],
                        #np.linspace(0., 0.7, 2)],
                        #np.concatenate((np.linspace(0,1.4,5)[:-1], [np.inf]))],
                        user = 'jpall_light', 
                        LEaxis = None,      
                        pol_degree = 1,     
                        weight_keys = ['weight_e', 'weight_mu'], 
                        extra_cuts = {},    
                        #break_energy = 800.,
                        Wcut = None,#10**0.2,
                        verbose = False)    


# loader =  dt.dataLoader(bin_edges   = [10**np.linspace(0.8,1.75,6),
#                                        np.arccos(np.linspace(-1,0.,5))[::-1],
#                                        np.array([ 0. ,  0.8, np.inf])],
#                         observables = ['reco_energy', 'reco_zenith', 'cascade_like'],
#                         user = 'jp_ms', 
#                         LEaxis = None,
#                         pol_degree = 1, 
#                         weight_keys = ['weight_e', 'weight_mu'], 
#                         verbose = False)

# loader =  dt.dataLoader(bin_edges   = [10**np.linspace(0.8,1.75,9),
#                                        np.arccos(np.linspace(-1,0.,9))[::-1],
#                                        np.array([ 0. ,  0.35,  0.7 ,  1.05, np.inf])],
#                         observables = ['reco_energy', 'reco_zenith', 'cascade_like'],
#                         user = 'jp_ms', 
#                         LEaxis = None,
#                         pol_degree = 1, 
#                         weight_keys = ['weight_ej12', 'weight_muj12'], 
#                         verbose = False)


# The local way is to get a fit for each bin in the histogram.
multiplier = 79988461.37448613 # This is only to get a feeling of Number of events
pol_degree  = 1
outdir = '/afs/ifh.de/user/y/yanezjua/scratch/analysis-results/systematic-checks/parameterizations_2015/'
histo_list = ['nue_histo', 'numu_histo', 'nutau_histo', 'nc_histo']


fmatrix_list = [loader.pmatrix_holeice, loader.pmatrix_domeff, loader.pmatrix_hadlight]
y_list       = [loader.y_holeice, loader.y_domeff, loader.y_hadlight]
yw2_list     = [loader.y_holeice_w2, loader.y_domeff_w2, loader.y_hadlight_w2]
x_list       = [loader.x_holeice, loader.x_domeff, loader.x_hadlight]
basevalue_list = [0.02, 1., 1.,]
xlabels      = [r'Scattering coefficient in hole ice (1/cm)',
                r'Optical efficiency variations',
                r'Hadronic light yield']
figname      = ['HoleIce', 'DomEff', 'HadLight']
norms        = [loader.holeice_norms, loader.domeff_norms, loader.hadlight_norms]

for index in range(len(fmatrix_list)):
    for myc in histo_list:
        fmatrix = fmatrix_list[index][myc]
        y       = y_list[index][myc]
        yw2     = yw2_list[index][myc]
        x       = x_list[index]
        basevalue = basevalue_list[index]


        for oneq in range(fmatrix.shape[2]):
            myfig2 = figure(figsize = (16,12))
            figaxes = []
            counter = 0

            for i in range(0, len(fmatrix)):
                for j in range(0, len(fmatrix[0])):
                    one_fcn = fmatrix[i,j,oneq,:]
                    figaxes.append(myfig2.add_subplot(len(fmatrix),len(fmatrix[0]),counter+1))
                    yvalues = y[i,j,oneq,:]*norms[index]
                    #yerror = yvalues*sqrt(yw2[i, j, oneq,:])/(one_fcn(basevalue)/norms[index])
                    yerror = yvalues*sqrt(yw2[i, j, oneq,:])/y[i,j,oneq,:]
                    figaxes[-1].errorbar(x, yvalues, yerror, linestyle = '', fmt = 'o', 
                                         ecolor = 'r', mec = 'r', mfc = 'r', capsize = 10)
                    # figaxes[-1].text(0.005, 1.25, r'$E=['+"%i" % loader.energy_axis[i]+'-'+"%i" % loader.energy_axis[i+1] +']$\n$' 
                    #                  + "%.2f" % cos(loader.zenith_axis[j]) +'< \cos\,\eta <'+"%.2f" % cos(loader.zenith_axis[j+1]) + '$', fontsize = 'large')
                    # figaxes[-1].text(0.005, 0.6, 'Counts: '+ "%i" % (fmatrix[i][j](basevalue)*multiplier))

                    newx = linspace(x.min(), x.max(), 100)
                    figaxes[-1].plot(newx, one_fcn(newx)/one_fcn(basevalue), 'b', linewidth = 3)
                    counter += 1
                    if i < (len(fmatrix)-1):
                        setp( figaxes[-1].get_xticklabels(), visible=False)
                    else:
                        figaxes[-1].axes.get_xaxis().set_ticks(floor(x*100.)/100.)
                    ylim([0., 2.])
                    grid(True)

            myfig2.text(0.5, 0.03, xlabels[index], #r'Scattering coefficient in hole ice (1/cm)',
                        horizontalalignment = 'center',
                        verticalalignment='center',
                        rotation = 'horizontal',
                        fontsize = 'large')

            myfig2.text(0.5, 0.95,myc + ' Qbin ' + "%i" % oneq,
                        horizontalalignment = 'center',
                        verticalalignment='center',
                        rotation = 'horizontal',
                        fontsize = 'large')
            myfig2.text(0.03, 0.5, 'Change in number of events from baseline',
                           horizontalalignment='center',
                           verticalalignment='center',
                           rotation = 'vertical',
                           fontsize = 'large')

            myfig2.subplots_adjust(hspace = 0.2, wspace = 0.3, right = 0.99, bottom = 0.1, left = 0.1)
            myfig2.savefig(outdir + figname[index] + '_' + myc  + '_Q' +"%i" % oneq + '_'+"%i" % pol_degree +'deg.png')
            show()
    print '\n\nFinished with ', figname[index]

# print 'DOM efficiency - Doing the figure'

# for myc in histo_list:
#     fmatrix = loader.pmatrix_domeff[myc]
#     y       = loader.y_domeff[myc]
#     yw2     = loader.y_domeff_w2[myc]
#     x       = loader.x_domeff
#     basevalue = 1.

#     for oneq in range(fmatrix.shape[2]):
#         myfig2 = figure(figsize = (16,12))
#         figaxes = []
#         counter = 0

#         for i in range(0, len(fmatrix)):
#             for j in range(0, len(fmatrix[0])):
#                 one_fcn = np.poly1d(fmatrix[i,j,oneq,:])
#                 figaxes.append(myfig2.add_subplot(len(fmatrix),len(fmatrix[0]),counter+1))
#                 yvalues = y[i,j,oneq,:]/one_fcn(basevalue)
#                 yerror = yvalues*sqrt(yw2[i, j, oneq,:])/one_fcn(basevalue)
#                 figaxes[-1].errorbar(x, yvalues, yerror, linestyle = '', fmt = 'o', ecolor = 'r', mec = 'r', mfc = 'r', capsize = 10)
#                 #figaxes[-1].text(0.92, 1.25, r'$E=['+"%i" % loader.energy_axis[i]+'-'+"%i" % loader.energy_axis[i+1] +']$\n$' 
#                 #                 + "%.2f" % cos(loader.zenith_axis[j]) +'< \cos\,\eta <'+"%.2f" % cos(loader.zenith_axis[j+1]) + '$', fontsize = 'large')
#                 #figaxes[-1].text(0.92, 0.6, 'Counts: '+ "%i" % (fmatrix[i][j](basevalue)*multiplier))

#                 newx = linspace(x.min(), x.max(), 100)
#                 figaxes[-1].plot(newx, one_fcn(newx)/one_fcn(basevalue), 'b', linewidth = 3)
#                 counter += 1
#                 if i < (len(fmatrix)-1):
#                     setp( figaxes[-1].get_xticklabels(), visible=False)
#                 else:
#                     figaxes[-1].axes.get_xaxis().set_ticks(x)
#                 ylim([0., 2.])
#                 xlim([x.min(), x.max()])
#                 grid(True)

#         title('Quality bin ' + "%i" % oneq)
#         myfig2.text(0.5, 0.03,r'Optical efficiency variations',
#                     horizontalalignment = 'center',
#                     verticalalignment='center',
#                     rotation = 'horizontal',
#                     fontsize = 'large')

#         myfig2.text(0.5, 0.95,myc,
#                     horizontalalignment = 'center',
#                     verticalalignment='center',
#                     rotation = 'horizontal',
#                     fontsize = 'large')

#         myfig2.text(0.03, 0.5, 'Change in number of events from baseline',
#                        horizontalalignment='center',
#                        verticalalignment='center',
#                        rotation = 'vertical',
#                        fontsize = 'large')

#         myfig2.subplots_adjust(hspace = 0.2, wspace = 0.3, right = 0.99, bottom = 0.1, left = 0.1)
#         myfig2.savefig(outdir + 'DomEff_' + myc  + '_Q' +"%i" % oneq + '_'+"%i" % pol_degree +'deg.png')

#         show()
