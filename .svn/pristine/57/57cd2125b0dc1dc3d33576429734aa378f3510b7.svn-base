import os, sys, pickle
import matplotlib as mpl
#mpl.use('Agg')
import matplotlib.pyplot as plt

import jp_mpl as jplot
import numpy as np
from matplotlib import gridspec
from copy import deepcopy
from scipy import interpolate
import oscfit_default_values

sys.path.append(os.path.dirname(os.path.realpath(__file__)))
import dataLoader

label_fontsize = 18
ticks_fontsize = 14
ticks_fontsize_small = 13
mpl.rc('xtick', labelsize=ticks_fontsize)#, labelweight = 'bold')
mpl.rc('ytick', labelsize=ticks_fontsize)
mpl.rc('ytick.minor', size=6)
mpl.rc('xtick.major', size=6)
mpl.rc("axes", linewidth=2.0, grid=True)#, labelweight = 'bold')
mpl.rc("lines", markeredgewidth=2.0, linewidth = 3)


minos_color  = 'g'
t2k_color    = 'm'
sk_color     = 'r'
ic_color     = ['#00688B', '#42C0FB', '#0276FD']
ic_color_old = 'b'


t2k_style = '--'
minos_style = '-.'
sk_style = ':'


# We have to assume a dm21 value. This is the one we have been using for long:
dm21 = 7.54E-5
oscillations_results_path = os.path.dirname(os.path.realpath(__file__)).split('/modules')[0] + \
    '/resources/oscillation_results'
sys.path.append(oscillations_results_path)


def plotNuisance(scan_results, parameters):
    return True

def plotScanResults(scan_results_list  = [], 
                    iccolor            = ['#283A90','g','b','k','m'],
                    iclabels           = [], 
                    compare_results    = False, 
                    outdir             = None, 
                    same_scale         = False,
                    side_boxes         = True):
    '''
    Give a list of results and it will produce contours
    This function still has to be fully updated
    '''

    contour_lists = []

    # Get the absolute min LLH
    llh_list = np.zeros(len(scan_results_list))
    for i, scan_results in enumerate(scan_results_list):
        llh_list[i] = scan_results['llh'].min()
    llh_min = np.min(llh_list)

    # Confidence contours in sin(theta)**2
    sin_theta_lim = [0.25, 0.75]
    dm32_lim      = [2., 4.]
    contour_list  = []

    if compare_results:
        delta_llh = np.array([4.6])*0.5
        linestyles = ['-']
        suffix = '_compare'
    else:
        delta_llh = np.array([2.3, 4.6])*0.5
        linestyles = ['--','-']
        suffix = '_single'

    ########################################
    # Confidence contours in sin(theta)**2 #
    ########################################
    if scan_results_list[0]['fit_settings']['oscMode'] != 'TwoNeutrino':

        for i in range(len(scan_results_list)):
            scan_results = scan_results_list[i]

            if scan_results['fit_settings']['oscMode'] == 'TwoNeutrino':
                print 'You cannot compare 2-neutrino fits like this!'
                exit()

            if scan_results['bestfit']['dm31'] > 0:
                #print 'plotTools: Normal hierarchy'
                normal_hierarchy = True
                dm32_bestfit = abs(scan_results['bestfit']['dm31']) - dm21
                dm32_values  = scan_results['dm31_list'] - dm21

            else:
                #print 'plotTools: Inverted hierarchy'
                normal_hierarchy = False
                dm32_bestfit = abs(scan_results['bestfit']['dm31']) + dm21
                dm32_values  = scan_results['dm31_list'] + dm21


            print '\n***Parameters at best fit (from the scan)***'
            for one_key in scan_results.keys():
                if '_bestFit' in one_key:
                    print one_key,'\t', scan_results[one_key]

            xstep = scan_results['theta23_list'][1] - scan_results['theta23_list'][0]
            xaxis = scan_results['theta23_list'] - xstep/2.
            ystep = scan_results['dm31_list'][1] - scan_results['dm31_list'][0]
            yaxis = scan_results['dm31_list'] - ystep/2.

            if normal_hierarchy:
                yaxis -= dm21
            else:
                yaxis += dm21

            llh_scan  = scan_results['llh'] - scan_results['llh'].min()
            if not same_scale:
                llh_scan_diff = 0
                suffix += '_fixScale'
            else:
                llh_scan_diff = scan_results['llh'].min() - llh_min


            if i == 0:
                fig_cont1 = plt.figure('Likelihood', figsize = (8,6))              
                gs = gridspec.GridSpec(2, 2, width_ratios=[5, 1], height_ratios=[1.33, 4.67]) 

                up_box = plt.subplot(gs[0])            
                up_box.axes.get_yaxis().set_ticks([0,1,2,3,4])
                #up_box.axes.get_yaxis().set_ticks([0,5, 10, 15, 20])
                up_box.yaxis.set_tick_params(labelsize =ticks_fontsize_small)

                plt.ylim([0,4])
                plt.xlim(sin_theta_lim)
                plt.ylabel(r'$-2\Delta\ln L$', fontsize = ticks_fontsize_small)

                side_box = plt.subplot(gs[3])
                side_box.axes.get_xaxis().set_ticks([0,1,2,3,4])
                #side_box.axes.get_yaxis().set_ticks([0,5, 10, 15, 20])
                side_box.xaxis.set_tick_params(labelsize=ticks_fontsize_small)

                plt.xlim([0,4])
                plt.ylim(dm32_lim)
                plt.xlabel(r'$-2\Delta\ln L$', fontsize = ticks_fontsize_small)

                plt.setp(up_box.get_xticklabels(), visible = False)
                plt.setp(side_box.get_yticklabels(), visible = False)


            #up_box = subplot(gs[0])
            print 'plotTools: One sigma errors for physics parameters'

            mixangle_projection = llh_scan.min(axis=1) + llh_scan_diff
            mixangle_values     = np.sin(scan_results['theta23_list'])**2
            yToFind = 1.
            yreduced = np.array(2*mixangle_projection) - yToFind
            freduced = interpolate.UnivariateSpline(mixangle_values, 
                                                    yreduced, s=0)
            mixangle_bestfit = np.sin(scan_results['bestfit']['theta23'])**2
            mixangle_errors  = mixangle_bestfit - freduced.roots()
            print 'Theta23', mixangle_bestfit, mixangle_errors, np.mean(np.abs(mixangle_errors))
        

            dm32_projection = llh_scan.min(axis=0) + llh_scan_diff
            yreduced = np.array(2*dm32_projection) - yToFind
            freduced = interpolate.UnivariateSpline(dm32_values, 
                                                    yreduced, s=0)
            dm32_errors  = dm32_bestfit - freduced.roots()
            print 'DM32', dm32_bestfit*10**3, dm32_errors*10**3, np.mean(np.abs(dm32_errors))*10**3

            if side_boxes:
                up_box.plot(np.sin(xaxis)**2, 2*mixangle_projection, color = iccolor[i])
                side_box.plot(2*dm32_projection,yaxis*10**3, color = iccolor[i])

            if i == 0:
                if side_boxes:
                    cont1_ax  = plt.subplot(gs[2], sharex = up_box, sharey=side_box)
                cont1_ax.axes.get_yaxis().set_ticks(np.arange(1., 4., 0.2))
                cont1_ax.axes.get_xaxis().set_ticks([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9])

            #cont1_ax  = subplot(gs[2])
            #pcolor(sin(xaxis)**2, yaxis*10**3, llh_scan.T)
            contour_list.append(plt.contour(np.sin(xaxis)**2, yaxis*10**3, llh_scan.T, delta_llh, 
                                        linewidths=(3,3), colors=(iccolor[i],iccolor[i]), 
                                        linestyles=linestyles, linewidth = 3))
            cont1_ax.plot(np.sin(scan_results['bestfit']['theta23'])**2, dm32_bestfit*10**3, 
                          markersize = 10, markeredgewidth = 3, markeredgecolor = iccolor[i], 
                          marker = 'x', zorder = 0) 

            if i == (len(scan_results_list)-1):
                if compare_results:
                    import sintheta2_results as st2
                    l1, =plt.plot(st2.t2k[:,0], st2.t2k[:,1], t2k_style,label = 'T2K 2014', 
                              color = t2k_color, zorder = -9)
                    # plot(st2.t2kbf[0], st2.t2kbf[1], marker = 'x', markeredgecolor = t2k_color, 
                    # markersize = 8, markeredgewidth = 4, zorder = -9)

                    l2, =plt.plot(st2.minos[:,0], st2.minos[:,1],minos_style, label = 'MINOS', 
                              color = minos_color, zorder = -9)
                    # plot(st2.minosbf[0], st2.minosbf[1], marker = 'x', markeredgecolor = minos_color, 
                    # markersize = 8, markeredgewidth = 4, zorder = -9)

                    l3, =plt.plot(st2.sk2015[:,0], st2.sk2015[:,1], sk_style,label = 'SK IV 2015', 
                              color = sk_color, zorder = -9)
                    # plot(st2.skbf[0], st2.skbf[1], marker = 'x', markeredgecolor = sk_color, 
                    #      markersize = 8, markeredgewidth = 4, zorder = -9)

                    # l4, =plot(st2.t2k_future[:,0], st2.t2k_future[:,1], '-',label = 'T2K 2020', 
                    #           color = t2k_color, zorder = -9)

                    l5, =plt.plot(st2.ic_x, st2.ic_y, '-', label = 'IC 2014 tracks [NH]', 
                              color = ic_color_old, zorder = -9)

                    plt.legend((contour_list[0].legend_elements()[0][0], l2, l1, l3, l5), 
                           (iclabels[i], 'MINOS w/atm [NH]', 'T2K 2014 [NH]', 'SK IV 2015 [NH]', 'IC2014 [NH]'), 
                           ncol = 2, loc = 9, numpoints = 2)

                    # legend((contour_list[0].legend_elements()[0][0], l2, l1, l3, l4), 
                    #        (iclabels[i], 'MINOS w/atm [NH]', 'T2K 2014 [NH]', 'SK IV [NH]', 'T2K 2020 [NH]'), 
                    #        ncol = 2, loc = 9, numpoints = 2)

                    # legend((contour_list[0].legend_elements()[0][0], l2, l1, l3), 
                    #        (iclabels[i], 'MINOS w/atm [NH]', 'T2K 2014 [NH]', 'SK IV [NH]'), 
                    #        ncol = 2, loc = 9, numpoints = 2)

                    # legend((c1.legend_elements()[0][0],c1.legend_elements()[0][1] ), 
                    #        ('IC 68% C.L.', 'IC 90% C.L.'), loc=0)

                    # legend((l1,l2,l3, l4, l5, c2.legend_elements()[0][0]), 
                    #        ('T2K 2014', 'MINOS', 'SK IV', 'IC79 (2012)', 'IC86 (2013)', 'IC86 (2014-3years)'),
                    #        loc=0, ncol = 2)

                    plt.text(0.44, 3.38,'90% CL contours',  color = '0.1', fontsize=ticks_fontsize)
                    #text(0.5, 3.7,'IceCube Preliminary',  color = '0.1', fontsize=16, style = 'italic')

            if not compare_results:
                contour_lists.append(contour_list[i].legend_elements()[0][-1])
                print contour_lists
                if i == (len(scan_results_list)-1):
                    plt.legend(contour_lists, (iclabels), ncol = 1, loc = 9)
                    plt.text(0.33, 3.4,'68% (dashed) and 90% (solid) CL contours',  
                         color = '0.1', fontsize=ticks_fontsize)
                    plt.text(0.26, 2.07,'IceCube Preliminary',  color = '0.1', 
                         fontsize=label_fontsize, style = 'italic')

            plt.ylim(dm32_lim)
            plt.xlim(sin_theta_lim)
            plt.xlabel(r'$\sin^{2}(\theta_{23})$', fontsize=label_fontsize)

            plt.ylabel(r'$\|\Delta m^{2}_{32}\|\,(10^{-3}$eV${}^2)$', fontsize=label_fontsize)

            fig_cont1.subplots_adjust(bottom = 0.12, left = 0.11, top = 0.97, 
                                      right = 0.98, wspace = 0, hspace = 0.)
            if outdir:
                fig_cont1.savefig(os.path.join(outdir,'Contours_v1'+suffix+'.pdf'))
                print 'plotTools: Saving!'
                fig_cont1.savefig(os.path.join(outdir,'Contours_v1'+suffix+'.png'), dpi = 160)
                if not compare_results:
                    #x: data_XXCL[:,0], y: data_XXCL[:,1]
                    data_68CL = contour_list[0].allsegs[0][0] # One sigma
                    data_90CL = contour_list[0].allsegs[1][0] # 90 %
                    pickle.dump({'68':data_68CL, '90':data_90CL}, 
                                open(os.path.join(outdir,'Contours_data.pckl'),'w'))
                    


    ##########################################
    # Confidence contours in sin(2*theta)**2 #
    ##########################################
    else:
        mixangle_lim = [0.6, 1.]
        suffix += '_s2t'
        fig_cont1 = plt.figure('Likelihood', figsize = (8,6)) 
        cont1_ax = fig_cont1.add_subplot(111)
        for i in range(len(scan_results_list)):
            scan_results = scan_results_list[i]

            dm32_bestfit = abs(scan_results['bestfit']['dm31'])

            print '***Parameters at best fit (from the scan)***'
            for one_key in scan_results.keys():
                if '_bestFit' in one_key:
                    print one_key, scan_results[one_key]

            xstep = scan_results['theta23_list'][1] - scan_results['theta23_list'][0]
            xaxis = scan_results['theta23_list']
            ystep = scan_results['dm31_list'][1] - scan_results['dm31_list'][0]
            yaxis = scan_results['dm31_list']

            llh_scan  = scan_results['llh'] - scan_results['llh'].min() #llh_min
            if not same_scale:
                #llh_min = scan_results['llh'].min()
                llh_scan_diff = 0
                suffix += '_fixScale'
            else:
                llh_scan_diff = scan_results['llh'].min() - llh_min

            print llh_scan.size
            contour_list.append(plt.contour(xaxis, yaxis*10**3, llh_scan.T, delta_llh, 
                                        linewidths=(3,3), 
                                        colors=(iccolor[i],iccolor[i]), 
                                        linestyles=linestyles, linewidth = 3))
            cont1_ax.plot(scan_results['bestfit']['mix_angle'], dm32_bestfit*10**3, 
                          markersize = 10, markeredgewidth = 3, 
                          markeredgecolor = iccolor[i], marker = 'x', zorder = 0) 

        plt.ylim(dm32_lim)
        plt.xlim(mixangle_lim)
        plt.xlabel(r'$\sin^{2}(\theta_{23})$', fontsize=label_fontsize)

        plt.ylabel(r'$\|\Delta m^{2}_{32}\|\,(10^{-3}$eV${}^2)$', fontsize=label_fontsize)
        fig_cont1.subplots_adjust(bottom = 0.12, left = 0.11, top = 0.97, right = 0.98, wspace = 0, hspace = 0.)
        plt.show()


    ########################
    # Nuisance parameters  #
    ########################

    physics_paramters  = ['dm31','theta23']
    atm_params         = {'norm_nu':[np.inf, 0.8, 1.2],
                          'norm_e': [1., -0.1, 0.1],
                          'gamma':  [ 0.05, -0.2, +0.2],
                          'atmmu_f':[np.inf, 0.8, 1.2],
                          'numubar_e2':[0., -2, 2],
                          'numubar_zw':[0., -2, 2]}
    xsec_params        = {'norm_nc':[1., -0.4, +0.4],
                          'axm_qe': [ 0., -2., +2],
                          'axm_res':[0., -2, +2]}
    detector_params    = {'pid_bias':[0, -0.4, 0.4],
                          'hole_ice':[0.02, -0.02, 0.02],
                          'domeff':[1., -0.2, 0.2],
                          'had_escale':[1., -0.2, 0.2]}
                          



    # # Doing the rest of the contours in sin(theta)**2
    # # FLUX figures

    # syst_name  = ['norm_nu', 'norm_e', 'gamma', 'norm_atmmu']
    # syst_title = ['Neutrino flux normalization', 'Electron neutrino deviation', 'Spectral index change', 'Atm. muon normalization']
    # syst_std   = [sqrt(scan_results['bestfit']['hor_events']), 0.2, 0.05, 1.]
    # syst_mean  = [0., 1., 0., 0.]
    # vmin_max   = [[-1., 1.],[-1., 1.],[-1., 1.], [0, 0.1]]#[0., 0.2]]

    # fsyst1    = figure( figsize=(2*fig_xsize, 2*fig_ysize))
    # fsyst1_subplots = []
    # fsyst1_axes = []

    # for i in range(len(syst_name)):
    #     fsyst1_subplots.append(fsyst1.add_subplot(2, 2, i+1))
    #     fsyst1_subplots[-1].set_title(syst_title[i])
    #     fsyst1_axes.append(fsyst1_subplots[-1].pcolor(sin(scan_results['theta23_list'])**2, yaxis*10**3, 
    #                                                   (scan_results[syst_name[i]].T - syst_mean[i])/syst_std[i], zorder = -9,
    #                                                   vmin = vmin_max[i][0], vmax = vmin_max[i][1]))
    #     fsyst1_subplots[-1].axis('tight')
    #     colorbar(fsyst1_axes[-1])
    #     fsyst1_subplots[-1].plot(sin(scan_results['bestfit']['theta23'])**2, dm32_bestfit*10**3, markersize = 12, markeredgewidth = 4, markeredgecolor = 'black', marker = 'x') 
    # if outdir:
    #     fsyst1.savefig(outdir + 'Systematics_flux.png')

    # syst_name  = ['domeff', 'hole_ice']
    # syst_title = ['Light collection eff.', 'Scattering coeff. hole']
    # syst_std   = [0.1, 0.01]
    # syst_mean  = [1., 0.02]

    # fsyst2    = figure( figsize=(2*fig_xsize, 2*fig_ysize))
    # fsyst2_subplots = []
    # fsyst2_axes = []
    # fsyst2 = figure( figsize=(2*fig_xsize, fig_ysize))

    # for i in range(len(syst_name)):
    #     fsyst2_subplots.append(fsyst2.add_subplot(1, 2, i+1))
    #     fsyst2_subplots[-1].set_title(syst_title[i])
    #     fsyst2_axes.append(fsyst2_subplots[-1].pcolor(sin(scan_results['theta23_list'])**2, yaxis*10**3, 
    #                                                   (scan_results[syst_name[i]].T-syst_mean[i])/syst_std[i], zorder = -9,
    #                                                   vmin=-1, vmax = 1.))
    #     fsyst2_subplots[-1].axis('tight')
    #     colorbar(fsyst2_axes[-1])
    #     plot(sin(scan_results['bestfit']['theta23'])**2, dm32_bestfit*10**3, 
    #          markersize = 12, markeredgewidth = 4, markeredgecolor = 'black', marker = 'x') 
    # if outdir:
    #     fsyst2.savefig(outdir + 'Systematics_detector.png')

    plt.show()


def plotCorrelationMatrix(scan_results = None, outdir = None):
    return True


def plotFitDeviations(scan_results=None, outdir = None):

    if not scan_results['result'].has_key('errors'):
    	print "Error estimates by the fit are not saved by default. If you want to calcuate the errors"
    	print " and see the deviations from the fit, rerun the fit or scan and store the fit details by" 
    	print " for a single_fit, rerun after setting: store_fit_details = True"
    	print " for scan using -r (refit) optioneg. ./plot_scan_results.py my_scan_name -r"
	print "please press any key to continue. All other plots will still be made." 
	raw_input()
    	return True



    nominal_values = oscfit_default_values.default_data_settings
    if scan_results['data_settings'] != []:
        nominal_values.update(scan_results['data_settings'])    

    # I don't want to see these labels. Remove them.
    skip_keys = ['norm','baseline_llh', 'theta23','dm31','dm32','theta24','theta13','theta34', 'dm41']
    result_error_dict = {}
    for one_key in scan_results['result']['errors'].keys():
        if one_key in skip_keys:
            continue
        result_error_dict[one_key] = \
            (scan_results['result'][one_key]- nominal_values[one_key])/scan_results['result']['errors'][one_key] 
        
    myfig = plt.figure(figsize=(10,8))
    ax = myfig.add_subplot(111)

    pos = np.arange(len(result_error_dict)) + 0.5

    plt.barh(pos,result_error_dict.values(), align='center')
    plt.yticks(pos, result_error_dict.keys())
    plt.xlabel('Standard deviations from nominal')
    plt.title('(Fit - Nominal) / Fit Error')
    plt.xlim([-2,2])
    plt.vlines([-1, -0.5, 0.5,1],pos.min(), pos.max(), 
               colors = ['r','b','b','r'], linestyle = ['-','--','--','-'])
    plt.grid(True)

    myfig.subplots_adjust(bottom = 0.10, left = 0.15, top = 0.95, 
                          right = 0.95, wspace = 0, hspace = 0.)
    if outdir:
        myfig.savefig(os.path.join(outdir, 'FitDeviations.png'), dpi=160)
        


def plotBestFit3D(scan_results = None, outdir = None):
    '''
    Start by plotting the fit deviations
    '''
    plotFitDeviations(scan_results, outdir)


    '''
    Produe Data / MC comparisons for a best fit point. 
    '''

    rylim = [0.6, 1.4]

    nofit_settings = deepcopy(scan_results['result'])
    nofit_settings['dm31'] = nofit_settings['theta23'] = 0.


    # Print the result
    skip_keys = ['fit_settings','corr_matrix','covariance', 'parameters']
    print '\n\n **** RESULTS ****'
    result_keys = scan_results['result'].keys()
    result_keys.sort()
    for one_key in result_keys:
        if one_key in skip_keys:
            continue
        print one_key, '\t',  scan_results['result'][one_key]
    print '\n\n'

    for sindex, loader_settings in enumerate(scan_results['loader_settings']):
        
        # Checking the number of PID bins
        if len(loader_settings['observables']) < 3:
            print 'plotTools: It seems like there is no PID in use.'
            pid_nbins = 1
            pid_used  = False
        else:
            pid_nbins = len(loader_settings['bin_edges'][-1])-1
            pid_used  = True


        #####
        # Doing the 2D histograms
        #####

        loader_settings['LEaxis'] = []

        data_histo     = scan_results['data_histograms'][sindex]
        one_loader     = dataLoader.dataLoader(**loader_settings)
        bestfit_histo  = one_loader.loadMCasData(scan_results['result'])
        nofit_histo    = one_loader.loadMCasData(nofit_settings)

        for pid in range(pid_nbins):

            #################
            # Doing it in energy bands, paper-like
            #################
            myfig = plt.figure(figsize=(9,8))
            plot_axis = np.cos(one_loader.bin_edges[1])
            myxaxis   = (plot_axis[:-1]+plot_axis[1:])/2.
            gs = gridspec.GridSpec((one_loader.bin_edges[0].size-1)/2,2)
            nestedgs = []
            figaxes = []

            # This is done 
            for i in range(0, one_loader.bin_edges[0].size-1):

                # Low i = low energy
                if pid_used:
                    data_1d    = data_histo[i,:,pid]
                    bestfit_1d = bestfit_histo[i,:,pid]
                    nofit_1d   = nofit_histo[i,:,pid]
                else:
                    data_1d    = data_histo[i,:]
                    bestfit_1d = bestfit_histo[i,:]
                    nofit_1d   = nofit_histo[i,:]                    

                dmc_1d     = data_1d/bestfit_1d
                derrors = dmc_1d*np.sqrt(data_1d)/data_1d

                these_axes = []
                nestedgs.append(gridspec.GridSpecFromSubplotSpec(2,1, subplot_spec=gs[i/2, i%2]))

                # Here plot all the contributions
                figaxes.append(myfig.add_subplot(nestedgs[-1][0,0],rasterization_zorder = -9))
                plt.errorbar(myxaxis, data_1d, yerr = [np.sqrt(data_1d), np.sqrt(data_1d)], 
                             fmt = 'o', ecolor = 'k', mfc = 'k', label = 'Data', zorder= 99, markersize = 6)
                jplot.unfilledBar(plot_axis, bestfit_1d, 
                                  color = 'b', label = 'MC with osc.', zorder = -9)
                jplot.unfilledBar(plot_axis, nofit_1d, 
                                  color = 'r', linestyle = '--', label = 'MC no osc.', zorder = -9)
                figaxes[-1].invert_xaxis()
                plt.ylim([0, 250])
                figaxes[-1].axes.get_xaxis().set_ticklabels([])
                figaxes[-1].axes.get_yaxis().set_ticks([70,140, 210])
                if (i % 2) != 0:
                    figaxes[-1].axes.get_yaxis().set_ticklabels([])
                if i in [0,1,2,3]:
                    xtext = -0.55
                    ytext = 150
                else:
                    xtext = -0.55
                    ytext = 30

                figaxes[-1].text(xtext,ytext, r'$E_\mathrm{reco} = ['+"%i" % 
                                 one_loader.bin_edges[0][i]  + '-' "%i" % one_loader.bin_edges[0][i+1]+
                                 ' ]\,\mathrm{GeV}$', fontsize = 14, backgroundcolor = 'white')

                figaxes.append(myfig.add_subplot(nestedgs[-1][1,0],rasterization_zorder = -9, 
                                                 sharex = figaxes[-1]))
                plt.errorbar(myxaxis, dmc_1d, yerr = [derrors, derrors], 
                         fmt = 'o', ecolor = 'k', mfc = 'k', label = 'Data', zorder= 99, markersize = 6)
                plt.plot(myxaxis, [1]*len(myxaxis), '--k')
                # print pid, i, 
                # print data_1d 
                # print bestfit_1d

                if i <= (one_loader.bin_edges[0].size-2):
                    figaxes[-1].axes.get_xaxis().set_ticklabels([])
                figaxes[-1].axes.get_yaxis().set_ticks([0.2,1, 1.8])
                plt.ylim(rylim)
                plt.xlim([plot_axis.min(), plot_axis.max()])
                if (i % 2) != 0:
                    figaxes[-1].axes.get_yaxis().set_ticklabels([])

                plt.subplots_adjust(hspace = 0.)
                nestedgs[-1].set_height_ratios([0.6, 0.4])

            gs.update(hspace = 0.18, left =0.1, wspace = 0.15, bottom = 0.1, top = 0.93)
            myfig.text(0.03, 0.5, r'$\mathrm{Events\,and\,Data/MC\,ratio\,per\,energy\,band}$',
                       horizontalalignment='center',
                       verticalalignment='center',
                       rotation = 'vertical',
                       fontsize = 'x-large')

            myfig.text(0.5, 0.04, r'$\cos\,(\theta_\mathrm{reco})$',
                       horizontalalignment = 'center',
                       verticalalignment='center',
                       rotation = 'horizontal',
                       fontsize = 'x-large')   
            #myfig.legend((l1,l2,l3), ('Data','Expectation: best fit', 'Expectation: no osc.'), 
            #             ncol=3, labelspacing=1.5, loc = 9, numpoints = 1, borderaxespad=.5)
            if outdir:
                myfig.savefig(os.path.join(outdir,
                                           'Data_MC_Full_PID_'+ "%i" % pid +'_'+loader_settings['user']+'.png'))
            plt.show()


            #############
            # Doing it in energy bands
            #############

            myfig     = plt.figure(figsize = (4+0.5, 1.*one_loader.bin_edges[0].size+0.5))
            plot_axis = np.cos(one_loader.bin_edges[1])
            myxaxis   = (plot_axis[:-1]+plot_axis[1:])/2.
            figaxes = []
            columns = 1
            rows    = one_loader.bin_edges[0].size-1

            for i in range(0, one_loader.bin_edges[0].size-1):
                if i == 0:
                    figaxes.append(myfig.add_subplot(rows,columns,i+1,rasterization_zorder = -9))
                else:
                    figaxes.append(myfig.add_subplot(rows,columns,i+1,rasterization_zorder = -9, 
                                                     sharex=figaxes[0], sharey=figaxes[0]))
                if i < (one_loader.bin_edges[0].size-2):
                    plt.setp(figaxes[-1].get_xticklabels(), visible=False)

                if pid_used:
                    indices_tuple = (i,slice(None),pid)
                else:
                    indices_tuple = (i,slice(None))
                jplot.unfilledBar(plot_axis, bestfit_histo[indices_tuple], color = 'b', zorder = -9)
                jplot.unfilledBar(plot_axis, nofit_histo[indices_tuple],   color = 'r', zorder = -9)
                plt.errorbar(myxaxis, data_histo[indices_tuple], 
                             yerr = [np.sqrt(data_histo[indices_tuple]),np.sqrt(data_histo[indices_tuple])],
                             fmt = 'o', ecolor = 'k', mfc = 'k', markersize = 6, zorder = 99)
            figaxes[0].invert_xaxis()
            figaxes[0].set_ylim([0, nofit_histo.max()*1.2])
            plt.xlim([plot_axis.min(), plot_axis.max()])


            myfig.text(0.06, 0.5, r'$\mathrm{Events\,per\,energy\,band}$',
                       horizontalalignment='center',
                       verticalalignment='center',
                       rotation = 'vertical',
                       fontsize = 'x-large')

            myfig.text(0.5, 0.04, r'$\cos\,(\theta_\mathrm{reco})$',
                       horizontalalignment = 'center',
                       verticalalignment='center',
                       rotation = 'horizontal',
                       fontsize = 'x-large')   

            plt.subplots_adjust(top = 0.98, right = 0.92, wspace=0., hspace =0., left = 0.2, bottom = 0.1)
            if outdir:
                myfig.savefig(os.path.join(outdir, 
                                           'Data_MC_Ebands_PID_'+ "%i" % pid +'_'+loader_settings['user']+'.png'))
            plt.show()

            #################
            # Doing it in zenith bands
            #################
            myfig     = plt.figure(figsize = (2.5*one_loader.bin_edges[1].size+0.5, 1.5+0.5))
            plot_axis = np.log10(one_loader.bin_edges[0])
            myxaxis   = (plot_axis[:-1]+plot_axis[1:])/2.
            figaxes = []
            columns = one_loader.bin_edges[1].size-1
            rows    = 1
            for i in range(0, one_loader.bin_edges[1].size-1):
                if i == 0:
                    figaxes.append(myfig.add_subplot(rows,columns,one_loader.bin_edges[1].size-i-1,
                                                     rasterization_zorder = -9))
                else:
                    figaxes.append(myfig.add_subplot(rows,columns,one_loader.bin_edges[1].size-i-1,
                                                     rasterization_zorder = -9, sharey = figaxes[0]))
                if i != one_loader.bin_edges[1].size-2:
                    plt.setp(figaxes[-1].get_yticklabels(), visible=False)

                if pid_used:
                    indices_tuple = (slice(None),i,pid)
                else:
                    indices_tuple = (slice(None),i)

                jplot.unfilledBar(plot_axis, bestfit_histo[indices_tuple], color = 'b', zorder = -9)
                jplot.unfilledBar(plot_axis, nofit_histo[indices_tuple],   color = 'r', zorder = -9)
                plt.errorbar(myxaxis, data_histo[indices_tuple], 
                             yerr = [np.sqrt(data_histo[indices_tuple]), np.sqrt(data_histo[indices_tuple])],
                             fmt = 'o', ecolor = 'k', mfc = 'k', markersize = 6, zorder = 99)
            figaxes[0].set_ylim([0, nofit_histo.max()*1.2])
            plt.xlim([plot_axis.min(), plot_axis.max()])

            myfig.text(0.02, 0.5, r'$\mathrm{Events}$',
                       horizontalalignment='center',
                       verticalalignment='center',
                       rotation = 'vertical',
                       fontsize = 'x-large')

            myfig.text(0.5, 0.1, r'$\log_{10}\,(E_\mathrm{reco}/\mathrm{GeV})$',
                       horizontalalignment = 'center',
                       verticalalignment='center',
                       rotation = 'horizontal',
                       fontsize = 'x-large')   

            plt.subplots_adjust(top = 0.9, right = 0.99, wspace=0., hspace =0, left = 0.06, bottom = 0.35)
            if outdir:
                myfig.savefig(os.path.join(outdir,
                                           'Data_MC_Zbands_PID_'+ "%i" % pid +'_'+loader_settings['user']+'.png'))
            plt.show()

## Separate this function - doing the L/E band implies re-doing all the systematic sets
def plotBestFitLE(scan_results = None, outdir = None):

    rylim = [0.6, 1.4]

    nofit_settings = deepcopy(scan_results['result'])
    nofit_settings['dm31'] = nofit_settings['theta23'] = 0.

    for sindex, loader_settings in enumerate(scan_results['loader_settings']):
        #################
        # Doing the LE band - this will not include any KDE stuff
        #################
        loader_settings = original_scan_results['loader_settings'][sindex]
        loader_settings['LEaxis'] = np.linspace(1., 3.3, 41)
        loader_settings['use_kde_bg'] = False
        loader_settings['use_kde_sys'] = False

        one_loader     = dataLoader.dataLoader(**loader_settings)
        bestfit_histo  = one_loader.loadMCasData(scan_results['result'])
        nofit_histo    = one_loader.loadMCasData(nofit_settings)

        if scan_results['data_settings'].has_key('year_data'):
            print 'plotTools: Comparing with data'
            data_histo     = one_loader.loadData(year = scan_results['data_settings']['year_data'])
        else:
            print 'plotTools: Comparing with MC'
            data_histo    = one_loader.loadMCasData(scan_results['data_settings'])


        for pid in range(pid_nbins):

            myfig  = plt.figure(figsize = (7,6))
            plot_axis = one_loader.LEaxis
            myxaxis = (plot_axis[:-1]+plot_axis[1:])/2

            ax3 = myfig.add_subplot(211)

            if pid_used:
                bestfit_1d = bestfit_histo[:,pid]
                nofit_1d   = nofit_histo[:,pid]
                data_1d    = data_histo[:,pid]
            else:
                bestfit_1d = bestfit_histo
                nofit_1d   = nofit_histo
                data_1d    = data_histo

            jplot.unfilledBar(plot_axis, bestfit_1d, color = 'b', label = 'MC best fit', zorder = -9)
            jplot.unfilledBar(plot_axis, nofit_1d, color = 'r', linestyle = '--', 
                              label = 'MC expectation', zorder = -9) 
            plt.errorbar(myxaxis, data_1d, yerr = [np.sqrt(data_1d), np.sqrt(data_1d)], 
                         fmt = 'o', ecolor = 'k', mfc = 'k', label = 'Data', zorder= 99)
            jplot.errorMarkHor(plot_axis, data_1d, color= 'k', label = None, zorder = 1)
            plt.ylabel(r'$\mathrm{Events}$', fontsize = 'x-large')
            #ax3.set_xscale('log')
            plt.setp(ax3.get_xticklabels(), visible=False)

            ax4 = myfig.add_subplot(212, sharex=ax3)
            dmc     = data_1d/bestfit_1d
            derrors = dmc*np.sqrt(data_1d)/data_1d
            plt.errorbar(myxaxis, dmc, yerr = [derrors, derrors], fmt = 'o', ecolor = 'k', 
                         mfc = 'k', label = 'Data', zorder= 99)
            plt.plot(myxaxis, [1]*len(myxaxis), '--k')
            plt.ylim(rylim)
            plt.ylabel(r'$\mathrm{Data / Simulation}$', fontsize = label_fontsize)    
            plt.xlabel(r'$\log_{10}\left(L_\mathrm{reco}/E_\mathrm{reco}\,\mathrm{km}/\mathrm{GeV} \right)$', 
                       fontsize = 'x-large')
            plt.subplots_adjust(wspace=0.2, hspace =0, left = 0.13, bottom = 0.13, top = 0.95, right = 0.97)
            plt.xlim([plot_axis.min(), plot_axis.max()])

            if outdir:
                myfig.savefig(os.path.join(outdir,
                                           'Data_MC_LE_PID_'+ "%i" % pid +'_'+loader_settings['user']+'.png'))
            plt.show()



def plotContributions(loader_settings= None, mc_settings=None, outdir = None):
    import dataLoader
    loader_settings['energy_axis'] = np.linspace(loader_settings['energy_axis'].min(), 
                                                 loader_settings['energy_axis'].max(), 21)
    loader_settings['zenith_axis'] = [0, np.pi]
    loader_settings['energy'] = 'reco_energy'
    loader_settings['zenith'] = 'reco_zenith'

    loader = dataLoader.dataLoader(**loader_settings)
    nu_histograms = loader.getNeutrinoHistograms(mc_settings)
    mc_settings['dm31'] = mc_settings['theta23'] = 0.
    nosc_histograms = loader.getNeutrinoHistograms(mc_settings)

    sumaxis = 1
    plot_axis = loader.energy_axis

    livetime = mc_settings['norm_nu']*3600.*365*24
    numu_osc  = np.sum(nu_histograms['numu_histo'], axis = sumaxis)*livetime
    numu_nosc = np.sum(nosc_histograms['numu_histo'], axis = sumaxis)*livetime
    nutau_osc = np.sum(nu_histograms['nutau_histo'], axis = sumaxis)*livetime
    background = np.sum(nu_histograms['nue_histo'] + nu_histograms['nc_histo'], axis = sumaxis)*livetime

    fig = plt.figure('Final level', figsize = (7, 3.3))
    color_order = ['r', '#b49d0d', '#852FAD', '#50AB46', '#4855E8']
    patterns = [ '+','XX', '\\\\', '//', '-', '*', 'o', 'O', '.']
    scale_factor = 1.85

    ax_e = fig.add_subplot(111, rasterization_zorder=-9)
    l1, a1, r1, = stackedBarPlot(plot_axis, [background, numu_nosc], 
                                 colors = [color_order[1], color_order[-1]], zorder = -12, 
                                 outer_line = True, hatch = [patterns[0], patterns[3]])
    l2, a2, r2, = stackedBarPlot(plot_axis, [background, nutau_osc, numu_osc], 
                                 colors = color_order[1:], outer_line = True, zorder = -9, hatch = patterns)


    rectangle = mpl.patches.Rectangle((0,0), 1, 1, fc = None, linewidth=1, zorder = -15,
                                  edgecolor = 'k')

    plt.legend([r1[1], r2[2], r2[1],  r2[0]], (r'$\nu_\mu\,\mathrm{disappeared}$', 
                                           r'$\nu_\mu$', r'$\nu_\tau\,\mathrm{from\,}\nu_\mu$', 
                                           r'$\nu_e$ + $\nu_\mathrm{NC}$'), 
           loc = 0, ncol = 1, columnspacing = 0.7,fontsize = label_fontsize, numpoints=2 ).set_zorder(-10)

    ltext = gca().get_legend().get_texts()
    for one_text in ltext:
        plt.setp(one_text, zorder = -8)
    
        plt.xlabel(r'$E_\nu\,\mathrm{(GeV)}}$', fontsize = label_fontsize)

    plt.legend(loc = 0)
    #ylim([0, 900])
    plt.ylabel(r'$\mathrm{Events}$', fontsize = label_fontsize)    
    plt.subplots_adjust(wspace=0.2, hspace =.35, left = 0.14, bottom = 0.2, top = 0.95, right = 0.96)
    #ax_e.axes.get_yaxis().set_ticks([0,200, 400, 600, 800])
    if outdir:
        fig.savefig(outdir + 'Sample_energy.png')
    plt.show()

def plotContributions2(hist_list = None, components = ['nue', 'numu', 'nutau', 'nc'], 
                      energy_axis = None, zenith_axis = None, ratios = False, hist_labels = None):
    myfig = plt.figure(figsize=(9,6))
    figaxes = []
    plot_axis = cos(zenith_axis)
    myxaxis = (plot_axis[:-1]+plot_axis[1:])/2
    multiplier = 79988461.37448613

    if ratios:
        text_x = -0.6
        text_y =  1.25
        myticks = [0.5, 1, 1.5]
    else:
        text_x = -0.98
        text_y = 90
        myticks = [0, 50, 100, 150]

    colors = ['#140000', '#660000','#CC0000','#FF4719','#ffc7ba']
    lstyles = ['-', '--', '--', ':']
    mylegends = []
    myartists = []

    all_histograms = []
    for i in range(0,8):
        figaxes.append(myfig.add_subplot(4,2,i+1,rasterization_zorder = -9))
        for hi, h in enumerate(hist_list):
            for ci, hist_key in enumerate(components):
                if ratios:
                    all_histograms.append(h[hist_key+'_histo'][i,:]/hist_list[hist_labels.index('baseline')][hist_key+'_histo'][i,:])
                else:
                    all_histograms.append(h[hist_key+'_histo'][i,:]*multiplier)
                if i == 0 and hist_labels != None:
                    myartists.append(unfilledBar(plot_axis, all_histograms[-1], color = colors[hi], label = hist_key, linestyle = lstyles[ci])[0])
                    mylegends.append(hist_labels[hi])# + '_' + hist_key)
                else:
                    unfilledBar(plot_axis, all_histograms[-1], color = colors[hi], label = hist_key, linestyle = lstyles[ci])

        # Some formatting
        figaxes[-1].text(text_x, text_y, r'$E_\mathrm{reco} = ['+"%i" % energy_axis[i]  + '-' "%i" % energy_axis[i+1]+' ]\,\mathrm{GeV}$', fontsize = 14)#, backgroundcolor = 'white')
        figaxes[-1].axes.get_yaxis().set_ticks(myticks)
        if i == 0 and hist_labels != None:
            figaxes[-1].legend(myartists, mylegends, ncol = 4, bbox_to_anchor=(0., 1.02, 1., .102), loc=3) 
        if (i % 2) != 0:
            figaxes[-1].axes.get_yaxis().set_ticklabels([])
        if i < 6:
            figaxes[-1].axes.get_xaxis().set_ticklabels([])

        #ylim([1, 200])
        #figaxes[-1].set_yscale('log')
        figaxes[-1].invert_xaxis()
    myfig.text(0.5, 0.03, r'$\cos\,(\theta_\mathrm{reco})$',
               horizontalalignment = 'center',
               verticalalignment='center',
               rotation = 'horizontal',
               fontsize = 'x-large')
    myfig.text(0.03, 0.5, 'Ratio to baseline',
               horizontalalignment='center',
               verticalalignment='center',
               rotation = 'vertical',
               fontsize = 'large')
    plt.show()


def plotBands(hist_name = 'allmc', energy_axis = None, zenith_axis = None, ratios = False, **kwargs):
    myfig = plt.figure(figsize=(9,6))
    figaxes = []
    plot_axis = cos(zenith_axis)
    myxaxis = (plot_axis[:-1]+plot_axis[1:])/2
    multiplier = 1.#79988461.37448613

    if ratios:
        text_x = -0.6
        text_y =  1.25
        myticks = [0.5, 1, 1.5]
    else:
        text_x = -0.98
        text_y = 90
        myticks = [0, 50, 100, 150]

    colors = ['#140000', '#660000','#CC0000','#FF4719','#ffc7ba']
    lstyles = ['-', '--', '--', ':']
    mylegends = []
    myartists = []

    all_histograms = []
    for i in range(0,8):
        figaxes.append(myfig.add_subplot(4,2,i+1,rasterization_zorder = -9))
        for hi, h in enumerate(hist_list):
            for ci, hist_key in enumerate(components):
                if ratios:
                    all_histograms.append(h[hist_key+'_histo'][i,:]/hist_list[hist_labels.index('baseline')][hist_key+'_histo'][i,:])
                else:
                    all_histograms.append(h[hist_key+'_histo'][i,:]*multiplier)
                if i == 0 and hist_labels != None:
                    myartists.append(unfilledBar(plot_axis, all_histograms[-1], color = colors[hi], label = hist_key, linestyle = lstyles[ci])[0])
                    mylegends.append(hist_labels[hi])# + '_' + hist_key)
                else:
                    unfilledBar(plot_axis, all_histograms[-1], color = colors[hi], label = hist_key, linestyle = lstyles[ci])

        # Some formatting
        figaxes[-1].text(text_x, text_y, r'$E_\mathrm{reco} = ['+"%i" % energy_axis[i]  + '-' "%i" % energy_axis[i+1]+' ]\,\mathrm{GeV}$', fontsize = 14)#, backgroundcolor = 'white')
        figaxes[-1].axes.get_yaxis().set_ticks(myticks)
        if i == 0 and hist_labels != None:
            figaxes[-1].legend(myartists, mylegends, ncol = 4, bbox_to_anchor=(0., 1.02, 1., .102), loc=3) 
        if (i % 2) != 0:
            figaxes[-1].axes.get_yaxis().set_ticklabels([])
        if i < 6:
            figaxes[-1].axes.get_xaxis().set_ticklabels([])

        #ylim([1, 200])
        #figaxes[-1].set_yscale('log')
        figaxes[-1].invert_xaxis()
    myfig.text(0.5, 0.03, r'$\cos\,(\theta_\mathrm{reco})$',
               horizontalalignment = 'center',
               verticalalignment='center',
               rotation = 'horizontal',
               fontsize = 'x-large')
    myfig.text(0.03, 0.5, 'Ratio to baseline',
               horizontalalignment='center',
               verticalalignment='center',
               rotation = 'vertical',
               fontsize = 'large')
    plt.show()



#     ####
#     ## Parameter scan plots
#     ####
#     mpl.rc('xtick', labelsize=10)  
#     mpl.rc('ytick', labelsize=10)
#     mpl.rc('ytick.minor', size=4)
#     mpl.rc('xtick.major', size=4)

#     fscan = figure('Parameter_scan', figsize=(2*fig_xsize, 2*fig_ysize))

#     fscan.text(0.03, 0.5, r'$\|\Delta m^{2}\|(10^{-3}$eV${}^2)$',
#                horizontalalignment='center',
#                verticalalignment='center',
#                rotation = 'vertical',
#                fontsize = 25)

#     fscan.text(0.5, 0.03, r'$\sin^{2}(2\theta)$',
#            horizontalalignment = 'center',
#            verticalalignment='center',
#            rotation = 'horizontal',
#            fontsize = 25)   

#     xstep = scan_results['theta23_list'][1] - scan_results['theta23_list'][0]
#     xaxis = concatenate((scan_results['theta23_list'] - xstep/2, [scan_results['theta23_list'][-1] + xstep/2]))
#     ystep = scan_results['dm31_list'][1] - scan_results['dm31_list'][0]
#     yaxis = concatenate((scan_results['dm31_list'] - ystep/2, [scan_results['dm31_list'][-1] + ystep/2]))


#     # LLH
#     likelihood = fscan.add_subplot(331)
#     lax = likelihood.pcolor(xaxis, yaxis , scan_results['llh_scan'].T - min_llh)
#     likelihood.set_title('Likelihood scan')
#     likelihood.axis('tight')
#     colorbar(lax)

#     # GAMMA
#     fgamma = fscan.add_subplot(332)
#     gamma_ticks = linspace(-0.05, 0.05, 11)
#     gax = fgamma.pcolor(xaxis, yaxis, scan_results['gamma_scan'].T)#, vmin = -0.05, vmax = 0.05)
#     fgamma.set_title('Gamma change')
#     fgamma.axis('tight')
#     colorbar(gax)#, ticks = gamma_ticks)

#     # NUMU NORM
#     numu_n = fscan.add_subplot(334)
#     numu_ticks = linspace(0.2, 1., 11)
#     muax = numu_n.pcolor(xaxis, yaxis, scan_results['norm_mu_scan'].T)#, vmin = 0.2, vmax = 1.)
#     numu_n.set_title('NuMu/Tau norm')
#     numu_n.axis('tight')
#     colorbar(muax)

#     # NUE NORM
#     nue_n = fscan.add_subplot(335)
#     #nue_ticks = linspace(-0.2, 0.2, 11)
#     eax = nue_n.pcolor(xaxis, yaxis, scan_results['norm_e_scan'].T)#, vmin = -0.2, vmax = 0.2)
#     nue_n.set_title('NuE norm')
#     nue_n.axis('tight')
#     #colorbar(eax, ticks = nue_ticks)
#     colorbar(eax)

#     # CORSIKA
#     cor_n = fscan.add_subplot(336)
#     #cor_ticks = linspace(0., 0.5, 11)
#     cax = cor_n.pcolor(xaxis, yaxis, scan_results['norm_c_scan'].T)#, vmin = 0.)#, vmax = 0.10)
#     cor_n.set_title('CORSIKA norm')
#     cor_n.axis('tight')
#     #colorbar(cax, ticks = cor_ticks)
#     colorbar(cax)
    

#     # DOMEFF
#     if scan_results.has_key('rpik_scan'):
#         # There should also be a scan for different DOM efficiencies
#         rpik = fscan.add_subplot(333)
#         rax = rpik.pcolor(xaxis, yaxis, scan_results['rpik_scan'].T)#, vmin = -10, vmax = 10)
#         rpik.set_title('rPiK')
#         rpik.axis('tight')
#         #cbar = colorbar(deax, ticks = dom_ticks)      
#         colorbar(rax)
#         #cbar.set_clim(-10,10)

#     # DOMEFF
#     if scan_results.has_key('domeff_scan'):
#         # There should also be a scan for different DOM efficiencies
#         dom_eff = fscan.add_subplot(337)
#         dom_ticks = linspace(-10, 10, 11)
#         deax = dom_eff.pcolor(xaxis, yaxis, scan_results['domeff_scan'].T)#, vmin = -10, vmax = 10)
#         dom_eff.set_title('DOM eff')
#         dom_eff.axis('tight')
#         #cbar = colorbar(deax, ticks = dom_ticks)      
#         colorbar(deax)
#         #cbar.set_clim(-10,10)

#     # HOLEICE
#     if scan_results.has_key('holeice_scan'):
#         holeice = fscan.add_subplot(338)
#         hax = holeice.pcolor(xaxis, yaxis , scan_results['holeice_scan'].T)#, vmin = 30., vmax = 100.)
#         holeice.set_title('HoleIce scan')
#         holeice.axis('tight')
#         colorbar(hax)

#     # RQE
#     if scan_results.has_key('rqe_scan'):
#         relqe = fscan.add_subplot(339)
#         rax = relqe.pcolor(xaxis, yaxis , scan_results['rqe_scan'].T)#, vmin = -3., vmax = 3.)
#         relqe.set_title('RQE diff')
#         relqe.axis('tight')
#         colorbar(rax)

#     if outname_base != 'DontSave':
#         fscan.savefig(outname_base + 'scans.png')

#     mpl.rc('xtick', labelsize=18)  
#     mpl.rc('ytick', labelsize=18)
#     mpl.rc('ytick.minor', size=10)
#     mpl.rc('xtick.major', size=10)

#     contfig = figure('Confidence', figsize = (fig_xsize, fig_ysize))
#     llh_scan2 = scan_results['llh_scan'] - min_llh
#     #delta_llh = array([ 1., 2.71]) * 0.5
#     delta_llh = array([ 2.30,4.61]) * 0.5
#     #delta_llh = linspace(0, 6., 30)
#     c1 = contour(scan_results['theta23_list'], scan_results['dm31_list']*10**3, llh_scan2.T, delta_llh, 
#                  linewidths=(2,2), colors=(iccolor,iccolor), linestyles=['--','-'] )
    
#     #clabel(c1, inline=1, fontsize=10)
#     from pylab import ylim, xlim
#     #ylim([1.6, 3])
#     #xlim([0.93, 1])
#     xlabel(r'$\sin^{2}(2\theta)$', fontsize=25)
#     ylabel(r'$\|\Delta m^{2}\|(10^{-3}$eV${}^2)$', fontsize=25)

#     plot(theta23_fit, dm31_fit*10**3, 'xr', markersize=12, markeredgewidth=5)
#     if not scan_results['exp_data'] == 'data': 
#         plot(scan_results['in_mixing_angle'], scan_results['in_mass_diff']*10**3, '+g', markersize=12, markeredgewidth=5)
#     gcf().subplots_adjust(bottom=0.15)

#     # Save 1 plot, 1 scan
#     if outname_base != 'DontSave':
#         contfig.savefig(outname_base + 'contours.png')
#         show()
#         fscan.clf()
#         contfig.clf()
#     show()
#     return fscan, contfig

#     #raw_input()

# '''
# import pickle
# from scanTools import showScanResults, getBestFits
# from numpy import *
# import LLHminimizer
# def bestFitStats(infile_name):
    
#     # Load the stuff inside that file
#     f = open(infile_name)
#     scan_results = pickle.load(f)
#     f.close()

#     fits = getBestFits(scan_results)

#     minimizer = LLHminimizer.fitOscParams(livetime = scan_results['livetime'], exp_data= scan_results['exp_data'] , event_selection= scan_results['event_selection'], 
#                              corsika_template= scan_results['corsika_template'], baseline_mc= scan_results['baseline_mc'], systematic=scan_results['systematic'], llh= scan_results['llh'])

#     hist = minimizer.bestFitHistograms(fits['dm2'], fits['theta23'], fits['corsika'], fits['numu'], fits['nue'], fits['gamma'])

#     return fits, hist
# '''


# if __name__=='__main__':
#     print 'Not defined yet'
    

   
