#################
###
### WARNING! this file is an example of how to produce the pickle files
### It depends very much on the user, and it has been stripped from many features
### Use it as guide to write your own. Consider writing an icetray-module as well.
### Icetray-modules might be easier to write, but they can be much slower than reading ROOT files.
###
#################



from optparse import OptionParser
import sys, os
from numpy import *
from numpy.polynomial.polynomial import polyfit

from icecube import icetray, dataclasses 
import ROOT

# This is necessary to calculate fluxes. If you have a different calculator, substitute this one
ROOT.gSystem.Load('/afs/ifh.de/user/y/yanezjua/scratch/icetray/FluxSteven/libInterpolater')
from ROOT import Flux
steven_flux = Flux() 

# I calculate the earth absorption probability and include it in the weight.
# That part has been removed it from this example script

# Number of reweight GENIE points
rw_points   = 4
# Values of the reweight GENIE points
rw_xvalues  = array([-2,-1,0,1,2])

# Function that receives the y_values of the reweight points
# and fits a quadratic polynomial a + bx + cx^2.
# The function returns only b and c (as a should be one).
def fitMaReweight(in_yvalues):
    if sum(in_yvalues) == 0:
        return array([0,0])
    yvalues = concatenate((in_yvalues[:2], [1.], in_yvalues[2:]))
    fitcoeff = polyfit(rw_xvalues, yvalues, deg = 2)[::-1]
    return fitcoeff[:2]

# Return an id for the scattering type.
# It will be used in future oscFit versions.
def getScatteringType(genie_interaction):
    if genie_interaction.GetLeaf('dis').GetValue(0) > 0:
        return 1
    elif genie_interaction.GetLeaf('res').GetValue(0) > 0:
        return 2
    elif genie_interaction.GetLeaf('qel').GetValue(0) > 0:
        return 3
    elif genie_interaction.GetLeaf('coh').GetValue(0) > 0:
        return 4
    else:
        return 0

# Get the root files
def getROOTfiles(options):
    indir = base_dir + options.level + '/'+ options.dtype +'/'
    filenames = [f for f in os.listdir(indir) if f.endswith('.root')]
    out_filenames = []
    for i in range(0, len(filenames)):
        if os.stat(indir + filenames[i]).st_size > 2000:
            out_filenames.append(indir + filenames[i])
        else:
            print 'Removing file - Too small: ', filenames[i]
            os.remove(indir + filenames[i])


    if options.nfiles > 0:
        return random.sample(out_filenames, options.nfiles)
    else:
        return out_filenames

# Check if the ROOT files are OK
def checkROOTfiles(filenames_full):
    total_files = 0
    filename_list = []
    print 'Checking the files ... '
    print sys.argv

    if True: #len(fnmatch.filter(sys.argv[1:], 'IC86_?*')) > 0:
        print 'Data!'
        for one_file in filenames_full:
            if not one_file[-5:] == '.root':
                print 'noRoot: Skipping file ' + one_file
                continue            
            myfile = TFile(one_file)
            try:
                merged_files = myfile.Get('I3EventHeader')
                merged_files.GetEntry(1)
                filename_list.append(one_file)
                total_files += 1
                myfile.Close()
            except:
                print 'damaged: Skipping file ' + one_file
                print os.stat(one_file).st_size
                myfile.Close()
                #print 'File will be deleted! Cancel to abort'
                #raw_input()
                #os.system('rm ' + one_file)
                
                continue
        print 'Done!\n'    
        return filename_list, total_files

# Go over the files
# 

def runEventSelection(options):

    from icecube import icetray

    data_livetime = 0
    aux = 1
    prefix = ''


    ########################################## PREAMBLE  ##########################################

    # Setting the directories
    filenames_full         = getROOTfiles(options)
    if len(filenames_full) == 0:
        print 'No files!'
        exit()

    # Book keeping
    total_events_kept  = 0
    total_events_range = 0
    tc = 0
    print '\nMax files used: ', options.nfiles

    # Making the file list
    filename_list, total_files = checkROOTfiles(filenames_full)
    n_files = len(filename_list)
    print 'Using ' + str(n_files), '  (', total_files, ') Factor 10 coming from last merger!'
    file_counter    = 0


    ########################################## TABLES MODE  ##########################################
    import pickle

    arr_size = 600000
    array_reco_energy = zeros(arr_size)
    array_reco_cascade = zeros(arr_size)
    array_reco_track = zeros(arr_size)
    array_reco_zenith= zeros(arr_size)
    array_reco_pid = zeros(arr_size)


    if not 'IC86' in options.dtype:
        array_mc_weight= zeros(arr_size)

    if not 'IC86' in options.dtype and not 'corsika' in options.dtype:
        array_mc_energy= zeros(arr_size)
        array_mc_zenith= zeros(arr_size)
        array_mc_weight= zeros(arr_size)
        array_mc_weight_e= zeros(arr_size)
        array_mc_weight_mu= zeros(arr_size)
        array_mc_oneweight = zeros(arr_size)
        array_mc_interaction= zeros(arr_size)
        array_mc_scattering = zeros(arr_size)
        array_mc_genie_y  = zeros(arr_size)
        array_mc_genie_x  = zeros(arr_size)
        array_mc_genie_Q2 = zeros(arr_size)
        array_mc_genie_W  = zeros(arr_size)
        array_mc_particle_type= zeros(arr_size)
        array_mc_maxial_res_weight=zeros([arr_size,2])
        array_mc_maxial_qe_weight=zeros([arr_size,2])

    tables_outdir = '/afs/ifh.de/user/y/yanezjua/scratch/analysis-results/finalLevelTables/2015_new'
    if not os.path.isdir(tables_outdir):
        print 'Creating out directory ', tables_outdir
        os.mkdir(tables_outdir)
    tables_outfile = os.path.join(tables_outdir,options.dtype + '.pckl')
    if os.path.isfile(tables_outfile):
        print 'Tables already exist. Press a key to redo.\n', tables_outfile
        #exit()
        #raw_input()


    ########################################## START WITH SELECTION ##########################################
    for one_file in filename_list:
        file_counter += 1
        if file_counter % aux == 0:
            print 'File : +' "%i" % file_counter
        if 'IC86' in options.dtype:
            run_start = 10e100
            run_stop  = 0.

        myfile = ROOT.TFile(one_file)

        ########################################## RETRIEVE TREES ##########################################

        # Header
        header = myfile.Get('I3EventHeader')
        total_events = int(header.GetEntries())

        # SANTA
        santa_muon    = myfile.Get(sg.santa_muon)
        santa_cascade = myfile.Get(sg.santa_cascade)
        santa_PID     = myfile.Get(sg.santa_pid)


        # MC trees
        if not 'IC86' in options.dtype:
            mcweight     = myfile.Get(sg.weight_name)
            trueNeutrino = myfile.Get(sg.neutrino_name)
            trueMuon     = myfile.Get(sg.muon_name)
	    if not options.nugen:
                # This is only for GENIE simulation
                genie_int    = myfile.Get('GENIE_InteractionType')
		genie_rw     = myfile.Get('GENIE_SystematicsReweight')
		genie_info   = myfile.Get('GENIE_InteractionInfo')

        if 'corsika' in options.dtype:
            mcweight     = myfile.Get('CorsikaWeightMap')

        ########################################## GO EVENT BY EVENT ##########################################
        for i in range(0, total_events):
            
            # DATA LIVETIME
            if 'IC86' in options.dtype:
                header.GetEntry(i)
                event_start_time = header.GetLeaf('time_start_utc_daq').GetValue()
                event_stop_time  = header.GetLeaf('time_end_utc_daq').GetValue()
                if event_start_time < run_start:
                    run_start = event_start_time
                if event_stop_time > run_stop:
                    run_stop  = event_stop_time
            else:
                mcweight.GetEntry(i)
                if not 'corsika' in options.dtype:
                    trueNeutrino.GetEntry(i)
                    trueMuon.GetEntry(i)
		    if not options.nugen:
			try:
			    genie_int.GetEntry(i)
			except:
			    print one_file
			    os.remove(one_file)
			    break
			genie_rw.GetEntry(i)
			genie_info.GetEntry(i)


            # Getting the observables
            santa_muon.GetEntry(i)     
            santa_cascade.GetEntry(i)     
            santa_PID.GetEntry(i)  
            

            ########################################## OBSERVABLES SELECTION #################################

            tracke_estimator = santa_muon.GetLeaf('energy').GetValue(0)
            cascade_estimator = santa_cascade.GetLeaf('energy').GetValue(0)
            energy_estimator = tracke_estimator + cascade_estimator
            zenith_fit       = santa_muon.GetLeaf('zenith').GetValue(0)

            pid_value = santa_PID.GetLeaf('value').GetValue(0)

            ###################################### END OF SELECTION / WEIGHTING! #############################

            if 'IC86' in options.dtype:
                weight_osc = 1.
                weight_nosc = 1
            elif 'corsika' in options.dtype:
                weight_osc = (mcweight.GetLeaf('Weight').GetValue(0)*mcweight.GetLeaf('Polygonato').GetValue(0)/
                                  (mcweight.GetLeaf('TimeScale').GetValue(0) * total_files*100000))
                weight_nosc = weight_osc
            else:
                # Retrieve interaction and particle
                interaction_type = mcweight.GetLeaf('InteractionType').GetValue(0)
                particle_type   = trueNeutrino.GetLeaf('pdg_encoding').GetValue(0)
                nu_energy = trueNeutrino.GetLeaf('energy').GetValue(0)
                nu_zenith = trueNeutrino.GetLeaf('zenith').GetValue(0)

                # Calculate Earth absorption probability
                # earth_probability = earth_abs.SurvivalProb(nu_energy, nu_zenith, particle_type)

                nu_theta = pi - nu_zenith
                nu_azimuth = trueNeutrino.GetLeaf('azimuth').GetValue(0)

                # Calculate the atm. weight (HONDA)
                norm = (earth_probability*10**9 * mcweight.GetLeaf('OneWeight').GetValue(0)*\
                        icetray.I3Units.GeV*icetray.I3Units.cm2*
                        icetray.I3Units.sr)/(mcweight.GetLeaf('NEvents').GetValue(0)* total_files)

                # Calculate the axial mass function
                axmass_res = zeros(rw_points)
                axmass_qe  = zeros(rw_points)
		if options.nugen:
		    scattering_type = 1
		    axial_mass_res = zeros(2)
		    axial_mass_qe  = zeros(2)
		else:
		    for rw_index in range(rw_points):
                        axmass_res[rw_index] = genie_rw.GetLeaf('MaCCRES').GetValue(rw_index)
			axmass_qe[rw_index]  = genie_rw.GetLeaf('MaCCQE').GetValue(rw_index)
		    axial_mass_res = fitMaReweight(axmass_res)
		    axial_mass_qe  = fitMaReweight(axmass_qe)
		    scattering_type = getScatteringType(genie_int)


                if particle_type > 0:
                    # Factor from new simulation
		    if options.nugen:
			norm /= 0.5 
		    else:
			norm /= 0.7
                    numu_weight = (10**-8)*norm*steven_flux.MuFlux(nu_energy, cos(nu_theta))/\
                        (icetray.I3Units.GeV*icetray.I3Units.cm2*icetray.I3Units.s*icetray.I3Units.sr)
                    nue_weight  = (10**-8)*norm*steven_flux.EFlux(nu_energy, cos(nu_theta))/\
                        (icetray.I3Units.GeV*icetray.I3Units.cm2*icetray.I3Units.s*icetray.I3Units.sr)
                else: 
                    # Factor from new simulation
		    if options.nugen:
			norm /= 0.5 
		    else:
			norm /= 0.3
                    numu_weight = (10**-8)*norm*steven_flux.MuBarFlux(nu_energy, cos(nu_theta))/\
                        (icetray.I3Units.GeV*icetray.I3Units.cm2*icetray.I3Units.s*icetray.I3Units.sr)
                    nue_weight  = (10**-8)*norm*steven_flux.EBarFlux(nu_energy, cos(nu_theta))/\
                        (icetray.I3Units.GeV*icetray.I3Units.cm2*icetray.I3Units.s*icetray.I3Units.sr)


            ###################################### END OF WEIGHTING ##########################################


            
            array_reco_energy[tc]     = energy_estimator
            array_reco_cascade[tc]    = cascade_estimator
            array_reco_track[tc]      = tracke_estimator 
            array_reco_zenith[tc]     = zenith_fit
            array_reco_pid[tc]        = pid_value

            if not 'IC86' in options.dtype:
                array_mc_weight[tc] = weight_osc

            if not 'IC86' in options.dtype and not 'corsika' in options.dtype:
                array_mc_energy[tc]   = nu_energy
                array_mc_zenith[tc]   = nu_zenith
                array_mc_oneweight[tc] = norm
                array_mc_weight_e[tc] = nue_weight
                array_mc_weight_mu[tc]     = numu_weight
                array_mc_interaction[tc]   = interaction_type
                array_mc_particle_type[tc] = particle_type
                array_mc_maxial_res_weight[tc,:] = axial_mass_res
                array_mc_maxial_qe_weight[tc,:]  = axial_mass_qe
                array_mc_scattering[tc] = scattering_type
		if options.nugen:
		    if trueMuon.GetLeaf('exists').GetValue(0):
		        muon_energy = trueMuon.GetLeaf('energy').GetValue(0)
		    else:
			muon_energy = 0.    
		    array_mc_genie_y[tc]   = 1 - muon_energy/nu_energy
		    array_mc_genie_x[tc]   = 0.
		    array_mc_genie_Q2[tc]  = 0.
		    array_mc_genie_W[tc]   = 0.

		else:
		    array_mc_genie_y[tc]   = genie_info.GetLeaf('y').GetValue(0)
		    array_mc_genie_x[tc]   = genie_info.GetLeaf('x').GetValue(0)
		    array_mc_genie_Q2[tc]  = genie_info.GetLeaf('Q2').GetValue(0)	
		    array_mc_genie_W[tc]   = genie_info.GetLeaf('W').GetValue(0)

            total_events_kept += weight_nosc
            tc += 1
            if energy_estimator != None and energy_estimator > 7. and energy_estimator < 60. and \
                    cos(zenith_fit) < 0:
                total_events_range += weight_nosc


        ########################################## DATA LIVETIME ##########################################
        if 'IC86' in options.dtype:
            data_livetime += (run_stop - run_start)*10**(-10)
        
        # Closing the file
        myfile.Close()

    ###################################### ALL FILES VISITED, FINISHING ######################################

    array_dict = {}
    array_dict = {'reco_energy':array_reco_energy[:tc],
                  'reco_cascade':array_reco_cascade[:tc],
                  'reco_track':array_reco_track[:tc],
                  'reco_zenith':array_reco_zenith[:tc],
                  'pid':array_reco_pid[:tc]}

    if not 'IC86' in options.dtype and not 'corsika' in options.dtype:
        # CC and NC are not separated any longer
        array_dict['oneweight']   = array_mc_oneweight[:tc]
        array_dict['interaction'] = array_mc_interaction[:tc]
        array_dict['weight']      = array_mc_weight[:tc]
        array_dict['weight_e']    = array_mc_weight_e[:tc]
        array_dict['weight_mu']   = array_mc_weight_mu[:tc]
        array_dict['energy']      = array_mc_energy[:tc]
        array_dict['zenith']      = array_mc_zenith[:tc]
        array_dict['ptype']       = array_mc_particle_type[:tc]
        array_dict['ma_qe']       = array_mc_maxial_qe_weight[:tc]
        array_dict['ma_res']      = array_mc_maxial_res_weight[:tc]
        array_dict['scattering']  = array_mc_scattering[:tc]
        array_dict['GENIE_x']     = array_mc_genie_x[:tc]
        array_dict['GENIE_y']     = array_mc_genie_y[:tc]
        array_dict['GENIE_Q2']    = array_mc_genie_Q2[:tc]
        array_dict['GENIE_W']     = array_mc_genie_W[:tc]

    if 'corsika' in options.dtype:
        array_dict['weight'] = array_mc_weight[:tc]

    pickle.dump(array_dict, open(tables_outfile, 'w'))


    print 'TOTAL EVENTS KEPT (weighted): ', total_events_kept
    print 'TOTAL EVENTS KEPT (unweighted): ', tc
    print 'Total events in desired range: ', total_events_range

    if 'IC86' in options.dtype:
        print 'DATA LIVETIME: ', data_livetime


if __name__ == '__main__':
    from optparse import OptionParser


    usage = "usage: %prog [options] inputfile"
    parser = OptionParser(usage)
    parser.add_option("-t", "--type", default="IC86_1",
                      dest="dtype", help = "Select the type of data you are analyzing (IC86_1, 1450)")
    parser.add_option("-n", "--nfiles", type="int", default=-1,
                      dest="nfiles", help = "Pick an upper limit for the number of files to analyze")
    parser.add_option("-v", "--level", default="finalLevel",
                      dest="level", help = "Pick the level at which the checks are done")  
    parser.add_option("-g", "--nugen", action="store_true", dest='nugen')



    (options,args) = parser.parse_args()
    if len(args) != 0:
        crap = "Got undefined options:"
        for a in args:
            crap += a
            crap += " "
            parser.error(crap)


    runEventSelection(options)
