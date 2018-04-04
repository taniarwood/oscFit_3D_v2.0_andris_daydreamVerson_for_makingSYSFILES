import sys, os, pickle
sys.path.append('/afs/ifh.de/user/y/yanezjua/scratch/oscFitBak/backup_oscFitScripts/modules/')
sys.path.append('/afs/ifh.de/user/y/yanezjua/i3scripts/systematics/crossSections')

import numpy as np
from icecube import dataio, icetray, dataclasses 
from findIntType import getMaVariations
import cut_values as cv


### Definition of the observables
muon_track       = 'SANTA_Fit_Muon'
hadronic_cascade = 'SANTA_Fit_CascadeHad'
fhlc_name        = 'FirstHLCvertex'
true_nu_name     = 'trueNeutrino'  # You need to get this information somehow
### End of observables definition

i3_basedir    = '/afs/ifh.de/user/y/yanezjua/scratch/datasets/finalLevel'
tables_outdir = '/afs/ifh.de/user/y/yanezjua/scratch/analysis-results/finalLevelTables/PRD_final_2'
if not os.path.isdir(tables_outdir):
    print 'Creating out directory ', tables_outdir
    os.mkdir(tables_outdir)

def checkFiles(flist):
    out_list = []
    faulty_list = []
    for one_file in flist:
        infile = dataio.I3File(one_file)
        faulty_file = False
        while infile.more():
            frame = infile.pop_frame()
            if frame.Stop == frame.Physics:
                continue
            frame.keys()
            try: 
                dummy = frame.Has('I3EventHeader')
            except: 
                faulty_file = True
                break
        if not faulty_file:
            out_list.append(one_file)
        else:
            faulty_list.append(one_file)
            
    print 'Files tested. Good files ', len(out_list), '/', len(flist)
    print 'Files with errors: '
    print faulty_list
    return out_list

def createPickleFile(data_type):

    # Setting the directory, listing the files
    indir = os.path.join(i3_basedir, data_type)
    filenames_all = [os.path.join(indir,f) for f in os.listdir(indir) if f.endswith('.i3.bz2')]

    # Go over the files once to see if they can be used
    # This step can be very slow. Skip it if you have checked this once!
    filename_list  = checkFiles(filenames_all)
    total_files    = len(filename_list)*1.0

    file_counter    = 0
    aux = 1

    # Counter of events in the table
    tc = 0

    if not 'IC86' in data_type and not 'corsika' in data_type:
        from icecube import  NewNuFlux
        import vacuumOscillations as vacOsc
        probCalculator = vacOsc.OscProb()
        flux_model = 'honda2006'
        atmFlux_jason06 = NewNuFlux.makeFlux(flux_model)
        atmFlux_jason12 = NewNuFlux.makeFlux('honda2012_spl_solmin')

        spectral_index = 0.05

    # Information stored for both data and MC
    arr_size = 300000 #Max array size defined. Make sure it's enough for you.
    array_reco_energy  = np.zeros(arr_size)
    array_reco_cascade = np.zeros(arr_size)
    array_reco_track   = np.zeros(arr_size)
    array_reco_zenith  = np.zeros(arr_size)
    array_reco_hlcz    = np.zeros(arr_size)

    if 'IC86' in data_type:
        array_eventid = np.zeros(arr_size)

    if not 'IC86' in data_type and not 'corsika' in data_type:
        array_mc_energy= np.zeros(arr_size)
        array_mc_zenith= np.zeros(arr_size)
        array_mc_weight= np.zeros(arr_size)
        array_mc_weight_ej06= np.zeros(arr_size)
        array_mc_weight_muj06= np.zeros(arr_size)
        array_mc_weight_ej12= np.zeros(arr_size)
        array_mc_weight_muj12= np.zeros(arr_size)
        array_mc_oneweight = np.zeros(arr_size)
        array_mc_interaction= np.zeros(arr_size)
        array_mc_particle_type= np.zeros(arr_size)
        array_mc_maxial_res_weight=np.zeros([arr_size,2])
        array_mc_maxial_qe_weight=np.zeros([arr_size,2])


        interactions_directory = \
           '/afs/ifh.de/user/y/yanezjua/scratch/analysis-results/systematic-checks/crossSections/interaction_dict/'
        if '12' in data_type:
            interactions_lib = pickle.load(open(interactions_directory + 'NuEIntDict3.pckl'))
        elif '14' in data_type:
            interactions_lib = pickle.load(open(interactions_directory + 'NuMuIntDict3.pckl'))
        elif '16' in data_type:
            interactions_lib = pickle.load(open(interactions_directory + 'NuTauIntDict3.pckl'))

    tables_outfile = os.path.join(tables_outdir, data_type.split('/')[-1] + '.pckl')
    print 'Saving to ', tables_outfile
    if os.path.isfile(tables_outfile):
        print 'Tables already exist. Press a key to redo.\n', tables_outfile
        raw_input()

    for one_file in filename_list:
        print one_file
        file_counter += 1
        if file_counter % aux == 0:
            print 'File : +' "%i" % file_counter

        infile = dataio.I3File(one_file)
        for frame in infile:
            if frame.Stop != frame.Physics:
                continue

            # Observables and additional info
            muon     = frame[muon_track]
            hadrons  = frame[hadronic_cascade]
            hlcz     = frame[fhlc_name]

            # Weighting
            if 'IC86' in data_type:
                weight_osc  = 1.
                weight_nosc = 1
            elif 'corsika' in data_type:
                print 'Weighting in CORSIKA is way too complicated now. Did not implement'
                exit()
            else:
                true_neutrino = frame[true_nu_name]

                # Retrieve interaction and particle
                interaction_type = frame['I3MCWeightDict']['InteractionType']
                particle_type    = true_neutrino.pdg_encoding
                nu_energy = true_neutrino.energy
                nu_zenith = true_neutrino.dir.zenith
                nu_theta = pi - nu_zenith
                nu_azimuth = true_neutrino.dir.azimuth

                # Calculate the atm. weight (HONDA)
                norm = (10**9 *frame['I3MCWeightDict']['OneWeight']*\
                            icetray.I3Units.GeV*icetray.I3Units.cm2 * icetray.I3Units.sr)/ \
                            (frame['I3MCWeightDict']['NEvents']* total_files/2.)

                if particle_type > 0: 
                    my_nuetype  = dataclasses.I3Particle.ParticleType.NuE
                    my_numutype = dataclasses.I3Particle.ParticleType.NuMu
                else: 
                    my_nuetype  = dataclasses.I3Particle.ParticleType.NuEBar
                    my_numutype = dataclasses.I3Particle.ParticleType.NuMuBar


                nue_jason06 = nu_energy**(spectral_index)* norm *\
                    atmFlux_jason06.getFlux(my_nuetype, nu_energy, cos(nu_theta))/\
                    (icetray.I3Units.GeV*icetray.I3Units.cm2*icetray.I3Units.s*icetray.I3Units.sr)
                nue_jason12 = nu_energy**(spectral_index)* norm *\
                    atmFlux_jason12.getFlux(my_nuetype, nu_energy, cos(nu_theta))/\
                    (icetray.I3Units.GeV*icetray.I3Units.cm2*icetray.I3Units.s*icetray.I3Units.sr)
                if nue_jason12 == 0.:
                    nue_jason12 = nue_jason06

                numu_jason06 = nu_energy**(spectral_index)* norm *\
                    atmFlux_jason06.getFlux(my_numutype, nu_energy, cos(nu_theta))/\
                    (icetray.I3Units.GeV*icetray.I3Units.cm2*icetray.I3Units.s*icetray.I3Units.sr)
                numu_jason12 = nu_energy**(spectral_index)* norm *\
                    atmFlux_jason12.getFlux(my_numutype, nu_energy, cos(nu_theta))/\
                    (icetray.I3Units.GeV*icetray.I3Units.cm2*icetray.I3Units.s*icetray.I3Units.sr)

                # numu_steven = nu_energy**(spectral_index)*norm*\
                #               steven_flux.MuFlux(nu_energy, cos(nu_theta))/\
                #               (nu_energy*\
                #                icetray.I3Units.GeV*icetray.I3Units.m2*icetray.I3Units.s*icetray.I3Units.sr)
                # nue_steven  = nu_energy**(spectral_index)*norm*\
                #               steven_flux.EFlux(nu_energy, cos(nu_theta))/\
                #               (nu_energy*\
                #                icetray.I3Units.GeV*icetray.I3Units.m2*icetray.I3Units.s*icetray.I3Units.sr)

                if numu_jason12 == 0.:
                    numu_jason12 = numu_jason06

                

                # Also add Ma variations only for CC interactions (NC interactions didn't play a role)
                if interaction_type == 1 and 'nu' in data_type:  # Cannot do this for nugen
                    axial_mass_qe, axial_mass_res = \
                        getMaVariations(idict = interactions_lib,\
                                            energy = nu_energy,\
                                            pdg    = int(particle_type),\
                                            xsecE  = frame['I3MCWeightDict']['Crosssection']/nu_energy)
                else:
                    axial_mass_res = np.zeros(2)
                    axial_mass_qe  = np.zeros(2)
            
            array_reco_energy[tc]  = muon.energy + hadrons.energy
            array_reco_cascade[tc] = hadrons.energy
            array_reco_track[tc]   = muon.energy
            array_reco_zenith[tc]  = muon.dir.zenith
            array_reco_hlcz[tc]    = hlcz.pos.z
                    
            if 'IC86' in data_type:
                array_eventid[tc]  = frame['I3EventHeader'].event_id

            if not 'IC86' in data_type and not 'corsika' in data_type:
                array_mc_energy[tc] = nu_energy
                array_mc_zenith[tc] = nu_zenith
                array_mc_oneweight[tc] = norm
                array_mc_weight_ej06[tc]  = nue_jason06
                array_mc_weight_muj06[tc] = numu_jason06
                array_mc_weight_ej12[tc]  = nue_jason12
                array_mc_weight_muj12[tc] = numu_jason12
                array_mc_interaction[tc]  = interaction_type
                array_mc_particle_type[tc] = particle_type
                array_mc_maxial_res_weight[tc,:] = axial_mass_res
                array_mc_maxial_qe_weight[tc,:]  = axial_mass_qe

            tc += 1

    # All files are done
    if not 'IC86' in data_type and not 'corsika' in data_type:
        cc_bool = array_mc_interaction == 1
        nc_bool = array_mc_interaction == 2
        array_dict = {'CC':{}, 'NC':{}}

        # Charged current
        array_dict['CC']['reco_energy']  = array_reco_energy[cc_bool]
        array_dict['CC']['reco_cascade'] = array_reco_cascade[cc_bool]
        array_dict['CC']['reco_track']  = array_reco_track[cc_bool]
        array_dict['CC']['reco_zenith'] = array_reco_zenith[cc_bool]
        array_dict['CC']['oneweight']   = array_mc_oneweight[cc_bool]
        array_dict['CC']['weight_ej06'] = array_mc_weight_ej06[cc_bool]
        array_dict['CC']['weight_ej12'] = array_mc_weight_ej12[cc_bool] 
        array_dict['CC']['weight_muj06']= array_mc_weight_muj06[cc_bool]
        array_dict['CC']['weight_muj12']= array_mc_weight_muj12[cc_bool]
        array_dict['CC']['energy']      = array_mc_energy[cc_bool]
        array_dict['CC']['zenith']      = array_mc_zenith[cc_bool]
        array_dict['CC']['ptype']       = array_mc_particle_type[cc_bool]
        array_dict['CC']['ma_qe']       = array_mc_maxial_qe_weight[cc_bool]
        array_dict['CC']['ma_res']      = array_mc_maxial_res_weight[cc_bool]
        array_dict['CC']['hlcz']        = array_reco_hlcz[cc_bool]

        # Neutral current
        array_dict['NC']['reco_energy']  = array_reco_energy[nc_bool]
        array_dict['NC']['reco_cascade'] = array_reco_cascade[nc_bool]
        array_dict['NC']['reco_track']  = array_reco_track[nc_bool]
        array_dict['NC']['reco_zenith'] = array_reco_zenith[nc_bool]
        array_dict['NC']['oneweight']   = array_mc_oneweight[nc_bool]
        array_dict['NC']['weight_ej06'] = array_mc_weight_ej06[nc_bool]
        array_dict['NC']['weight_ej12'] = array_mc_weight_ej12[nc_bool] 
        array_dict['NC']['weight_muj06']= array_mc_weight_muj06[nc_bool]
        array_dict['NC']['weight_muj12']= array_mc_weight_muj12[nc_bool]
        array_dict['NC']['energy']      = array_mc_energy[nc_bool]
        array_dict['NC']['zenith']      = array_mc_zenith[nc_bool]
        array_dict['NC']['ptype']       = array_mc_particle_type[nc_bool]
        array_dict['NC']['hlcz']        = array_reco_hlcz[nc_bool]


    else:
        array_dict = {'reco_energy':array_reco_energy[:tc],
                      'reco_cascade':array_reco_cascade[:tc],
                      'reco_track':array_reco_track[:tc],
                      'reco_zenith':array_reco_zenith[:tc],
                      'hlcz':array_reco_hlcz[:tc]} 
    if 'IC86' in data_type:
        array_dict['eventid'] = array_eventid[:tc]

    pickle.dump(array_dict, open(tables_outfile, 'w'))
    
    print 'TOTAL EVENTS KEPT (unweighted): ', tc


if __name__ == '__main__':
    data_type = sys.argv[1]
    createPickleFile(data_type)
