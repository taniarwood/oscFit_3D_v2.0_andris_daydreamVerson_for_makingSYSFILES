# Each sample that you will fit needs a file like this
########################################################
mc_sets = {'baseline':{'baseline':'600'},
           'domeff':{0.88:'601',
                     0.94:'603',
                     0.97:'604',
                     1.:'600',
                     1.03:'605',
                     1.06:'606',
                     1.12:'608'},
           'hole_ice':{#0:'562', # Note that these sets are not really corresponding to the R value of dima
                      0.010:'610', # I changed them to keep the code THE SAME
                      0.015:'611',
                      0.020:'600',
                      0.025:'612',
                      0.030:'613'},
           'hi_fwd':{  -5.:'621',
                      -3.:'622',
                      -1.:'624',
                      0. :'612',
                      1. :'623',
                      2. :'620'},
 } 
# Declare as True if you have systematic sets. False otherwise
########################################################
systematic_mc= True 
# Directory where your pickle files are located
########################################################
data_dir = '/gs/project/ngw-282-ac/trwood/jasper_home/MSU_contain_removed/data/'
sim_dir='/gs/project/ngw-282-ac/trwood/jasper_home/MSU_contain_removed/flat_tania3_DPM_interm_BAKOct12017_jaspert_BKBK/model_tests_auto/GH_qgsjet-II-04/'
########################################################
genie_p1 = {'nue':sim_dir+'Level6.nue.12',
            'numu':sim_dir+'Level6.numu.14',
            'nutau':sim_dir+'Level6.nutau.16'}
#genie_p2 = '.09232015.pckl'
genie_p2 = '.tania.pckl'
##genie_p2 = '.muflux.pckl'
### in JP's version this was .nuNu.pckl
nugen_nue  = None 
nugen_numu = None 
atmmu_data_files = [data_dir + 'Level6.0000.data_bkg1.IC86_2.tania.pckl', 
                     data_dir + 'Level6.0000.data_bkg1.IC86_3.tania.pckl',
                     data_dir + 'Level6.0000.data_bkg1.IC86_4.tania.pckl']
atmmu_data_files_aux = [data_dir + 'Level6.0000.data_bkg2.IC86_2.tania.pckl',
                         data_dir + 'Level6.0000.data_bkg2.IC86_3.tania.pckl',
                         data_dir + 'Level6.0000.data_bkg2.IC86_4.tania.pckl']
atmmu_sets = {}
pure_noise_files = [ ]
data = {}
