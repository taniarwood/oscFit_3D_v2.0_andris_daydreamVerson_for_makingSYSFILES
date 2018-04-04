######
# Each sample that you will fit needs a file like this
# It tells the dataLoader and fitter where to find the data files
######

mc_sets = {'baseline':{'baseline':'600'}, 
#           'domeff':{0.88:'601',
#                     0.94:'603',
#                     0.97:'604',
#                     1.:'600',
#                     1.03:'605',
#                     1.06:'606',
#                     1.12:'608'},
#           'hole_ice':{#0:'562', # Note that these sets are not really corresponding to the R value of dima
#                      0.010:'610', # I changed them to keep the code THE SAME
#                      0.015:'611',
#                      0.020:'600',
#                      0.025:'612',
#                      0.030:'613'},
#           'hifwd':{  -5.:'621',
#                      -3.:'622',
#                      -1.:'624',
#                      0. :'612',
#                      1. :'623',
#                      2. :'620'},
           # 'hadlight':{1.:'550',
           #              0.95:'568',
           #              1.05:'567'},
           # 'alternative':{'newhi':'570',
           #                #'alt1':'574',
           #                #'alt2':'575',
           #                #'alt3':'576',
           #                'mie':'800',
           #                'new_baseline':'750',
           #                'hadron_raw':'650',
           #                },
           }

# Declare as True if you have systematic sets. False otherwise
systematic_mc= False

# Directory where your pickle files are located
# This set now has a cut on the reduced chi2 of the cascade-like reco
#data_dir = '/lustre/home/jpa14/oscfit_files/02182016/'
data_dir = '/home/trwood/MSU_sample/MSU_sample_sept2016/data/'
#sim_dir = '/Users/trwood/MSU_sample_sept2016/oscfit/'

#sim_dir ='/home/trwood/Baseline_set/'
#sim_dir = '/gs/project/ngw-282-ac/trwood/MSU_sample/MSU_sample_sept2016/oscfit/MSU_tania_repickle/oscfitv2_repickle/'
#sim_dir = '/gs/project/ngw-282-ac/trwood/MSU_sample/MSU_sample_sept2016/oscfit/MSU_tania_repickle/oscfitv2_repickle_analysis_keys_only/'
sim_dir = '/home/trwood/Baseline_set/'
### This has to be adjusted by the user, depending on the naming convention of the pickle files.
### It helps data_loader to do everything in a loop instead of writing it for each systematic
genie_p1 = {'nue':sim_dir+'Level6.nue.12',
            'numu':sim_dir+'Level6.numu.14',
            'nutau':sim_dir+'Level6.nutau.16'}
#genie_p2 = '.09232015.pckl'
genie_p2 = '.10082016.pckl'
#TANIA NOTE>> I forget what this 'genie_p2', what it refers to, and if i ned to change it
#it has a date different/one year newer then the other files. .. (in origininal msc script, 
#see msu_20122014_hole_ice.py

### The old GENIE MC does not cover high energies. You HAVE TO patch it with NuGen stuff
### Define the location of the files and the scaling factors below
### Deciding the crossover energy is up to you. The code will just MERGE THEM with the baseline GENIE. 
### There is no clever interpolation or anything similar.
nugen_nue  = None #data_dir + '10601.pckl'
nugen_numu = None #data_dir + '10602.pckl'
# Fudge factors
# nugen_nuecc  = 130.
# nugen_nuenc  = 130.
# nugen_numucc = 140.
# nugen_numunc = 120.


# Atm. muons (data and corsika)
### bkg1 : invert CC without padding (2+ CC hits) and no VICH cut in inverted sample
### bkg2 : invert CC with    padding (3+ CC hits) and no VICH cut in inverted sample
#atmmu_sets
#atmmu_new_data_files
#atmmu_sets = [ data_dir + 'Level6_IC86.2_data_bkg1.11082016.pckl', 
#                     data_dir + 'Level6_IC86.3_data_bkg1.11082016.pckl', 
#                     data_dir + 'Level6_IC86.4_data_bkg1.11082016.pckl']

#atmmu_histo = {data_dir + 'Level6_IC86.2_data_bkg1.11082016.pckl', 
#                     data_dir + 'Level6_IC86.3_data_bkg1.11082016.pckl', 
#                     data_dir + 'Level6_IC86.4_data_bkg1.11082016.pckl'}

atmmu_data_files = [data_dir + 'Level6.0000.data_bkg1.IC86_2.11082016.pckl', 
                     data_dir + 'Level6.0000.data_bkg1.IC86_3.11082016.pckl', 
                     data_dir + 'Level6.0000.data_bkg1.IC86_4.11082016.pckl']


atmmu_data_files_aux = [data_dir + 'Level6.0000.data_bkg2.IC86_2.11082016.pckl',
                     data_dir + 'Level6.0000.data_bkg2.IC86_3.11082016.pckl',
                     data_dir + 'Level6.0000.data_bkg2.IC86_4.11082016.pckl']

#Level6.0000.data.IC86_4.11082016.pckl		Level6.0005.data.IC86_4.11082016.pckl
#Level6.0000.data_bkg1.IC86_2.11082016.pckl
#

#atmmu_data_files = []
#old atmmu_data_files = [ data_dir + 'Level6.0000.data_bkg2.IC86_2.11082016.pckl', 
#                     data_dir + 'Level6.0000.data_bkg2.IC86_3.11082016.pckl', 
#                     data_dir + 'Level6_IC86.4_data_bkg2.IC86_4.11082016.pckl']

atmmu_sets = {}
#pure_noise_files = {}
pure_noise_files = [ ]
#atmmu_corsika_files = {}
# Data
data = {'2012':data_dir + 'Level6_IC86.2012_data.08102015.pckl',
        '2013':data_dir + 'Level6_IC86.2013_data.08102015.pckl',
        '2014':data_dir + 'Level6_IC86.2014_data.08102015.pckl'}
