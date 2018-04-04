######
# Each sample that you will fit needs a file like this
# It tells the dataLoader and fitter where to find the data files
######
data_dir = '/lustre/fs19/group/icecube/terliuk/fit_files/data_PRD_full_selection/'
# in Madison: /data/user/terliuk/PRD_extend_finalLevel/
mc_sets = {'baseline':{'baseline':'550'}, 
           'domeff':{0.85:'551',
                     0.9:'552',
                     0.95:'553',
                     1.:'550',
                     1.05:'554',
                     1.10:'555',
                     1.15:'556'},
           'hole_ice':{1/100.:'560',
                      1/50.:'550',
                      1/30.:'561',
                      1/40.:'563',
                      0.022:'564',
                      0.018:'565',
                      0.015:'566',
                      0.03:'571',
                      0.0275:'572',
                      0.0125:'573'
                       }, 
           'newhole_ice': {0.25: '570'}
             }

# Declare as True if you have systematic sets. False otherwise
systematic_mc= True

# Directory where your pickle files are located


### This has to be adjusted by the user, depending on the naming convention of the pickle files.
### It helps data_loader to do everything in a loop instead of writing it for each systematic
genie_p1 = {'nue':data_dir+'12',
            'numu':data_dir+'14',
            'nutau':data_dir+'16'}
genie_p2 = '.pckl'

### The old GENIE MC does not cover high energies. You HAVE TO patch it with NuGen stuff
### Define the location of the files and the scaling factors below
### Deciding the crossover energy is up to you. The code will just MERGE THEM with the baseline GENIE. 
### There is no clever interpolation or anything similar.
nugen_nue  = data_dir + '10601.pckl'
nugen_numu = data_dir + '10602.pckl'
# Fudge factors
nugen_nuecc  = 130.
nugen_nuenc  = 130.
nugen_numucc = 140.
nugen_numunc = 120.


# Atm. muons (data and corsika)
atmmu_data_files = [ data_dir + 'IC86_1_muons.pckl',
                     data_dir + 'IC86_2_muons.pckl',
                     data_dir + 'IC86_3_muons.pckl']
atmmu_corsika_files = []
# Data
data = {        } # No data in this user to avoid accidental unblinding



atmmu_sets = {}
pure_noise_files = []