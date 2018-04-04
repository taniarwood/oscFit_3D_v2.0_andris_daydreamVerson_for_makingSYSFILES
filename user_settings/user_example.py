######
# Each sample that you will fit needs a file like this
# It tells the dataLoader and fitter where to find the data files
######

mc_sets = {'baseline':{'baseline':'500'}, 
           # # If you don't have systematic sets, comment them out.
           # 'domeff':{0.85:'501',
           #           0.9:'502',
           #           0.95:'503',
           #           1.:'500',
           #           1.05:'504',
           #           1.10:'505',
           #           1.15:'506'},
           # 'hole_ice':{0:'530',
           #            1/100.:'510',
           #            1/50.:'500',
           #            1/30.:'520'},
           # 'hadlight':{1.:'500',
           #             0.95:'596',
           #             1.05:'597'}
           }



# Declare as True if you have systematic sets. False otherwise
systematic_mc= True

# Directory where your MC and data pickle files are located
data_dir = '/afs/ifh.de/user/y/yanezjua/scratch/analysis-results/finalLevelTables/2015_new/'

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

# Fudge factors. My fudge factors
nugen_nuecc  = 1.3
nugen_nuenc  = 1.2
nugen_numucc = 1.4
nugen_numunc = 1.1


# Atm. muons (data and corsika)
atmmu_data_files = [ data_dir + 'IC86_1_muons.pckl',
                     data_dir + 'IC86_2_muons.pckl',
                     data_dir + 'IC86_3_muons.pckl']
atmmu_corsika_files = []
# Data
data = {'2011':data_dir + 'IC86_1.pckl',
        '2012':data_dir + 'IC86_2.pckl',
        '2013':data_dir + 'IC86_3.pckl'}
# Pure noise files 
pure_noise_files = ['MyFileOfPureNoise.pckl']

atmmu_sets = {}