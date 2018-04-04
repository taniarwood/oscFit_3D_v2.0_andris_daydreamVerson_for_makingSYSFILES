######
# User file for QUICKSTART
######

# Directory where you have extracted the MC files
data_dir = '/home/jp/projects/icecube/oscFit_data/tables/quickstart_Tania/'

# Leave everything else the same
mc_sets = {'baseline':{'baseline':'550'},
           # 'domeff':{0.9:'552',
           #           1.10:'555',
           #           1.05:'554',
           #           0.95:'553',
           #           1.:'550'},
           # 'hole_ice':{0.022:'564',
           #             1/30.:'561',
           #             1/100.:'560',
           #             1/50.:'550'}
           }

# Declaring the naming convention. Useful in case the user has variations of a selection.
genie_p1 = {'nue':data_dir+'12',
            'numu':data_dir+'14',
            'nutau':data_dir+'16'}
genie_p2 = '.pckl'

atmmu_sets = {}
atmmu_sets = {'baseline': {'baseline':'fakemu_100'},
              'domeff': {1.04:'fakemu_104',
                         1.2:'fakemu_104',
                         1.3:'fakemu_104',
                         0.85:'fakemu_85',
                         1.:'fakemu_100'}
              }

# Atm. muons (data and corsika)
atmmu_data_files    = [data_dir + 'IC86_1_muons.pckl',
                       data_dir + 'IC86_2_muons.pckl',
                       data_dir + 'IC86_3_muons.pckl']
atmmu_data_files_aux    = []

# Data
data = {}

# Pure noise files 
pure_noise_files = []
