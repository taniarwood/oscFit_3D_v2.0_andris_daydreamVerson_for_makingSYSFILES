#####
# These settings depend on the system of the user.
# You have to set the paths where your oscillation calculators are located
# You can also set the paths for storing your results
#####
import sys,os
# nuCraft_dir      = '/afs/ifh.de/user/y/yanezjua/scratch/icetray/ICEREC_2013/src/nuCraft/python/nuCraft'
# prob3_directory  = '/afs/ifh.de/group/amanda/software_test/RHEL_6.0_amd64/python-tools/lib/python2.7/site-packages/pisa/pisa/oscillations/prob3'
#prob3_earthmodel = '/afs/ifh.de/group/amanda/software_test/RHEL_6.0_amd64/python-tools/lib/python2.7/site-packages/pisa/pisa/resources/oscillations/PREM_10layer.dat'
# KDE_path = "/afs/ifh.de/user/t/terliuk/scratch/software/KDE_and_similar" # XXX: should include folder SS_KDE or modify import line in modules/KDE_tools.py

resources_dir = os.path.dirname(os.path.realpath(__file__)).split('/user_settings')[0] + '/resources'
# prob3_earthmodel = resources_dir + 'PREM_12layer.dat'

fits_dir  = '/afs/ifh.de/user/y/yanezjua/scratch/analysis-results/oscFit/bestfits'
scans_dir = '/afs/ifh.de/user/y/yanezjua/scratch/analysis-results/oscFit/scans'
kde_database = '/home/jp/projects/icecube/oscFit_results/kdedb'
detsys_database = '/home/jp/projects/icecube/oscFit_results/histdb'

# globes_wrapper="/afs/ifh.de/user/t/terliuk/scratch/software/GLoBES/wrapper" # GLoBES wrapper
