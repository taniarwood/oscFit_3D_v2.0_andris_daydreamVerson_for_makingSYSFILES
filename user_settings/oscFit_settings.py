#####
# These settings depend on the system of the user.
# You have to set the paths where your oscillation calculators are located
# You can also set the paths for storing your results
#####
import sys,os

#Tania note to tania:  I have pisa installed on gullimen.. i am not sure of the path to where it went though. i installed it like this:
#also note, it contains prob3++ and PREM_10layer.dat
#pip install --user "git+https://github.com/tarlen5/pisa@master#egg=home

# nuCraft_dir      = '/afs/ifh.de/user/y/yanezjua/scratch/icetray/ICEREC_2013/src/nuCraft/python/nuCraft'
# prob3_directory  = '/afs/ifh.de/group/amanda/software_test/RHEL_6.0_amd64/python-tools/lib/python2.7/site-packages/pisa/pisa/oscillations/prob3'
#prob3_earthmodel = '/afs/ifh.de/group/amanda/software_test/RHEL_6.0_amd64/python-tools/lib/python2.7/site-packages/pisa/pisa/resources/oscillations/PREM_10layer.dat'
# KDE_path = "/afs/ifh.de/user/t/terliuk/scratch/software/KDE_and_similar" # XXX: should include folder SS_KDE or modify import line in modules/KDE_tools.py

resources_dir = os.path.dirname(os.path.realpath(__file__)).split('/user_settings')[0] + '/resources'
# prob3_earthmodel = resources_dir + 'PREM_12layer.dat'
#/gs/project/ngw-282-ac/trwood/oscFit3D_v2.0_Tania_Jan8_Drangon_chi/oscFit3D_v2.0_Tania/
kde_database = '/home/trwood/oscfit_output/kdedb'
detsys_database = '/home/trwood/oscfit_output/histdb'

oscillograms_dir = None #'/afs/ifh.de/user/y/yanezjua/i3scripts/oscFit/oscillograms/'
#nuCraft_dir      = '/afs/ifh.de/user/y/yanezjua/scratch/icetray/ICEREC_2013/src/nuCraft/python/nuCraft'
fits_dir  = '/home/trwood/oscfit_output/msu_repickled_tania/ChiFits/'
scans_dir = '/home/trwood/oscfit_output/msu_repickled_tania/ChiScans/'

# globes_wrapper="/afs/ifh.de/user/t/terliuk/scratch/software/GLoBES/wrapper" # GLoBES wrapper
