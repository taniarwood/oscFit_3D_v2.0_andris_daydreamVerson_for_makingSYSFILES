#!/usr/bin/env python

# Interactive setup


import sys, os, shutil, fileinput, pickle
import numpy as np

def get_answer(message):
    valid_answers = ['y','n']

    ans = raw_input(message)
    ans = ans.lower()
    
    while len(ans) > 1 or (ans not in valid_answers):
        ans = raw_input('Answer with a single character (y/n): ')

    if ans.lower() == 'y':
        return True
    else:
        return False
   
def validPath(message):
    valid_path = False
    counter = 0
    while not valid_path:
        path = raw_input(message)
        if counter > 0:
            print 'The path you provided ', path, 'is not valid. Try again.'
        if os.path.isdir(path):
            valid_path = True
    return path

def createPath(path):
    if not os.path.isdir(path):
        try:
            os.mkdir(path)
        except:
            print 'The path ', path, ' could not be created. Verify why!'
            sys.exit()
    return 

def verifySettings(oscFit_path):
    # Go over the user files, and verify the pickle files inside
    user_dir = os.path.join(oscFit_path, 'user_settings')
    sys.path.append(user_dir)

    settings = __import__('oscFit_settings')

    # Create a new directory for the KDE, systematic databases
    if not hasattr(settings, 'kde_database'):
        if not hasattr(settings, 'detsys_database'):
            print '\nThis version of oscFit stores intermediate steps to load data faster'
            message =  'Create a directory in a suitable location and type it down here\n'
            dbdir  = validPath(message)
        else:
            dbdir = settings.detsys_database.split('histdb')[0]
    else:
        dbdir = settings.detsys_database.split('histdb')[0]

    # KDE path/database
    add_kdedb_path = False
    if hasattr(settings, 'KDE_path') and not hasattr(settings, 'kde_database'):
        kde_db = os.path.join(dbdir, 'kdedb')
        createPath(kde_db)
        add_kdedb_path = True

    # Detector systematic histograms database
    add_detsys_path = False
    if not hasattr(settings, 'detsys_database'):
        detsys_path = os.path.join(dbdir, 'histdb')
        createPath(detsys_path)
        add_detsys_path = True
    
    # If no need to create a path, nothing to do here!
    if not add_detsys_path and not add_kdedb_path:
        print 'Settings file seems valid'
        return 

    # Otherwise I need to open the file
    with open(os.path.join(user_dir,'oscFit_settings.py'), 'ab') as f:
        if add_kdedb_path:
            f.write('kde_database = ' + '"'+kde_db+'"')
        if add_detsys_path:
            f.write('detsys_database = ' + '"'+detsys_path+'"')
        f.close()


# Run this script from its own directory
def copyPreviousSettings(oscFit_path):
    
    ## Define the path of your previous oscFit version:
    message =  'Write down your previous oscFit path:\n'
    valid_path = False
    counter = 0
    while not valid_path:
        if counter > 0:
            print 'You provided the path ', old_settings_path, ' but no settings file was found there. Try again!'
        old_oscfit_path = raw_input(message)
        counter += 1
        old_settings_path = os.path.join(old_oscfit_path, 'user_settings')
        settings_file     = os.path.join(old_settings_path, 'oscFit_settings.py')

        if os.path.isfile(settings_file):
                valid_path = True

    new_settings_path = os.path.join(oscFit_path,'user_settings')
    new_settings      = os.path.join(new_settings_path, 'oscFit_settings.py')

    sys.path.append(new_settings_path)

    # Copying the file
    shutil.copyfile(settings_file, new_settings)
    print 'Configuration file successfully copied!\n'

    print 'Proceeding to copy old users. You will be asked one by one if you want to import them.'

    # Figuring out if the users should be copied as well
    user_required_dicts = ['mc_sets','genie_p1','genie_p2','atmmu_sets','data']
    user_required_lists = ['atmmu_data_files','pure_noise_files']

    skip_names = ['oscFit_settings.py','oscFit_settings_template.py']

    user_names = os.listdir(old_settings_path)
    for one_user in user_names:
        if ((one_user in skip_names) or 
            ('~' in one_user) or 
            ('.pyc' in one_user)):
            continue

        message= 'Copy user ' + one_user + ' to new oscFit? (y/n) : '
        if get_answer(message):
            new_user = os.path.join(new_settings_path, one_user)
            shutil.copyfile(os.path.join(old_settings_path, one_user),
                            new_user + '.bak')
            
            # Change holeice for hole_ice in the systematics definition
            with open(new_user, 'wt') as fout:
                with open(new_user + '.bak', 'rt') as fin:
                    for line in fin:
                        fout.write(line.replace('holeice','hole_ice'))
                        if 'atmmu_corsika_files' in line:
                            print 'Not copying the atmmu_corsika_files key. Needs strucutre. See example'
                            continue
            fout.close()
            os.remove(new_user + '.bak')

            # Add information to the new user in case it wasn't there yet
            user_module = __import__(one_user.split('.py')[0])
            with open(os.path.join(new_settings_path, one_user), 'ab') as f:
                for one_dict in user_required_dicts:
                    if not hasattr(user_module, one_dict):
                        f.write('\n'+one_dict + ' = {}')
                for one_list in user_required_lists:
                    if not hasattr(user_module, one_list):
                        f.write('\n'+one_list + ' = []')
                f.close()


def convertPickleFiles(path):
    file_list = os.listdir(path)
    for one_file in file_list:
        if not '.pckl' in one_file:
            print 'Skipping ', one_file
            continue

        infile = open(os.path.join(path, one_file))
        data = pickle.load(infile)
        infile.close()

        # Do not touch files that arent neutrino MC
        if not (data.has_key('energy') and data.has_key('ptype')):
            print 'Skipping ', one_file
            continue
            
        # This is a neutrino MC file. Add the "scattering" field
        if data.has_key('scattering'):
            print 'File was already compatible ', one_file
            continue

        print 'Updating ', one_file
        data['scattering'] = np.zeros_like(data['energy'])

        qel_bool = np.sum(data['ma_qe'], axis=1) > 0
        res_bool = np.sum(data['ma_res'], axis=1) > 0
        dis_bool = ~(qel_bool+res_bool)

        data['scattering'][dis_bool] = 1
        data['scattering'][res_bool] = 2
        data['scattering'][qel_bool] = 3

        pickle.dump(data, open(os.path.join(path, one_file), 'w'))
        
    return

if __name__ == '__main__':

    oscFit_path = os.path.dirname(os.path.realpath(__file__)).split('setup.py')[0]
    print '\n** oscFit3D installer **\n'
    print 'Path of the new oscFit is', oscFit_path
    print 'If you are in a different path, go to that path and re-run setup.py from there\n'

    print 'Pickle files requirements have changed.'
    print 'A new key is needed: "scattering" -> DIS=1, RES=2, QEL=3, COH=4.'
    print 'A converter in this script will solve this.'
    print 'Additional interaction information (Q^2 and x) are optional for making use of new features.'
    print 'You need to re-produce pickle files with that information yourself'


    message = '\nWould you like to let the converter try to get the "scattering" key? '
    if get_answer(message):
        pickle_path = validPath('Type the path where your pickle files are stored\n')
        print 'Starting conversion ...'
        convertPickleFiles(pickle_path)
        print 'Finished!'
    else:
        print 'You will have to do this conversion on your own!\n'

    message= 'Would you like to import your settings from a previous oscFit version? (y/n) : '
    if get_answer( message):
        copyPreviousSettings(oscFit_path)
        verifySettings(oscFit_path)

        print '\noscFit setup done. Now go fit those neutrinos! FIT.THEM.ALL.\n'

    else:
        print '\nInteractive creation of a configuration file coming. You will have to do it manually for now'
        print 'Have a look at the examples and documentation: https://sites.ualberta.ca/~yaezgarz/oscfit/ \n'


