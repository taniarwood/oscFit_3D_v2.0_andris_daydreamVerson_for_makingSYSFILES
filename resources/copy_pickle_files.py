import pickle
import numpy as np
import sys, os, shutil

indir_base = '/home/jp/projects/icecube/oscFit_data/tables'
indir  = os.path.join(indir_base, 'oscFit_quickstart_data')
outdir = os.path.join(indir_base, 'quickstart_Tania')

infiles = os.listdir(indir)
for one_file in infiles:
    data = pickle.load(open(os.path.join(indir, one_file)))

    # Is it neutrino MC?
    if not data.has_key('weight_mu'):
        shutil.copy(os.path.join(indir, one_file),
                    os.path.join(outdir, one_file))
        continue

    print 'Doing ', one_file
    # It is, now get the new key
    #data['weight_mu_pi'] = np.zeros_like(data['weight_mu'])
    data['weight_mu_pi'] = data['weight_mu']

    pickle.dump(data, open(os.path.join(outdir, one_file),'w'))
