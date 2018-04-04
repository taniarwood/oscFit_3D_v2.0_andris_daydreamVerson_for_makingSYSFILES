#!/usr/bin/env python
import os, sys

modules_dir2 = os.path.dirname(os.path.realpath(__file__)).split('/scripts')[0] + '/user_settings'
sys.path.append(modules_dir2)

import oscFit_settings as user
import numpy as np

import pickle, matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import jp_mpl as jplot


fit_name = sys.argv[1]

data = pickle.load(open(os.path.join(user.fits_dir, fit_name +'.pckl')))


cm = np.array(data['result']['corr_matrix'])
params = np.array(data['result']['parameters'])


# oscFit v1.+
cm = np.array(data['result']['corr_matrix'])
params_in = np.array(data['result']['parameters'])
params = []
for one_param in params_in:
    if one_param in np.unique(data['result']['covariance'].keys()):
        params.append(one_param)


fig = plt.figure(figsize=(14,10))
ax = fig.add_subplot(111)
colormap = plt.get_cmap('RdBu',21)

p, xp, yp = jplot.apcolor(data=cm, cmap = colormap, vmin=-1, vmax=1)

plt.colorbar()

ax.set_xticks(xp)
ax.set_xticklabels(params)
plt.xticks(rotation=90)
ax.set_yticks(yp)
ax.set_yticklabels(params)
fig.subplots_adjust(bottom = 0.15, right=0.95)
plt.show()

fig.savefig(os.path.join(user.fits_dir, fit_name + '.png'))





