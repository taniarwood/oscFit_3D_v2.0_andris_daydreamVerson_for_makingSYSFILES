#!/usr/bin/env python

import matplotlib
matplotlib.use('Agg')
import os, sys, pickle
modules_dir = os.path.dirname(os.path.realpath(__file__)).split('/scripts')[0] + '/modules'
sys.path.append(modules_dir)
settings_dir = os.path.dirname(os.path.realpath(__file__)).split('/scripts')[0] + '/user_settings'
sys.path.append(settings_dir)
import oscFit_settings as settings
import plotTools as pt


fitname = sys.argv[1].replace('.pckl','')
fitdata = pickle.load(open(os.path.join(settings.fits_dir, fitname + '.pckl')))

newdir = os.path.join(settings.fits_dir, fitname)
if not os.path.isdir(newdir):
    os.mkdir(newdir)
pt.plotBestFit3D(fitdata, outdir = newdir)#, result_key = 'result')
