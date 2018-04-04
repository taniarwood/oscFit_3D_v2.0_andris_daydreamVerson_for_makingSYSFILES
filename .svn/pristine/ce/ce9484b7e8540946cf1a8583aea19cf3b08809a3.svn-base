#!/usr/bin/env python

import os, sys
modules_dir = os.path.dirname(os.path.realpath(__file__)).split('/scripts')[0] + '/modules'
sys.path.append(modules_dir)
modules_dir2 = os.path.dirname(os.path.realpath(__file__)).split('/scripts')[0] + '/user_settings'
sys.path.append(modules_dir2)

import numpy as np
import pickle
import oscFit_settings as user

def doOscFit(indir, task_nr):
    import os
    job_script = '/afs/ifh.de/user/y/yanezjua/oscFit_trunk/scripts/scans/oscFit_oneJob_mctest.py'
    qsub_string = '  '.join(['python',job_script,"%i" % task_nr,'1',indir])
    print qsub_string
    os.system(qsub_string)

def runScan(options):
    import os, sys, pickle
    import numpy as np

    scan_settings = pickle.load(open(os.path.join(user.scans_dir, options.NAME, 
                                               'OscFit_ScanSettings.pckl')))
    print '\nLLH scan: Running with the following fit settings'
    for one_key in scan_settings['fit_settings']:
        print '\t', one_key ,'\t', scan_settings['fit_settings'][one_key]

    total_jobs = scan_settings['dm31_map'].size
    print '\nLLH scan: Total jobs ', total_jobs
    job_script = os.path.join(os.path.dirname(os.path.realpath(__file__)),'oscFit_oneJob_mctest.py')

    if len(options.JOBS) == 0:
        job_range = range(0, total_jobs)
    else:
        job_range = np.array(options.JOBS,dtype=int)-1

    if options.TEST:
        print 'LLH scan: Running first point as a test!' 
        print '\tYou will be asked if you wish to continue before moving on'
        os.system(' '.join(['python',job_script,'1','1',options.NAME]))

        print 'LLH scan: Continue execution? ... '
        raw_input()

    if options.MODE == 'uge_farm':
        #print job_script
        farm_multiplicity = 1

        if len(options.JOBS) == 0:
            total_jobs = np.ceil(total_jobs*1./farm_multiplicity)
            job_array =  '1-'+"%i" % total_jobs
            qsub_line = '  '.join(['qsub -t',job_array,'./oscFit_farmScript.sh',
                                   job_script,"%i" % farm_multiplicity,options.NAME])
            print qsub_line
            os.system(qsub_line)
        else:
            for job_array in options.JOBS:
                qsub_line = '  '.join(['qsub -t',job_array,'./oscFit_farmScript.sh',
                                       job_script,"1",options.NAME])
                print qsub_line
                os.system(qsub_line)

    elif options.MODE == 'iparallel':
        from IPython.parallel import Client, interactive
        import iparallel
        rc = Client(profile = 'sge')
        lview = rc.load_balanced_view()
        result = lview.map_async(doOscFit, [options.NAME]*len(job_range), job_range)
        iparallel.waitOn(result)

    elif options.MODE == 'local':
        import os
        for i in job_range:
            os.system('  '.join(['python',job_script,"%i" % i,'1', options.NAME]))
    print 'Finished!'


def foo_callback(option, opt, value, parser):
  setattr(parser.values, option.dest, value.split(','))

if __name__ == '__main__':
    from optparse import OptionParser


    usage = "usage: %prog [options] inputfile"
    parser = OptionParser(usage)
    parser.add_option("-m", "--mode", default="local",
                      dest="MODE", help = "Select where the scan will be ran (add your own!)")
    parser.add_option("-n", "--name", default=None,
                      dest="NAME", help = "Name of the scan created")
    parser.add_option("-t", "--test", action="store_true", dest='TEST',
                      help='Toggle to run first scan point locally (for testing purposes)')
    parser.add_option('-j', '--jobs',type='string',action='callback', callback=foo_callback,
                      default=[], dest="JOBS")
    (options,args) = parser.parse_args()


    runScan(options)
