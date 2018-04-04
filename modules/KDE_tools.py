import os, sys
import numpy as np
import pickle
from copy import deepcopy
import time

from miscFunctions import compare_dictionaries, keepKeys
np.set_printoptions(precision=2, linewidth=120)
settings_dir = os.path.dirname(os.path.realpath(__file__)).split('/modules')[0] + '/user_settings'
sys.path.append(settings_dir)
try:
    KDE_path = __import__("oscFit_settings").KDE_path
except:
    print "KDE_path is not specified in settings"
    exit()

try:
    sys.path.append(KDE_path)
    from SS_KDE import *
    from SS_KDE.classes import KDE
    print "KDE_tools: KDE tool imported successfully."
except:
    print "KDE_tools: Could not import KDE tool - not found or broken!"
    print "KDE_tools: If you want to use KDEs, you can get the tool in: \n\t\thttp://code.icecube.wisc.edu/svn/sandbox/schoenen/kde" 
    print "KDE_tools: Then adapt the corresponding path to oscFit_settings"
    exit()
    



def smoothByKDE(observables, bin_edges, weights=None, saveFile=None, truth={}, matchKeys=[], fileSuffix="", bootstrap=False, mirroring=True, verbose=False):
        '''
        use Kernel Density Estimation (based on tool in S. Schoenen's sandbox) to make pdfs smooth, store these to database
        and read them from database again if needed. Also allows for bootstrapping the resulting KDEs.
        '''
        if not saveFile == None:
            if len(fileSuffix) == 0:    kde_file = saveFile.split(".pckl")[0]+"_KDE_v2.pckl"
            else:                       kde_file = saveFile.split(".pckl")[0]+"_KDE"+str(fileSuffix)+"_v2.pckl"
        truth = deepcopy(  truth  )
        print 'SMOOTH by KDE ', saveFile
        # reading kde file from database, if available for this set of parameters #
        was_read = False
        print "Reading kdes from file and compare truth dicts ..."
        if not bootstrap and not saveFile == None and os.path.exists(kde_file):
            inFile = open(kde_file)
            kde_in = pickle.load(inFile)
            inFile.close()
            
            for t_i, tIn in enumerate(kde_in["kde"]):
                if verbose: print "\t ... testing dict No. ", (t_i+1)
                test_1 = deepcopy(  truth  )
                test_2 = deepcopy(  tIn["truth"]    )
                if len(matchKeys) == 0:
                    print "\t -> use matchKeys from truth.keys() ..."
                    matchKeys = test_1.keys()
                keepKeys(  test_1, matchKeys, 0  )
                keepKeys(  test_2, matchKeys, 0  )
                if compare_dictionaries(test_1, test_2, verbose=verbose):
                    if verbose: print "\t -> reading kde-smoothed histogram from file "+kde_file+" ..."
                    smooth, edges = tIn["hist"], tIn["edges"]
                    was_read = True
                    print "... done."
                    break


        if not was_read: print "... found none matching =/ ... done."
            
        # ... or do kde now #
        if not was_read:
            starttime = time.time()
            if weights == None:
                weights = np.ones(len(observables[0]))
            if not len(bin_edges) == 3:
                print "KDE smoothing only implemented for 3D histograms with PID dimension not being KDE-ed. This is not the case here -> Skip KDE!"
                atmmu_histo = np.histogramdd(observables, bin_edges)[0]
                return atmmu_histo, bin_edges

            N_average = 10
            print "running KDE on "+str(len(bin_edges[2])-1)+" PID bins separately for "+str(len(observables[0]))+" data points ..."
            smooth  = np.zeros(  [ len(bin_edges[2])-1, len(bin_edges[0])-1, len(bin_edges[1])-1 ]  ) # PID first

            # if lists are empty, skip kde and return empty array #
            if len(observables[0]) < 2:
                print "Skip KDE, since less than 2 values for KDE calculation."
                return np.transpose(smooth, [1,2,0]), bin_edges

             # bootstrap kde if requested #
            if bootstrap:
                size        = len(observables[0])
                print "Bootstrap input values to absolute size: "+str(size)
                
                indices     = np.random.choice(size, size=size, replace=True)
                bs_obs_0    = np.array(observables[0])[indices]
                bs_obs_1    = np.array(observables[1])[indices]
                bs_obs_2    = np.array(observables[2])[indices]
                bs_w        = np.array(weights)[indices]
            else:
                bs_obs_0    = deepcopy(observables[0])
                bs_obs_1    = deepcopy(observables[1])
                bs_obs_2    = deepcopy(observables[2])
                bs_w        = deepcopy(weights)
                
            # loop PID bins #
            for i in range(0, len(bin_edges[2])-1):
                if verbose: print "\t ... initialize ..."
                inittime = time.time()

                mask    = ( np.array(bs_obs_2 > bin_edges[2][i]) ) & ( np.array(bs_obs_2 <= bin_edges[2][i+1]) )
                obs0    = bs_obs_0[mask]
                obs1    = bs_obs_1[mask]
                w       = bs_w[mask]
                myKde_c = KDE([obs0, obs1], weights=w, use_cuda=False, method='silverman')
                myKde_c.calcLambdas(weights=False, weightedCov=False)
                calctime = time.time()

                # calculate kde bin by bin #
                if verbose: print "\t ... evaluate ..."
                for ver in (-1,0,1):
                    for hor in (-1,0,1):
                        if not (  (ver == 0 and hor == 0) or mirroring  ): continue
                        for x in range(0,len(smooth[i])):
                            for y in range(0,len(smooth[i][x])):

                                if hor == 0:        x_min, x_max = bin_edges[0][x], bin_edges[0][x+1]
                                elif hor == -1:     x_min, x_max = 2*bin_edges[0][0] - bin_edges[0][x+1], 2*bin_edges[0][0] - bin_edges[0][x]
                                elif hor == +1:     x_min, x_max = 2*bin_edges[0][-1] - bin_edges[0][x+1], 2*bin_edges[0][-1] - bin_edges[0][x]
                                
                                if ver == 0:        y_min, y_max = bin_edges[1][y], bin_edges[1][y+1]
                                elif ver == -1:     y_min, y_max = 2*bin_edges[1][0] - bin_edges[1][y+1], 2*bin_edges[1][0] - bin_edges[1][y]
                                elif ver == +1:     y_min, y_max = 2*bin_edges[1][-1] - bin_edges[1][y+1], 2*bin_edges[1][-1] - bin_edges[1][y]
                                
                                x_vals              = np.arange(  x_min, x_max, (bin_edges[0][x+1]- bin_edges[0][x]) * 1.0 / (N_average+1)  )[1:]
                                y_vals              = np.arange(  y_min, y_max, (bin_edges[1][y+1]- bin_edges[1][y]) * 1.0 / (N_average+1)  )[1:]
                                
                                eval_grid = np.meshgrid(x_vals, y_vals)
                                
                                myKde_c.kde(  [ eval_grid[0].flatten(),
                                                eval_grid[1].flatten() ], weights=True  )
                                phase_space_elem    = ( bin_edges[0][x+1]-bin_edges[0][x] ) * ( bin_edges[1][y+1]-bin_edges[1][y] )
                                smooth[i][x][y]     += np.sum(  myKde_c.values  ) / (N_average**2)  * phase_space_elem
        
                if verbose: print "Normalization: "+str(np.sum(w) *1.0 / np.sum(smooth[i]))
                smooth[i]   *= np.sum(w) *1.0 / np.sum(smooth[i])

                del myKde_c
                if verbose: print "\t -> "+str(i+1)+". PID bin done. (time initializing: "+str(calctime - inittime)+"s, time evaluating: "+str(time.time() - calctime)+")."
                
            print "... done (time needed: "+str(time.time() - starttime)+"s)."
            smooth = np.transpose(smooth, [1,2,0])

            # saving kde to database #
            if not bootstrap:
                if not saveFile == None:
                    if not was_read:
                        if os.path.exists(kde_file):
                            print "\t => Writting kde-smoothed histogram to existing(!) file "+kde_file+" ... "
                            out         = open(kde_file, 'rb')
                            kde_all     = pickle.load(out)
                            out.close()

                            kde         = {     "hist"          : smooth,
                                                "edges"         : bin_edges,
                                                "truth"         : truth
                                           }
                            kde_all["kde"].append(kde)
                            
                            out         = open(kde_file, 'wb')
                            pickle.dump(kde_all, out)
                            out.close()
                        else:
                            print "\t => Writting kde-smoothed histogram to file "+kde_file+" ..."
                            out     = open(kde_file, 'wb')
                            kde     = { "kde" :  [{     "hist"          : smooth,
                                                        "edges"         : bin_edges,
                                                        "truth"         : truth
                                                   }]   }
                            pickle.dump(kde, out)
                            out.close()
                        print "... done."
                    else:   print "\t => Already have this (similar) truth dict in kde file, so don't save this kde."
                else:       print "\t => Can't write kde-smoothed histograms to file, since no file to save was given!"
            else:           print "\t => Don't save kde to file, since bootstrapped pdfs are not meant to be used more than once."

        return smooth, bin_edges
