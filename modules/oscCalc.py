import os, sys
from numpy import *
import numpy as np
import sys, os
import pickle
# from miscFunctions import find_nearest
settings_dir = os.path.dirname(os.path.realpath(__file__)).split('/modules')[0] + '/user_settings'
sys.path.append(settings_dir)
import oscFit_settings as settings
#from scipy.ndimage import filters

def propagationDistance(zenith):
    L1 = 19.
    R = 6378.2 + L1
    phi = np.arcsin((1-L1/R)*np.sin(zenith))
    psi = zenith - phi
    return np.sqrt( (R-L1)**2 + R**2 - (2*(R-L1)*R*np.cos(psi)))

class OscCalc(object):
    def __init__(self,
                 oscMode   = 'Prob3', 
                 doTables  = True,
                 nbins     = 150,
             ):


        self.dm21=self.dm31=self.theta23=self.theta12=self.theta13=self.dcp=self.mix_angle= None
        self.dm41=self.theta24=self.theta34=None #sterile neutrino paramters

        self.oscMode   = oscMode
        self.tables    = doTables
        self.ebins     = nbins
        self.zbins     = nbins+10

        if self.oscMode.lower() == 'twoneutrino':
            self.GetProb = self.TwoNeutrino

        elif self.oscMode.lower() == 'vacuum':
            self.GetProb = self.Vacuum
            c_speed = 299792458.           # m/s
            hbar   =  6.58211899*10**(-16) # eV*s
            Nav    =  6.02214179*10**(23)  # 1/mol
            Gf     =  1.166364*10**(-5)    # GeV^-2
            Gfreal = (Gf*(hbar*c_speed/10**7)**3)*10**9# eV*cm
            self.vacuum_factor =  2*10**9*1.6j*10**(-19)*10**(-15)*(4*1.05*10**(-34)*c_speed)**(-1)


        elif self.oscMode.lower() == 'nucraft':
            self.GetProb = self.nuCraft
            sys.path.append(settings.nuCraft_dir)
            global NuCraft
            import NuCraft
            print 'Doing nucraft'
            self.numPrec = 1e-3

        elif self.oscMode.lower() == 'prob3':
            self.GetProb = self.Prob3
            sys.path.append(settings.prob3_directory)
            import BargerPropagator
            detector_depth = 2.
            self.prop_height   = 20.
            self.barger_prop = BargerPropagator.BargerPropagator(settings.prob3_earthmodel, detector_depth)
            self.barger_prop.UseMassEigenstates(False)
            self.barger_prop.SetOneMassScaleMode(False)
            self.barger_prop.SetWarningSuppression(True)
        elif self.oscMode.lower() == "globes_standard":
            sys.path.append(settings.globes_wrapper)
            import GLoBES
            self.GetProb = self.GLoBES_standard
            print "Loading GLoBES"
            curdir = os.getcwd()
            os.chdir(settings.globes_wrapper)
            self.globes_calc =  GLoBES.GLoBESCalculator("test")
            self.globes_calc.InitSteriles(0)
            rad = []; dens = []
            premfile = file(settings.prob3_earthmodel)
            for l in premfile.readlines():
                l_sp = l.split()
                if len(l_sp)==0: continue
                rad.append(float(l_sp[0]))
                dens.append(float(l_sp[1]))
            self.globes_calc.SetEarthModel(rad, dens)
            os.chdir(curdir)  
        elif self.oscMode.lower() == "globes_sterile":
            sys.path.append(settings.globes_wrapper)
            import GLoBES
            self.GetProb = self.GLoBES_sterile
            self.GetProbToNuS = self.GLoBES_sterile_to_nus
            print "Loading GLoBES"
            curdir = os.getcwd()
            os.chdir(settings.globes_wrapper)
            self.globes_calc =  GLoBES.GLoBESCalculator("test")
            self.globes_calc.InitSteriles(2)
            rad = []; dens = []
            premfile = file(settings.prob3_earthmodel)
            for l in premfile.readlines():
                l_sp = l.split()
                if len(l_sp)==0: continue
                rad.append(float(l_sp[0]))
                dens.append(float(l_sp[1]))
            self.globes_calc.SetEarthModel(rad, dens)
            os.chdir(curdir)  
        else: 
            print "Error! Unknown option !"
            exit()
        if self.tables:
            self.GetProbDir = self.GetProb
            self.GetProb = self.TablesProb
            self.prob_maps = {'nue':{},
                              'numu':{}
                              }
            for flavor in self.prob_maps.keys():
                for pdg_code in [12, 14, 16]:
                    for nubar in [1, -1]:
                        self.prob_maps[flavor][nubar*pdg_code] = np.zeros([self.ebins, self.zbins])
                        self.prob_maps[flavor][nubar*pdg_code] = np.zeros([self.ebins, self.zbins])
            self.prob_maps['eedges'] = np.linspace(0., 3., self.ebins+1)
            self.prob_maps['zedges'] = np.linspace(-1, 0.2, self.zbins+1)
            self.prob_maps['ebins'] = (self.prob_maps['eedges'][1:] + self.prob_maps['eedges'][:-1])/2.
            self.prob_maps['zbins'] = (self.prob_maps['zedges'][1:] + self.prob_maps['zedges'][:-1])/2.

            # Taking care of overflow by hand
            self.prob_maps['eedges'][0]  = 0.
            self.prob_maps['eedges'][-1] = np.inf
            self.prob_maps['zedges'][-1] = np.inf

            self.prob_maps['ezbins'] = np.meshgrid(self.prob_maps['ebins'], self.prob_maps['zbins'])

    def setParameters(self,
                      dm21    = 7.53E-5, # PDG2014
                      dm31    = 2.51E-3, # PDG2014 (check the indices)
                      theta23 = np.arcsin(sqrt(1.))/2.,    # PDG 2014
                      theta12 = np.arcsin(sqrt(0.846))/2., # PDG 2014
                      theta13 = np.arcsin(sqrt(0.093))/2., # PDG 2014
                      deltacp = 0.,
                      mix_angle = 1.,
                      #### Sterile neutrinos ###
                      dm41 = 1.0, 
                      theta24 = 0.0, 
                      theta34 = 0.0
                      ):
        self.dm21 = dm21
        self.dm31 = dm31
        self.theta23 = theta23
        self.theta13 = theta13
        self.theta12 = theta12
        self.deltacp = deltacp
        self.mix_angle = mix_angle
        # sterile parameters 
        self.dm41 = dm41
        self.theta24 = theta24
        self.theta34 = theta34
        if self.tables:
            # In table mode tables have to be calculated at the point of changing parameters
            # Calculate each oscillation probability of NuX (orig. flux) to NuX' (oscillated flux)

            tot_bins = len(self.prob_maps['ezbins'][0].flatten())
            #print 'Doing the tables with nbins: ', tot_bins

            # To access any of these probabilities
            # map[from][to][zenith_bin, energy_bin]
            # zenith_bin from self.prob_maps['zbins']
            # energy_bin from self.prob_maps['ebins']

            ebins = 10**self.prob_maps['ezbins'][0].flatten()
            zbins = np.arccos(self.prob_maps['ezbins'][1].flatten())

            nux_to_nue      = self.GetProbDir(nu_energy = ebins,
                                              nu_zenith = zbins,
                                              pdg_encoding = np.array([12]*tot_bins))
            self.prob_maps['nue'][12]  = nux_to_nue[:,0].reshape((len(self.prob_maps['zbins']),
                                                                  len(self.prob_maps['ebins'])))
            self.prob_maps['numu'][12] = nux_to_nue[:,1].reshape((len(self.prob_maps['zbins']),
                                                                  len(self.prob_maps['ebins'])))

            nux_to_nuebar   = self.GetProbDir(nu_energy = ebins,
                                              nu_zenith = zbins,
                                              pdg_encoding = np.array([-12]*tot_bins))
            self.prob_maps['nue'][-12]  = nux_to_nuebar[:,0].reshape((len(self.prob_maps['zbins']),
                                                                      len(self.prob_maps['ebins'])))
            self.prob_maps['numu'][-12] = nux_to_nuebar[:,1].reshape((len(self.prob_maps['zbins']),
                                                                      len(self.prob_maps['ebins'])))
            
            nux_to_numu     = self.GetProbDir(nu_energy = ebins,
                                              nu_zenith = zbins,
                                              pdg_encoding = np.array([14]*tot_bins))
            self.prob_maps['nue'][14]  = nux_to_numu[:,0].reshape((len(self.prob_maps['zbins']),
                                                                   len(self.prob_maps['ebins'])))
            self.prob_maps['numu'][14] = nux_to_numu[:,1].reshape((len(self.prob_maps['zbins']),
                                                                   len(self.prob_maps['ebins'])))
            
            nux_to_numubar  = self.GetProbDir(nu_energy = ebins,
                                              nu_zenith = zbins,
                                              pdg_encoding = np.array([-14]*tot_bins))
            self.prob_maps['nue'][-14]  = nux_to_numubar[:,0].reshape((len(self.prob_maps['zbins']),
                                                                       len(self.prob_maps['ebins'])))
            self.prob_maps['numu'][-14] = nux_to_numubar[:,1].reshape((len(self.prob_maps['zbins']),
                                                                       len(self.prob_maps['ebins'])))

            nux_to_nutau     = self.GetProbDir(nu_energy = ebins,
                                               nu_zenith = zbins,
                                               pdg_encoding = np.array([16]*tot_bins))
            self.prob_maps['nue'][16]  = nux_to_nutau[:,0].reshape((len(self.prob_maps['zbins']),
                                                                    len(self.prob_maps['ebins'])))
            self.prob_maps['numu'][16] = nux_to_nutau[:,1].reshape((len(self.prob_maps['zbins']),
                                                                    len(self.prob_maps['ebins'])))

            nux_to_nutaubar  = self.GetProbDir(nu_energy = ebins,
                                               nu_zenith = zbins,
                                               pdg_encoding = np.array([-16]*tot_bins))
            self.prob_maps['nue'][-16]  = nux_to_nutaubar[:,0].reshape((len(self.prob_maps['zbins']),
                                                                        len(self.prob_maps['ebins'])))
            self.prob_maps['numu'][-16] = nux_to_nutaubar[:,1].reshape((len(self.prob_maps['zbins']),
                                                                        len(self.prob_maps['ebins'])))

            

    def TwoNeutrino(self, 
                    nu_energy   = np.array([24.]), # GeV
                    nu_zenith   = np.array([pi]), # km
                    pdg_encoding = np.array([14])):

        # In the 2 neutrino approach, NuE's are left untouched
        if abs(pdg_encoding[0]) == 12:
            return np.vstack((np.ones(len(nu_energy)), np.zeros_like(nu_energy))).T

        # Only numu goes to nutau
        numu_to_nutau = self.mix_angle*np.sin(1.267*self.dm31*propagationDistance(nu_zenith)/nu_energy)**2

        if abs(pdg_encoding[0]) == 14:
            return np.vstack((np.zeros_like(nu_energy), 1.-numu_to_nutau)).T
        elif abs(pdg_encoding[0]) == 16:
            return np.vstack((np.zeros_like(nu_energy), numu_to_nutau)).T
        print 'Something went wrong, you should not end up here'
        exit()

    def Vacuum(self, 
               nu_energy   = np.array([24.]), # GeV
               nu_zenith   = np.array([pi]), # km
               pdg_encoding = np.array([14])):
        
        s12 = np.sin(self.theta12)
        s13 = np.sin(self.theta13)
        s23 = np.sin(self.theta23)

        c12 = np.cos(self.theta12)
        c13 = np.cos(self.theta13)
        c23 = np.cos(self.theta23)

        nu_baseline = propagationDistance(nu_zenith)

        dm21exp = np.exp(-self.vacuum_factor*self.dm21*nu_baseline/nu_energy)
        dm31exp = np.exp(-self.vacuum_factor*self.dm31*nu_baseline/nu_energy)

        # Nue to x
        A_nue_to_nue = (c12**2*c13**2 + s12**2*c13**2*dm21exp+dm31exp*s13**2)
        P_nue_to_nue = np.real(A_nue_to_nue*conjugate(A_nue_to_nue))
        A_nue_to_numu = (-s12*c23-c12*s23*s13)*c12*c13 + \
                        (c12*c23 - s12*s23*s13)*dm21exp*s12*c13 + s23*c13*dm31exp*s13
        P_nue_to_numu = np.real(A_nue_to_numu*conjugate(A_nue_to_numu))
        A_nue_to_nutau = (s12*s23 - c12*c23*s13)*c12*c13 + (-c12*s23-s12*c23*s13)*dm21exp*s12*c13 + \
                         c23*c13*dm31exp*s13
        P_nue_to_nutau = np.real(A_nue_to_nutau*conjugate(A_nue_to_nutau))

        # Numu to x
        A_numu_to_nue = (-s12*c23-c12*s23*s13)*c12*c13 + \
                        (c12*c23 - s12*s23*s13)*dm21exp*s12*c13 + s23*c13*dm31exp*s13
        P_numu_to_nue = np.real(A_numu_to_nue*conjugate(A_numu_to_nue))
        A_numu_to_numu = (-s12*c23 -c12*s23*s13)**2 + dm21exp*(c12*c23 - s12*s23*s13)**2 + s23**2*c13**2*dm31exp
        P_numu_to_numu = np.real(A_numu_to_numu*conjugate(A_numu_to_numu))
        A_numu_to_nutau = (s12*s23 - c12*c23*s13)*(-s12*c23 - c12*s23*s13) + \
                          (-c12*s23 - s12*c23*s13)*dm21exp*(c12*c23-s12*s23*s13) + c23*c13**2*dm31exp*s23
        P_numu_to_nutau = np.real(A_numu_to_nutau*conjugate(A_numu_to_nutau))


        if abs(pdg_encoding[0])  == 12:
            return np.vstack([P_nue_to_nue, P_numu_to_nue]).T
        elif abs(pdg_encoding[0])  == 14:
            return np.vstack([P_nue_to_numu, P_numu_to_numu]).T
        elif abs(pdg_encoding[0])  == 16:
            return np.vstack([P_nue_to_nutau, P_numu_to_nutau]).T
        print 'Something went wrong!'
        exit()

    def nuCraft(self,
               nu_energy    = np.array([23., 32.]), # GeV
               nu_zenith    = np.array([pi,  7*pi/8.]), # km
               pdg_encoding = np.array([ -14, -12])
               ):

        self.nucraft = NuCraft.NuCraft((1., self.dm21, self.dm31), 
                                       [(1,2,self.theta12*180./np.pi),(1,3,self.theta13*180./np.pi, 
                                                                       self.deltacp*180./np.pi),
                                        (2,3,self.theta23*180./np.pi)])
        
        # Calculate the transition probability of an (anti-)particle with give energy and zenith
        NuE_to_X   = np.array(self.nucraft.CalcWeights((np.sign(pdg_encoding)*12, 
                                                     nu_energy, nu_zenith), numPrec = self.numPrec))
        NuMu_to_X  = np.array(self.nucraft.CalcWeights((np.sign(pdg_encoding)*14, 
                                                     nu_energy, nu_zenith), numPrec = self.numPrec))

        # Return the resulting flux for the desired particle
        if abs(pdg_encoding[0])  == 12:
            return np.vstack((NuE_to_X[:,0], NuMu_to_X[:,0])).T
        elif abs(pdg_encoding[0])  == 14:
            return np.vstack((NuE_to_X[:,1], NuMu_to_X[:,1])).T
        elif abs(pdg_encoding[0])  == 16:
            return np.vstack((NuE_to_X[:,2], NuMu_to_X[:,2])).T


    def Prob3(self,
              nu_energy    = np.array([23., 32.]), # GeV
              nu_zenith    = np.array([pi,  7*pi/8.]), # km
              pdg_encoding = np.array([ -14, -12])
              ):

        NuE = 1; NuMu = 2; NuTau = 3;
        out_nu = (  NuE * (np.abs(pdg_encoding) == 12) +
                   NuMu * (np.abs(pdg_encoding) == 14) +
                  NuTau * (np.abs(pdg_encoding) == 16))

        sin_sq_theta12 = np.float64( np.sin(self.theta12)**2 ) 
        sin_sq_theta13 = np.float64( np.sin(self.theta13)**2 ) 
        sin_sq_theta23 = np.float64( np.sin(self.theta23)**2 ) 
        dm32 = self.dm31 - self.dm21

        kSquared = True
        kNuType = np.array(np.sign(pdg_encoding), dtype=np.int)
        cos_nu_zenith = np.cos(nu_zenith)

        osc_prob = np.zeros([len(nu_energy), 2])
        for index in range(len(nu_energy)):
            self.barger_prop.SetMNS(sin_sq_theta12, sin_sq_theta13, sin_sq_theta23,
                                    self.dm21, dm32, self.deltacp,
                                    nu_energy[index], kSquared, int(kNuType[index]))
            self.barger_prop.DefinePath(cos_nu_zenith[index], self.prop_height)
            self.barger_prop.propagate(int(kNuType[index]))

            osc_prob[index, :] = [self.barger_prop.GetProb(NuE,int(out_nu[index])), 
                                  self.barger_prop.GetProb(NuMu,int(out_nu[index]))]

        return osc_prob

    def GLoBES_standard(self,
              nu_energy    = np.array([23., 32.]), # GeV
              nu_zenith    = np.array([pi,  7*pi/8.]), # km
              pdg_encoding = np.array([ -14, -12])
              ):

        cur_params = np.array( [self.theta12, self.theta13, self.theta23, 
                            self.deltacp, self.dm21, self.dm31], dtype=float )
        self.globes_calc.SetParametersArr(cur_params)

        ini_type = np.array(np.sign(pdg_encoding)*np.ones( len(pdg_encoding) ) , dtype = int)
        
        nue_to_nuX = self.globes_calc.MatterProbPDGArr( 12*ini_type, np.array(pdg_encoding, dtype = int),  nu_energy, np.cos(nu_zenith) )
        numu_to_nuX = self.globes_calc.MatterProbPDGArr( 14*ini_type, np.array(pdg_encoding, dtype = int),  nu_energy, np.cos(nu_zenith) )
        
        return np.vstack((nue_to_nuX, numu_to_nuX)).T

    def TablesProb(self,
                   nu_energy    = np.array([12, 21]),   # Ebin
                   nu_zenith    = np.array([12,  45]),  # Zbin
                   pdg_encoding = np.array([ -14, -12])
               ):

        # Return the probability of going into that flavor.
        # Make use of the fact that I only do one flavor at a time, but nu/nubar are mixed

        nubar = pdg_encoding < 0
        ptype = abs(pdg_encoding[0])

        osc_prob = np.zeros([len(nu_energy), 2])

        # Neutrinos
        osc_prob[~nubar,0] = self.prob_maps['nue'][ptype][nu_zenith[~nubar], nu_energy[~nubar]]
        osc_prob[~nubar,1] = self.prob_maps['numu'][ptype][nu_zenith[~nubar], nu_energy[~nubar]]

        # Antineutrinos
        osc_prob[nubar,0] = self.prob_maps['nue'][(-1)*ptype][nu_zenith[nubar], nu_energy[nubar]]
        osc_prob[nubar,1] = self.prob_maps['numu'][(-1)*ptype][nu_zenith[nubar], nu_energy[nubar]]

        return osc_prob

    def GLoBES_sterile(self,
              nu_energy    = np.array([23., 32.]), # GeV
              nu_zenith    = np.array([pi,  7*pi/8.]), # km
              pdg_encoding = np.array([ -14, -12])
              ):
        params = np.array( [self.theta12, self.theta13, self.theta23, 
                                       self.deltacp, self.dm21, self.dm31, 
                                       self.dm41, 0.0,  self.theta24, self.theta34, 0.0, 0.0 ], dtype=float )
        self.globes_calc.SetParametersArr(params ) 
        
        ini_type = np.array(np.sign(pdg_encoding)*np.ones( len(pdg_encoding) ) , dtype = int)
        nue_to_nuX = self.globes_calc.MatterProbPDGArr( 12*ini_type,  np.array(pdg_encoding, dtype = int), nu_energy, np.cos(nu_zenith) )
        numu_to_nuX = self.globes_calc.MatterProbPDGArr( 14*ini_type,  np.array(pdg_encoding, dtype = int),  nu_energy, np.cos(nu_zenith) )
        return np.vstack((nue_to_nuX, numu_to_nuX)).T

    def GLoBES_sterile_to_nus(self,
              nu_energy    = np.array([23., 32.]), # GeV
              nu_zenith    = np.array([pi,  7*pi/8.]), # km
              pdg_encoding = np.array([ -14, -12])
              ):
        params = np.array( [self.theta12, self.theta13, self.theta23, 
                            self.deltacp, self.dm21, self.dm31, 
                            self.dm41, 0.0,  self.theta24, self.theta34, 0.0, 0.0 ], dtype=float )
        self.globes_calc.SetParametersArr(params ) 
        
        ini_type = np.array(np.sign(pdg_encoding)*np.ones( len(pdg_encoding) ) , dtype = int)
        nuX_to_nue= self.globes_calc.MatterProbPDGArr( np.array(pdg_encoding, dtype = int), 12*ini_type, nu_energy, np.cos(nu_zenith) )
        nuX_to_numu= self.globes_calc.MatterProbPDGArr( np.array(pdg_encoding, dtype = int), 14*ini_type, nu_energy, np.cos(nu_zenith) )
        nuX_to_nutau= self.globes_calc.MatterProbPDGArr( np.array(pdg_encoding, dtype = int), 16*ini_type, nu_energy, np.cos(nu_zenith) )
        return np.vstack((1.0 - nuX_to_nue - nuX_to_numu - nuX_to_nutau)).T   
