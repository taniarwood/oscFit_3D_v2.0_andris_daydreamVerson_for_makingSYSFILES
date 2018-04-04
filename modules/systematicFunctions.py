import numpy as np
from copy import deepcopy
from miscFunctions import cartesian
from numpy import polyfit as polyfit
#from scipy.interpolate import interp1d
import scipy
from scipy.optimize import curve_fit

############
# Auxiliary, since there is no exponential class predefined in numpy like there is for polynomials
# M. Leuermann
#############

class exp_func:
    def __init__(self,a, b, c, x0):
        self.a, self.b, self.c, self.x0 = a, b, c, x0
        
    def __call__(self,x):
        return self.eval(x, self.a, self.b, self.c, self.x0)
        
    def derive(self,x):
        return np.array([np.exp(self.b*(x-self.x0)/self.x0), 
                         self.a*self.b/self.x0*np.exp(self.b*(x-self.x0)/self.x0),
                         0])
        
    @staticmethod
    def eval(x, a, b, c, x0):
        return c + a*np.exp(b*(x-x0)/x0)

    @property
    def coeffs(self):
        return [self.a, self.b, self.c]

    @coeffs.setter
    def coeffs(self, val):
        self.a, self.b, self.c = val

class dummy_func:
    def __init__(self, c):
        self.c = c

    def __call__(self, dummy):
        return self.c

############
# Detection uncertainties
# J.P. Yanez
#############


def getSystematicFunctions(y_values, y_valuesw2, xpoints, xi_nominal = None,
                           poly_degree = 2):
    '''
    Rewrite the description

    Parameters
    ----------
    y_values (w2): 

    xpoints: point in x at which the histograms are evaluated

    poly_degree: degree of the polinomial used in the fit

    free_norm: removes normalization, fit only the shapes (set to False to reproduce MC)
    '''

    # Reverting to the code used in the PRL result. Fixing the number of dimensions
    # Keeping the function general enough to do Atmospheric muons as well as neutrinos
    # The bins of y_values are [xpoints, energy, zentih, pid]

    # # Initialize the functions
    if poly_degree <= 0:
        one_fit = exp_func(1.0, 0.0,0.0, x0=xpoints[xi_nominal])
        cov     = [1.0]
        ncoeffs = 3
    else:
        one_fit = np.poly1d([1.0])
        cov     = [1.0]
        ncoeffs = poly_degree + 1

    poly_matrix = np.zeros(list(y_values.shape[1:]) + [ncoeffs])
    cov_matrix  = np.zeros(list(y_values.shape[1:]) + [ncoeffs, ncoeffs])

    # Parameterizing as deviations from the nominal MC
    nominal_y  = y_values[xi_nominal,:]
    werrors    = np.sqrt(y_valuesw2)/nominal_y
    werrors[werrors==0] = 1.
    y_values   /= nominal_y

    # This is my error (sigma) - removing zero errors
        

    for i in range(y_values.shape[1]): # Looping over energy
        for j in range(y_values.shape[2]): # Looping over the zenith
            for k in range(y_values.shape[3]): # Looping over the particle ID

                # Polynomial fits
                if poly_degree > 0:
                    one_fit, cov = polyfit(xpoints,
                                           y_values[:,i,j,k],
                                           deg = poly_degree,
                                           rcond = None,
                                           full = False,
                                           cov  = True,
                                           w = 1./werrors[:,i,j,k])
                # Exponential fits
		#ran out of trys.. got error about max iterations being reached with happens in bins with
		# too few events.  at least this way you can run witha  default value until a better way is 
		#found
                else:
                    one_fit = np.array([1,0,0])
		    cov =  np.identity(3)
		    try:
			f = lambda x,a,b,c: exp_func.eval(x,a,b,c,x0=xpoints[xi_nominal])
                    	one_fit, cov = curve_fit(f, xpoints, y_values[:,i,j,k],
                                             sigma = werrors[:,i,j,k],
                                             p0 = (0.3, -15.0, 0.9))
                    except: pass
                    # one_fit = exp_func(*opt, x0=xpoints[xi_nominal])
                poly_matrix[i,j,k,:] = one_fit
                cov_matrix[i,j,k] = cov

                if np.inf in cov or np.sum(cov) == np.inf or np.sum(cov) == 0:
                    print 'This is COV', cov
                    print 'These are the errors', werrors[:,i,j,k]

    return poly_matrix, cov_matrix

############
# Cross section uncertainties
# J.P. Yanez & Teppei Katori
#############

def axialMassVar(coeff = np.zeros(2), Ma = 0.):
    # The 1 is summed because I simplified the formula to save space
    return 1 + coeff[:,0]*Ma**2 + coeff[:,1]*Ma

# From Teppei Katori (p.11 in https://drive.google.com/file/d/0B8TQi1F3KxYmQkwwU0VHOTdnNEU/view)
def tkDISreweight(a = 0., b = 1., bjorken_x = np.zeros(1)):
    return b*bjorken_x**(-a)


############
# Sub-leading atmospheric flux uncertainties
# J.P. Yanez & A. Terliuk
#############

# See the ipython notebook for the details of how these functions were derived
# There are two basic functions: a gaussian bell and a Log-Log parameterization


def norm_fcn(x, A, sigma = 0.3):
    #x *= 0.9
    return A/np.sqrt(2*np.pi*sigma**2) * np.exp(-x**2/(2*sigma**2))
def LogLogParam(energy = 1., y1 = 1., y2 = 1.,
                x1=0.5, x2=3.,
                cutoff_value = False):
    nu_nubar = np.sign(y2)
    if nu_nubar == 0.0:
        nu_nubar = 1.
    y1 = np.sign(y1)*np.log10(np.abs(y1)+0.0001)
    y2 = np.log10(np.abs(y2+0.0001))
    modification = nu_nubar*10**(((y2-y1)/(x2-x1))*(np.log10(energy)-x1)+y1-2.)
    if cutoff_value:
        modification *= np.exp(-1.*energy/cutoff_value)
    return modification

# These parameters are obtained from fits to the paper of Barr
# E dependent ratios, max differences per flavor (Fig.7)
e1max_mu = 3.
e2max_mu = 43
e1max_e  = 2.5
e2max_e  = 10
e1max_mu_e = 0.62
e2max_mu_e = 11.45
# Evaluated at
x1e = 0.5
x2e = 3.

# Zenith dependent amplitude, max differences per flavor (Fig. 9)
z1max_mu = 0.6
z2max_mu = 5.
z1max_e  = 0.3
z2max_e  = 5.
nue_cutoff  = 650.
numu_cutoff = 1000.
# Evaluated at
x1z = 0.5
x2z = 2.


# This is the neutrino/antineutrino ratio modification for NuMu
def ModNuMuFlux(energy, czenith,
                e1=1., e2=1., z1=1., z2=1.):
    A_ave = LogLogParam(energy=energy, 
                        y1=e1max_mu*e1, 
                        y2=e2max_mu*e2,
                        x1=x1e, x2=x2e)
    A_shape = 2.5*LogLogParam(energy=energy, 
                              y1=z1max_mu*z1, 
                              y2=z2max_mu*z2,
                              x1=x1z, x2=x2z, 
                              cutoff_value = numu_cutoff)
    return A_ave - (norm_fcn(czenith, A_shape, 0.32) - 0.75*A_shape)

# The NuE flux modification requires that you know the NuMu parameters
# The uncertainty is added to that of NuMu. Assuming they are correlated
def ModNuEFlux(energy, czenith,
               e1mu=1., e2mu=1., z1mu=1., z2mu=1.,
               e1e=1., e2e=1., z1e=1., z2e=1.):

    A_ave = LogLogParam(energy=energy, 
                        y1=e1max_mu*e1mu + e1max_e*e1e, 
                        y2=e2max_mu*e2mu + e2max_e*e2e,
                        x1=x1e, x2=x2e)
    A_shape = 1.*LogLogParam(energy=energy, 
                             y1=z1max_mu*z1mu + z1max_e*z1e, 
                             y2=z2max_mu*z2mu + z2max_e*z2e,
                             x1=x1z, x2=x2z,
                             cutoff_value = nue_cutoff)
    return A_ave - (1.5*norm_fcn(czenith, A_shape, 0.4) - 0.7*A_shape)


def modRatioNuMuBar(particle_list, nu_nubar, nubar_mu):
    modfactor = np.ones(len(particle_list['energy']))
    modfactor = nubar_mu * ModNuMuFlux(particle_list['energy'],
                            np.cos(particle_list['zenith']),
                            1.0, 1.0,1.0,1.0)


    weightmod = 1. + modfactor*nu_nubar
    antineutrinos = particle_list['ptype']<0
    weightmod[antineutrinos] = 1./(1+(1-nu_nubar)*modfactor[antineutrinos])

    return weightmod

def modRatioNuEBar(particle_list, nu_nubar, nubar_e):
    modfactor = np.ones(len(particle_list['energy']))
    modfactor = nubar_e*ModNuEFlux(particle_list['energy'],
                           np.cos(particle_list['zenith']),
                           1.0, 1.0, 1.0, 1.0,
                           1.0, 1.0, 1.0, 1.0)

    weightmod = 1. + modfactor*nu_nubar
    antineutrinos = particle_list['ptype']<0
    weightmod[antineutrinos] = 1./(1+(1-nu_nubar)*modfactor[antineutrinos])

    return weightmod

# This I put for completeness. One could do without it.
# Don't use. Has no effect.
#def modRatioMuE(energy, elow, ehigh):
#    return 1.+LogLogParam(energy,
#                          elow*e1max_mu_e, ehigh*e2max_mu_e,
#                          x1e, x2e, 1000.)

def modRatioUpHor_NuMu(particle_list, numu_uphor):
    A_shape   = 1.*np.abs(numu_uphor)*LogLogParam(energy=particle_list['energy'], 
                               y1=z1max_mu,
                               y2=z2max_mu,
                               x1=x1z, x2=x2z,
                               cutoff_value = numu_cutoff)


    return 1-3.5*np.sign(numu_uphor)*norm_fcn(np.cos(particle_list['zenith']), A_shape, 0.35)

def modRatioUpHor_NuE(particle_list, nue_uphor):
    A_shape   = 1.*np.abs(nue_uphor)*LogLogParam(energy=particle_list['energy'], 
                               y1=(z1max_e+z1max_mu),
                               y2=(z2max_e+z2max_mu),
                               x1=x1z, x2=x2z,
                               cutoff_value = nue_cutoff)


    return 1-3.5*np.sign(nue_uphor)*norm_fcn(np.cos(particle_list['zenith']), A_shape, 0.35)


############
# Correction to the Honda flux (based on 2016 results, to be applied to pre-2016 fluxes)
# Taken from Honda's talk at ANW'16 (http://indico.universe-cluster.de/indico/contributionDisplay.py?contribId=0&confId=3533)
#############

# Correction function
coeff = [0.02812,-0.05848,1.045]
honda_fcn = np.poly1d(coeff)
ecritical = -coeff[1]/(coeff[0]*2)

# This factor is introduced so that the modification at low-energies is simply 1.
flux_renorm = honda_fcn(ecritical)-1.

# Formula that returns the correction factor
def honda_fix(energy):
    energy = np.log10(np.array(energy))
    ebool  = energy > ecritical
    output = np.ones_like(energy)
    output[ebool]  = honda_fcn(energy[ebool]) - flux_renorm
    return output
