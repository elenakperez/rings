
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

 Utilities for 12.805 Data Analysis class
 
    ) import_PIES                 : import PIES data for final project 
    ) band_avg_autospec_dens      : compute band averaged autospectral density
    ) band_avg_autospec_dens_plot : plot band averaged autospectral density
    ) e_dag                       : compute E dagger matrix
    ) x_tilde                     : compute x tilde 
    ) C_nn                        : compute noise covariance matrix 
    ) C_xx                        : compute solution covariance matrix
    ) x_ML                        : compute maximum likelihood solution
    ) uncertainty_P               : compute uncertainty 
    ) eof_covar                   : EOF using covariance approach
    ) eof_svd                     : EOF using singular value decomposition approach
    ) compute_var                 : compute variance of a time series
    ) deg_of_freedom              : compute degrees of freedom
    ) 
    
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

# import Tom's functions
from utils.functions_12805 import * # script of functions Tom made for this pset

# import necessary packages and functions
import numpy as np
import xarray as xr
import pandas as pd
import netCDF4 as nc
import scipy.io
from scipy.interpolate import griddata
import pickle
import math
import matplotlib
import statistics
import matplotlib.pyplot as plt
from scipy import signal

# turn off warnings
import warnings
warnings.filterwarnings("ignore")


#-------------------------------------------------------------------------------------------------------------------------------------------
# ) this function import data from Magdalena's PIES and returns the variables as Numpy arrays

def import_PIES(filepath, site):
    """"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    Input:
        filepath (string)    : path to the data (e.g. '/Users/elenaperez/Desktop/chatts/data/classes/PIES_data/L2/P1/P1.mat')
        site (string)        : which site to import data for (e.g. 'P1')
    
    Output:
     
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    
    'READ IN FILE'
    P = scipy.io.loadmat(filepath)

    'CONVERT DATA TO COLUMN VECTORS'
    P_amp = P[site]['amplitude'][0][0][0][0]
    P_comm = P[site]['comments'][0][0][0]
    P_drift = P[site]['drift'][0][0][0]
    P_drift_coeff = P[site]['driftcoef'][0][0][0]
    P_paros = P[site]['paros'][0][0][0]
    P_phase = P[site]['phase'][0][0][0]
    P_prs = P[site]['prs'][0][0][0]
    P_prsave = P[site]['prsave'][0][0][0]
    P_prsdd = P[site]['prsdd'][0][0][0]
    P_refyr = P[site]['refyr'][0][0][0]
    P_sn = P[site]['sn'][0][0][0]
    P_tau = P[site]['tau'][0][0][0] 
    P_tau_tide = P[site]['tau_tide'][0][0][0]
    P_taudd = P[site]['taudd'][0][0][0]
    P_tide = P[site]['tide'][0][0][0]
    P_tmp = P[site]['tmp'][0][0][0]
    P_tmpdd = P[site]['tmpdd'][0][0][0]
    
    return P_amp, P_comm, P_drift, P_drift_coeff, P_paros, P_phase, P_prs, P_prsave, P_prsdd, P_refyr, P_sn, P_tau, P_tau_tide, P_taudd, P_tide, P_tmp, P_tmpdd
    

#-------------------------------------------------------------------------------------------------------------------------------------------

def band_avg_autospec_dens(x, dt, M, taper):
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    This function produces a one-sided, band-averaged estimate of the autospectral density of a given 
    time series x, with time increment dt, and number of bands M. 
    
    Input:
        x (array)                   : vector to be transformed (real & uniformly spaced)
        dt (float)                  : time increment
        M (float)                   : number of bands to average 
        taper (bool)                : if True, taper the time series before the FFT

    Output: 
        spec_dens_avg (array)       : band-averaged autospectral density of the time series x
        fm_avg (array)              : band-averaged frequency of the time series x
        lower_ci (float)            : lower confidence limit for a 95% confidence interval
        upper_ci (float)            : upper confidence limit for a 95% confidence interval
    
     
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    
    'DEFINE VARIABLES'
    N = x.size # number of observations
    T = N * dt # total time of observations 
    
    '1. COMPUTE THE MEAN OF THE TIME SERIES x(t) AND SUBTRACT IT FROM THE TIME SERIES'
    x_anom = x - x.mean()
    
    if taper: 
        '2. TAPER THE RECORD USING A TAPER WINDOW xtaper(t) = w(t)x(t)'
        w = scipy.signal.windows.tukey(len(x_anom)) # Tukey window
        xtaper = w*x
    
    '3. COMPUTE THE FFT OF THE ENTIRE TAPERED RECORD'
    x_fft, fm = centeredFFT(x_anom, dt)
    
    '4. GENERATE A FREQUENCY VECTOR fm = m/(NΔt) WITH m RANGING FROM −N/2 TO N/2 − 1 (FOR EVEN N)'
    # generated from the centeredFFT function, fm
    
    '5. COMPUTE THE RAW, ONE-SIDED SPECTRAL DENSITY AND DISCARDED FREQUENCIES ≤ 0'
    # equation 4.27 from the notes, psi_k = autospectral density
    spec_dens = ((2*T)/(N**2))*(x_fft*(np.conj(x_fft)))
    
    # only keep positive frequencies, negative freq don't have physical meaning
    spec_dens = spec_dens[np.where(fm > 0)]
    fm = fm[np.where(fm > 0)]
        
    '6. BAND AVERAGE THE RAW SPECTRUM OVER nd FREQUENCY BANDS'
    spec_dens_avg = band_avg(spec_dens, M)
    fm_avg = band_avg(fm, M)
        
    '7. ESTIMATE A CONFIDENCE INTERVAL (TYPICALLY A 95% CONFIDENCE INTERVAL)'
    lower_ci, upper_ci = confid(0.05,2*M)
    
    return spec_dens_avg, fm_avg, lower_ci, upper_ci

#-------------------------------------------------------------------------------------------------------------------------------------------

def band_avg_autospec_dens_plot(fm_avg, spec_dens_avg, M, title, fig_quality, lineColor, confidColor, ifColor):
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    This function produces a plot of the autospectral density of a time sereis spec_dens_avg plotted against
    the frequency averaged fm_avg with a 95% confidence interval.
    
    Input:
        spec_dens_avg (array)       : band-averaged autospectral density of the time series x
        fm_avg (array)              : band-averaged frequency of the time series x
        title (string)              : title of the figure 
        fig_quality (int)           : quality of the figure, e.g. 100 dpi
        lineColor (string)          : color of the line plot
        confidColor (string)        : color of the confidence interval
        ifColor (bool)              : if True use custom colors, else False use default color scheme

    Output: 
    * returns plot of the band-averaged autospectral density of a time-series plotted against the frequencies
    
     
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    
    # create figure
    fig,ax = plt.subplots(figsize=(10,7))
    fig.set_dpi(fig_quality)
    fig.suptitle(title, fontsize=14, y=0.925, fontweight='bold')

    # plot
    if ifColor:
        ax.loglog(fm_avg, spec_dens_avg, color=lineColor, alpha=0.8);
        confidence_interval(0.05,M*2,confidColor)
    else:
        ax.loglog(fm_avg, spec_dens_avg, color='#2c7fb8', alpha=0.8);
        confidence_interval(0.05,M*2,'#253494')

    # axes formatting
    ax.set_xlabel('Frequency [1/s]',fontweight='bold');
    ax.set_ylabel('Autospectral Density [(m/s)^/(1/s)]',fontweight='bold');
    
    return fig, ax


#-------------------------------------------------------------------------------------------------------------------------------------------

"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Input:
    E (array)       : numpy array of the E matrix

Output: 
    E_dag (array)   : numpy array of the E dagger matrix


"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

def e_dag(E):
    return np.matmul(np.linalg.inv((np.matmul(E.T,E))),E.T)

#-------------------------------------------------------------------------------------------------------------------------------------------


"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Input:
    E_dag (array)   : numpy array of the E dagger matrix
    y (array)       : numpy array of y

Output: 
    x_tilde (array) : numpy array of x tilde


"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

def x_tilde(E_dag, y):
    return E_dag@y

#-------------------------------------------------------------------------------------------------------------------------------------------


"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Input:
    E_dag (array)   : numpy array of the E dagger matrix
    C_nn (array)    : numpy array of the noise covariance matrix C_nn

Output: 
    C_xx (array)    : numpy array solution covariance matrix C_xx


"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

def C_xx(E_dag, C_nn):
    return E_dag@C_nn@E_dag.T

#-------------------------------------------------------------------------------------------------------------------------------------------

"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Input:
    E (array)       : numpy array of the E dagger matrix
    C_nn (array)    : numpy array of the noise covariance matrix C_nn
    y (array)       : numpy array of y, the knowns

Output: 
    x_ML (array)    : numpy array of maximum likelihood estimation


"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

def x_ML(E, C_nn, y):
    return (np.linalg.inv(E.T@np.linalg.inv(C_nn)@E))@E.T@np.linalg.inv(C_nn)@y


#-------------------------------------------------------------------------------------------------------------------------------------------

"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Input:
    E (array)       : numpy array of the E dagger matrix
    C_nn (array)    : numpy array of the noise covariance matrix C_nn

Output: 
    x_ML (array)    : numpy array of uncertainity P


"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

def uncertainty_P(E, C_nn):
    return np.sqrt(np.linalg.inv(E.T@np.linalg.inv(C_nn)@(E)))

#-------------------------------------------------------------------------------------------------------------------------------------------




"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Input:
    y (array)           : numpy array y, the knowns

Output: 
    cov_eigval (array)  : eigenvalues for covariance approach to EOF analysis
    cov_eigvec  (array) : eigenvectors for covariance approach to EOF analysis


"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

def eof_covar(y):
    
    y = y - y.mean(axis = 1, keepdims = True) # we must remove the mean so that <y>=0, but do it for each column (axis=1)
    C_yy = (y@y.T)/len(y[1]) # calculate the covariance matrix <y*y.T>
    cov_eigval, cov_eigvec = np.linalg.eig(C_yy) # we want to find the eigenvectors of C_yy, which are the EOFs of U
    
    return cov_eigval, cov_eigvec
    

#-------------------------------------------------------------------------------------------------------------------------------------------


"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Input:
    y (array)     : numpy array y, the knowns

Output: 
    U (array)     : numpy array of U matrix for SVD
    Lamda (array) : numpy array of Lambdas for SVD
    Vh (array)    : numpy array of V matrix for SVD
    b (array)     : expansion coefficients for the EOF modes of the SVD approach
    


"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

def eof_svd(y):
    (U,lambdas,Vh)  = np.linalg.svd(y) # use numpy linalg to take SVD 
    Lambda = np.diag(lambdas)
    
    b = U.T@y
    
    return U, Lambda, Vh, b


#-------------------------------------------------------------------------------------------------------------------------------------------

    
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

I believe I found this function online and modified it for my homework. I can't remeber where I grabbed it from but I want to acknowledge that this is not completely my own work. Although it is useful, hence why it's included in my software toolbox.

Input:
    data (array)     : numpy array of the data

Output: 
    variance (array) : numpy array of the variance of the data
    


"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
def compute_var(data):
    
    n = len(data)
    mean = sum(data) / n     # Mean of the data
    deviations = [(x - mean) ** 2 for x in data]     # Square deviations
    variance = sum(deviations) / n # Variance
    
    return variance

#-------------------------------------------------------------------------------------------------------------------------------------------

    
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Input:
    M (int)          : number of observations
    N (int)          : number of constraints
Output: 
    degOfF (int)     : degrees of freedom for a problem
    


"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

def deg_of_freedom(M, N):
    return M-N # degrees of fredom = number of obs - number of constraints


#-------------------------------------------------------------------------------------------------------------------------------------------

