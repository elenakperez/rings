import datetime

def centeredFFT(x, dt):
    """
    Computes FFT, with zero frequency in the center, and returns 
    dimensional frequency vector.
    X, freq = centeredFFT(x, dt)
    
    Parameters
    ----------
    x : numpy.ndarray of shape (N,) or (N,1)
        1D array to be transformed by FFT
    dt : numeric
        Time increment (used to compute dimensional freuency array)

    Returns (tuple)
    -------
    X: FFT of input x, with zero frequency in the center
    freq: Dimensional frequency vector corresponding to X

    #function [X,freq]=centeredFFT(x,dt)
    #
    # Adapted from a matlab function written by Quan Quach of blinkdagger.com 
    # Tom Farrar, 2016, 2020 jfarrar@whoi.edu
    # converted from matlab 2020
    # This code was written for MIT 12.805 class
    """
    import numpy as np
    from scipy import fft

    N = len(x)
    x = x.reshape(N,)
    
    #Generate frequency index
    if N % 2 == 0:
        m= np.arange(-N/2,N/2,1) # N even; this includes start (-N/2) and does not include stop (+N/2)
        #m = mslice[-N / 2:N / 2 - 1]    # N even (Matlab syntax)
    else:
        m= np.arange(-(N-1)/2,(N-1)/2+1,1) # N odd
        #m = mslice[-(N - 1) / 2:(N - 1) / 2]    # N odd (Matlab syntax)
    
    freq = m / (N * dt) #the dimensional frequency scale
    X = fft.fft(x)
    X = fft.fftshift(X) #swaps the halves of the FFT vector so that the zero frequency is in the center\\
    #If you are going to compute an IFFT, first use X=ifftshift(X) to undo the shift}
    return (X, freq) # Return tuple; could instead do this as dictionary or list

def band_avg(yy, num):
    '''
    Compute block averages for band averaging.

    Parameters
    ----------
    yy : np.array
        1D array to be averaged.
    num : numeric
        number of adjacent data points to average.

    Returns
    -------
    Bin-averaged version of input data, subsampled by the factor num

    # Tom Farrar, 2016, 2020 jfarrar@whoi.edu
    # This code was written for MIT 12.805 class

    '''
    #MATLAB code:
    #yyi=0;
    #for n=1:num
    # yyi=yy(n:num:[end-(num-n)])+yyi;
    #end
    import numpy as np
    
    yyi = 0
    for n in np.arange(0, num): # 1:num
        yyi=yy[n:-(num-n):num]+yyi;
        
    yy_avg=yyi/num
    
    return yy_avg


def confid(alpha,nu):
    """
    Computes the upper and lower 100(1-alpha)% confidence limits for 
    a chi-square variate (e.g., alpha=0.05 gives a 95% confidence interval).
    Check values (e.g., Jenkins and Watts, 1968, p. 81) are $\nu/\chi^2_{19;0.025}=0.58$
    and $\nu/\chi^2_{19;0.975}=2.11$ (I get 0.5783 and 2.1333 in MATLAB).
    
   
    Parameters
    ----------
    alpha : numeric
        Number of degrees of freedom
    nu : numeric
        Number of degrees of freedom

    Returns (tuple)
    -------
    lower: lower bound of confidence interval
    upper: upper bound of confidence interval

    # Tom Farrar, 2020, jfarrar@whoi.edu
    # converted from matlab 2020
    # This code was written for MIT 12.805 class
    """
    
    # requires:
    from scipy import stats
    
    upperv=stats.chi2.isf(1-alpha/2,nu)
    lowerv=stats.chi2.isf(alpha/2,nu)
    lower=nu / lowerv
    upper=nu / upperv
    
    return (lower, upper) # Return tuple; could instead do this as dictionary or list
    
def confidence_interval(alpha,nu,cstr,yspot=None,xspot=None,width=None,ax=None):
    """
    Plot (1-alpha)*100% spectral confidence interval on a log-log scale
    
    Parameters
    ----------
    alpha: numeric, between 0 and 1
        100*alpha is the percentage point of the chi-square distribution
        For example, use alpha=0.05 for a 95% confidence interval
    nu: numeric
        number of degrees of freedom
    cstr:
        color for confidence interval 
        For example, cstr = 'r' or cstr = h[-1].get_color()
    xspot: (1,) numpy.ndarray
        horizontal coordinate for confidence interval (e.g., xspot=freq(3);)
    yspot: (1,) numpy.ndarray
        vertical coordinate for confidence interval
    width: numeric
        width (in points) for top and bottom horizontal bars


    Returns
    -------
 
    # Tom Farrar, 2020, jfarrar@whoi.edu
    # converted from matlab 2020
    # This code was written for MIT 12.805 class
    """
    import numpy as np
    import matplotlib.pyplot as plt

    if ax is None:
      ax = plt.gca()
    if width is None:
      plt.sca(ax)
      fig = plt.gcf()
      # Get size of figure in pixels (would be preferable to use axis size instead)
      size = fig.get_size_inches()*fig.dpi # size in pixels
      width = np.round(0.01*size[0]) # make cap width 1.5% of figure width
      # Get width of axis in pixels
      # ax_width = np.diff(ax.get_xlim())*fig.dpi
      # width = np.round(0.015*ax_width) # make cap width 1.5% of figure width
    if yspot is None:
      plt.sca(ax)
      yax = ax.get_ylim() # Get xlim of current axis
      yax_midpoint = 0.75*np.diff(np.log10(yax)) # find point 75% from left in log space
      yspot = 10**(yax_midpoint + np.log10(yax[0])) # set default xspot to 75% from left
    if xspot is None:
      plt.sca(ax)
      xax = ax.get_xlim() # Get xlim of current axis
      xax_midpoint = 0.75*np.diff(np.log10(xax)) # find point 75% from left in log space
      xspot = 10**(xax_midpoint + np.log10(xax[0])) # set default xspot to 75% from left

    lower, upper = confid(alpha, nu)

    # Definition of lowz and upz differs from matlab version because 
    # plt.errorbar plots log10(yspot-lowz) and log10(yspot+upz), whereas, in 
    # matlab version I was plotting log10(yspot*lower) and log10(yspot*upper)
    lowz = yspot*(1-lower)
    upz = yspot*(upper-1)
    err = [lowz, upz]

    # plot confidence interval
    plt.errorbar(xspot, yspot, yerr=err, fmt='', capsize=width, ecolor=cstr)
    plt.text(xspot,yspot,'  ' + str(100*(1-alpha)) + '%',horizontalalignment='left');
    # plt.show()
    return(ax)


# From a Stackoverflow post by Rich Signell
# https://stackoverflow.com/questions/13965740/converting-matlabs-datenum-format-to-python
def matlab2datetime(matlab_datenum):
    day = datetime.datetime.fromordinal(int(matlab_datenum))
    dayfrac = datetime.timedelta(days=matlab_datenum%1) - datetime.timedelta(days = 366)
    return day + dayfrac
    
    # This doesn't work, but it would be nice
    # matlab2datetime_vec = np.vectorize(tt.matlab2datetime)
    # time = matlab2datetime_vec(mat['mday'])