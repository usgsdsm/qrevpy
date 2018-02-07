'''
Created on Jul 24, 2017

@author: gpetrochenkov
'''
"""
This module implements the Lowess function for nonparametric regression.
Functions:
lowess Fit a smooth nonparametric regression curve to a scatterplot.
For more information, see
William S. Cleveland: "Robust locally weighted regression and smoothing
scatterplots", Journal of the American Statistical Association, December 1979,
volume 74, number 368, pp. 829-836.
William S. Cleveland and Susan J. Devlin: "Locally weighted regression: An
approach to regression analysis by local fitting", Journal of the American
Statistical Association, September 1988, volume 83, number 403, pp. 596-610.
"""

# Authors: Alexandre Gramfort <alexandre.gramfort@telecom-paristech.fr>
#
# License: BSD (3-clause)

from math import ceil
import numpy as np
from scipy import linalg
from scipy.signal import butter, filtfilt




def lowess(x, y, f=2. / 3., _iter=3):
    """lowess(x, y, f=2./3., iter=3) -> yest
    Lowess smoother: Robust locally weighted regression.
    The lowess function fits a nonparametric regression curve to a scatterplot.
    The arrays x and y contain an equal number of elements; each pair
    (x[i], y[i]) defines a data point in the scatterplot. The function returns
    the estimated (smooth) values of y.
    The smoothing span is given by f. A larger value for f will result in a
    smoother curve. The number of robustifying iterations is given by iter. The
    function will run faster with a smaller number of iterations.
    """
    n = len(x)
    r = int(ceil(f * n))
    h = [np.sort(np.abs(x - x[i]))[r] for i in range(n)]
    w = np.clip(np.abs((x[:, None] - x[None, :]) / h), 0.0, 1.0)
    w = (1 - w ** 3) ** 3
    yest = np.zeros(n)
    delta = np.ones(n)
    for iteration in range(iter):
        for i in range(n):
            weights = delta * w[:, i]
            b = np.array([np.sum(weights * y), np.sum(weights * y * x)])
            A = np.array([[np.sum(weights), np.sum(weights * x)],
                          [np.sum(weights * x), np.sum(weights * x * x)]])
            beta = linalg.solve(A, b)
            yest[i] = beta[0] + beta[1] * x[i]

        residuals = y - yest
        s = np.median(np.abs(residuals))
        delta = np.clip(residuals / (6.0 * s), -1, 1)
        delta = (1 - delta ** 2) ** 2

    return yest


def lowpass_filt(data, fs, cutoff, order):
#     cutoff = .001
#     cutoff = .0001
    
    lowcut = cutoff / (.5 * fs)
        
    b, a = butter(order, [lowcut], btype='highpass')

    filtered_data = filtfilt(b, a, data)
   
    return filtered_data

if __name__ == '__main__':
    import math
    from scipy.io import loadmat
#     a = loadmat('02DB005_20170427_000_depthbeams.mat')
#     a = loadmat('20130705145405_depthbeams.mat')
    a = loadmat('20130705145744_depthbeams.mat')
    mat = a['dt'][0]
    x = np.arange(mat.shape[0])
    y = mat
    no_nans =  np.where(~np.isnan(mat))
    new_y = np.interp(x, x[no_nans], y[no_nans])
    y = new_y
#     n = 100
#     x = np.linspace(0, 2 * math.pi, n)
#     y = np.sin(x) + 0.3 * np.random.randn(n)
# 
    f = 1./20
    yest = lowess(x, new_y, f=f, iter=1)
    

    yest2 = lowpass_filt(y, 1, .25, 2)
  
    window_size = 20
    b = 0 - window_size/2
    e = 0 + window_size/2
    idx = np.repeat(0, len(yest2))
    w_sizes = []
    while e < len(yest2) + window_size/2:
        
        if b < 0:
            begin = 0
            end = b + window_size
            
        elif e > mat.shape[0] - 1:
            begin = e - window_size
            end = mat.shape[0] - 1
        
        else:
            begin = b
            end = e
             
        print (begin, end)
        
        w_sizes.append(end-begin)
        
        q75, q25 = np.percentile(yest2[begin:end], [75 ,25])
        iqr = q75 - q25
        
        
        idx[begin:end][yest2[begin:end] > q75 + (1.5 * iqr)] += 1
        idx[begin:end][yest2[begin:end] < q25 - (1.5 * iqr)] += 1
        
        b += 1
        e += 1
        
    yest2_final = yest2.copy()
    yest2_final[np.where(idx >= w_sizes)] = np.nan
    
    window_size = 40
    b = 0 - window_size/2
    e = 0 + window_size/2
    idx = np.repeat(0, len(yest))
    w_sizes = []
    while e < len(yest) + window_size/2:
        
        if b < 0:
            begin = 0
            end = b + window_size
            
        elif e > mat.shape[0] - 1:
            begin = e - window_size
            end = mat.shape[0] - 1
        
        else:
            begin = b
            end = e
             
        print (begin, end)
        
        w_sizes.append(end-begin)
        
        q75, q25 = np.percentile(yest[begin:end], [75 ,25])
        iqr = q75 - q25
        
        
        idx[begin:end][yest[begin:end] > q75 + (1.5 * iqr)] += 1
        idx[begin:end][yest[begin:end] < q25 - (1.5 * iqr)] += 1
        
        b += 1
        e += 1
        
    yest_final = yest.copy()
    yest_final[np.where(idx > 10)] = np.nan
    import pylab as pl
    pl.clf()
    ax = pl.subplot(4,1,1)
    ax.plot(x, y, label='Raw Signal')
    ax.plot(x, yest, label='Loess Smoothed')
    
    ax.legend()
     
    ax2 = pl.subplot(4,1,2)
#     ax2.plot(x, y, label='y noisy')
    ax2.plot(x, yest_final, label='Butterworth')
    ax2.scatter(x[np.isnan(yest_final)], yest[np.isnan(yest_final)], label='Loess Flag Error', color='red')
    ax2.legend()
    
    ax3 = pl.subplot(4,1,3)
#     ax2.plot(x, y, label='y noisy')
    ax3.plot(x, yest2, label='Butterworth')
    ax3.legend()
    
    ax4 = pl.subplot(4,1,4)
    ax4.plot(yest2_final)
    ax4.scatter(x[np.isnan(yest2_final)], yest2[np.isnan(yest2_final)], label='Buttworth Flag Error', color='red')
    ax4.legend()
   
    pl.show()