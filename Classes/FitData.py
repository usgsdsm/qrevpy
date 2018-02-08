'''
Created on Sep 26, 2017

@author: gpetrochenkov
'''
import numpy as np
from numpy import poly1d
from scipy.optimize.minpack import curve_fit
from scipy.stats import t


class FitData(object):
    '''  Class adapted from the extrap program.
    
     Class definition for profile extrapolation fit properties. This class
     consists of a constructor method and a method to plot the fit as a
     solid line. 
     Data required for the constructor method include data of class
     NormData, threshold for the minimum number of points for a valid
     median, top extrapolation method, bottom extrapolation method, type
     of fit, and if a manual fit, the exponent.
     David S. Mueller, 2/18/2011
    
     Modified 6/17/2011, dsm
     1) Added fit statistics
    
     Modified 10/17/2011, dsm
     2) Moved application of threshold criteria to NormData as property
     validData.
    
     Last modificaitons / validation 5/15/2012'''
    
    def __init__(self):
        self.__file_name = None #Name of transect file
        self.__top_method = None #Top extrapolation method
        self.__bot_method = None #Bottom extrapolation method
        self.__coef = None #Power fit coefficient
        self.__exponent = None #Power fit exponent
        self.__u = None #Fit values of the variable
        self.__u_auto = None #fit values from automatic fit
        self.__z_auto = None #z values for automtic fit
        self.__z = None #Distance from the streambed for fit variable
        self.__exp_method = None #Method to determine exponent (default, optimize, or manual)
        self.__data_type = None #Type of data (velocity or unit discharge)
        self.__exponent_95_ci = None #95% confidence intervals for optimized exponent
        self.__residuals = None #Residuals from fit
        self.__rsqr = None #adjusted r^2 for optimized exponent
        self.__fit_func = None #method to use in curve fit
        self.__bounds = None #Bounds for curve fitting coefficients (None if not necessary)
        self.__p0 = None #Initial guess in curve fit (None if not necessary)
        self.__r_squared = None #R squared of model
        
    def populate_data(self, norm_data, top, bot, method, kargs = None):
        
        #If no arguments just create object
        if kargs is None:
            unit_norm_no = norm_data._NormData__unit_normalized_no
            avg_z = norm_data._NormData__unit_normalized_z
            y = norm_data._NormData__unit_normalized_med
            idxz = np.squeeze(norm_data._NormData__valid_data)
            zc = np.nan
            
            lower_bound = [-np.inf, 0.01]
            upper_bound = [np.inf, 1]
            
            #Process data if available
            if idxz is not None:
                idx_power = idxz
                
                #Create arrays for data fitting
                #Select median values to use in extrapolation methods selected and create
                #methods selected and create fir output data arrays
                
                #If bottom is No Slip, Power at top is not allowed
                if bot == 'No Slip':
                    if top == 'Power':
                        top = 'Constant'
                        
                fit_combo = ''.join([top, bot])
                if fit_combo == 'PowerPower':
                    self.__z = np.arange(0,1.01,.01)
                    self.__z = self.__z.T
                    zc = np.nan
                    uc = np.nan
                elif fit_combo == 'ConstantPower':
                    self.__z = np.arange(0, np.max(avg_z[idxz])+0.01, 0.01)
                    self.__z = np.vstack([self.__z, np.nan])
                    zc = np.arange(np.max(avg_z[idxz] + 0.01), 1.01, 0.01)
                    zc = zc.T
                    uc = np.tile(y[idxz[0]], zc.shape)
                elif fit_combo == '3-PointPower':
                    self.__z = np.arange(0, np.max(avg_z[idxz]) + 0.01, 0.01)
                    self.__z = np.vstack([self.__z, np.nan])
                    #If less than 6 bins use contatnt at the top
                    if len(idxz) < 6:
                        zc = np.arange(np.max(idxz) + 0.01, 1.01, 0.01)
                        zc = zc.T
                        uc = np.tile(y[idxz[0]], zc.shape)
                    else:
                        p = poly1d(avg_z[0:3], y[0:3])
                        zc = np.max(avg_z[idxz] + 0.01, 1.01, 0.01)
                        zc = zc.T
                        uc = zc * p[0] + p[1]
                        
                elif fit_combo == 'ConstantNo Slip':
                    #Optimize constant / no slip if sufficient cells are available
                    if method == 'optimize':
                        idx = idxz[1+len(idxz)- np.floor(len(avg_z[idxz]) / 3):-1];
                        if len(idx) < 4:
                            method = 'default'
                            
                    #Compute Constant / No Sli using WinRiver II and
                    #RiverSurveyor Live defaault cells
                    else:
                        idx = np.where(avg_z[idxz] <= .2)
                        if len(idx) < 1:
                            idx = idxz[-1]
                        else:
                            idx = idxz[idx[0]]
                            
                    #Configures u and z arrays
                    idxns = idx
                    self.__z = np.arange(0, avg_z[idxns[0]] + 0.01, 0.01)
                    self.__z = np.hstack([self.__z, [np.nan]])
                    idx_power = idx
                    
                    #If less than 6 bins use constatnt at the top
                    if len(idxz) < 6:
                        zc = np.arange(np.max(avg_z) +0.01, 1.01, 0.01) 
                        zc = zc.T
                        uc = np.tile(y[idxz[0]], zc.shape)
                    else:
                        p = np.polyfit(avg_z[idxz[0:3]], y[idxz[0:3]], deg=1)
                        
                        zc = np.arange(np.max(avg_z[idxz]), 1.01, .01)
                        zc = zc.T
                        uc = zc * p[0] + p[1]
                        
                
                #Compute exponent
                zfit = avg_z[idx_power]
                yfit = y[idx_power]
                
                #Check data validity
                ok1 = [np.isfinite(z) for z in zfit]
                ok2 = [np.isfinite(y) for y in yfit]
                ok_ =  np.array([z == 1 and y == 1 for z,y in zip(ok1,ok2)])
                if np.all(ok_) == False:
                    pass
                    #Add warning
                    
                self.__exponent = np.nan
                self.__exponent_95_ci = np.nan
                self.__rsqr = np.nan
                
                lower_method = method.lower()
                
                if lower_method == 'manual':
                    self.__fit_func = 'linear'
                    self.__exponent = kargs[0]
                    self.__bounds = None
                    self.__p0 = None
                    
                elif lower_method == 'default':
                    self.__fit_func = 'linear'
                    self.__exponent = 1./6.
                    self.__bounds = None
                    self.__p0 = None
                    
                elif lower_method == 'optimize':
                    self.__fit_func = 'power'
                    self.__bounds = [lower_bound, upper_bound]
                    strt = yfit[ok_]
                    self.__p0 = [strt[-1], 1./6]
                    
                fit_funcs = {
                    'linear': lambda x, a: a * x**(self.exponent),
                    'power': lambda x, a, b: a * x**b
                }
                
                if len(ok_) > 1:
                    popt, pcov = curve_fit(fit_funcs[self.__fit_func], 
                                       zfit, yfit, p0 = self.__p0, bounds = self.__bounds)
                    
                    #Extract exponent and confidence intervals from fit
                    if lower_method == 'optimize':
                        self.exponent = popt[1]
                        if self.__exponent < 0.05:
                            self.__exponent = 0.05
                            
                        if len(zfit[ok_]) > 2:
                            n = len(y)    # number of data points
                            p = len(popt) # number of parameters
                            
                            t_val = t.cdf(.025, n-1)

                            #get 95% confidence intervals
                            upper, lower = [], []
                            for j in range(len(popt)):
                                if self.__bounds[0][j] == -np.inf and self.__bounds[1][j] == np.inf:
                                    lower.append(popt[j] - t_val * np.sqrt(np.diag(pcov)[j])) 
                                    upper.append(popt[j] + t_val * np.sqrt(np.diag(pcov)[j]))
                                else:
                                    lower.append(np.nan)
                                    upper.append(np.nan)
                             
                            #Stack the confidence intervals
                            self.__exponent_95_ci  = np.vstack([lower, upper])       
                            if method == 'optimize':
                                self.__exponent_95_ci  = self.__exponent_95_ci[:,1]
                                
                            #Get the rsquared for the model
                            ss_tot = np.sum((y - np.mean(yfit))**2)
                            ss_res = np.sum((y[idx_power] - fit_funcs[self.__fit_func](zfit, *popt))**2)
                            self.__r_squared = 1 - (ss_res/ss_tot) 
                        else:     
                            self.__exponent_95_ci = np.nan
                            self.__r_squared = np.nan
                
                #Fit power curve to appropriate data
                self.__coef = ((self.__exponent + 1) * .05 * np.nansum(y[idx_power])) / np.nansum(((avg_z[idx_power] + .5 * .05)**(self.__exponent - 1) - ((avg_z[idx_power] - .5 * .05)**(self.__exponent+1))))
                   
                #Compute residuals
                self.__residuals = y[idx_power] - self.__coef * avg_z[idx_power]**self.__exponent
                
                #Compute values (velocity or discharge) based on exponent and compute coefficient
                self.__u = self.__coef * self.__z**self.__exponent
                if type(zc) is float and zc == np.nan:
                    self.__u = np.hstack([self.__u, [uc]])
                    self.__z = np.hstack([self.__z, [zc]])
                    
                #Assign variables to object properties
                self.__file_name = norm_data._NormData__file_name
                self.__top_method = top
                self.__bot_method = bot
                self.__exp_method = method
                self.__data_type = norm_data._NormData__data_type
                
            else:
                #If not data are valid simply apply methods
                self.__exponent = np.nan
                self.__exponent_95_ci = [np.nan, np.nan]
                self.__r_squared = np.nan
                self.__file_name = norm_data._NormData__file_name
                self.__top_method = top
                self.__bot_method = bot
                self.__exp_method = method
                self.__data_type = norm_data._NormData.__data_type
                
                    
                    
                    
                                   
        
        