'''
Created on Sep 26, 2017

@author: gpetrochenkov
'''
import numpy as np
import statsmodels.api as sm
from Classes.FitData import FitData

class SelectFit(object):
    '''Class definition for data class that contains all of the automated
    extrapolation method selection information. This inherits all the
    properties from FitData and addes the automated selection data to it.
    David S. Mueller 6/17/2011
    
    
    Last modificaitons / validation 5/15/2012
    7/24/2013
    DSM modified to catch C shaped profile lines 63-66, 123
    DSM modified to fix manual fits and how it displays line 172
    7/25/2013
    DSM fixed bug in line 172, changed 0.1667 to obj.exponent
    
    4/4/2014
    DSM modified for use in QRev'''
    
    def __init__(self):
        self.__fit_method = None #User selected method Automatic or Manual
        self.__top_method = None
        self.__bot_method = None
        self.__exponent = None
        self.__exp_method = None
        self.__u = None
        self.__u_auto = None
        self.__z = None
        self.__z_auto = None
        self.__residuals = None
        self.__coef = None
        self.__bot_method_auto = None #Selected extrapolation for top
        self.__top_method_auto = None #Selected extrapolation for bottom
        self.__exponent_auto = None #Selected exponent
        self.__top_fit_r2 = None #Top fit custom r^2
        self.__top_max_diff = None #Maximum difference between power and 3-pt at top
        self.__bot_diff = None #Difference between power and no slop at z = 0.1
        self.__bot_r2 = None #Bottom fit r^2
        self.__fit_r2 = None #Selected fit of selected power/no slip fit
        self.__ns_exponent = None # No slip optimized exponent
        self.__pp_exponent = None #Power Power optimized exponent
        self.__top_r2 = None
        self.__rsqr = None
        self.__exponent_95_ci = None
        
    def populate_data(self, normalized, fit_method, kargs = None):
        
        valid_data = normalized._NormData__valid_data
        
        #Store data in properties to object
        self.__fit_method = fit_method
        
        #Compute power fit with optimized exponent as reference to determine
        #if constant no slip will be more appropriate
        ppobj = FitData()
        ppobj.populate_data(normalized, 'Power', 'Power', 'optimize')
        
        #Store results in obhect
        self.__pp_exponent = ppobj._FitData__exponent
        self.__residuals = ppobj._FitData__residuals
        self.__rsqr = ppobj.__rsqr
        self.__exponent_95_ci = ppobj.__exponent_95_ci
        
        #Begin automatic fit
        
        #More than 6 cells are required to compute an optimized fit.  For fewer
        #than 7 cells the default power/power fit is selected due to lack of sufficient
        #data for a good analysis
        if len(self.__residuals) > 6:
            #Compute the difference between the top two cells of data and the optimized power fit
            top2 = np.sum(normalized._NormData__unit_normalized_med[valid_data[-2:]] \
                          - ppobj._FitData__coef * normalized._NormData__unit_normalized_z[valid_data[-2:]]) \
                          ** ppobj._FitData__exponent
                          
            #Compute the difference between the bottom two cells of data and the optimized power fit
            bot2 = np.sum(normalized._NormData__unit_normalized_med[valid_data[:2]] \
                          - ppobj._FitData__coef * normalized._NormData__unit_normalized_z[valid_data[:2]]) \
                          ** ppobj._FitData__exponent
                          
            #Compute the difference between the middle two cells of data and the optimized power fit
            mid1 = np.floor(len(np.isnan(valid_data) == False) / 2)
            mid2 = np.sum(normalized._NormData__unit_normalized_med[valid_data[mid1:mid1+2]] \
                          - ppobj.coef * normalized._NormData__unit_normalized_z[valid_data[mid1:mid1+2]]) \
                          ** ppobj._FitData__coef
                          
            self.__top_method_auto = 'Power'
            self.__bot_method_auto = 'Power'
            
            #Evaluate difference in data and power fit at water surface using a linear fit throught the top 4
            #median cells and save results
            x = normalized._NormData__unit_normalized_med[valid_data[:4]]
            x = sm.add_constant(x)
            y = normalized._NormData__unit_normalized_z[valid_data[:4]]
            lin_fit = sm.OLS(x,y)
            result = lin_fit.fit()
            dsmfitr2 = 1 - (np.sum(result.resid ** 2) / np.mean(np.abs(result.resid)))
            self.__top_r2 = dsmfitr2
            self.__top_r2 = result.rsquared
            
            #Evaluate overall fit
            #If the optimized power fit does not have an r^2 better than 0.8 or if the optimized
            #exponent if 0.1667 falls within the 95% confidence interval of the optimized fit,
            #there is insufficient justification to change the exponent from 0.1667
            if ppobj.__rsqr < 0.8 or (0.1667 > self.__exponent_95_ci[0] and 0.1667 < self.__exponent_95_ci):
                #If an optimized exponent cannot be justified the linear fit is used to determine if a constant
                #fit at the top is a better alternative than a power fit.  If the power fit is the better
                #alternative the exponent is set to the default 0.1667 and the data is refit
                if np.abs(self.__top_fit_r2 < 0.8 or self.__top_r2 < 0.9):
                    fd = FitData()
                    ppobj = fd(normalized, 'Power', 'Power', 'Manual', 0.1667)
                    
            #Evaluate fit of top and bottom portions of the profile
            #Set save selected exponent and associated fit statistics
            self.__exponent_auto = ppobj._FitData__exponent
            self.__fit_r2 = ppobj._FitData__r_squared
            
            #Compute the difference at the water surface between a linear fit of the top 4 measured cells
            #and the best selected power fit of the whole profile
            self.__top_max_diff = ppobj._FitData__u[-1] - np.sum(result.params)
            
            #Evaluate the difference at the bottom between power using the whole profile and power using
            #only the bottom third
            ns_fd = FitData()
            ns_fd.populate_data(normalized, 'Constant', 'No Slip', 'Optimize')
            self.__ns_exponent = ns_fd._FitData__exponent
            self.__bot_r2 = ns_fd._FitData__rsquared
            self.__bot_diff = ppobj._FitData__u[np.round(ppobj._FitData__z,2) == 0.1] \
            - ns_fd._FitData__u[np.round(ns_fd._FitData__z, 2) == 0.1]
            
            #Begin automatic selection logic
            #-----------------------------------
            
            # A constant no slip fit condition is selected if:
            # 
            # 1)The top of the power fit doesn't fit the data well.
            # This is determined to be the situation when
            # (a) the difference at the water surface between the
            # linear fit and the power fit is greater than 10% and
            # (b) the difference is either positive or the difference
            # of the top measured cell differs from the best
            # selected power fit by more than 5%.
            # OR
            # 2) The bottom of the power fit doesn't fit the data 
            # well. This is determined to be the situation when (a)  
            # the difference between and optimized no slip fit
            # and the selected best power fit of the whole profile
            # is greater than 10% and (b) the optimized on slip fit has
            # and r^2 greater than 0.6.
            # OR
            # 3) Flow is bidirectional. The sign of the top of the
            # profile is different from the sign of the bottom of
            # the profile.
            # OR
            # 4) The profile is C-shaped. This is determined by
            # (a) the sign of the top and bottom difference from
            # the best selected power fit being different than the
            # sign of the middle difference from the best selected
            # power fit and (b) the combined difference of the top
            # and bottom difference from the best selected power
            # fit being greater than 10%.
            
            if (np.abs(self.__top_max_diff > 0.1) \
                and (self.__top_max_diff > 0 or np.abs(normalized._NormData__unit_normalized_med[valid_data[0]] - ppobj._FitData__u[-1]) > 0.05) \
                or (( np.abs(self.__bot_diff) > 0.1) and self.__bot_r2 > 0.6 ) \
                or (np.sign(normalized._NormData__unit_normalized_med)[valid_data[0]] != np.sign(normalized._NormData__unit_normalized_med[valid_data[-1]]))) \
                or np.sign(bot2) * np.sign(top2) == np.sign(mid2) and np.abs(bot2+top2) > 0.1:
                
                #Set the bottom to no slip
                self.__bot_method_auto = 'No Slip'
                #If the no slip fit with an optimized exponent does not have r^2 better than 0.8 use the default 0.1667 for the no slip exponent
                if ns_fd.__r_squared > 0.8:
                    self.__exponent_auto = ns_fd._FitData__exponent
                    self.__fit_r2 = ns_fd._FitData__rsquared
                else:
                    self.__exponent_auto = 0.1667
                    self.__fit_r2 = np.nan
                    
                #Use the no slip 95% confidence intervals if they are available
                if ns_fd.__exponent_95_ci is not None and np.isnan(ns_fd.exponent95confint) == False:
                    self.__exponent_95_ci[0] = ns_fd._FitData__exponent_95_ci[0]
                    self.__exponent_95_ci[1] = ns_fd._FitData__exponent_95_ci[1]
                else:
                    self.__exponent_95_ci[0] = np.nan
                    self.__exponent_95_ci[1] = np.nan
                    
                #Set the top method to constant
                self.__top_method_auto = 'Constant'
                
            else:
                
                #Leave the fit to power/power and set the best selected optimized exponent as the automatic fit exponent
                self.__exponent_auto = ppobj._FitData__exponent
                
        else:
            
            #if the data are insufficient for a valid analysis use the power/power fit with the default 0.1667 exponent
            self.__top_method_auto = 'Power'
            self.__bot_method_auto = 'Power'
            self.__exponent_auto = 0.1667
            self.__ns_exponent = 0.1667
            
        #Update the fit uysing the automatically selected methods
        update_fd = FitData()
        update_fd.populate_data(normalized, self.__top_method_auto, self.__bot_method_auto, 'Manual', self.__exponent_auto)
        update_auto = update_fd
        
        if fit_method == 'Manual':
            
            if len(kargs) == 1:
                trans_data = kargs[0]
                update_fd = FitData()
                update_fd.populate_data(normalized, trans_data.extrap._ExtrapData__top_method, trans_data.extrap._ExtrapData__bot_method, 'Manual', trans_data.extrap._ExtrapData__exponent)
            else:
                update = FitData()
                update_fd.populate_data(normalized, kargs[1], kargs[2], 'Manual', kargs[3])
                
        #Store fit data in object
        self.__top_method = update._FitData__top_method
        self.__bot_method = update._FitData__bot_method
        self.__exponent = update._FitData__exponent
        self.__coef = update._FitData__coef
        self.__u = update._FitData__u
        self.__u_auto = update_auto._FitData__u_auto
        self.__z_auto = update_auto._FitData__z_auto
        self.__z = update._FitData__z
        self.__exp_method = update._FitData__exp_method
        self.__residuals = update._FitData__residuals
                
              
            
            
        
