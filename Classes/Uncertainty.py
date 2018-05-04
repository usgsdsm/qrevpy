"""
Created on Sep 26, 2017

@author: gpetrochenkov
<<<<<<< HEAD
"""
=======
'''
import numpy as np
from scipy.stats import t

class Uncertainty(object):
    
    def __init__(self, qa_data):
        self.__cov = None #Coefficient of variation for all used transect discharges
        self.__cov_95 = None #COV inflated by the 95# coverage factor
        self.__invalid_95  = None   # Estimated 95# uncertainty for dicharge in invalid bins and ensembles
        self.__edges_95 = None # Estimated 95# uncertainty for the computed edge discharges
        self.__extrapolation_95 = None # Estimated 95# uncertainty in discharge due to top and bottom extrapolations
        self.__moving_bed_95 = None  # Estimated 95# uncertainty due to moving-bed tests and conditions
        self.__systematic = None # Systematic error estimated at 1.5#
        self.__total_95 = None   # Estimated 95# uncertainty in discharge using automated values
        self.__cov_95_user = None  # User provided value for random uncertainty
        self.__invalid_95_user = None # User provided estimate of uncertainty for invalid data
        self.__edges_95_user  = None # User provided estimate of uncertainty for edges
        self.__extrapolation_95_user = None # User provided estimate of uncertainty for top and bottom extrapolation
        self.__moving_bed_95_user = None     # User provided estimate of uncertainty due to moving-bed conditions
        self.__systematic_user = None      # User provided estimate of systematic uncertainty
        self.__total_95_user = None         # Estimated 95# uncertainty in discharge using user provide values to override automated values
        self.__qa_data = qa_data
        
        
    def populate_data(self, meas):
        
        discharge = [x for x in meas.discharge]
        
        #Determine transects used to compute discharge
        checked = [x.checked == 1 for x in meas.transects]
        
        #Assign automatically generated uncertainties to properties
        self.ran_uncert_q(discharge[checked], 'total')
        self.uncert_invalid_data(discharge, checked)
        self.uncert_edge(discharge, checked)
        self.uncert_extrap(meas, discharge, checked)
        self.moving_bed_uncert(meas)
        self.__systematic = 1.5
        
        #Set user define uncertainties to empty
        self.__cov_95_user = [];
        self.__invalid_95_user = [];
        self.__edges_95_user = [];
        self.__extrapolation_95_user = [];
        self.__moving_bed_95_user = [];
        self.__systematic_user = [];
        self.__total_95_user = [];
        
        self.est_uncert()
        
    def est_uncert(self):
        '''Compute the uncertainty of the measurement using the automatically computed
        uncertainties and user overrides
        '''
        
        self.__total_95 = 2 * np.sqrt( (self.__cov95 / 2)**2 + (self.__invalid_95 /2)**2 \
                                       + (self.__edges_95 / 2)**2 + (self.__extrapolation_95 /2)**2 \
                                       + (self.__moving_bed_95 / 2) **2 + self.__systematic**2)
        
        
        #Override random uncertainty if user input provided
        if self.__cov_95_user is None:
            random = self.__cov_95
        else:
            random = self.__cov_95_user
            
            
        #Override uncertainty for invalid data if user input is provided
        if self.__invalid_95_user is None:
            invalid = self.__invalid_95
        else:
            invalid = self.__invalid_95_user
            
        #Override uncertainty edges if user input provided
        if self.__edges_95_user is None:
            edges = self.__edges_95
        else:
            edges = self.__edges_95_user
            
        #Override uncertainty for top and bottom extrapolation if user input provided
        if self.__extrapolation_95_user is None:
            extrapolation = self.__extrapolation_95
        else:
            extrapolation = self.__extrapolation_95_user
        
        #Override uncertainty for moving bed conditions if user input provided
        if self.__moving_bed_95_user is None:
            moving_bed = self.__moving_bed_95
        else:
            moving_bed = self.__moving_bed_95_user
            
        #Override systematic uncertainty if user input is provided
        if self.__systematic_user is None:
            system = self.__systematic
        else:
            system = self.__systematic_user
            
        #Compute total uncertainty using available user overrides
        self.__total_95_user = 2 * np.sqrt((random/2)**2 + (invalid/2)**2 + (edges/2)**2 \
                                           (extrapolation/2)**2 + (moving_bed/2)**2 + system**2)
        
    def add_user_input(self, prop, value):
        '''Sets specified property ro provided user input and recomputes the total uncertainty'''
        
        if np.isnan(value):
            setattr(self,prop,None)
        else:
            setattr(self,prop,value)
            
        self.est_uncert()
        
    def remove_user_input(self, prop):
        '''Removes user input from specified property'''
        
        setattr(prop, None)
        self.est_uncert()
        
                
    def ran_uncert_q(self, discharge, prop):
        '''Compute 95% random uncertainty for property of discharge object
        use simplified method for 2 allocate_transects
        Input:
        discharge: array of discharge measurements
        prop: name of property in object
        '''
        data = np.array([])
        #Determine number of objects instances
        n_max = len(discharge)
        if n_max > 0:
            #combine data into single array
            for n in range(0, n_max):
                val = getattr(discharge[n], prop)
                np.append(data, val)
                
                #Compute coefficient of variation
                avg = np.nanmean(data)
                std = np.nanstd(data)
                cov = np.abs(std/avg) * 100.
                
                #Inflate COV to 95% level and report as percent
                cov_95 = np.abs(t.ppf(.025, n_max -1)) * cov / np.sqrt(n_max)
                
                #Use approximate method as taught in class to reduce the high
                #coverage factor for 2 transects and account for prior knowledge related to_affine_scalar
                #720 second duration analysis
                if n_max == 2:
                    cov_95 = cov * 3.3
        else:
            cov = np.nan
            cov_95 = np.nan
            
        self.__cov = cov
        self.__cov_95 = cov_95
        
    def uncert_edge(self, discharge, checked):
        '''Compute uncertainty of edge discharge.  Currently assuming ransom plus bias in edge is within 
        30% of actual value
        
        Input:
        discharge: object of QComp
        checked: logical vecotr of transects to be used to compute total
        '''
        
        q_data = self.__qa_data
        
        #Compute total mean discharge
        mean_q = q_data.mean_q(discharge[checked], 'total')
        
        #Compute QDAta
        mean_left = q_data.mean_q(discharge[checked], 'left')
        per_left = (mean_left / mean_q) * 100
        
        #Compute percent of total discharge in right edge
        mean_right = q_data.mean_q(discharge[checked], 'right')
        per_right = (mean_right / mean_q) * 100
        #???Why is per left and per right still in the code???
        
        #Compute combined edge uncertainty
        per_edge = ((mean_left + mean_right) / mean_q) * 100
        self.__edges_95 = per_edge * 0.3
        
    def extrap_uncert(self, meas, discharge, checked):
        '''Compute uncertainty of top and bottom extrapolations
        
        Input:
        meas: object of Measurement
        discharge: object of QComp
        checked: logical vector of transects used to compute final
        '''
        
        #Compute mean total discharge
        q_select = np.nanmean([x.total_uncorrected for x in discharge[checked]])
                              
        #Create array of discharges from various extrapolation methods
        q_possible = [meas.extrap_fit.q_sensitivity.__q_pp_mean,
                      meas.extrap_fit.q_sensitivity.__q_pp_optmean,
                      meas.extrap_fit.q_sensitivity.__q_cns_mean,
                      meas.extrap_fit.q_sensitivity.__q_cns_optmean,
                      meas.extrap_fit.q_sensitivity.__q_3p_NS_mean,
                      meas.extrap_fit.q_sensitivity.__q_3p_NS_optmean]
        
        #Compute difference in discharges from selected method
        diff = np.abs(q_possible - q_select)
        
        #Sort differences
        diff = np.sort(diff)
        
        #Find the 4 smallest differences
        idx = np.where(diff > 0)[0][:4]
        
        #Estimate the uncertainty as the average of the wo smallest differences and report as a percent
        self.__extrapolation_95 = np.mean(diff[idx] / q_select) * 100
        
    def invalid_uncert(self, discharge, checked):
        '''Computes an estimate of the uncertainty for the discharge computed
        for invalid bins and ensembles
        
        Input:
        discharge: object of QComp
        checked: logical vector of transects used to compute final discharge
        '''
        
        q_data = self.__qa_data
        #compute mean total discharge
        mean_q = q_data.mean_q(discharge[checked], 'total')
        
        #Compute the percent discharge in invalid cells
        per_cells = (q_data.mean_q(discharge[checked], 'int_cells') / mean_q) * 100
        
        #Compute the percent discharge in invalid ensembles
        per_ens = (q_data.mean_q(discharge[checked], 'int_ens') / mean_q) * 100
        
        #Compute the uncertainty for combined invalid cells and ensembles
        self.__invalid_95 = (per_cells + per_ens) * 0.2
        
    def moving_bed_uncert(self, meas):
        '''Estimates the 95% uncertainty of the discharge due to the moving-bed tests,
        moving-bed conditions, and navigation reference'''
        
        #Compute quality code from moving-bed tests
        code = meas.qa.__moving_bed.code
        
        if code == 1:
            mb_uncert = 1
        elif code == 2:
            mb_uncert = 1.5
        elif code == 3:
            mb_uncert = 3
            
        #If bottom track is not used then the moving-bed uncertainty is zero
        if meas.transects[0].boat_vel.selected == 'bt_vel':
            mb_uncert = 0
            
        self.__moving_bed_95 = mb_uncert
        
                              
        
    
        
        
>>>>>>> 6ca6c50c231afa610ed3a693864074d7104a5f20
