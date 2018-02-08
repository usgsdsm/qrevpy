'''
Created on Sep 26, 2017

@author: gpetrochenkov
'''
import numpy as np
from Classes.SelectFit import SelectFit
from Classes.ExtrapQSensitivity import ExtrapQSensitivity
from Classes.NormData import NormData

class ComputeExtrap(object):
    '''Class to compute the optimized or manually specified extrapolation methods'''
    
    def __init__(self):
        self.threshold = None #threshold as a percent for determining if a median is valid 
        self.subsection = None #percent of discharge, does not account for transect direction
        self.fit_method = None #method used to determine fit.  Automatic or manual
        self.norm_data = [] #object of class norm data
        self.sel_fit = [] #Object of class SelectFit
        self.q_sensitivity = None #Object of class ExtrapQSensitivity
        self.messages = [] #Variable for messages to UserWarning
        
    def populate_data(self, trans_data, kargs = None):
        
        self.threshold = 20
        self.subsection = [0, 100]
        self.fit_method = 'Automatic'
        self.process_profiles(trans_data, 'q')
        #Compute the sensitivity of the final discharge to changes in extrapolation methods
        if kargs is None:
            self.q_sensitivity = ExtrapQSensitivity(trans_data,self.sel_fit)
            
    def process_profiles(self, trans_data, data_type):
        '''Function that serves and the main control for other classes and functions'''
        
        for n in trans_data:
        #Compute normalized data
            nd = NormData()
            nd.populate_data(n, data_type, self.threshold, self.subsection)
            self.norm_data.append(nd)
            
        #Get composite norm data
        comp_nd = NormData()
        comp_nd.get_composite_data(trans_data, self.norm_data)
        self.norm_data.append(comp_nd)
        
        for n in range(len(trans_data)):
        #Compute the fit for the selected  method
            if self.fit_method == 'Manual':
                sf = SelectFit()
                sf.populate_data(self.norm_data[n], self.fit_method, trans_data[n])
                self.sel_fit.append(sf)
            else:
                sf = SelectFit()
                sf.populate_data(self.norm_data[n], self.fit_method)
                self.sel_fit.append(sf)
                
            if sf._SelectFit__top_fit_r2 is not None:
                #Evaluate if there is a potential that a 3-point top method may be appropriate
                if sf.__top_fit_r2 > 0.9 or sf.__top_r2 > 0.9 and np.abs(sf.__top_max_diff) > 0.2:
                    self.messages.append('The measurement profile may warrant a 3-point fit at the top')
                
    def update_q_sensitivity(self, trans_data):
        self.q_sensitivity = ExtrapQSensitivity()
        self.q_sensitivity.populate_data(trans_data, self.sel_fit)
        
    def change_fit_method(self, trans_data, new_fit_method, n, kargs = None):
        '''Function to change the extrapolation methods and update the discharge sensitivity computations'''
        self.fit_method = new_fit_method
        self.sel_fit = SelectFit()
        self.sel_fit.populate_data(self.norm_data, new_fit_method, kargs)
        self.q_sensitivity = ExtrapQSensitivity()
        self.q_sensitivity.populate_data(trans_data, self.sel_fit)
        
    def change_threshold(self, trans_data, data_type, threshold):
        '''Function to change the threshold for accepting the increment median as valid.  The threshold
        is in percent of the median number of points in all increments'''
        
        self.threshold = threshold
        self.process_profiles(trans_data, data_type)
        self.q_sensitivity = ExtrapQSensitivity()
        self.q_sensitivity.populate_data(trans_data, self.sel_fit)
        
        
    def change_extents(self, trans_data, data_type, extents):
        '''Function allows the data to be subsection by specifying the percent cumulative discharge
        for the start and end points.  Currently this function does not consider transect direction'''
        
        self.subsection = extents
        self.process_profiles(trans_data, data_type)
        self.q_sensitivity = ExtrapQSensitivity()
        self.q_sensitivity.populate_data(trans_data, self.sel_fit)
        
    def change_data_type(self, trans_data, data_type):
        self.process_profiles(trans_data, data_type)
        self.q_sensitivity = ExtrapQSensitivity(trans_data, self.selfit)
        
        
        
            
            
        