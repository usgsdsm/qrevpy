'''
Created on Sep 26, 2017

@author: gpetrochenkov
'''
from Classes import ExtrapQSensitivity

class ComputeExtrap(object):
    '''Class to compute the optimized or manually specified extrapolation methods'''
    
    def __init__(self):
        self.threshold = None #threshold as a percent for determining if a median is valid 
        self.subsection = None #percent of discharge, does not account for transect direction
        self.fit_method = None #method used to determine fit.  Automatic or manual
        self.norm_data = None #object of class norm data
        self.sel_fit = None #Object of class SelectFit
        self.q_sensitivity = None #Object of class ExtrapQSensitivity
        self.messages = None #Variable for messages to UserWarning
        
    def populate_data(self, trans_data, kargs = None):
        
        self.threshold = 20
        self.subsection = [0, 100]
        self.fit_method = 'Automatic'
        self.process_profiles(trans_data, 'q')
        #Compute the sensitivity of the final discharge to changes in extrapolation methods
        if kargs is not None:
            self.q_sensitivity = ExtrapQSensitivity(trans_data,self.sel_fit)
            
    def process_profiles(self, trans_data, data_type):
        '''Function that serves and the main control for other classes and functions'''
        
        #Determine number of transects
        n_transects = len(trans_data)
        
        #Compute normalized data
        self.norm_data = NormData(trans_data, data_type, self.threshold, self.subsection)
            
        