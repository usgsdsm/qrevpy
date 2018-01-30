"""
Created on Sep 14, 2017

@author: gpetrochenkov
"""
import numpy as np
from _operator import xor

class HeadingData(object):
    """This class stores and manipulates self.data data"""
    
    def __init__(self):
        
        self.data  = None #Corrected self.data data
        self.original_data = None #original uncorrected self.data data
        self.source = None #Source of self.data data (internal, external)
        self.mag_var_deg = None #Magnetic variation for these self.data data
        self.mag_var_orig_deg = None #Original magnetic variation
        self.align_correction_deg = None #Alignment correction to align compass with instrument
        self.mag_error = None #Percent change in mean magnetic field from calibration`
        self.pitch_limit = None
        self.roll_limit = None
        
    def populate_data(self, data_in, source_in, kargs = None):
        arg = kargs[0]
        self.original_data = data_in
        self.source = source_in
        self.mag_var_deg = arg
        self.mag_var_orig_deg = arg
        self.align_correction_deg = kargs[1]
        #Correct the original data for the magvar and alignment
        self.data = self.original_data + self.mag_var_deg + self.align_correction_deg
        self.fix_upper_limit()
        self.interp_heading()
        if kargs is not None and len(kargs) > 2:
            self.mag_error = kargs[2]
            
    def set_mag_var(self, mag_var, h_source):
        """Applies a new magvar to the object"""
        self.mag_var_deg = mag_var
        if h_source == 'internal':
            self.data = self.original_data + self.mag_var_deg
            self.fix_upper_limit()
            
    def set_align_correction(self,align_correction, h_source):
        """Applies a new alignment correction to the object"""
        self.align_correction_deg = align_correction
        if h_source == 'external':
            self.data = self.original_data + self.align_correction_deg
            self.fix_upper_limit()
            
    def set_PR_Limit(self, type_prop, limits):
        setattr(self, type_prop, limits)
        
    def fix_upper_limit(self):
        idx = np.where(self.data > 360)[0]
        if len(idx) > 0:
            self.data[idx] = self.data[idx] - 360   
            
    def interp_heading(self):
        """Interpolate invalid self.data.  Use linear interpolation if there are
        valid values on either side of the invalid self.data.  If the invalid self.data
        occurs at the beginning of the time series back fill using the 1st valid.
        If the invalid self.data occurs at the end of the time series forward fill
        with the last valid self.data"""
        
        idx = np.where(np.isnan(self.data))[0]
        
        if len(idx) > 0:
            
            first_idx = np.where(np.isnan(self.data) == False)[0][0]
            last_idx = np.where(np.isnan(self.data) == False)[0][-1]
        
            #Process each invalid self.data
            for n in range(len(idx)):
                before_idx = np.where(np.isnan(self.data[1:idx[n]]) == False)[0]
                after_idx = np.where(np.isnan(self.data[idx[n]:]) == False)[0]
                
                #If invalid self.data is beginning back fill
                if len(before_idx) < 1:
                    self.data[idx[n]] = self.data[first_idx]
                #If invalid self.data is at end forward fill
                elif len(after_idx) < 1:
                    self.data[idx[n]] = self.data[last_idx]
                #If invalid self.data is in middle interpolate
                else:
                    before_idx = before_idx[-1]
                    after_idx = after_idx[0] + idx[n] - 1
                    
                    test1 = self.data[before_idx] > 180
                    test2 = self.data[after_idx] > 180
                    if not xor(test1, test2):
                        c = 0
                    elif test1 == True:
                        c = 360
                    elif test2 == True:
                        c = -360
                    self.data[idx[n]] = (((self.data[after_idx] - self.data[before_idx] + c) / (before_idx - after_idx)) * (before_idx - idx[n])) + self.data[before_idx]
                    if self.data[idx[n]] > 360:
                        self.data[idx[n]] - 360
            
            
    
            
    
            