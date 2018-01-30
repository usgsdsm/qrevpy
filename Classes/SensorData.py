"""
Created on Sep 14, 2017

@author: gpetrochenkov
"""

class SensorData(object):
    """the class stores typically time series data for (pitch, roll,
    temperature, salinity, and speed of sound and its source"""
    
    def __init__(self):
        self.data = None
        self.data_orig = None
        self.source = None
        
        
    def populate_data(self,data_in, source_in):
        
        self.data = data_in
        self.data_orig = data_in
        self.source = source_in
        
    def change_data(self, data_in):
        self.data = data_in
        
    def set_source(self, source_in):
        self.source = source_in
