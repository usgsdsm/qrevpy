'''
Created on Sep 14, 2017

@author: gpetrochenkov
'''

class SensorData(object):
    '''the class stores typically time series data for (pitch, roll,
    temperature, salinity, and speed of sound and its source'''
    
    def __init__(self):
        self.__data = None
        self.__data_orig = None
        self.__source = None
        
        
    def populate_data(self,data_in, source_in):
        
        self.__data = data_in
        self.__data_orig = data_in
        self.__source = source_in
        
    def change_data(self, data_in):
        self.__data = data_in
        
    def set_source(self, source_in):
        self.__source = source_in
