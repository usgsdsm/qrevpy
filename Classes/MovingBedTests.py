'''
Created on Sep 26, 2017

@author: gpetrochenkov
'''
import numpy as np
from Classes.TransectData import TransectData, allocate_transects

class MovingBedTests(object):
    
    def __init__(self):
        
        self.__type = None #Loop or Stationary
        self.__transect = None #Object of TransectData
        self.__duration_sec = None #Duration of test in secs
        self.__percent_invlaid_BT = None #percent of invalid bottom track
        self.__compass_diff_deg = None #Difference in heading for out and back of loop
        self.__flow_dir = None #Mean flow direction from loop test
        self.__mb_dir = None #Moving bed or closure error direction 
        self.__dist_us_m = None #Distance moved upstream in m
        self.__flow_spd_mps = None # Magnitude of water velocity in mps
        self.__mb_spd_mps = None #Magnitude of moving=bed velocity in mps
        self.__percent_mb = None #Potential error due to moving bed in percent
        self.__moving_bed = None #Moving-bed determined 'Yes' 'No'
        self.__user_valid = None #Logical to allow user to determine if test should be considered a valid test
        self.__test_quality = None # Quality of test 'Valid' 'Warnings' 'Errors'
        self.__use_2_correct = None #Use this test to correct discharge
        self.__selected = None #Selected ad valid moving-bed test to user for correction or determine moving-bed condition
        self.__messages = None #Cell array of warning and error messages based on data processing
        self.__near_bed_speed_mps = None # Mean near-bed water speed for test in mps
        self.__stationary_us_track = None #Upstream component of the bottom track referenced ship track
        self.__stationary_cs_track = None #Cross=stream component of the bottom track referenced ship track
        self.__stationary_mb_vel = None #Moving-bed velocity by ensemble
        
    def populate_data(self, source, kargs = None):
        
        if source == 'TRDI':
            self.mb_TRDI(kargs[0], kargs[1])
        else:
            self.mb_SonTek(kargs)
            
    def mb_TRDI(self, mmt, transect):
        '''Function to create object properties for TRDI moving-bed tests'''
        
        self.__transect = transect
        self.__user_valid = True
        self.__type = transect.moving_bed_type
        
        
         
    
        
        