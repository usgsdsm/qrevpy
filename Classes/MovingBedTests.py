'''
Created on Sep 26, 2017

@author: gpetrochenkov
'''
import numpy as np
from Classes.TransectData import TransectData, allocate_transects, adjusted_ensemble_duration

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
        
        #Convert to earth coordinates and set the navigation reference to BT
        #for both boat and water data    
        self.__transect.change_coord_sys('Earth')
        self.__transect.change_nav_reference(1,'BT')
            
        #Adjust data for default manufacturer specific handling of invalid data
        delta_t = adjusted_ensemble_duration(self.__transect, 'mbt')
        
        if self.__type == 'Loop':
            if source == 'TRDI':
                self.loop_test(delta_t)
            else:
                self.loop_test()
        elif self.type == 'Stationary':
            self.stationary_test()
        else:
            pass
        
        
            
    def mb_TRDI(self, transect, mmt_transect):
        '''Function to create object properties for TRDI moving-bed tests'''
        
        self.__transect = transect
        self.__user_valid = True
        self.__type = mmt_transect.moving_bed_type
        
    def loop_test(self, kargs = None):
        '''Process loop moving bed test'''
        #Assign data from transect to local variables
        self.__transect.boat_interpolations(False, 'BT', kargs=['Linear'])
        self.__transect.boat_interpolations(False, 'GPS', kargs=['Linear'])
        trans_data = self.__transect
        in_transect_idx = trans_data.in_transect_idx
        n_ensembles = len(in_transect_idx)
        bt_valid = trans_data.boat_vel.bt_vel._BoatData__valid_data[0,in_transect_idx]
        #Check that there is some valid BT data
        if np.nansum(bt_valid) > 0:
            wt_U = trans_data.w_vel._WaterData__u_processed_mps[:,in_transect_idx]
            wt_V = trans_data.w_vel._WaterData__v_processed_mps[:,in_transect_idx]
            if kargs is None:
                ens_duration = trans_data.date_time.ens_duration_sec[in_transect_idx]
            else:
                ens_duration = kargs[:]
            bt_U = trans_data.boat_vel.bt_vel._BoatData__u_processed_mps[in_transect_idx]
            bt_V = trans_data.boat_vel.bt_vel._BoatData__v_processed_mps[in_transect_idx]
            bin_size = trans_data.depths.bt_depths.depth_cell_size_m[:,in_transect_idx]
            
            #Compute flow speed and direction
            #Compute discharge weighted mean velocity components for the
            #purposed of computing the mean flow direction
            
        
         
    
        
        