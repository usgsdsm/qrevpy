'''
Created on Sep 26, 2017

@author: gpetrochenkov
'''
import numpy as np
from Classes.TransectData import TransectData, allocate_transects, adjusted_ensemble_duration
from Classes.QComp import QComp
from MiscLibs.convenience import cart2pol

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
            qcomp = QComp()
            xprod = qcomp.cross_product(kargs=[trans_data])
            q = qcomp.discharge_middle_cells(xprod, trans_data, ens_duration)
            wght = np.abs(q)
            se = np.nansum(np.nansum(wt_U * wght)) / np.nansum(np.nansum(wght))
            sn = np.nansum(np.nansum(wt_V * wght)) / np.nansum(np.nansum(wght))
            dir, flow_speed_q = cart2pol(se,sn)
            self.__flow_dir = np.rad2deg(dir)
            
            #compute the area weighted mean velocity components for the
            #purposed of computing the mean flow speed
            #SEEMS LIKE THE FLOW SPEED AND DIRECTION SHOULD BE HANDLED BY
            #THE SAME NOT DIFFERENTLY
            wght_area = np.multiply(np.multiply(np.sqrt(bt_U**2 + bt_V**2),  bin_size), ens_duration)
            idx = np.where(np.isnan(wt_U) == False)
            se = np.nansum(np.nansum(wt_U[idx] * wght_area[idx])) / np.nansum(np.nansum(wght_area[idx]))
            sn = np.nansum(np.nansum(wt_V[idx] * wght_area[idx])) / np.nansum(np.nansum(wght_area[idx]))
            dir_a, self.__flow_spd_mps = cart2pol(se,sn)
            flow_dir_a = np.rad2deg(dir_a)
            
            #Compute moving bed velocity and potential erro in discharge
            #compute closure distance and direction
            bt_X = np.nancumsum(bt_U * ens_duration)
            bt_Y = np.nancumsum(bt_V * ens_duration)
            dir, self.__dist_us_m = cart2pol(bt_X[-1], bt_Y[-1])
            self.mb_dir = np.rad2deg(dir)
            
            #compute duration of test
            self.duration_sec = np.nansum(ens_duration)
            
            #Compute the moving-bed velocity
            self.__mb_spd_mps = self.__dist_us_m /  self.duration_sec
            
            #Compute potential error in BT referenced discharge
            self.__percent_mb = (self.__mb_spd_mps / (self.__flow_spd_mps + self.__mb_spd_mps)) * 100
            
            #Assess invalid bottom track
            #Compute percent invalid bottom track
            self.__percent_invlaid_BT = (np.nansum(bt_valid == False) / len(bt_valid)) * 100
            
            #Determine if more than 9 consecutive seconds of invalid BT occurred
            consect_BT_time = np.zeros(n_ensembles)
            for n in range(1, n_ensembles):
                if bt_valid[n] == False:
                    consect_BT_time = consect_BT_time(n-1) + ens_duration[n]
                else:
                    consect_BT_time[n] = 0
                    
            max_consect_BT_time = np.nanmax(consect_BT_time)
            
            #Evaluate compass calibration based on flow direction
            
            #Find apex of loop
            #adapted from
            #http://www.mathworks.de/matlabcentral/newsreader/view_thread/164048
            L1 = np.array([bt_X[0], bt_Y[0], 0])
            L2 = np.array([bt_X[-1], bt_Y[-1], 0])
            
            distance = np.zeros(n_ensembles)
            for n in range(n_ensembles):
                P = np.array([bt_X[n], bt_Y[n], 0])
                distance[n] = np.linalg.norm(np.cross(L2-L1,P-L1)) /  np.linalg.norm(L2-L1)
                
            dmg_idx = np.where(distance == np.nanmax(distance))[0][0]
            
            #Compute flow direction on outgoing part of loop
            u_out = wt_U[:,:dmg_idx]    
            v_out = wt_V[:,:dmg_idx]
            wght = np.abs(q[:,:dmg_idx])
            se = np.nansum(u_out * wght) / np.nansum(wght)
            sn = np.nansum(v_out * wght) / np.nansum(wght)  
            dir, _ = cart2pol(se, sn)
            flow_dir1 = np.rad2deg(dir)
            
            #Compute unweighted flow direction in each cell
            dir, _ = cart2pol(u_out, v_out)
            flow_dir_cell = np.rad2deg(dir)
            
            #compute difference from mean and correct to +/- 180
            v_dir_corr = flow_dir_cell - flow_dir1
            v_dir_idx = v_dir_corr > 180
            v_dir_corr[v_dir_idx] = 360-v_dir_corr[v_dir_idx]
            v_dir_idx = v_dir_corr < -180
            v_dir_corr[v_dir_idx] = 360 + v_dir_corr[v_dir_idx]
            
            #number of invalid weights
            idx2 = np.where(np.isnan(wght) == False)     
            nwght = len(idx2)             
            
            #Compute 95% uncertainty using wieghted standard deviation
            uncert1 = 2 * np.sqrt(np.nansum(np.nansum(wght * v_dir_corr**2)) / (((nwght - 1) * np.nansum(np.nansum(wght))) / wght)) /  np.sqrt(nwght)
            
            #Compute flow direction on returning part of loop
            u_ret = wt_U[:, dmg_idx:]
            v_ret = wt_V[:, dmg_idx:]
            wght = np.abs(q[:, dmg_idx:])
            se = np.nansum(u_ret * wght) /  np.nansum(wght)

            
            
        
         
    
        
        