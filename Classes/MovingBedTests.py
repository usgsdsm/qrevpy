'''
Created on Sep 26, 2017

@author: gpetrochenkov
'''
import numpy as np
from Classes.TransectData import adjusted_ensemble_duration
from Classes.QComp import QComp
from MiscLibs.convenience import cart2pol, sind, pol2cart, rad2azdeg

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
        elif self.__type == 'Stationary':
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
        
        self.__messages = []
        vel_criteria = 0.012
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
            direct, flow_speed_q = cart2pol(se,sn)
            self.__flow_dir = rad2azdeg(direct)
            
            #compute the area weighted mean velocity components for the
            #purposed of computing the mean flow speed
            #SEEMS LIKE THE FLOW SPEED AND DIRECTION SHOULD BE HANDLED BY
            #THE SAME NOT DIFFERENTLY
            wght_area = np.multiply(np.multiply(np.sqrt(bt_U**2 + bt_V**2),  bin_size), ens_duration)
            idx = np.where(np.isnan(wt_U) == False)
            se = np.nansum(np.nansum(wt_U[idx] * wght_area[idx])) / np.nansum(np.nansum(wght_area[idx]))
            sn = np.nansum(np.nansum(wt_V[idx] * wght_area[idx])) / np.nansum(np.nansum(wght_area[idx]))
            dir_a, self.__flow_spd_mps = cart2pol(se,sn)
            flow_dir_a = rad2azdeg(dir_a)
            
            #Compute moving bed velocity and potential error in discharge
            #compute closure distance and direction
            bt_X = np.nancumsum(bt_U * ens_duration)
            bt_Y = np.nancumsum(bt_V * ens_duration)
            direct, self.__dist_us_m = cart2pol(bt_X[-1], bt_Y[-1])
            self.__mb_dir = rad2azdeg(direct)
            
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
                    consect_BT_time[n] = consect_BT_time[n-1] + ens_duration[n]
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
            u_out = wt_U[:,:dmg_idx+1]    
            v_out = wt_V[:,:dmg_idx+1]
            wght = np.abs(q[:,:dmg_idx+1])
            se = np.nansum(u_out * wght) / np.nansum(wght)
            sn = np.nansum(v_out * wght) / np.nansum(wght)  
            direct, _ = cart2pol(se, sn)
            flow_dir1 = rad2azdeg(direct)
            
            #Compute unweighted flow direction in each cell
            direct, _ = cart2pol(u_out, v_out)
            flow_dir_cell = rad2azdeg(direct)
            
            #compute difference from mean and correct to +/- 180
            v_dir_corr = flow_dir_cell - flow_dir1
            v_dir_idx = v_dir_corr > 180
            v_dir_corr[v_dir_idx] = 360-v_dir_corr[v_dir_idx]
            v_dir_idx = v_dir_corr < -180
            v_dir_corr[v_dir_idx] = 360 + v_dir_corr[v_dir_idx]
            
            #number of invalid weights
            idx2 = np.where(np.isnan(wght) == False)     
            nwght = len(idx2[0])             
            
            #Compute 95% uncertainty using wieghted standard deviation
            uncert1 = 2. * np.sqrt(np.nansum(np.nansum(wght * v_dir_corr**2)) / (((nwght - 1) * np.nansum(np.nansum(wght))) / nwght)) /  np.sqrt(nwght)
            
            #Compute flow direction on returning part of loop
            u_ret = wt_U[:, dmg_idx+1:]
            v_ret = wt_V[:, dmg_idx+1:]
            wght = np.abs(q[:, dmg_idx+1:])
            se = np.nansum(u_ret * wght) /  np.nansum(wght)
            sn = np.nansum(v_ret * wght) / np.nansum(wght)
            direct, _ = cart2pol(se, sn)
            flow_dir2 = rad2azdeg(direct)
            
            #compute unwieghted flow direction in each cell
            direct, _ = cart2pol(u_ret, v_ret)
            flow_dir_cell = rad2azdeg(direct)
            
            #Compute difference from mean and correct to +/- 180
            v_dir_corr = flow_dir_cell - flow_dir2
            v_dir_idx = v_dir_corr > 180
            v_dir_corr[v_dir_idx] = 360 - v_dir_corr[v_dir_idx]
            v_dir_idx = v_dir_corr < -180
            v_dir_corr[v_dir_idx] = 360 + v_dir_corr[v_dir_idx]
            
            #Number of valid weights
            idx2 = np.where(np.isnan(wght) == False)
            nwght = len(idx2[0])
            
            #Compute 95% uncertainty using weighted standard deviation
            uncert2 = 2.*np.sqrt(np.nansum(np.nansum(wght * (v_dir_corr)**2)) / (((nwght-1)*np.nansum(np.nansum(wght))) / nwght)) / np.sqrt(nwght)
            
            #Compute and report difference in flow direction
            diff_dir = np.abs(flow_dir1 - flow_dir2)
            if diff_dir > 180:
                diff_dir = diff_dir - 360
            self.__compass_diff_deg = diff_dir
            uncert = uncert1 + uncert2
            
            #Compute potential compass error
            idx = np.where(np.isnan(bt_X)==False)
            if len(idx[0]) > 0:
                idx = idx[0][-1]
            width = np.sqrt((bt_X[dmg_idx] - bt_X[idx] /2)**2 + (bt_Y[dmg_idx] - bt_Y[idx] / 2)**2)
            compass_error = (2*width * sind(diff_dir / 2) * 100) / (self.duration_sec * self.__flow_spd_mps)
            
            #Initialize message counter
            self.__test_quality = 'Good'
            
            
            #Low water velocity
            if self.__flow_spd_mps < 0.25:
                self.__messages.append('WARNING: The water velocity is less than recommended minimum for' \
                'this test and could cause the loop method to be inaccurate.  ' \
                'CONSIDER USING A STATIONARY TEST TO CHECK MOVING-BED CONDITIONS')
                self.__test_quality = 'Warnings'
                
            #Percent invalid bottom track
            if self.__percent_invlaid_BT > 20:
                self.__messages.append('ERROR: Percent invalid bottom track exceeds 20 percent. THE LOOP IS NOT ACCURATE. TRY A STATIONARY MOVING-BED TEST.')
                self.__test_quality = 'Errors'
            elif self.__percent_invlaid_BT > 5:
                self.__messages('WARNING: Percent invalid bottom track exceeds 5 percent. Loop may not be accurate. PLEASE REVIEW DATA.')
                self.__test_quality = 'Warnings'
                
            #More than 9 consecutive seconds of invalid BT
            if max_consect_BT_time > 9:
                self.__messages.append('ERROR: Bottom track is invalid for more than 9 consecutive seconds. THE LOOP IS NOT ACCURATE. TRY A STATIONARY MOVING-BED TEST.')
                self.__test_quality = 'Errors'
                
            if np.abs(compass_error) > 5 and np.abs(diff_dir) > 3 and np.abs(diff_dir) > uncert:
                self.__messages.append('ERROR: Difference in flow direction between out and back sections of loop could result in a 5 percent or greater error in final discharge. REPEAT LOOP AFTER COMPASS CAL. OR USE A STATIONARY MOVING-BED TEST.')
                self.__test_quality = 'Errors'
        
        else:
            self.__messages.append('ERROR: Loop has no valid bottom track data. REPEAT OR USE A STATIONARY MOVING-BED TEST.')  
            self.__test_quality = 'Errors'
            
        #If loop is valid then evaluate moving-bed condition
        if self.__test_quality != 'Errors':
            
            #Check minimum moving-bed velocity criteria
            if self.__mb_spd_mps > vel_criteria:
                #Check that closure error is in upstream direction
                if np.abs(self.__flow_dir - self.__mb_dir) > 135 and np.abs(self.__flow_dir - self.__mb_dir) < 225:
                    #Check if moving-bed is greater than 1% of the mean flow speed
                    if self.__percent_mb > 1:
                        self.__messages.append('Loop Indicates a Moving Bed -- Use GPS as reference. If GPS is unavailable or invalid use the loop method to correct the final discharge.')
                        self.__moving_bed = 'Yes'
                    else:
                        self.__messages.append('Moving Bed Velocity < 1% of Mean Velocity -- No Correction Recommended')
                        self.__moving_bed = 'No'
                else:
                    self.__messages.append('ERROR: Loop closure error not in upstream direction. REPEAT LOOP or USE STATIONARY TEST') 
                    self.__test_quality = 'Errors'
                    self.__moving_bed = 'Uknown'
            else:
                self.__messages.append('Moving-bed velocity < Minimum moving-bed velocity criteria -- No correction recommended')
                self.__moving_bed = 'No'
        else:
            self.__messages.append('ERROR: Due to ERRORS noted above this loop is NOT VALID. Please consider suggestions.')
            self.__moving_bed = 'Uknown'
            
            
    def stationary_test(self):
        '''Processed the stationary moving-bed tests'''
        #Assign data from treansect to local variables
        trans_data = self.__transect
        in_transect_idx = self.__transect.in_transect_idx
        n_ensembles = len(in_transect_idx)
        bt_valid = trans_data.boat_vel.bt_vel.valid_data[0,in_transect_idx]
        #Check to see that there is some valid bottom track data
        self.__messages = []
        if np.nansum(bt_valid) > 0:
            wt_U = trans_data.w_vel.__u_processed_mps[:, in_transect_idx]
            wt_V = trans_data.w_vel.__v_processed_mps[:, in_transect_idx]
            ens_duration = trans_data.datetime.ens_duration_sec[in_transect_idx]
            bt_U = trans_data.boat_vel.bt_vel.__u_processed_mps[in_transect_idx]
            bt_V = trans_data.boat_vel.bt_vel.__v_processed_mps[in_transect_idx]
            
            bin_depth = trans_data.depths.bt_depths.depth_cell_depth_m[:, in_transect_idx]
            trans_select = getattr(trans_data.depths, trans_data.depths.selected)
            depth_ens = trans_select.depth_processed_m[in_transect_idx]
            
            nb_U, nb_V, unit_NBU, unit_NBV = self.near_bed_velocity(wt_U, wt_V, depth_ens, bin_depth)
            
            #Compute bottom track parallel to water velocity
            unit_NB_vel = [[unit_NBU],[unit_NBV]]
            bt_vel = [[bt_U],[bt_V]]
            bt_vel_up_strm = -1 * np.dot(bt_vel, unit_NB_vel)
            bt_up_strm_dist = bt_vel_up_strm * ens_duration
            bt_up_strm_dist_cum = np.nancumsum(bt_up_strm_dist)
            
            #Compute bottom track perpendicular to water velocity
            nb_vel_ang, _ = cart2pol(unit_NBU, unit_NBV)
            nb_vel_unit_cs1, nb_vel_unit_cs2 = pol2cart(nb_vel_ang + np.pi / 2, np.ones(nb_vel_ang.shape))
            nb_vel_unit_cs = np.hstack([nb_vel_unit_cs1, nb_vel_unit_cs2])
            bt_vel_cs = np.dot(bt_vel, nb_vel_unit_cs.T,1)
            bt_cs_strm_dist = bt_vel_cs * ens_duration
            bt_cs_strm_dist_cum = np.nancumsum(bt_cs_strm_dist)
            
            #Compute cumulative mean moving bed velocity
            valid_bt_vel_up_strm = np.isnan(bt_vel_up_strm) == False
            mb_vel = np.nancumsum(bt_vel_up_strm) / np.nancumsum(valid_bt_vel_up_strm)
            
            #Compute the average ensemble velocities corrected for moving bed
            if mb_vel[-1] > 0:
                u_corrected = np.add(wt_U, (unit_NB_vel[0,:]) * bt_vel_up_strm)
                v_corrected = np.add(wt_V, (unit_NB_vel[1,:]) * bt_vel_up_strm)
            else:
                u_corrected = wt_U
                v_corrected = wt_V
                
            ''' Compute the mean of the ensemble magnitudes
                Mean is computed using magnitudes because if a Streampro with no
                 compass is the data source the change in direction could be
                either real change in water direction or an uncompensated turn of
                 the floating platform. This approach is the best compromise when
                 there is no compass or the compass is unreliable, which is often
                 why the stationary method is used. A weighted average is used
                 to account for the possible change in cell size within and
                 ensemble for the RiverRay and RiverPro.'''
            
            mag = np.sqrt(u_corrected**2 + v_corrected**2)
            depth_cell_size = trans_data.depths.bt_depths.depth_cell_size_m[:, in_transect_idx]
            depth_cell_size[np.isnan(mag)] = np.nan
            mag_w = mag * depth_cell_size
            avg_vel = np.nansum(mag_w) / np.nansum(depth_cell_size)
            pot_error_per = (mb_vel[-1] / avg_vel) * 100
            if pot_error_per < 0:
                pot_error_per = 0   
                
            #Compute percent invalid bottom track
            self.__percent_invlaid_BT = (np.nansum(bt_valid == False) / len(bt_valid)) * 100
            self.__dist_us_m = bt_up_strm_dist_cum[-1]
            self.duration_sec = np.nansum(ens_duration)
            self.__compass_diff_deg = []
            self.__flow_dir = []
            self.__mb_dir = []
            self.__flow_spd_mps = avg_vel
            self.__mb_spd_mps = mb_vel[-1]
            self.percent_mb = pot_error_per
            self.__near_bed_speed_mps = np.sqrt(np.nanmean(nb_U)**2 + np.nanmean(nb_V)**2)
            self.__stationary_us_track = bt_up_strm_dist_cum
            self.__stationary_cs_track = bt_cs_strm_dist_cum
            self.__stationary_mb_vel = mb_vel
            
            #Quality check
            self.test_quality = 'Good'
            #check duration
            if self.__duration_sec < 300:
                self.__messages.append('WARNING - Duration of stationary test is less than 5 minutes')
                self.__test_quality = 'Warnings'
                
            #Check validity of mean moving-bed velocity
            if self.duration_sec > 60:
                mb_vel_std = np.nanstd(mb_vel[-30:])
                cov = mb_vel_std / mb_vel[-1]
                if cov > 0.25 and mb_vel_std > 0.03:
                    self.__messages.append('WARNING - Moving-bed velocity may not be consistent. Average maybe inaccurate.')
                    self.__test_quality = 'Warnings'
                else:
                    cov = np.nan
                    
            #Check percentage of invalid BT data
            #if sum(bt_valid) <= 150
            if np.nansum(ens_duration[valid_bt_vel_up_strm]) <= 120:
                
                self.__messages.append('ERROR - Total duration of valid BT data is insufficient for a valid test.')
                self.__test_quality = 'Errors'
                self.__moving_bed = 'Unknown'
            elif self.__percent_invlaid_BT > 10:
                self.__messages.append('WARNING - Number of ensembles with invalid bottom track exceeds 10%')
                self.__test_quality = 'Warnings'
                
            #Determine if the test indicates a moving bed
            if self.__test_quality != 'Errors':
                if self.__percent_mb > 1:
                    self.__moving_bed = 'Yes'
                else:
                    self.__moving_bed = 'No'
                    
        else:
            self.__messages.append('ERROR - Stationary moving-bed test has no valid bottom track data.')
            self.__test_quality = 'Errors'
            self.__moving_bed = 'Unknown'
            
        
    def near_bed_velocity(self, u, v, depth, bin_depth):
        '''Compute near bed velocities'''
        #Compute z near bed as 10% of depth
        z_near_bed = depth * 0.1
        
        #Begin computing near-bed velocities
        
        n_ensembles = u.shape[1]
        nb_U = np.tile([np.nan], (1,n_ensembles))
        nb_V = np.tile([np.nan], (1,n_ensembles))
        unit_NBU = np.tile([np.nan], (1,n_ensembles))
        unit_NBV = np.tile([np.nan], (1,n_ensembles))
        z_depth = np.tile([np.nan], (1,n_ensembles)) 
        u_mean = np.tile([np.nan], (1,n_ensembles))
        v_mean = np.tile([np.nan], (1,n_ensembles))               
        speed_near_bed = np.tile([np.nan], (1, n_ensembles))
        for n in range(n_ensembles):
            idx = np.where(np.isnan(u[:,n])==False)
            if len(idx) > 0:
                idx = idx[1][-1]
                #compute near-bed velocity
                z_depth[n] = depth[n] - np.nanmean(bin_depth[idx,n])
                u_mean[n] = np.nanmean(u[idx,n])
                v_mean[n] = np.nanmean(v[idx,n])
                nb_U[n] = (u_mean[n] / z_depth[n]**(1./6.)) * (z_near_bed[n]**(1./6.))
                nb_V[n] = (v_mean[n] / z_depth[n]**(1./6.)) * (z_near_bed[n]**(1./6.))
                speed_near_bed[n] = np.sqrt(nb_U**2 + nb_V[n]**2)
                unit_NBU[n] = nb_U[n] / speed_near_bed[n]
                unit_NBV[n] = nb_V[n] / speed_near_bed[n]

        return (nb_U, nb_V, unit_NBU, unit_NBV)    
        
        
         
    
        
        