"""
Created on Sep 26, 2017

@author: gpetrochenkov
"""
import numpy as np
from Classes.TransectData import adjusted_ensemble_duration
from Classes.TransectData import TransectData
from Classes.QComp import QComp
from MiscLibs.convenience import cart2pol, sind, pol2cart, rad2azdeg
from Classes.MatSonTek import MatSonTek


class MovingBedTests(object):
    """Stores and processes moving-bed tests.

    Attributes
    ----------
        type: str
            Loop or Stationary
        transect: object
            Object of TransectData
        duration_sec: float
            Duration of test, in secs
        percent_invalid_BT: float
            Percent of invalid bottom track
        compass_diff_deg: float
            Difference in heading for out and back of loop
        flow_dir: float
            Mean flow direction from loop test
        mb_dir: float
            Moving bed or closure error direction
        dist_us_m: float
            Distance moved upstream, in m
        flow_spd_mps: float
            Magnitude of water velocity, in mps
        mb_spd_mps: float
            Magnitude of moving=bed velocity, in mps
        percent_mb: float
            Potential error due to moving bed, in percent
        moving_bed: bool
            Moving-bed determined (True or False)
        user_valid: bool
            Boolean to allow user to determine if test should be considered a valid test (True or False)
        test_quality: str
            Quality of test, 'Valid' 'Warnings' 'Errors'
        use_2_correct: bool
            Use this test to correct discharge (True or False)
        selected: bool
            Selected as valid moving-bed test to use for correction or determine moving-bed condition
        messages: list
            List of strings for warning and error messages based on data processing
        near_bed_speed_mps: float
            Mean near-bed water speed for test, in mps
        stationary_us_track: np.array(float)
            Upstream component of the bottom track referenced ship track
        stationary_cs_track: np.array(float)
            Cross=stream component of the bottom track referenced ship track
        stationary_mb_vel: np.array(float)
            Moving-bed velocity by ensemble, m/s
    """
    
    def __init__(self):
        """Initialize class and instance variables."""

        self.type = None  # Loop or Stationary
        self.transect = None  # Object of TransectData
        self.duration_sec = None  # Duration of test in secs
        self.percent_invalid_BT = None  # Percent of invalid bottom track
        self.compass_diff_deg = None  # Difference in heading for out and back of loop
        self.flow_dir = None  # Mean flow direction from loop test
        self.mb_dir = None  # Moving bed or closure error direction
        self.dist_us_m = None  # Distance moved upstream in m
        self.flow_spd_mps = None  # Magnitude of water velocity in mps
        self.mb_spd_mps = None  # Magnitude of moving=bed velocity in mps
        self.percent_mb = None  # Potential error due to moving bed in percent
        self.moving_bed = None  # Moving-bed determined 'Yes' 'No'
        self.user_valid = True  # Logical to allow user to determine if test should be considered a valid test
        self.test_quality = None  # Quality of test 'Valid' 'Warnings' 'Errors'
        self.use_2_correct = None  # Use this test to correct discharge
        self.selected = None  # Selected ad valid moving-bed test to user for correction or determine moving-bed condition
        self.messages = None  # Cell array of warning and error messages based on data processing
        self.near_bed_speed_mps = None  # Mean near-bed water speed for test in mps
        self.stationary_us_track = None  # Upstream component of the bottom track referenced ship track
        self.stationary_cs_track = None  # Cross=stream component of the bottom track referenced ship track
        self.stationary_mb_vel = None  # Moving-bed velocity by ensemble
        
    def populate_data(self, source, file=None, type=None):
        """Process and store moving-bed test data.

        Parameters
        ----------
        source: str
            Manufacturer of ADCP, SonTek or TRDI
        file: object or str
            Object of TransectData for TRDI and str of filename for SonTek
        type: str
            Type of moving-bed test (Loop or Stationary)
        """

        if source == 'TRDI':
            self.mb_TRDI(file, type)
        else:
            self.mb_sontek(file, type)
        
        # Convert to earth coordinates and set the navigation reference to BT
        # for both boat and water data
        self.transect.change_coord_sys(new_coord_sys='Earth')
        self.transect.change_nav_reference(update=True, new_nav_ref='BT')
            
        # Adjust data for default manufacturer specific handling of invalid data
        delta_t = adjusted_ensemble_duration(self.transect, 'mbt')
        
        if self.type == 'Loop':
            if source == 'TRDI':
                self.loop_test(delta_t)
            else:
                self.loop_test()
        elif self.type == 'Stationary':
            self.stationary_test()
        else:
            raise ValueError('Invalid moving-bed test identifier specified.')
        


    def mb_TRDI(self, transect, type):
        """Function to create object properties for TRDI moving-bed tests

        Parameters
        ----------
        transect: object
            Object of TransectData
        type: str
            Type of moving-bed test."""
        
        self.transect = transect
        self.user_valid = True
        self.type = type

    def mb_sontek(self, file_name, type):
        self.type = type

        # Read Matlab file for moving-bed test
        rsdata = MatSonTek(file_name)

        # Create transect objects for each discharge transect
        self.transect = TransectData()
        self.transect.SonTek(rsdata, file_name)
        
    def loop_test(self, ens_duration=None):
        """Process loop moving bed test.

        Parameters
        ----------
        ens_duration: np.array(float)
            Duration of each ensemble, in sec
        """
        # Assign data from transect to local variables
        # TODO edit BoatData to remove use of kargs
        self.transect.boat_interpolations(update=False, target='BT', method='Linear')
        self.transect.boat_interpolations(update=False, target='GPS', method='Linear')
        trans_data = self.transect
        in_transect_idx = trans_data.in_transect_idx
        n_ensembles = len(in_transect_idx)
        bt_valid = trans_data.boat_vel.bt_vel.valid_data[0, in_transect_idx]

        # Set variables to defaults
        self.messages = []
        vel_criteria = 0.012

        # Check that there is some valid BT data
        if np.nansum(bt_valid) > 0:
            wt_u = trans_data.w_vel.u_processed_mps[:, in_transect_idx]
            wt_v = trans_data.w_vel.v_processed_mps[:, in_transect_idx]
            if ens_duration is None:
                ens_duration = trans_data.date_time.ens_duration_sec[in_transect_idx]
            # else:
            #     ens_duration = kargs[:]
            bt_u = trans_data.boat_vel.bt_vel.u_processed_mps[in_transect_idx]
            bt_v = trans_data.boat_vel.bt_vel.v_processed_mps[in_transect_idx]
            bin_size = trans_data.depths.bt_depths.depth_cell_size_m[:, in_transect_idx]
            
            # Compute flow speed and direction
            # Compute discharge weighted mean velocity components for the
            # purposed of computing the mean flow direction
            # qcomp = QComp()
            xprod = QComp.cross_product(transect=trans_data)
            q = QComp.discharge_middle_cells(xprod, trans_data, ens_duration)
            wght = np.abs(q)
            se = np.nansum(np.nansum(wt_u * wght)) / np.nansum(np.nansum(wght))
            sn = np.nansum(np.nansum(wt_v * wght)) / np.nansum(np.nansum(wght))
            direct, flow_speed_q = cart2pol(se,sn)
            self.flow_dir = rad2azdeg(direct)
            
            #compute the area weighted mean velocity components for the
            #purposed of computing the mean flow speed
            #SEEMS LIKE THE FLOW SPEED AND DIRECTION SHOULD BE HANDLED BY
            #THE SAME NOT DIFFERENTLY
            wght_area = np.multiply(np.multiply(np.sqrt(bt_u ** 2 + bt_v ** 2), bin_size), ens_duration)
            idx = np.where(np.isnan(wt_u) == False)
            se = np.nansum(np.nansum(wt_u[idx] * wght_area[idx])) / np.nansum(np.nansum(wght_area[idx]))
            sn = np.nansum(np.nansum(wt_v[idx] * wght_area[idx])) / np.nansum(np.nansum(wght_area[idx]))
            dir_a, self.flow_spd_mps = cart2pol(se,sn)
            flow_dir_a = rad2azdeg(dir_a)
            
            #Compute moving bed velocity and potential error in discharge
            #compute closure distance and direction
            bt_X = np.nancumsum(bt_u * ens_duration)
            bt_Y = np.nancumsum(bt_v * ens_duration)
            direct, self.dist_us_m = cart2pol(bt_X[-1], bt_Y[-1])
            self.mb_dir = rad2azdeg(direct)
            
            #compute duration of test
            self.duration_sec = np.nansum(ens_duration)
            
            #Compute the moving-bed velocity
            self.mb_spd_mps = self.dist_us_m / self.duration_sec
            
            #Compute potential error in BT referenced discharge
            self.percent_mb = (self.mb_spd_mps / (self.flow_spd_mps + self.mb_spd_mps)) * 100
            
            #Assess invalid bottom track
            #Compute percent invalid bottom track
            self.percent_invlaid_BT = (np.nansum(bt_valid == False) / len(bt_valid)) * 100
            
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
            u_out = wt_u[:, :dmg_idx + 1]
            v_out = wt_v[:, :dmg_idx + 1]
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
            u_ret = wt_u[:, dmg_idx + 1:]
            v_ret = wt_v[:, dmg_idx + 1:]
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
            self.compass_diff_deg = diff_dir
            uncert = uncert1 + uncert2
            
            #Compute potential compass error
            idx = np.where(np.isnan(bt_X)==False)
            if len(idx[0]) > 0:
                idx = idx[0][-1]
            width = np.sqrt((bt_X[dmg_idx] - bt_X[idx] /2)**2 + (bt_Y[dmg_idx] - bt_Y[idx] / 2)**2)
            compass_error = (2*width * sind(diff_dir / 2) * 100) / (self.duration_sec * self.flow_spd_mps)
            
            #Initialize message counter
            self.test_quality = 'Good'
            
            
            #Low water velocity
            if self.flow_spd_mps < 0.25:
                self.messages.append('WARNING: The water velocity is less than recommended minimum for' \
                'this test and could cause the loop method to be inaccurate.  ' \
                'CONSIDER USING A STATIONARY TEST TO CHECK MOVING-BED CONDITIONS')
                self.test_quality = 'Warnings'
                
            #Percent invalid bottom track
            if self.percent_invlaid_BT > 20:
                self.messages.append('ERROR: Percent invalid bottom track exceeds 20 percent. THE LOOP IS NOT ACCURATE. TRY A STATIONARY MOVING-BED TEST.')
                self.test_quality = 'Errors'
            elif self.percent_invlaid_BT > 5:
                self.messages('WARNING: Percent invalid bottom track exceeds 5 percent. Loop may not be accurate. PLEASE REVIEW DATA.')
                self.test_quality = 'Warnings'
                
            #More than 9 consecutive seconds of invalid BT
            if max_consect_BT_time > 9:
                self.messages.append('ERROR: Bottom track is invalid for more than 9 consecutive seconds. THE LOOP IS NOT ACCURATE. TRY A STATIONARY MOVING-BED TEST.')
                self.test_quality = 'Errors'
                
            if np.abs(compass_error) > 5 and np.abs(diff_dir) > 3 and np.abs(diff_dir) > uncert:
                self.messages.append('ERROR: Difference in flow direction between out and back sections of loop could result in a 5 percent or greater error in final discharge. REPEAT LOOP AFTER COMPASS CAL. OR USE A STATIONARY MOVING-BED TEST.')
                self.test_quality = 'Errors'
        
        else:
            self.messages.append('ERROR: Loop has no valid bottom track data. REPEAT OR USE A STATIONARY MOVING-BED TEST.')  
            self.test_quality = 'Errors'
            
        #If loop is valid then evaluate moving-bed condition
        if self.test_quality != 'Errors':
            
            #Check minimum moving-bed velocity criteria
            if self.mb_spd_mps > vel_criteria:
                #Check that closure error is in upstream direction
                if np.abs(self.flow_dir - self.mb_dir) > 135 and np.abs(self.flow_dir - self.mb_dir) < 225:
                    #Check if moving-bed is greater than 1% of the mean flow speed
                    if self.percent_mb > 1:
                        self.messages.append('Loop Indicates a Moving Bed -- Use GPS as reference. If GPS is unavailable or invalid use the loop method to correct the final discharge.')
                        self.moving_bed = 'Yes'
                    else:
                        self.messages.append('Moving Bed Velocity < 1% of Mean Velocity -- No Correction Recommended')
                        self.moving_bed = 'No'
                else:
                    self.messages.append('ERROR: Loop closure error not in upstream direction. REPEAT LOOP or USE STATIONARY TEST') 
                    self.test_quality = 'Errors'
                    self.moving_bed = 'Uknown'
            else:
                self.messages.append('Moving-bed velocity < Minimum moving-bed velocity criteria -- No correction recommended')
                self.moving_bed = 'No'
        else:
            self.messages.append('ERROR: Due to ERRORS noted above this loop is NOT VALID. Please consider suggestions.')
            self.moving_bed = 'Uknown'

    def stationary_test(self):
        """Processed the stationary moving-bed tests"""
        #Assign data from treansect to local variables
        trans_data = self.transect
        in_transect_idx = self.transect.in_transect_idx
        n_ensembles = len(in_transect_idx)
        bt_valid = trans_data.boat_vel.bt_vel.valid_data[0,in_transect_idx]
        #Check to see that there is some valid bottom track data
        self.messages = []
        if np.nansum(bt_valid) > 0:
            wt_U = trans_data.w_vel.u_processed_mps[:, in_transect_idx]
            wt_V = trans_data.w_vel.v_processed_mps[:, in_transect_idx]
            ens_duration = trans_data.date_time.ens_duration_sec[in_transect_idx]
            bt_U = trans_data.boat_vel.bt_vel.u_processed_mps[in_transect_idx]
            bt_V = trans_data.boat_vel.bt_vel.v_processed_mps[in_transect_idx]
            
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
                
            """ Compute the mean of the ensemble magnitudes
                Mean is computed using magnitudes because if a Streampro with no
                 compass is the data source the change in direction could be
                either real change in water direction or an uncompensated turn of
                 the floating platform. This approach is the best compromise when
                 there is no compass or the compass is unreliable, which is often
                 why the stationary method is used. A weighted average is used
                 to account for the possible change in cell size within and
                 ensemble for the RiverRay and RiverPro."""
            
            mag = np.sqrt(u_corrected**2 + v_corrected**2)
            depth_cell_size = trans_data.depths.bt_depths.depth_cell_size_m[:, in_transect_idx]
            depth_cell_size[np.isnan(mag)] = np.nan
            mag_w = mag * depth_cell_size
            avg_vel = np.nansum(mag_w) / np.nansum(depth_cell_size)
            pot_error_per = (mb_vel[-1] / avg_vel) * 100
            if pot_error_per < 0:
                pot_error_per = 0   
                
            #Compute percent invalid bottom track
            self.percent_invlaid_BT = (np.nansum(bt_valid == False) / len(bt_valid)) * 100
            self.dist_us_m = bt_up_strm_dist_cum[-1]
            self.duration_sec = np.nansum(ens_duration)
            self.compass_diff_deg = []
            self.flow_dir = []
            self.mb_dir = []
            self.flow_spd_mps = avg_vel
            self.mb_spd_mps = mb_vel[-1]
            self.percent_mb = pot_error_per
            self.near_bed_speed_mps = np.sqrt(np.nanmean(nb_U)**2 + np.nanmean(nb_V)**2)
            self.stationary_us_track = bt_up_strm_dist_cum
            self.stationary_cs_track = bt_cs_strm_dist_cum
            self.stationary_mb_vel = mb_vel
            
            #Quality check
            self.test_quality = 'Good'
            #check duration
            if self.duration_sec < 300:
                self.messages.append('WARNING - Duration of stationary test is less than 5 minutes')
                self.test_quality = 'Warnings'
                
            #Check validity of mean moving-bed velocity
            if self.duration_sec > 60:
                mb_vel_std = np.nanstd(mb_vel[-30:])
                cov = mb_vel_std / mb_vel[-1]
                if cov > 0.25 and mb_vel_std > 0.03:
                    self.messages.append('WARNING - Moving-bed velocity may not be consistent. Average maybe inaccurate.')
                    self.test_quality = 'Warnings'
                else:
                    cov = np.nan
                    
            #Check percentage of invalid BT data
            #if sum(bt_valid) <= 150
            if np.nansum(ens_duration[valid_bt_vel_up_strm]) <= 120:
                
                self.messages.append('ERROR - Total duration of valid BT data is insufficient for a valid test.')
                self.test_quality = 'Errors'
                self.moving_bed = 'Unknown'
            elif self.percent_invlaid_BT > 10:
                self.messages.append('WARNING - Number of ensembles with invalid bottom track exceeds 10%')
                self.test_quality = 'Warnings'
                
            #Determine if the test indicates a moving bed
            if self.test_quality != 'Errors':
                if self.percent_mb > 1:
                    self.moving_bed = 'Yes'
                else:
                    self.moving_bed = 'No'
                    
        else:
            self.messages.append('ERROR - Stationary moving-bed test has no valid bottom track data.')
            self.test_quality = 'Errors'
            self.moving_bed = 'Unknown'

    @staticmethod
    def auto_use_2_correct(moving_bed_tests, boat_ref=None):
        """Apply logic to determine which moving-bed tests should be used
        for correcting bottom track referenced discharges with moving-bed conditions.

        Parameters
        ----------
        moving_bed_tests: list
            List of MovingBedTests objects.
        boat_ref: str
            Boat velocity reference.

        Returns
        -----
        moving_bed_tests: list
            List of MovingBedTests objects.
        """
        
        # Initialize variables
        lidx_user = []
        lidx_no_errors = []
        test_type = []
        lidx_stationary = []
        lidx_loop = []
        flow_speed = []
        for test in moving_bed_tests:
            test.use_2_correct = False
            test.selected = False
            # Valid test according to user
            lidx_user.append(test.user_valid == True)
            # Valid test according to quality assessment
            lidx_no_errors.append(test.test_quality == 'Errors')
            # Identify type of test
            test_type.append(test.type)
            lidx_stationary.append(test_type == 'Stationary')
            lidx_loop.append(test_type == 'Loop')
            flow_speed.append(test.flow_spd_mps)

        # Combine
        lidx_valid_loop = np.all(np.vstack((lidx_user, lidx_no_errors, lidx_loop)))
        lidx_valid_stationary = np.all(np.vstack((lidx_user, lidx_no_errors, lidx_stationary)))
        
        # Check flow speed
        lidx_flow_speed = np.array(flow_speed) > 0.25
        
        # Determine if there are valid loop tests
        if np.any(lidx_valid_loop) and np.any(lidx_flow_speed):
            lidx_loops_2_select = np.all(np.vstack((lidx_flow_speed, lidx_valid_loop)), 0)

            # Select last loop
            idx_select = np.where(lidx_loops_2_select == True)[0]
            if len(idx_select) > 0:
                idx = np.where(lidx_loop == True)[0]
                idx_select = idx(idx_select[-1])
            else:
                idx_select = len(moving_bed_tests)
            test_select = moving_bed_tests[idx_select]
            test_select.selected = True
            
            if test_select.moving_bed == 'Yes':
                test_select.use_2_correct = True
                
        # If there are no valid loop loof for valid stationary tests
        elif np.any(lidx_valid_stationary):
            for lidx in lidx_valid_stationary:
                if lidx:
                    moving_bed_tests[lidx].selected = True
                    # Determine if any stationary test resulted in a moving bed
                    if moving_bed_tests[lidx].moving_bed == 'Yes':
                        moving_bed_tests[lidx].use_2_correct = True

        elif lidx_valid_loop:
            # Select last loop
            idx_select = np.where(lidx_valid_loop)[0][0]
            moving_bed_tests[idx_select].selected = True
            if moving_bed_tests[idx_select].moving_bed == 'Yes':
                moving_bed_tests[idx_select].use_2_correct = True
        
        # If the navigation reference for discharge computations is set
        # GPS then none of test should be used for correction. The
        # selected test should be used to determine if there is a valid
        # moving-bed and a moving-bed condition. 
        if boat_ref is None:
            ref = 'BT'
        else:
            ref = boat_ref

        if ref != 'BT':
            for test in moving_bed_tests:
                test.use_2_correct = False
        return moving_bed_tests