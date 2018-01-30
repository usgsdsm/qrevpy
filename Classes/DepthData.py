"""
Created on Jul 20, 2017

@author: gpetrochenkov

Modified DSM 1/24/2018
    - removed kargs
    - added docstrings
"""
import numpy as np
from numpy.matlib import repmat
from MiscLibs.lowess import lowess

class DepthData(object):
    """Process and store depth data.
    Supported sources include bottom track
    vertical beam, and external depth sounder.

    Attributes
    ----------
        depth_orig_m: np.array
            Original multi-beam depth data from transect file (includes draft_orig) in meters.
        depth_beams_m: np.array
            Depth data from transect file adjusted for any draft changes, in meters.
        depth_processed_m: np.array
            Depth data filtered and interpolated.
        depth_freq_kHz: float
            Defines acoustic frequency used to measure depth.
        depth_invalid_index:
            Index of depths marked invalid.
        depth_source: str
            Source of depth data ("BT", "VB", "DS").
        depth_source_ens: str
            Source of each depth value ("BT", "VB", "DS", "IN").
        draft_orig_m: float
            Original draft from data files, in meters.
        draft_use_m: float
            Draft used in computation of depth_beams_m and depth_cell_depths_m.
        depth_cell_depth_orig_m: np.array
            Depth to centerline of depth cells in raw data, in meters.
        depth_cell_depth_m: np.array
            Depth to centerline of depth cells adjusted for draft or speed of sound changes, in meters.
        depth_cell_size_orig_m: np.array
            Size of depth cells in meters from raw data, in meters.
        depth_cell_size_m:
            Size of depth cells adjusted for draft or speed of sound changes, in meters.
        smooth_depth: np.array
            Smoothed beam depth, in meters.
        smooth_upper_limit: np.array
            Smooth function upper limit of window, in meters.
        smooth_lower_limit: np.array
            Smooth function lower limit or window, in meters.
        avg_method:str
            Defines averaging method: "Simple", "IDW", only applicable to bottom track.
        avg_depth:??????????????????
        filter_type: str
            Type of filter: "None", "TRDI", "Smooth".
        interp_type: str
            Type of interpolation: "None", "Linear", "Smooth".
        valid_data_method: str
            QRev or TRDI.
        valid_data: np.array
            Logical array of valid mean depth for each ensemble.
        valid_beams: np.array
            Logical array, 1 row for each beam identifying valid data.
    """
    
    def __init__(self):
        self.depth_orig_m = None #Original multi-beam depth data from transect file (includes draft_orig) in meters
        self.depth_beams_m = None #Depth data from transect file adjusted for any draft changes, in meters
        self.depth_processed_m = None # Depth data filtered and interpolated
        self.depth_freq_hz = None #Defines ADCP frequency used of each raw data point
        self.depth_invalid_index = None # Index of depths marked invalid
        self.depth_source = None # Source of depth data ("BT", "VB", "DS")
        self.depth_source_ens = None # Source of each depth value ("BT", "VB", "DS", "IN")
        self.draft_orig_m = None # Original draft from data files, in meters
        self.draft_use_m = None # Draft used in computation of depth_beams_m and depth_cell_depths_m
        self.depth_cell_depth_orig_m = None # Depth cell range from the transducer, in meters
        self.depth_cell_depth_m = None # Depth to centerline of depth cells, in meters
        self.depth_cell_size_orig_m = None # Size of depth cells in meters from raw data
        self.depth_cell_size_m = None # Size of depth cells in meters
        self.smooth_depth = None # Smoothed beam depth 
        self.smooth_upper_limit = None # Smooth function upper limit of window 
        self.smooth_lower_limit = None # Smooth function lowerl limit or window 
        self.avg_method = None # Defines averaging method: "Simple", "IDW"
        self.avg_depth = None
        self.filter_type = None # Type of filter: "None", "TRDI", "Smooth"
        self.interp_type = None # Type of interpolation: "None", "Linear", "Smooth"
        self.valid_data_method = None # QRev or TRDI 
        self.valid_data = None # Logical array of valid mean depth for each ensemble
        self.valid_beams = None # Logical array, 1 row for each beam identifying valid data
        
    def populate_data(self, depth_in, source_in, freq_in, draft_in, cell_depth_in, cell_size_in):
        """Stores data in DepthData.

        Parameters
        ----------
        depth_in: np.array
            Raw depth data, in meters.
        source_in: str
            Source of raw depth data.
        freq_in: float
            Acoustic frequency used to measure depths, in kHz.
        draft_in: float
            Draft of transducer used to measure depths, in meters.
        cell_depth_in: np.array
            Depth to centerline of each depth cell, in meters. If source does not have depth cells the depth cell depth
            from bottom track should be used.
        cell_size_in: np.array
            Size of each depth cell, in meters. If source does not have depth cells the depth cell size
            from bottom track should be used.
        """

        self.depth_orig_m = depth_in
        self.depth_beams_m = depth_in
        self.depth_source = source_in
        self.depth_source_ens = np.array([source_in] * depth_in.shape[-1])
        self.depth_freq_hz = freq_in
        self.draft_orig_m = draft_in
        self.draft_use_m = draft_in
        self.filter_type = 'None'
        self.interp_type = 'None'
        self.valid_data_method = 'QRev'
        
        #For BT data assign additional properties
        if source_in == 'BT':
            self.avg_method = 'IDW'
        else:
            self.avg_method = 'None'
        # DSM changed 1/24/2018
        self.depth_cell_depth_orig_m = cell_depth_in
        self.depth_cell_size_orig_m = cell_size_in
        self.depth_cell_size_m = cell_size_in
        self.depth_cell_depth_m = cell_depth_in

# ===============================
# DSM stopped here 1/24/2018
        self.apply_filter('dummy', filter_type='None')
        

    def change_draft(self, draft):
        """Changes the draft for object
        
        draft: new draft for object
        """
        #Compute draft change
        draft_change = draft - self.draft_use_m
        self.draft_use = draft
        
        #Apply draft to ensemble depths if BT or VB
        if self.depth_source == 'DS':
            self.depth_beams_m = self.depth_beams_m + draft_change
            self.depth_processed_m = self.depth_processed_m + draft_change 
            
        #Apply draft to depth cell locations
        if len(self.depth_cell_depth_m) > 0:
            self.depth_cell_depth_m = self.depth_cell_depth_m + draft_change

    def add_cell_data(self, bt_depths):
        """Adds cell data to depth cources with no cell data.
        such as the vertical beam and depth sounder.  This allows
        a single object to contain all the required depth data
        
        bt_depths: object of DepthData with bottom track depths
        """
        
        self.depth_cell_depth_orig_m = bt_depths.depth_cell_depth_orig_m
        self.depth_cell_size_m = bt_depths.depth_cell_size_m
        self.depth_cell_depth_m = bt_depths.depth_cell_depth_m
        
    def set_avg_method(self, setting):
        self.avg_method = setting
        
    def set_valid_data_method(self, setting):
        self.valid_data_method = setting
        
    def compute_avg_BT_depth(self, kargs=None):
        """Computes average depth for BT_Depths
        
        kargs: averaging method (Simple or IDW)
        """
        
        if self.depth_source == 'BT':
            
            if kargs is not None:
                self.avg_method = kargs[0]
            
            #Get valid depths
            depth = self.depth_beams_m
            depth[self.valid_beams == False] = np.nan
            
            #Compute average depths
            self.depth_processed_m = self.average_depth(depth, self.draft_use_m, self.avg_method)
            self.depth_processed_m = self.avg_depth
            
            #Set depths to nan if depth are not valid beam depths
            self.depth_processed_m[self.valid_data == 0] = np.nan

    def apply_filter(self, transect, filter_type=[]):
        """Coordinate the application of depth filters.

        Parameters
        ----------
        transect: TransectData
            Object of transect data.
        filter_type: str
            Type of filter to apply (None, Smooth, TRDI).
        """
        # DSM changed 1/24/2018 if kargs is not None:
        #     setting = kargs[0]
        # else:
        #     setting = self.filter_type
        if filter_type:
            self.filter_type = filter_type

        # Compute selected filter
        if self.filter_type == 'None':
            # No filter
            self.filter_none()
        elif self.filter_type == 'Smooth':
            # Lowess smooth filter
            # DSM changed 1/24/2018 self.filter_smooth()
            self.filter_smooth(transect)
        elif self.filter_type == 'TRDI':
            # TRDI filter for multiple returns
            # DSM changed 1/24/2018  self.filter_trdi(transect)
            self.filter_trdi()
            self.filter_type = 'TRDI'
            
        self.valid_mean_data()
        # Update processed depth with filtered results
        if self.depth_source == 'BT':
            self.compute_avg_BT_depth()
        else:
            self.depth_processed_m = np.array(self.depth_beams_m)
            # DSM changed 1/24/2018 self.depth_processed_m[np.squeeze(self.valid_data) == 0] = np.nan
            self.depth_processed_m[self.valid_data == 0] = np.nan
            
    def apply_interpolation(self, transect, kargs=None):
        """Coordinates application of interpolations
        
        transect: object of TransectData
        kargs: interpolation setting
        """
        
        #Determine filters to apply
        if kargs is not None:
            setting = kargs[0]
        else:
            setting = self.interp_type
            
        #Apply selected interpolation
        
        #No filtering
        if setting == 'None':
            self.interpolate_none()
        #Hold last valid depth indefinitely
        elif setting == 'HoldLast':
            self.interpolate_hold_last()
        #Use values form a Loess smooth
        elif setting == 'Smooth':
            self.interpolate_smooth()
        #Linear interpolation
        else:
            self.interpolate_linear(transect)
            
        #Identify ensembles with interpolated depths
        idx = np.where(self.valid_data[:] == False)
        if len(idx) > 0:
            idx = idx[0]
            idx2 = np.where(np.isnan(self.depth_processed_m[idx]) == False)
            if len(idx2) > 0:
                idx2 = idx2[0]
                self.depth_source_ens[idx[idx2]] = 'IN'
        
    def apply_composite(self,comp_depth, comp_source):
        """Applies the data from CompDepth compted in DepthStructure
        to DepthData object
        
        comp_depth: composite depth computed in DepthStructure
        comp_source: source of composite depth (BT, VB, DS)
        """
        
        #Assign composite depth to property
        self.depth_processed_m = comp_depth
        
        #Assign appropriate composite source for each ensemble
        self.depth_source_ens[comp_source == 1] = 'BT'
        self.depth_source_ens[comp_source == 2] = 'VB'
        self.depth_source_ens[comp_source == 3] = 'DS'
        self.depth_source_ens[comp_source == 4] = 'IN'
        self.depth_source_ens[comp_source == 0] = 'NA'
        
    def sos_correction(self, ratio):
        """Correct depth for new speed of sound setting
        
        ratio: ration of new to old speed of sound value
        """
        
        # Correct unprocessed depths
        self.depth_beams_m = self.draft_use_m+np.multiply(self.depth_beams_m-self.draft_use_m,ratio)
        
        #Correct processed depths
        self.depth_processed_m = self.draft_use_m+np.multiply(self.depth_processed_m-self.draft_use_m,ratio)
        
        #Correct cell size and location
        self.depth_cell_size = np.multiply(self.depth_cell_size_m,ratio)
        self.depth_cell_depth_m = self.draft_use_m + np.multiply(self.depth_cell_depth_m - self.draft_use_m, ratio)
        
    #------------Private Methods-------------------------------------------
        
    def valid_mean_data(self):
        
        if self.depth_source == 'BT':
            self.valid_data = np.array([True for x in range(self.valid_beams.shape[1])])
            nvalid = np.sum(self.valid_beams, axis=0)
            
            if self.valid_data_method == 'TRDI':
                self.valid_data[nvalid<3] = False
            else:
                self.valid_data[nvalid<2] = False
        else:
            self.valid_data = self.valid_beams
            
    def filter_none(self):
        """Applies no filter to depth data"""
        
        #Set all ensembles to have valid data
        self.valid_beams = np.array([True for x in range(self.depth_beams_m.shape[0])])
        
        #Set ensembles with no depth data to invalid
        self.valid_beams[self.depth_beams_m == 0] = False
        self.valid_beams[np.isnan(self.depth_beams_m)] = False
        
        self.filter_type = 'None'
        
    def filter_smooth(self, transect):
        """Filters depths Measured by each beam
        This filter uses a moving InterQuartile Range filter on residuals from a robust
        Lowess smooth of the depths in each beam to identify unnatural spikes in the depth
        measurements from each beam.  Each beam is filtered independently.  The filter
        criteria are set to be the maximum of the IQR filter, 5% of the measured depth, or 0.1 meter
        
        Reccomended filter settings:
        filter_width = 20
        half_width = 10
        multiplier = 15
        
        # TODO: Make a user defined settings for each of the parameters
        
        Input Requirements:
        
            depth_raw - Bt.depth_beams_m - number of beams (row) by number of ensembles
                matrix (col) of measured depths in m
                
            filter width - number of points used in robust Lowess smooth
            
            half_width - number of points to each side of target point used in computing IQR.
                This is the raw number of points actual points used may be less if some are bad.
                
            multiplier - number multiplied times the IQR to determine the filter criteria
        
        """
        #If the smoothed depth has not been computed
        if self.smooth_depth is None:
            
            #Set filter charactersitics
            self.filter_type = 'Smooth'
            cycles = 3
            filter_width = 20
            half_width = 10
            multiplier = 15
            
            #Determine number of beams
            n_beams, n_ensembles = self.depth_orig_m.shape[0], self.depth_orig_m.shape[1]
            depth_raw = self.depth_orig_m
            
            #Set bad depths to nan
            depth = repmat([np.nan], n_beams, n_ensembles)
            
            #Arrays initialized
            depth_smooth = np.nan
            
            #Arrays initialized
            depth_smooth = repmat([np.nan], n_beams, n_ensembles)
            depth_res = repmat([np.nan], n_beams, n_ensembles)
            depth_filter = repmat([np.nan], n_beams, n_ensembles)
            upper_limit = repmat([np.nan], n_beams, n_ensembles)
            lower_limit = repmat([np.nan], n_beams, n_ensembles)
            bad_depth = repmat([np.nan], n_beams, n_ensembles)
            depth_filtered = depth
            
            #Create position array
            if transect.boat_vel[transect.boatvel.selected] is not None:
                boat_vel_x = transect.boat_vel[transect.boat_vel.selected]['u_processed_mps']
                boat_vel_y = transect.boat_vel[transect.boat_vel.selected]['v_processed_mps']
                track_x = boat_vel_x * transect.datetime['ens_duraction_sec']
                track_y = boat_vel_y * transect.datetime['ens_duraction_sec']
            else:
                track_x = np.nan
                track_y = np.nan
            
            #Create position array
            idx = np.where(np.isnan(track_x))
            if len(idx[0]) < 1:
                x = np.nancumsum(np.sqrt(track_x**2+track_y**2))
            else:
                x = np.nancumsum(transect.datetime['ens_duration_sec'])
                
            #Loop for each beam, smooth is applied to each beam
            for j in range(n_beams):
                if np.nansum(depth_filtered[j] / depth_filtered[j].shape[0]) > .5:
                    #Compute residuals based on robust loess smooth
                    if len(x) > 1:
                        depth_smooth[j] = lowess(x,depth_filtered[j],filter_width/len(depth_filtered[j]),1)
                    else:
                        depth_smooth[j] = depth_filtered[j]
                    
                    depth_res[j] = depth[j] - depth_smooth[j]
                    
                    #Run the filter multiple times
                    for n in range(cycles):
                        
                        #Compute inner quartile range
                        fill_array = self.runIQR(half_width,depth_res[j])
                        
                        #Compute filter criteria and apply appropriate
                        criteria = multiplier * fill_array
                        idx = np.where(criteria < np.max([depth[j] * .05,np.ones(depth.shape) / 10]))
                        criteria[idx] = np.max([depth[j] * .05,np.ones(depth.shape) / 10])
                        
                        #compute limits
                        upper_limit[j] = depth_smooth[j] + criteria
                        lower_limit[j] = depth_smooth[j] - criteria
                        
                else:
                    depth_smooth[j] = np.nan
                    upper_limit[j] = np.nan
                    lower_limit[j] = np.nan
                    
            #Save smooth results to avoid recomputing them if needed later
            self.smooth_depth = depth_smooth
            self.smooth_upper_limit = upper_limit
            self.smooth_lower_limit = lower_limit
        
        #Reset valid data
        self.filter_none()
        
        #Set filter type
        self.filter_type = 'Smooth'
        
        #Determine number of beams
        n_beams, n_ensembles = self.depth_orig_m.shape
        
        #Get depth
        depth_raw = np.array(self.depth_orig_m)
        
        #Set bad depths to nan
        depth = repmat([np.nan, depth_raw.shape[0], depth_raw.shape[1]])
        depth[depth_raw > 0] = depth_raw[depth_raw > 0];
        
        #Apply filter
        for j in range(n_beams):
            if np.nansum(self.smooth_upper_limit[j]) > 0:
                bad_idx = np.where(depth[j] > self.smooth_upper_limit[j] | depth[j] < self.smooth_lower_limit[j])
                #Update depth matrix
                depth_res[j,bad_idx] = np.nan
                #Update valid data matrix
                self.valid_beams[j, bad_idx] = False
            else:
                bad_idx=np.isnan(depth[j])
                self.valid_beams[j, bad_idx] = False
                
                
    def interpolate_none(self):
        """Applies no interpolation"""
        
        #Compute processed depth without interpolation
        if self.depth_source == 'BT':
            #Bottom track methods
            self.compute_avg_BT_depth()
        else:
            #Vertical beam or depth sounder depths
            self.depth_processed_m = self.depth_beams_m
            
        self.depth_processed_m[not self.valid_data] = np.nan
        
        #Set interpolation type
        self.interp_type = 'None'
        
    def interpolate_hold_last(self):
        """This function folds the last valid value until the next valid data point"""
        
        #Get number of ensembles
        n_ensembles = len(self.depth_processed_m)
        
        #Process data by ensemble
        for n in range(1,n_ensembles):
            
            #If current ensemble's depth is invalid assign depth from previous example
            if np.isnan(self.depth_processed_m[n]):
                self.depth_processed_m[n] = self.depth_processed_m[n-1]
                
                
    def interpolate_smooth(self):
        """Apply interpolation based on the lwoess smooth"""
        
        self.interp_type = 'Smooth'
        
        #Get depth data from object
        depth_new = self.depth_beams_m
        
        #Update depth data with interpolated depths
        depth_new[not self.valid_beams] = self.smooth_depth[not self.valid_beams]
        
        #Compute processed depths with interpolated values
        if self.depth_source == 'BT':
            
            #Temporarily change self.depth_beams_m to compute average
            #for bottom track based depths
            temp_save = self.depth_beams_m
            self.depth_beams_m = depth_new
            self.computeAvgBTDepth()
            self.depth_beams_m = temp_save
            
    def interpolate_linear(self, transect):
        """Apply linear interpolation"""
        
        #Set interpolation type
        self.interp_type='Linear'
        
        #Create position array
        if transect.boat_vel[transect.boat_vel.selected] is None:
            boat_vel_x = transect.boat_vel[transect.boat_vel.selected].u_processed_mps
            boat_vel_y = transect.boat_vel[transect.boat_vel.selected].v_processed_mps         
            track_x = boat_vel_x * transect.datetime.ens_duration_sec
            track_y = boat_vel_y * transect.datetime.ens_duration_sec
        else:
            sizeu = transect.boat_vel[transect.boat_vel.selected].u_processed_mps.shape
            sizev = transect.boat_vel[transect.boat_vel.selected].v_processed_mps.shape
            track_x = repmat([np.nan], sizeu[0], sizeu[1])
            track_y = repmat([np.nan], sizev[0], sizev[1])
              
        idx = np.where(np.isnan(track_x[1:]))
        
            #If the navigation reference has no gaps use it for interpolation, if not use time 
        if len(idx[0]) < 1:
            x = np.nancumsum(np.sqrt(track_x**2 + track_y**2))
        else:
            #Compute accumulated time
            x = np.nancumsum(transect.datetime.ens_duration_sec)
            
        #Determine number of beams
        n_beams = self.depth_beams_m.shape[0]
        
#             Create strict monotonic arrays for depth and track by idnetifying duplicate
#             track values.  The first track value is used and the remaining duplicates
#             are set to nan.  The depth assogned to that first track value is the average
#             of all duplicates.  The depths for the duplicates are then set to nan.  Only
#             valid strictly montonic track and depth data are used for the input in to linear
#             interpolation.   Only the interpolated data for invalid depths are added
#             to the valid depth data to create depth_new
        
        depth_mono = self.depth_beams_m
        x_mono = x
        
        idx0=np.where(np.diff(x) == 0)
        if len(idx0[0]) > 0:
            #split array in to subarrays in proper sequence e.g [[2,3,4],[7,8,9]] etc.
            idx1 = np.add(np.where(np.diff(idx0) != 1)[0], 1)
            group = np.split(idx0[0], idx1[0])
            n_group = len(group)
            for k in range(n_group):
                indices = group[k]
                indices = np.append(indices, indices[-1] + 1)
                first_idx = indices[0]
                depth_avg = np.nanmean(depth_mono[:,indices], axis=2) #maybe axis = 1
                depth_mono[:,indices[0]] = depth_avg
                depth_mono[:,indices[2:]] = np.nan
                x[indices[2:]] = np.nan
                
            #Interpolate each beam
            depth_new = self.depth_beams_m
            for n in range(n_beams):
                valid_depth_mono = not np.isnan(depth_mono[n])
                valid_x_mono = self.valid_beams[n]
                valid_data = self.valid_beams
                
                valid = np.vstack([valid_depth_mono, valid_x_mono, valid_data])
                valid = np.sum(valid.T, axis=1)
                valid = np.where(valid == 3)[0]
                
                if len(valid) > 1:
                    #compute interpolation function from all valid data
                    depth_int = np.interp(x_mono[valid], depth_mono[valid], x_mono)
                    depth_new[n, not self.valid_beams[n]] = depth_int(n,not self.valid_beams[n])
                
                else:
                    #No valid data
                    depth_int[n] = repmat([np.nan],1,len(valid))
                    
                if self.depth_source == 'BT':
                    self.depth_processed_m = self.average_depth(depth_new,self.draft_use_m, self.avg_method)
                
                
        
    def filter_butter(self, transect):
        """Filters depth of every beam using a Butterworth filter at a specified order and cutoff frequency"""
        pass  
    
    def average_depth(self, depth, draft, method):
        """Compute average depth from bottom track beam depths
        
        depth - array of beam depths
        
        """
        if method == 'Simple':
            self.avg_depth = np.nanmean(depth)
        else:
            #compute inverse weighted mean depth
            rng = depth - draft
            w = 1-rng / repmat(np.nansum(rng), 4, 1)
            self.avg_depth = np.array(draft+np.nansum((rng*w) / repmat(np.nansum(w,0), 4, 1), 0))
            self.avg_depth[self.avg_depth == draft] = np.nan
            
           
            
    def run_IQR(self, half_width, mydata):
        """Computes a running Innerquartile Range
        The routine accepts a column vector as input.  "halfWidth" number of data
        points for computing the Innerquartile Range are selected before and
        after the target data point, but no including the target data point.
        Near the ends of the series the number of points before or after are reduced.
        Nan in the data are counted as points.  The IQR is computed on the slected
        subset of points.  The process occurs for each point in the provided column vector.
        A column vector with the computed IQR at each point is returned.
        
        half_width - number of ensembles before and after current ensemble which are used to compute the IQR
        my_data - data for which the IQR is computed
        """
        npts = len(mydata)
        
        if npts<20:
            half_width = np.floor(npts/2)
        
        iqr_array = []
         
        #Compute IQR for each point
        for n in range(npts):
            
            #Sample selection for 1st point
            if n == 0:
                sample = mydata[1:half_width]
                
            #Sample selection a end of data set
            elif n + half_width > npts:
                sample = np.hstack([mydata[n-half_width:n],mydata[n+1:npts]])
                
            #Sample selection at beginning of data set
            elif half_width >= n:
                sample = np.hstack([mydata[1:n],mydata[n+1:n+half_width]])
                
            #Sample selection in body of data set
            else:
                sample = np.hstack([mydata[n-half_width:n],mydata[n+1:n+half_width]])
                
            iqr_array.append(np.nanpercentile(sample, 75)-np.nanpercentile(sample, 25))
            
            return np.array(iqr_array)
        
            
        
    def filter_none(self):
        """Applies no filter to depth data"""
        
        #Set all ensembles to have valid data
        if len(self.depth_beams_m.shape) > 1:
            self.valid_beams = repmat([True], self.depth_beams_m.shape[0], self.depth_beams_m.shape[1])
        else:
            self.valid_beams = repmat([True], self.depth_beams_m.shape[0], 1)
        
        #Set ensembles with no depth data to invalid
        self.valid_beams[self.depth_beams_m == 0] = False
        self.valid_beams[np.isnan(self.depth_beams_m)] = False
        
        #Set filter type
        self.filter_type = 'None'
    
    def filter_smooth(self):
        pass
    
    def filter_trdi(self):
        pass
    
    def interpolate_none(self):
        """Applies no interpolation"""
        
        if self.depth_source == 'BT':
            self.compute_avg_BT_depth()
        else:
            self.depth_processed_m = self.depth_beams_m
            
        self.depth_processed_m[self.valid_data == False] = np.nan
        
        self.interp_type = 'None'
        
        
        
        