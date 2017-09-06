'''
Created on Sep 5, 2017

@author: gpetrochenkov
'''
import numpy as np
from numpy.matlib import repmat

class BoatData(object):
    '''Class to process and store boat velocity data'''
    
    def __init__(self):
        
        #Variables passed to the constructor
        self.__raw_vel_mps = None   #contains the raw unfiltered velocity data in m/s.
                                    #Rows 1-4 are beams 1,2,3,3 if if beam or u,v,w,d if otherwise 
        self.__frequency_hz = None #Defines ADCP frequency used for velocity Measurement
        self.__orig_coord_sys = None #Defines the original raw data velocity Coordinate
                                    #"Beam", "Inst", "Ship", "Earth"
        self.__nav_ref = None #Defines the original raw data navigation reference
                                #"None", "BT", "GGA" "VTG"
                                
        #Coordinate transformed data
        self.__coord_sys = None #Defines the current coordinate system "Beam", "Inst", "Ship", "Earth"
                                #For u, v, w, and d
        self.__u_mps = None #Horizontal velocity in x-direction, in m/s
        self.__v_mps = None #Horizontal velocity in y-direction, in m/s
        self.__w_mps = None #Vertical velocity (+ up), m/s
        self.__d_mps = None #Difference in vertical velocities compute from opposing beam pairs in m/s
        self.__num_invalid = None #Number of ensembles with invalid velocity data
        self.__bottom_mode = None #BT mode for TRDI, 'Variable' for SonTek
        
        #Processed data
        self.__u_processed_mps = None # Horizontal velocity in x-direction filtered and interpolated
        self.__v_processed_mps = None # Horizontal velocity in y-direction filtered and interpolated
        self.__proecessed_source = None
        
        #filter and interpolation properties
        self.__d_filter = None #Difference velocity filter "Manual", "Off", "Auto"
        self.__d_filter_threshold = None #Threshold for difference velocity filter
        self.__w_filter = None #Vertical velocity filter "On", "Off"
        self.__w_filter_threshold = None #Threshold for vertical velocity filter
        self.__gps_diff_qual_filter = None #Differential correction quality (1,2,4)
        self.__gps_altitude_filter = None # Change in altitude filter "Auto", "Manual", "Off"
        self.__gps_altitude_filter_change = None # Threshold from mean for altitude filter
        self.__gps_HDOP_filter = None # HDOP filter "Auto", "Manual", "Off"
        self.__gps_HDOP_filter_max = None # Max acceptable value for HDOP
        self.__gps_HDOP_filter_change = None # Maximum change allowed from mean
        self.__smooth_filter = None #Filter based on smoothing function 
        self.__smooth_speed = None #Smoothed boat speed
        self.__smooth_upper_limit = None #Smooth function upper limit of window
        self.__smooth_lower_limit = None #Smooth function lower limit of window
        self.__interpolate = None #Type of interpolation: "None", "Linear", "Smooth" etc.
        self.__beam_filter = None #3 for 3-beam solutions, 4 for 4-beam SolutionStackDescription
        self.__valid_data = None #Logical array of identifying valid and invalid data for each filter applied
                                # Row 1 [0] - composite
                                # Row 2 [1] - original
                                # Row 3 [2] - d_filter of diff_qual
                                # Row 4 [3] - w_filter or altitude
                                # Row 5 [4] - smooth_filter
                                # Row 6 [5] - beam_filter or HDOP
                                
    def populate_data(self, source, vel_in, freq_in, coord_sys_in, nav_ref_in, kargs = None):
        '''Creates object and sets properties with data provided
        
        Input:
        source: manufacturer (TRDI, SonTek)
        vel_in: boat velocity array
        freq_in: acoustic frequency boat velocity
        coord_sys_in: coordinate system of boat velocity
        nav_ref_in: source of boat velocity (BT, GGA, VTG)
        kargs: [0]: for TRDI data and specifies whether 3 beam solutions were selected
                    in the mmt file and the bottom mode.  These are not specified for SonTek data.
               [1]: bottom mode for TRDI ADCP
        '''
        
        if source == 'SonTek':
            vel_in = self.__filter_sontek(vel_in)
            
        self.__raw_vel_mps = vel_in
        self.__frequency_hz = freq_in
        self.__coord_sys = coord_sys_in
        self.__orig_coord_sys = coord_sys_in
        self.__nav_ref = nav_ref_in
        self.__beam_filter = 3
        
        #Bottom mode is variable unless TRDI with BT reference
        self.__bottom_mode = 'Variable'
        if source == 'TRDI' and nav_ref_in == 'BT':
            self.bottom_mode = kargs[1]
            #Apply 3-beam setting from mmt file
            if kargs[0] < .5:
                self.__beam_filter = 4
                
        if nav_ref_in == 'BT':
            
            #Boat velocities are referenced to ADCP not the streambed and thus must be reversed
            self.__u_mps = -1 * vel_in[0,:]
            self.__v_mps = -1 * vel_in[1,:]
            self.__w_mps = vel_in[2,:]
            self.__d_mps = vel_in[3,:]
            
            #Default filtering applied during initial construction of object
            self.__d_filter = 'Off'
            self.__d_filter_threshold = 99
            self.__wfilter = 'Off'
            self.__w_filter_threshold = 99
            self.__smooth_filter = 'Off'
            self.interpolate = 'None'
            
        else:
            
            #GPS referenced boat velocity
            self.__u_mps = vel_in[0,:]
            self.__v_mps = vel_in[1,:]
            self.__w_mps = np.nan
            self.__d_mps = np.nan
            
            #Default filtering
            self.__gps_diff_qual_filter = 2
            self.__gps_altitude_filter = 'Off'
            self.__gps_altitude_filter_change = 3
            self.__gps_HDOP_filter = 'Off'
            self.__gps_HDOP_filter_max = 2.5
            self.__gps_HDOP_filter_change = 1
            self.__smooth_filter = 'Off'
            self.__interpolate = 'None'
            
        #Assign data to processed property
        self.__u_processed_mps = self.__u_mps
        self.__v_processed_mps = self.__v_mps
     
        #Preallocate arrays
        n_ensembles = vel_in.shape[1]
        self.__valid_data = repmat([True], 6, n_ensembles)
        self.__smooth_speed = repmat([np.nan], 1, n_ensembles)
        self.__smooth_upper_limit = repmat([np.nan], 1, n_ensembles)
        self.__smooth_lower_limit = repmat([np.nan], 1, n_ensembles)
        
        #Determine number of raw invalid
        #--------------------------------
        # Find invalid raw data
        valid_vel = np.tile([True], self.__raw_vel_mps.shape)
        valid_vel[np.isnan(self.__raw_vel_mps)] = False
        
        #Identify invalid ensembles
        if nav_ref_in == 'BT':
            self.__valid_data[1, np.sum(valid_vel) < 3] = False
        else:
            self.__valid_data[2, np.sum(valid_vel) < 2] = False
            
        #Combine all filter data to composite valid data
        self.__valid_data[0,:] = np.all(self.__valid_data[1:,:])
        self.__num_invalid = np.sum(self._valid_data[0,:] == False)
        self.__processed_source = self.__u_mps.shape
        self.__processed_source[self.__valid_data[0,:] == True] = nav_ref_in
        self.__processed_source[self.__valid_data[0,:] == False] = "INT"
        
        
    def interpolate_composite(self, transect):
        '''This function interpolates processed data flagged invalid using linear interpolation
        
        Input:
        transect: object of TransectData
        '''
        u = self.__u_processed_mps
        v = self.__v_processed_mps
        
        valid = np.isnan(u) == False
        
        #Check for valid data
        if np.sum(valid) > 1:
            
            #Compute ensTime
            ens_time = np.nancumsum(transect.datetime.ens_duration_sec)
        
            
        
            
        
            
        
                
        