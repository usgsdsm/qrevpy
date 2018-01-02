'''
Created on Sep 5, 2017

@author: gpetrochenkov
'''
import numpy as np
from numpy.matlib import repmat
from MiscLibs.convenience import cosd, sind, cart2pol, pol2cart, iqr
from MiscLibs.lowess import lowess

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
        self.__u_processed_mps = np.copy(self.__u_mps)
        self.__v_processed_mps = np.copy(self.__v_mps)
     
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
            self.__valid_data[1, np.sum(valid_vel,0) < 3] = False
        else:
            self.__valid_data[1, np.sum(valid_vel,0) < 2] = False
            
        #Combine all filter data to composite valid data
        self.__valid_data[0,:] = np.all(self.__valid_data[1:,:],0)
        self.__num_invalid = np.sum(self.__valid_data[0,:] == False)
        self.__processed_source = np.tile('',self.__u_mps.shape)
        self.__processed_source[np.where(self.__valid_data[0,:] == True)] = nav_ref_in
        self.__processed_source[np.where(self.__valid_data[0,:] == False)] = "INT"
        
    def change_coord_sys(self, new_coord_sys, sensors, adcp):
        '''This function allows the coordinate system to be changed.  Current implementation
        is only to allow change to a higher order coordinate system Beam - Inst - Ship - Earth
        
        Input:
        new_coord_sys: new coordinate_sys (Beam, Inst, Ship, Earth)
        sensors: object of Sensors
        adcp: object of InstrumentData
        '''
        
        #Remove any trailing spaces
        if isinstance(self.__orig_coord_sys, str):
            o_coord_sys = self.__orig_coord_sys.strip()
        else:
            o_coord_sys = self.__orig_coord_sys[0].strip()
        
        if self.__orig_coord_sys[0].strip() != new_coord_sys.strip():
            #Assign the transformation matrix and retrieve the sensor data
            t_matrix = adcp.t_matrix.matrix
            t_matrix_freq = adcp.frequency_hz
            
            pitch_select = getattr(sensors.pitch_deg, sensors.pitch_deg.selected)
            p =  getattr(pitch_select, '_SensorData__data')
            roll_select = getattr(sensors.roll_deg, sensors.roll_deg.selected)
            r = getattr(roll_select, '_SensorData__data')
            heading_select = getattr(sensors.heading_deg, sensors.heading_deg.selected)
            h = getattr(heading_select, 'data')
            
            #Modify the transformation matrix and heading, pitch, and roll values base on
            #the original coordinate system so that only the needed values are used in
            #computing the new coordinate system
            if o_coord_sys == 'Beam':
                orig_sys = 1
            elif o_coord_sys == 'Inst':
                orig_sys = 2
                t_matrix[:] = np.eye(t_matrix.shape[0])
            elif o_coord_sys == 'Ship':
                orig_sys = 3
                p = np.zeros(h.shape)
                r = np.zeros(h.shape)
                t_matrix[:] = np.eye(t_matrix.shape[0])
            elif o_coord_sys == 'Earth':
                orig_sys = 4
                
            #Assign a value to the new coordinate system
            if new_coord_sys == 'Beam':
                new_sys = 1
            elif new_coord_sys == 'Inst':
                new_sys = 2
            elif new_coord_sys == 'Ship':
                new_sys = 3
            elif new_coord_sys == 'Earth':
                new_sys = 4
                
            #Check to ensure the new coordinate system is a higher order than the original system
            if new_sys - orig_sys > 0:
                
                #Compute trig function for heaing, pitch and roll
                CH = cosd(h)
                SH = sind(h)
                CP = cosd(p)
                SP = sind(p)
                CR = cosd(r)
                SR = sind(r)
                
                vel_changed = np.tile([np.nan], self.__raw_vel_mps.shape)
                n_ens = self.__raw_vel_mps.shape[1]
                
                for ii in range(n_ens):
                    
                    #Compute matrix for heading, pitch, and roll
                    hpr_matrix = [[((CH[ii] * CR[ii]) + (SH[ii]*SP[ii]*SR[ii])),
                                (SH[ii] * CP[ii]),
                                ((CH[ii] * SR[ii]) - SH[ii]*SP[ii]*CR[ii])],
                                [(-1 * SH[ii] * CR[ii])+(CH[ii] * SP[ii] * SR[ii]),
                                CH[ii] * CP[ii], 
                                (-1 * SH[ii] * SR[ii])-(CH[ii] * SP[ii] * CR[ii])],
                                [(-1.*CP[ii] * SR[ii]),
                                SP[ii],
                                CP[ii] * CR[ii]]]
                    
                    #Transofm beam coordinates
                    if o_coord_sys == 'Beam':
                        
                        #Determine frequency index for transformation matrix
                        if len(t_matrix.shape) > 2:
                            idx_freq = np.where(t_matrix_freq==self.__frequency_hz[ii])
                            t_mult = np.copy(t_matrix[idx_freq])
                        else:
                            t_mult = np.copy(t_matrix)
                            
                        #Get velocity data
                        vel = np.squeeze(self.__raw_vel_mps[:,ii])
                        
                        #Check for invalid beams
                        idx_3_beam = np.where(np.isnan(vel))
                        
                        #3-beam solution
                        if len(idx_3_beam[0]) == 1:        
                            
                            #Special processing for RiverRay
                            if adcp.model == 'RiverRay':
                                
                                #Set beam pairing
                                beam_pair_1a = 0
                                beam_pair_1b = 1
                                beam_pair_2a = 2
                                beam_pair_2b = 3
                                
                                #Set speed of sound correction variables Note: Currently (2013-09-06) 
                                #WinRiver II does not use a variable correction and assumes the speed 
                                #of sound and the reference speed of sound are the same.
                                #sos = sensors.speed_ofs_sound_mps.selected.data[ii]
                                #sos_reference = 1536
                                #sos_correction = np.sqrt(((2 * sos_reference) / sos) **2 -1)
                                
                                sos_correction = np.sqrt(3)
                                
                                #Reconfigure transformation matrix based on which beam is invalid
                                
                                #Beam 1 invalid
                                if idx_3_beam[0][0] == beam_pair_1a:
                                    
                                    #Double valid beam in invalid pair
                                    t_mult[0:2, beam_pair_1b] *= 2
                                    
                                    #Eliminate invalid pair from vertical velocity computations
                                    t_mult[2,:] = [0, 0, 1/sos_correction, 1/sos_correction]
                                    
                                    #Reconstruct beam velocity matrix to use only valid beams
                                    t_mult = t_mult[0:3, [beam_pair_1b,beam_pair_2a,beam_pair_2b]]
                                    
                                    #Reconstruct beam velocity matrix to use only valid beams
                                    vel = vel[[beam_pair_1b, beam_pair_2a, beam_pair_2b]]
                                    
                                    #Apply transformation matrix
                                    temp_t = t_mult.dot(vel)
                                    
                                    #Correct horizontal velocity for invalid pair with the vertical velocity
                                    #and speed of sound correction
                                    temp_t[0] = temp_t[0] + temp_t[2] * sos_correction
                                
                                #Beam 2 invalid
                                if idx_3_beam[0][0] == beam_pair_1b:
                                    
                                    #Double valid beam in invalid pair
                                    t_mult[0:2, beam_pair_1a] = t_mult[0:2, beam_pair_1a] * 2
                                    
                                    #Eliminate invalid pair from vertical velocity computations
                                    t_mult[2,:] = [0, 0, 1/sos_correction, 1/sos_correction]
                                    
                                    #Reconstruct transformation matrix as a 3x3 matrix
                                    t_mult = t_mult[0:3, [beam_pair_1a, beam_pair_2a, beam_pair_2b]]
                                    
                                    #Reconstruct beam velocity matrix to use only valid beams
                                    vel = vel[[beam_pair_1a, beam_pair_2a, beam_pair_2b]]
                                    
                                    #Apply transformation matrix
                                    temp_t = t_mult.dot(vel)
                                    
                                    #Correct horizontal velocity for invalid pair with the vertical
                                    #velocity and speed of sound correction
                                    temp_t[0] = temp_t[0] - temp_t[2] * sos_correction
                                    
                                #Beam 3 invalid
                                if idx_3_beam[0][0] == beam_pair_2a:
                                    
                                    #Double valid beam in invalid pair
                                    t_mult[0:2, beam_pair_2b] = t_mult[:2, beam_pair_2b] * 2
                                    
                                    #Eliminate invalid pair from vertical velocity computations
                                    t_mult[2,:] = [1/sos_correction, 1/sos_correction, 0, 0]
                                    
                                    #Reconstruct transformation matrix as a 3x3 matrid
                                    t_mult = t_mult[:3, [beam_pair_1a, beam_pair_1b, beam_pair_2b]]
                                    
                                    #Reconstruct beam velocity matrix to use only valid beams
                                    vel = vel[[beam_pair_1a, beam_pair_1b, beam_pair_2b]]
                                    
                                    #Apply transformation matrix
                                    temp_t = t_mult.dot(vel)
                                    
                                    #Correct horizontal velocity for invalid pair with the vertical
                                    #velocity and speed of sound correction
                                    temp_t[1] = temp_t[1] - temp_t[2] * sos_correction
                                    
                                #Beam 4 invalid
                                if idx_3_beam[0][0] == beam_pair_2b:
                                    
                                    #Double valid beam in invalid pair
                                    t_mult[:2, beam_pair_2a] *= 2
                                    
                                    #Eliminate invalid pair from vertical velocity computations
                                    t_mult[2,:] = [1/sos_correction, 1/sos_correction, 0, 0]
                                    
                                    #Reconstruct transformations matrix as a 3x3 matrix
                                    t_mult = t_mult[:3, [beam_pair_1a, beam_pair_1b, beam_pair_2a]]
                                    
                                    #Reconstruct beam velocity matrix to use only valid beams
                                    vel = vel[[beam_pair_1a, beam_pair_1b, beam_pair_2a]]
                                    
                                    #Apply transformation matrix
                                    temp_t = t_mult.dot(vel)
                                    
                                    #Correct horiaontal velocity for invalid pair with the vertical
                                    #velocity and speed of sound correction
                                    temp_t[1] = temp_t[1] + temp_t[2] * sos_correction
                                    
                            else:
                                
                                #3 Beam solution for non-RiverRay
                                vel_3_beam_zero = vel
                                vel_3_beam_zero[np.isnan(vel)] = 0
                                vel_error = t_mult[3,:] * vel_3_beam_zero
                                vel[idx_3_beam] = -1 * vel_error / t_mult[3,idx_3_beam]
                                temp_t = t_mult.dot(vel)
                                
                            #apply transformation matrix for 3 beam solutions
                            temp_THPR = np.array(hpr_matrix).dot(temp_t[:3])
                            temp_THPR = np.hstack([temp_THPR, np.nan])
                            
                        else:
                            
                            #Apply transormation matrix for 4 beam solutions
                            temp_t = t_mult.dot(np.squeeze(self.__raw_vel_mps[:,ii]))
                            
                            #Apply hpr_matrix
                            temp_THPR = np.array(hpr_matrix).dot(temp_t[:3])
                            temp_THPR = np.hstack([temp_THPR, temp_t[3]])
                            
                    else:
                        
                        #Getvelocity data
                        vel = np.squeeze(self.__raw_vel_mps[:,ii])
                        
                        #Apply heading pitch roll for inst and ship coordinate data
                        temp_THPR = np.array(hpr_matrix).dot(vel[:3])
                        temp_THPR = np.hstack([temp_THPR, vel[3]])
                            
                    
                    vel_changed[:,ii] = temp_THPR.T
                
                #Assign results to object
                self.__u_mps = -1 * vel_changed[0,:]
                self.__v_mps = -1 * vel_changed[1,:]
                self.__w_mps = -1 * vel_changed[2,:]
                self.__d_mps = -1 * vel_changed[3,:]
                self.__coord_sys = new_coord_sys
                self.__u_processed_mps = np.copy(self.__u_mps)
                self.__v_processed_mps = np.copy(self.__v_mps)
                
                
                        
    def apply_interpolation(self, transect, kargs = None):
        '''Function to apply interpolations to navigation data
        
        Input:
        transect: object of TransectData
        kargs: specified interpolation method if different from that
        in saved object
        '''
        
        #Reset processed data
        self.__u_processed_mps = self.__u_mps
        self.__v_processed_mps = self.__v_mps               
        self.__u_processed_mps[self.__valid_data[0,:] == False] = np.nan
        self.__v_processed_mps[self.__valid_data[0,:] == False] = np.nan
        
        #Determine interpolation methods to apply
        interp = self.interpolate
        if kargs is not None:
            interp = kargs[0]
            
        #Apply specified interpolation method
        
        #Applied specified interpolation method
        self.interpolate = interp
        if interp == 'None':
            #sets invalid data to nan with no interpolation
            self.interpolate_none()
            
        elif interp == 'ExpandedT':
            #Set interpolate to none as the interpolation done is in the QComp
            self.interpolate_next()
            
        elif interp == 'Hold9': #Sontek Method
            #Interpolates using SonTeks method of holding last valid for up to 9 samples
            self.interpolate_hold_9()
            
        elif interp == 'HoldLast':
            #Interpolates by holding last valid indefinitely
            self.interpolate_hold_last()
            
        elif interp == 'Linear':
            #Interpolates using linear interpolation
            self.interpolate_linear(transect)
            
        elif interp == 'Smooth':
            #Interpolates using smooth interpolation
            self.interpolate_smooth(transect)
            
        elif interp == 'TRDI':
            #TRDI interpolation is done in discharge.
            # For TRDI the interpolation is done on discharge not on velocities
            self.interpolate_none()
            
            
            
    def __interpolate_hold_9(self):
        '''This function applies Sontek's  approach to maintaining the last valid boat speed
        for up to nine invalid samples
        '''
        
        #Initialize variables
        n_ensembles = self.__u_mps.shape[1]
       
        #Get data from object
        self.__u_processed_mps = self.__u_mps
        self.__v_processed_mps = self.__v_mps
        self.__u_processed_mps[self.__valid_data[0,:] == False] = np.nan
        self.__v_processed_mps[self.__valid_data[0,:] == False] = np.nan
        
        n_invalid = 0
        #Process data by ensembles
        for n in range(1,n_ensembles):
            #Check if ensemble is invalid and number of consecutive invalids is less than 9
            if self.__valid_data[0,n] == False and n_invalid < 9:
                self.__u_processed_mps = self.__u_processed_mps[n - 1]
                self.__v_processed_mps = self.__v_processed_mps[n - 1]
                n_invalid += 1
            else:
                n_invalid = 0
                
    def interpolate_none(self):
        '''This function removes any interpolation from the data and sets filtered data to nan'''
        
        #Reset processed data
        self.__u_processed_mps = self.__u_mps
        self.__v_processed_mps = self.__v_mps
        self.__u_processed_mps[self.__valid_data[0,:] == False] = np.nan
        self.__v_processed_mps[self.__valid_data[0,:] == False] = np.nan  
        
    def interpolate_hold_last(self):
        '''This function holds the last valid value until the next valid data point'''
        
        #Initialize variables
        n_ensembles = self.__u_mps.shape[1]
       
        #Get data from object
        self.__u_processed_mps = self.__u_mps
        self.__v_processed_mps = self.__v_mps
        self.__u_processed_mps[self.__valid_data[0,:] == False] = np.nan
        self.__v_processed_mps[self.__valid_data[0,:] == False] = np.nan
        
        n_invalid = 0
        #Process data by ensembles
        for n in range(1,n_ensembles):
            #Check if ensemble is invalid and number of consecutive invalids is less than 9
            if self.__valid_data[0,n] == False and n_invalid < 9:
                self.__u_processed_mps = self.__u_processed_mps[n - 1]
                self.__v_processed_mps = self.__v_processed_mps[n - 1]
               
               
    def interpolate_next(self):
        '''This function uses the next valid data to back fill for invalid'''
        
        #Get valid ensembles
        valid_ens = self.__valid_data[0,:]
        
        #Process ensembles
        n_ens = len(valid_ens)
        
        for n in np.arange(0,n_ens-1)[::-1]:
            if valid_ens[n] == False:
                self.__u_processed_mps[n] = self.__u_processed_mps[n+1]
                self.__v_processed_mps[n] = self.__v_processed_mps[n+1]
                
                
    def interpolate_smooth(self, transect):
        '''This function interpolates data flagged invalid using the smooth function'''
        
        #Get data from object
        
        u = self.__u_mps
        v = self.__v_mps
        u[self.__valid_data[0,:] == False] = np.nan
        v[self.__valid_data[0,:] == False] = np.nan
        
        #Compute ens_time
        ens_time = np.nancumsum(transect.datetime.ens_duration_sec)
        
        #Apply smooth to each component
        u_smooth = lowess(ens_time, u, 10/len(u))
        v_smooth = lowess(ens_time, v, 10/len(v))
        
        #Save data in object
        self.__u_processed_mps = u
        self.__v_processed_mps = v
        self.__u_processed_mps[np.isnan(u)] = u_smooth[np.isnan(u)]
        self.__v_processed_mps[np.isnan(v)] = v_smooth[np.isnan(v)]
            
            
    def interpolate_linear(self, transect):
        '''This function interpolates data flagged invalid using linear interpolation
        
        Input:
        transect: object of TransectData
        '''  
        
        u = self.__u_mps
        v = self.__v_mps            
        
        valid = np.isnan(u) == False
        
        #Check for valid data
        if sum(valid) > 1 and sum(self.__valid_data[0,:]) > 1:
            
            #Compute ens_time
            ens_time = np.nancumsum(transect.datetime.ens_duration_sec)
            
            #Apply linear interpolation
            self.__u_processed_mps = np.interp(ens_time,
                                               ens_time[self.__valid_data[0,:]],
                                               u[self.__valid_data[0,:]])
            #Apply linear interpolation
            self.__v_processed_mps = np.interp(ens_time,
                                               ens_time[self.__valid_data[0,:]],
                                               v[self.__valid_data[0,:]])
            
         
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
            
            
    def apply_filter(self, transect, kargs):
        '''Function to apply filters to navigation data
        
        Input:
        transect: object of TransectData
        kargs: specified filter method(s) and associated thresholds
        if different from that saved in object.  More than one filter
        can be applied during a single call
        '''
        
        if kargs is not None:
            nargs = len(kargs)
            n=1
            while n < nargs:
                
                #Filter based on number of valid beams
                if kargs[n] == 'Beam':
                    n += 1
                    beam_filter_setting = kargs[n]
                    self.__filter_beam(beam_filter_setting)
                
                #Filter based on difference velocity
                elif kargs[n] == 'Difference':
                    n += 1 
                    d_filter_setting = kargs[n]
                    if d_filter_setting == 'Manual':
                        n+=1
                        self.__filter_diff_vel(d_filter_setting, [kargs[n]])
                    else:
                        self.__filter_diff_vel(d_filter_setting)
                        
                #Filter based on vertical velocity
                elif kargs[n] == 'Vertical':
                    n+=1
                    w_filter_setting = kargs[n]
                    if w_filter_setting == 'Manual':
                        n+=1
                        setting = kargs[n]
                        if setting is not None:
                            setting = self.__w_filter_threshold
                        self.__filter_vert_vel(w_filter_setting, [setting])
                    else:
                        self.__filter_vert_vel(w_filter_setting)
                        
                #Filter based on lowess smooth
                elif kargs[n] == 'Other':
                    n += 1
                    self.__filter_smooth(transect, [kargs[n]])
                    
                n += 1
        else:
            
            #Apply all filters based on stored settings
            self.__filter_beam(self.__beam_filter)
            self.__filter_diff_vel(self.__d_filter, [self.__d_filter_threshold])
            self.__filter_vert_vel(self.__w_filter, [self.__w_filter_threshold])
            self.__filter_smooth(transect, self.__smooth_filter)
            
        #Apply previously specified interpolation method
        self.apply_interpolation(transect)
                    
    #-----------------------------------------------Private Methods--------------
    def __filter_beam(self, setting):
        '''The determination of invalid data depends on the whether 
        3-beam or 4-beam solutions are acceptable. This function can be
        applied by specifying 3 or 4 beam solutions are setting
        obj.beamFilter to -1 which will trigger an automatic mode. The
        automatic mode will find all 3 beam solutions and then compare
        the velocity of the 3 beam solutions to nearest 4 beam solution
        before and after the 3 beam solution. If the 3 beam solution is
        within 50% of the average of the neighboring 3 beam solutions the
        data are deemed valid if not invalid. Thus in automatic mode only
        those data from 3 beam solutions that appear sufficiently
        than the 4 beam solutions are marked invalid. The process happens
        for each ensemble. If the number of beams is specified manually
        it is applied uniformly for the whole transect.
        
        Input:
        setting: setting for beam filter (3, 4, -1)  
        '''
        
        self.__beam_filter = setting
        
        #In manual mode determine number of raw invalid and number of 3 beam solutions  
        #3 beam solutions if selected
        if self.__beam_filter > 0:
            
            #Find invalid raw data
            valid_vel = np.ones(self.__raw_vel_mps.shape)
            valid_vel[np.isnan(self.__raw_vel_mps)] = 0
            
            #Determine how many beams transformed coordinates are valid
            valid_vel_sum = np.sum(valid_vel, 0)
            valid = np.ones(valid_vel_sum.shape)
            
            #Compare number of valid beams or coordinates to filter value
            valid[valid_vel_sum < self.beam_filter] = False
            
            #Save logical of valid data to object
            self.__valid_data[5,:] = valid
            
        else:
            
            #Apply automatic filter
            #----------------------
            #Find all 3 beam solutions
            temp = np.copy(self)
            temp.__filter_beam(4)
            temp3 = np.copy(temp)
            temp3.__filter_beam(3)
            valid_3_beams = temp3.__valid_data[5,:] - temp.__valid_data[5,:]
            n_ens = len(temp.__valid_data[5,:])
            idx = np.where(valid_3_beams == True)[0]
            
            #If 3 beam solutions exist evaluate there validity
            if len(idx) > 0:
                
                #Identify 3 beam solutions that appear to be invalid
                n3_beam_ens = len(idx)
                
                #Check each three beam solution for validity
                
                for m in range(n3_beam_ens):
                    
                    #Check if 3 beam idx is first or last ensemble
                    if idx[m] > 1 and idx[m] < n_ens:
                        
                        #Find nearest 4 beam solutions before and after
                        #3 beam solution
                        ref_idx_before = np.where(temp.__valid_data[5,:idx[m]] == True)[0]
                        if len(ref_idx_before) > 0:
                            ref_idx_before = ref_idx_before[0][-1]
                        else:
                            ref_idx_before = None
                            
                        ref_idx_after = np.where(temp.__valid_data[5, :idx[m]] == True)[0]
                        if len(ref_idx_after) > 0:
                            ref_idx_after = idx[m] + ref_idx_after[0]
                        else:
                            ref_idx_after = None
                            
                        if ref_idx_after is not None and ref_idx_before is not None:
                            u_ratio = (temp.__u_mps[idx[m]]) / ((temp.__u_mps[ref_idx_before] + temp.__u_mps[ref_idx_after]) / 2.) - 1            
                            v_ratio = (temp.__v_mps[idx[m]]) / ((temp.__v_mps[ref_idx_before] + temp.__v_mps[ref_idx_after]) / 2.) - 1
                        else:
                            u_ratio = 1
                            v_ratio = 1
                            
                        #If 3-beam differs from 4-beam by more than 50% mark it invalid
                        if np.abs(u_ratio) > 0.5 or np.abs(v_ratio) > 0.5:
                            temp.__valid_data[5,idx[m]] = 0
                        else:
                            temp.__valid_data[5,idx[m]] = 1  
            
            self.__beam_filter = -1
        
        #Combine all filter data to composite valid data
        self.__valid_data[0,:] = np.all(self.__valid_data[1:,:])
        self.__num_invalid = np.sum(self.__valid_data, 0)
        
    def __filter_diff_vel(self, setting, kargs = None):
        '''Applies either manual or automatic filtering of the difference
        (error) velocity. The automatic mode is based on the following:
        This filter is based on the assumption that the water error velocity
        should follow a gaussian distribution. Therefore, 5 iqr
        should encompass all of the valid data. The standard deviation and
        limits (multiplier*standard deviation) are computed in an iterative 
        process until filtering out additional data does not change the computed 
        standard deviation. 
        
        Input:
        setting: difference velocity setting (Off, Manual, Auto)
        kargs: if manual, the user specified threshold
        '''
        
        self.__d_filter = setting
        if kargs is not None:
            self.__d_filter_threshold = kargs[0]
            
        #Set multiplier
        multiplier = 5
        minimum_window = 0.01
        
        d_vel = np.copy(self.__d_mps)
        
        #apply selected method
        if self.__d_filter == 'Manual':
            d_vel_max_ref = np.abs(self.__d_filter_threshold)
            d_vel_min_ref = -1 * d_vel_max_ref
        elif self.__d_filter == 'Off':
            d_vel_max_ref = np.nanmax(d_vel) + 99
            d_vel_min_ref = np.nanmin(d_vel) - 99
        elif self.__d_filter == 'Auto':
            #Initialize variables
            d_vel_filtered = np.copy(d_vel)
            
            #Fix to zeros in Sontek M9 data
            d_vel_filtered[d_vel_filtered == 0] = np.nan
            
            #Initialize variables
            std_diff = np.repeat(1, 1000)
            k = 0
            
            #Loop until no addictional data are removed
            while std_diff[k] != 0 and k < 1000 and np.isnan(std_diff[k]) == False:
                k += 1
                
                #Compute standard deviation
                d_vel_std = iqr(d_vel_filtered)
                threshold_window = multiplier * d_vel_std
                if threshold_window < minimum_window:
                    threshold_window = minimum_window
                    
                #Compute maximum and minimum thresholds
                d_vel_max_ref = np.nanmedian(d_vel_filtered) + threshold_window
                d_vel_min_ref = np.nanmedian(d_vel_filtered) - threshold_window
                
                #Identify valid and invalid data
                d_vel_bad_idx = np.where((d_vel_filtered > d_vel_max_ref) 
                                         or (d_vel_filtered < d_vel_min_ref))
                d_vel_good_idx = np.where((d_vel_filtered <= d_vel_max_ref)
                                          and (d_vel_filtered >= d_vel_min_ref))
                
                #Determine differences due to last filter iteration
                d_vel_std2 = iqr(d_vel_filtered)
                std_diff[k] = d_vel_std2 - d_vel_std
        
        #Set valid data row 3 for difference velocity filter results
        self.__valid_data[2,:] = False
        self.__valid_data[3, (d_vel <= d_vel_max_ref) and (d_vel >= d_vel_min_ref)] = True
        self.__valid_data[3, self.__valid_data[2,:] == False] = True
        self.__valid_data[3, np.isnan(self.__d_mps)] = True
        self.__d_filter_threshold = d_vel_max_ref
        
        #Combine all filter data to composite filter data
        self.__valid_data[0,:] = np.all(self.__valid_data[1:,:])
        self.__num_invalid = np.sum(self.__valid_data == False, 0)
    
    def __filter_vert_vel(self, setting, kargs=None):
        '''Applies either manual or automatic filtering of the vertical
        velocity.  Uses same assumptions as difference filter.
        
        Input:
        setting: filter setting (Off, Manual, Auto)
        kargs: if setting is manual, the user specified threshold
        '''
        
        #Set vertical velocity filter properties
        self.__w_filter = setting
        if kargs is not None:
            self.__w_filter_threshold = kargs[0]
            
        #Set multiplier
        multiplier = 5
        
        #Get vertical velocity data from object
        w_vel = np.copy(self.__w_mps)
        
        #Apply selected method
        if self.__w_filter == 'Manual':
            w_vel_max_ref = np.abs(self.__w_filter_threshold)
            w_vel_min_ref = -1 * w_vel_max_ref
            
        elif self.__w_filter == 'Off':
            w_vel_max_ref = np.nanmax(w_vel) + 1
            w_vel_min_ref = np.nanmin(w_vel) - 1 
              
        elif self.__w_filter == 'Auto':     
            
            #Initialize variables
            w_vel_filtered = np.copy(w_vel)
            std_diff = 1
            i = 0
            
            #Loop until no additional data are removed
            while std_diff != 0 and i < 1000 and np.isnan(std_diff) == False:
                i += 1
                
                #Compute inner quartile range
                w_vel_std = iqr(w_vel_filtered)
                
                #Compute maximum and minimum thresholds
                w_vel_max_ref = np.nanmedian(w_vel_filtered) + multiplier * w_vel_std
                w_vel_min_ref = np.nanmedian(w_vel_filtered) - multiplier * w_vel_std
                
                #Identify valid and invalid data
                w_vel_bad_idx = np.where((w_vel_filtered > w_vel_max_ref)
                                         or (w_vel_filtered < w_vel_min_ref))
                w_vel_good_idx = np.where((w_vel_filtered <= w_vel_max_ref)
                                          and w_vel_filtered >= w_vel_min_ref)
                
                #Determine differences due to last filter iteration
                w_vel_std2 = np.nanstd(w_vel_filtered)
                std_diff = w_vel_std2 - w_vel_std
        
        #Set valid data row 4 for difference velocity filter results     
        self.__valid_data[3,:] = False
        self.__valid_data[3, np.where((w_vel <= w_vel_max_ref) & (w_vel >= w_vel_min_ref))] = True
        self.__valid_data[3, self.__valid_data[1,:] == False] = True
        self.__w_filter_threshold = w_vel_max_ref
        
        #Combine all filter data to composite valid data
        self.__valid_data[0,:] = np.all(self.__valid_data[1:,:])
        self.__num_invalid = np.sum(self.__valid_data[0,:])
        
    def __filter_smooth(self, transect, setting):
        '''This filter employs a running trimmed standard deviation filter to
        identify and mark spikes in the boat speed. First a robust Loess 
        smooth is fitted to the boat speed time series and residuals between
        the raw data and the smoothed line are computed. The trimmed standard
        deviation is computed by selecting the number of residuals specified by
        "halfwidth" before the target point and after the target point, but not
        including the target point. These values are then sorted, and the points
        with the highest and lowest values are removed from the subset, and the 
        standard deviation of the trimmed subset is computed. The filter
        criteria are determined by multiplying the standard deviation by a user
        specified multiplier. This criteria defines a maximum and minimum
        acceptable residual. Data falling outside the criteria are set to nan.
          
        Recommended filter setting are:
        filterWidth=10;
        halfWidth=10;
        multiplier=9;
        
        David S. Mueller, USGS, OSW
        9/8/2005
        
        Input:
        transect: object of clsTransectData
        setting: filter setting (On, Off)
        '''
        
        #Set property
        self.__smooth_filter = setting
        
        #Compute ens_time
        ens_time = np.nancumsum(transect.date_time.ens_duration_sec)
        
        #Determine if smooth filter should be applied
        if self.__smooth_filter == 'On':
            
            #Boat velocity components
            b_vele = np.copy(self.__u_mps)
            b_veln = np.copy(self.__v_mps)
            
            #Set filter parameters
            filter_width = 10
            half_width = 10
            multiplier = 9
            cycles = 3
            
            #Initialize variables
            direct = np.tile([np.nan], b_vele.shape)
            speed = np.tile([np.nan], b_vele.shape)
            speed_smooth = np.tile([np.nan], b_vele.shape)
            speed_red = np.tile([np.nan], b_vele.shape)
            speed_filtered = np.tile([np.nan], b_vele.shape)
            b_vele_filtered = np.copy(b_vele)
            b_veln_filtered = np.copy(b_veln)
            bt_bad = np.tile([np.nan], b_vele.shape)
            
            #Compute speed and direction of boat
            direct, speed = cart2pol(b_vele, b_veln)
            direct = 90 - (direct * 180 / np.pi)
            
            #Compute residuals from a robust Loess smooth
            speed_smooth = lowess(ens_time, speed, filter_width / len(speed))
            speed_res = speed - speed_smooth
            
            #Apply a trimmed standard deviation filter multiple times
            speed_filtered = np.copy(speed)
            
            for i in range(cycles):
                fil_array = run_std_trim(half_width, speed_res.T)
                
                #Compute filter bounds
                upper_limit = speed_smooth + multiplier * fil_array
                lower_limit = speed_smooth + multiplier * fil_array
                
                #Apply filter to residuals
                bt_bad_idx = np.where((speed > upper_limit) or (speed < lower_limit))
                speed_res[bt_bad_idx] = np.nan
                
            #Update valid_data properth
            self.__valid_data[4,:] = True
            self.__valid_data[4, bt_bad_idx] = False
            self.__valid_data[4, self.__valid_data[1,:] == False] = True
            self.__smooth_upper_limit = upper_limit
            self.__smooth_lower_limit = lower_limit
            self.__smooth_speed = speed_smooth
        
        else:
            
            #Not filter applied all data assumed valid
            self.__valid_data[4,:] = True
            self.__smooth_upper_limit = np.nan
            self.__smooth_lower_limit = np.nan
            self.__smooth_speed = np.nan
            
        #Combine all filter data to composite valid data
        self.__valid_data[0,:] = np.all(self.__valid_data[1:,])
        self.__num_invalid = np.sum(self.__valid_data[0,:] == False, 0)
        
    def __filter_diff_qual(self, gps_data, kargs = None):
        '''Filters GPS data based on the minimum acceptable differential correction quality
        
        Input:
        gps_data: object of GPSData
        kargs: new setting for filter.
        '''
        
        #New filter setting if provided
        if kargs is not None:
            self.__gps_diff_qual_filter = kargs[0]
            
        #Reset valid_data property
        self.__valid_data[2,:] = True
        self.__valid_data[5,:] = True
        
        #Determine and apply appropriate filter type
        self.__valid_data[2, np.isnan(gps_data.diff_qual_ens)] = False
        if self.__gps_diff_qual_filter is not None:
            
            if self.__gps_diff_qual_filter == 1: #autonomous
                self.__valid_data[2, gps_data.diff_qual_ens < 1] = False
            elif self.__gps_diff_qual_filter == 2: #differential correction
                self.__valid_data[2, gps_data.diff_qual_ens < 2] = False
            elif self.__gps_diff_qual_filter == 4: #RTK
                self.__valid_data[2, gps_data.diff_qual_ens < 4] = False
                
            #If there is no indication of the quality assume 1 fot vtg
            if self.__nav_ref == 'VTG':
                self.__valid_data[2, np.isnan(gps_data.diff_qual_ens)] = True
                
        #Combine all filter data to composite valid data
        self.__valid_data[0,:] = np.all(self.__valid_data[1:,:])
        self.__num_invalid = np.sum(self.__valid_data[0,:] == False)
        
        
    def filter_altitude(self, gps_data, kargs):
        '''Filter GPS data based on a change in altitude. Assuming the data
        are collected on the river the altitude should not change
        substantially during the transect. Since vertical resolution is
        about 3 x worse that horizontal resolution the automatic filter
        threshold is set to 3 m, which should ensure submeter horizontal 
        accuracy.
        
        Input:
        gps_data: GPSData
        kargs:
        kargs[0]: new setting for filter (Off, Manual, Auto)
        kargs[1]: change threshold
        '''
        
        #New filter settings if provided
        if kargs is not None:
            self.__gps_altitude_filter = kargs[0]
            if len(kargs) > 1:
                self.__gps_altitude_filter_change = kargs[1]
                
        #Set threshold for Auto
        if self.__gps_altitude_filter == 'Auto':
            self.__gps_altitude_filter_change = 3
            
        #Set all data to valid
        self.__valid_data[3,:] = True
        self.__valid_data[5,:] = True
        
        #Manual or Auto is selected, apply filter
        if self.__gps_altitude_filter == 'Off':
            
            #Initialize variables
            num_valid_old = np.sum(self.__valid_data[3,:])
            k = 0
            change = 1
            
            #Loop until no change in the number of valid ensembles
            while k < 100 and change > 0.1:
                
                #Compute mean using valid ensembles
                alt_mean = np.nanmean(gps_data.altitude_ens_m[self.valid_data[1,:]])
                
                #compute difference for each ensemble
                diff = np.abs(gps_data.altitude_ens_m - alt_mean)
                
                #Mark invalid those ensembles with differences greater than the change threshold
                self.__valid_data[3, diff > self.__gps_altitude_filter_change] = False
                k += 1
                num_valid = np.sum(self.__valid_data[3,:])
                change = num_valid_old - num_valid
                
        
        #Combine all filter data to composite valid data
        self.__valid_data[0,:] = np.all(self.__valid_data[1:,:])
        self.__num_invalid = np.sum(self.__valid_data[0,:] == False)
        
    def __filter_HDOP(self, gps_data, kargs = None):
        '''Filter GPS data based on both a maxumym HDOP and a change in HDOP
        over the transect
        
        Input:
        gps_data: GPS_Data
        
        kargs:
        kargs[0]: filter setting (On, off, Auto)
        kargs[1]: maximum threshold
        kargs[2]: change threshold
        '''
        
        if gps_data.hdop_ens is None:
            self.__valid_data[5,:self.__valid_data.shape[1]] = True
        else:
            #New settings if provided
            self.__gps_HDOP_filter = kargs[0]
            if len(kargs) > 1:
                self.__gps_HDOP_filter_max = kargs[1]
                self.__gps_HDOP_filter_change = kargs[2]
                
            #Settings for auto mode
            if self.__gps_HDOP_filter == 'Auto':
                self.__gps_HDOP_filter_change = 3
                self.__gps_HDOP_filter_max = 4
                
            #Set all ensembles to valid
            self.__valid_data[5,:] = True
            
            #Apply filter for manual or auto
            if self.__gps_HDOP_filter == 'Off':
                
                #Initialize variables
                num_valid_old = np.sum(self.__valid_data[5,:])
                k = 0
                change = 1
                
                #Loop until the number of valid ensembles does not change
                while k < 100 and change > 0.1:
                    
                    #Compute mean HDOP for all valid ensembles
                    hdop_mean = np.nanmean(gps_data.hdop_ens[self.__valid_data[5,:]])
                    
                    #Compute the difference in HDOP and the mean for all ensembles
                    diff = np.abs(gps_data.hdop_ens - hdop_mean)
                    
                    #If the change is HDOP or the value of HDOP is greater
                    #than the threshold setting mark the data invalid
                    self.__valid_data[5, diff > self.__gps_HDOP_filter_change] = False
                    
                    k+=1
                    num_valid = np.sum(self.__valid_data[5,:])
                    change = num_valid_old - num_valid
                    num_valid_old = num_valid
                    
            #combine all filter data to composite data
            self.__valid_data[0,:] = np.all(self.__valid_data[1:,:])
            self.__num_invalid = np.sum(self.__valid_data[0,:] == False)
    
                
def run_std_trim(half_width, my_data):
    ''' The routine accepts a column vector as input. "halfWidth" number of data
          points for computing the standard deviation are selected before and
          after the target data point, but not including the target data point.
          Near the ends of the series the number of points before or after are
          reduced. nan in the data are counted as points. The selected subset of
          points are sorted and the points with the highest and lowest values are
          removed from the subset and the standard deviation computed on the
          remaining points in the subset. The process occurs for each point in the
          provided column vector. A column vector with the computed standard
          deviation at each point is returned.
          
        Input:
        half_width: number of ensembles on each side of target ensemble to used
            for compuring trimmed standard deviation
        my_data: data to be processed
        
        Output:
        fill_array: column vector with computed standard
    '''
    #Determine number of points to process
    n_pts = my_data.shape[0]
    if n_pts < 20:
        half_width = np.floor(n_pts/2)
        
    fill_array = []
    #Compute standard deviation for each point
    for i in range(n_pts):
        
        #Sample selection for 1st point
        if i == 0:
            sample = my_data[1:1+half_width]
            
        #Sample selection at end of data set
        elif i+half_width > n_pts:
            sample = np.hstack([my_data[i-half_width:i-1]], my_data[i+1:n_pts])
            
        #Sample selection at beginning of data set
        elif half_width >= i:
            sample = np.hstack([my_data[:i-1], my_data[i+1:i+half_width]])
            
        #Samples selection in body of data set
        else:
            sample = np.hstack([my_data[i-half_width:i-1], my_data[i+1:i+half_width]])
            
        #Sort and ompute trummed standard deviation
        sample = np.sort(sample)
        fill_array.append(np.nanstd(sample[1:sample.shape[0] - 1]))
        
    return np.array(fill_array)
            
            
          
            
            
    
        
            
        
            
        
            
        
                
        