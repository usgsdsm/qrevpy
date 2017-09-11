'''
Created on Sep 6, 2017

@author: gpetrochenkov
'''

import numpy as np

class GPSData(object):
    '''Class containing the raw GPS data and algorithms to convert 
    that raw data to boat velocity'''
    
    def __init__(self):
        
        #raw_ properties
        self.raw_GGA_lat_deg = None # self.raw_ latitude in degress, [n,ensemble]
        self.raw_GGA_lon_deg = None # self.raw_ longitude in degrees, [n,ensemble]
        self.raw_GGAA_altitude_m = None # self.raw_ altitude in meters, [n,ensemble]
        self.raw_GGA_differential = None # Differential correction indicatore, [n,ensemble]
        self.raw_GGA_hdop = None # Horizontal dilution of precision, [n,ensemble]
        self.raw_GGA_utc = None # UTC time, hhmmss.ss, [n, ensemble]
        self.raw_GGA_serial_time = None # UTC time of gga data in seconds past midnight, [n,ensemble]
        self.raw_GGA_num_sats = None # Number of satellites reported in GGA sentence, [n,ensemble]
        self.raw_VTG_course_deg = None # Course in degress, [n, ensemble]
        self.raw_VTG_speed_mps = None # Speed in m/s, [n, ensemble]
        self.raw_VTG_delta_time = None # VTG delta time (sec)
        self.raw_VTG_mode_indicator = None # VTG mode indicator
        self.raw_GGA_delta_time = None # GGA delta time (sec)
        
        # Manufacturer assigned ensemble values
        self.ext_GGA_lat_deg = None # Raw latitude in degress, [1,ensemble]
        self.ext_GGA_lon_deg = None # Raw longitude in degrees, [1,ensemble]
        self.ext_GGA_altitude_m = None # Raw altitude in meters, [1,ensemble]
        self.ext_GGA_differential = None # Differential correction indicatore, [1,ensemble]
        self.ext_GGA_hdop = None # Horizontal dilution of precision, [1,ensemble]
        self.ext_GGA_utc = None # UTC time, hhmmss.ss, [1, ensemble]
        self.ext_GGA_serial_time = None # UTC time of gga data in seconds past midnight, [1,ensemble]
        self.ext_GGA_num_sats = None # Number of satellites reported by software [1,ensemble]
        self.ext_VTG_course_deg = None # Course in degress, [1, ensemble]
        self.ext_VTG_speed_mps = None # Speed in m/s, [1, ensemble]
       
        # User specification
        self.gga_position_method = None # Method used to process gga data for position ('End', 'Average' 'External')
        self.gga_velocity_method = None # Method used to process gga data for velocity ('End','Average' 'External')
        self.vtg_velocity_method = None # Method used to process vtg data for velocity ('Average' 'External)
        
        # Computed properties for ensembles
        self.gga_lat_ens_deg = None # Processed latitude in degrees, [1,ensemble]
        self.gga_lon_ens_deg = None # Processed longitude in degrees, [1,ensemble]
        self.utm_ens_m = None # UTM position from processed gga data, [2,ensemble]
        self.gga_velocity_ens_mps = None # Boat velocity computed from gga data [2,ensemble]
        self.gga_serial_time_ens = None # UTC time of gga data in seconds past midnight, [1,ensemble]
        self.vtg_velocity_ens_mps = None # Boat velocity computed from vtg data [2,ensemble]
        self.per_good_ens = None # Percentage of available data used to compute ensemble value
        self.hdop_ens = None # HDOP for each ensemble using velocity method
        self.num_sats_ens = None # Number of satellites for each ensemble, using velocity method
        self.altitude_ens_m = None # Altitude for each ensemble, using velocity method
        self.diff_qual_ens = None # Differential quality for each ensemble, using velocity method
        
    def populate_data(self, raw_GGA_utc,raw_GGA_lat,raw_GGA_lon,raw_GGA_alt,raw_GGA_diff,
                raw_GGA_hdop,raw_GGA_num_sats,raw_GGA_delta_time,raw_VTG_course,raw_VTG_speed,raw_VTG_delta_time,
                raw_VTG_mode_indicator,ext_GGA_utc,ext_GGA_lat,ext_GGA_lon,ext_GGA_alt,ext_GGA_diff,
                ext_GGA_hdop,ext_GGA_num_sats,ext_VTG_course,ext_VTG_speed,
                GGA_p_method,GGA_v_method,VTG_method):
        
        #assign input to raw properties
        if raw_GGA_utc is None:
            self.raw_GGA_utc = np.tile([np.nan], raw_GGA_lat.shape)
            self.raw_GGA_serial_time = np.tile([np.nan], raw_GGA_lat.shape)
        else:
            self.raw_GGA_utc = raw_GGA_utc
            self.raw_GGA_serial_time = np.floor(raw_GGA_utc / 10000) * 3600 + \
             np.floor(np.mod(raw_GGA_utc,10000) / 100) * 60 + np.mod(raw_GGA_utc,100)
            
       
        self.raw_GGA_lat_deg = raw_GGA_lat
        self.raw_GGA_lon_deg = raw_GGA_lon
        self.raw_GGA_lat_deg[self.raw_GGA_lat_deg == 0 & self.raw_GGA_lon_deg == 0] = np.nan
        self.raw_GGA_lat_deg[raw_GGA_diff < 1] = np.nan
        self.raw_GGA_lon_deg[np.isnan(self.raw_GGA_lat_deg)] = np.nan
        self.raw_GGA_altitude_m = raw_GGA_alt
        self.raw_GGA_altitude_m[np.isnan(self.raw_GGA_lat_deg)] = np.nan
        self.raw_GGA_differential = raw_GGA_diff
        self.raw_GGA_differential[np.isnan(self.raw_GGA_lat_deg)] = np.nan
        self.raw_GGA_hdop = raw_GGA_hdop
        self.raw_GGA_hdop[np.isnan(self.raw_GGA_lat_deg)] = np.nan
        self.raw_GGA_num_sats = raw_GGA_num_sats
        self.raw_GGA_num_sats[np.isnan(self.raw_GGA_lat_deg)] = np.nan
        self.raw_GGA_serial_time[np.isnan(self.raw_GGA_lat_deg)] = np.nan
        self.raw_GGA_delta_time = raw_GGA_delta_time
        self.raw_VTG_course_deg = raw_VTG_course
        self.raw_VTG_speed_mps = raw_VTG_speed
        self.raw_VTG_course[self.raw_VTG_course_deg == 0 & self.raw_VTG_speed_mps == 0] = np.nan
        self.raw_VTG_speed_mps[np.isnan(self.raw_VTG_course_deg)] = np.nan
        self.raw_VTG_delta_time = raw_VTG_delta_time
        self.raw_VTG_mode_indicator = raw_VTG_mode_indicator
        
        #Assign input data to ensemble values computed by other software
        self.ext_GGA_utc = ext_GGA_utc
        self.ext_GGA_lat_deg = ext_GGA_lat
        self.ext_GGA_lon_deg = ext_GGA_lon
        self.ext_GGA_altitude_m = ext_GGA_alt
        self.ext_GGA_differential = ext_GGA_diff
        self.ext_GGA_hdop = ext_GGA_hdop
        self.ext_GGA_num_sats = ext_GGA_num_sats
        self.ext_GGA_serial_time = np.floor(ext_GGA_utc / 10000) * 3600 + \
            np.floor(np.mod(ext_GGA_utc,10000) / 100) * 60 + np.mod(ext_GGA_utc, 100)
        self.ext_VTG_course_deg = ext_VTG_course
        self.ext_VTG_speed = ext_VTG_speed
        
        #Assign input data to method properties
        self.gga_position_method = GGA_p_method
        self.gga_velocity_method = GGA_v_method
        self.vtg_velocity_method = VTG_method
        
        #If GGA data exist compute position and velocity
        if np.sum(np.sum(np.isnan(raw_GGA_lat) == False)) > 0:
            self.process_GGA()
        
        #If VTG data exist compute velocity
        if np.sum(np.sum(np.isnan(raw_VTG_speed) == False)) > 0:
            self.process_VTG()
        
    def process_GGA(self, kargs = None):
        '''Computes boat velocity from GGA DATA.  Kargs used to specify the methods
        for computing the position and velocity'''
        
        if kargs is not None:
            p_setting = kargs[0]
            v_setting = kargs[1]
        else:
            p_setting = self.gga_position_method
            v_setting = self.gga_velocity_method
            
        #Use only valid GGA data
        valid = self.raw_GGA_num_sats
        valid[np.isnan(valid)] = 0
        valid[valid > 0] = 1
        GGA_lat_deg = self.raw_GGA_lat_deg
        GGA_lat_deg[valid == False] = np.nan
        GGA_lon_deg = self.raw_GGA_lon_deg
        GGA_lon_deg[valid == False] = np.nan
        GGA_serial_time = self.raw_GGA_serial_time
        GGA_serial_time[valid == False] = np.nan
        GGA_delta_time = self.raw_GGA_delta_time
        GGA_delta_time[valid == False] = np.nan
        GGA_hdop = self.raw_GGA_hdop
        GGA_hdop[valid == False] = np.nan
        GGA_altitude = self.raw_GGA_altitude_m
        GGA_altitude[valid == False] = np.nan
        GGA_differential = self.raw_GGA_differential
        GGA_differential[valid == False] = np.nan
        
        #Apply method for computing position of ensemble
        
        #Use ensemble data from other software
        if p_setting == 'External':
            self.gga_lat_ens_deg = self.ext_GGA_lat_deg
            self.gga_lon_ens_deg = self.ext_GGA_lon_deg
        #Uses last valid data for each ensemble
        elif p_setting == 'End':
            n_ensembles = GGA_lat_deg.shape[0]
            for n in range(n_ensembles):
                idx = np.where(np.isnan(GGA_lat_deg[n,:]))[0][-1]
                if len(idx) < 1:
                    idx = 0
                    
                self.gga_lat_ens_deg[n] = GGA_lat_deg[n,idx]
                self.gga_lon_ens_deg[n] = GGA_lat_deg[n,idx]
        elif p_setting == 'First':
            n_ensembles = GGA_lat_deg.shape[0]
            for n in range(n_ensembles):
                idx = 0
                self.gga_lat_ens_deg[n] = GGA_lat_deg[n,idx]
                self.gga_lon_ens_deg[n] = GGA_lat_deg[n,idx]
        elif p_setting == 'Mindt':
            d_time = np.abs(self.raw_GGA_delta_time)
            d_time[valid == False] = np.nan
            d_time_min = np.nanmin(d_time.T).T
            
            use = []
            for n in range(len(d_time_min)):
                use.append(np.abs(d_time[n,:]) == d_time_min[n])
            self.gga_lat_ens_deg = np.tile([np.nan], (1, len(d_time_min))) 
            self.gga_lon_ens_deg = np.tile([np.nan], (1, len(d_time_min)))
            for n in range(len(d_time_min)):
                idx = np.where(use[n,:] == True)[0][0]
                if len(idx) > 0:
                    self.gga_lat_ens_deg[n] = GGA_lat_deg[n, idx]
                    self.gga_lon_ens_deg[n] = GGA_lon_deg[n, idx]
                    
        y_UTM, x_UTM = self.compute_UTM(self.gga_lat_ens_deg, self.gga_lon_ens_deg)
        self.utm_ens_m = (x_UTM, y_UTM)
        
         
            
            
        
        