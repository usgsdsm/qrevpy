'''
Created on Sep 6, 2017

@author: gpetrochenkov
'''

import numpy as np
import utm
from MiscLibs.convenience import azdeg2rad, pol2cart

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
        self.raw_GGA_lat_deg[np.where((self.raw_GGA_lat_deg == 0) & (self.raw_GGA_lon_deg == 0))] = np.nan
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
        self.raw_VTG_course_deg[np.where((self.raw_VTG_course_deg == 0) & (self.raw_VTG_speed_mps == 0))] = np.nan
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
        self.ext_GGA_serial_time = np.floor(np.array(ext_GGA_utc) / 10000) * 3600 + \
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
        GGA_num_sats = self.raw_GGA_num_sats
        GGA_num_sats[valid == False] = np.nan
        GGA_altitude_m = self.raw_GGA_altitude_m
        GGA_altitude_m[valid == False] = np.nan
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
            d_time_min = np.nanmin(d_time.T, 0).T
            
            use = []
            for n in range(len(d_time_min)):
                use.append(np.abs(d_time[n,:]) == d_time_min[n])
                
            use = np.array(use)
            self.gga_lat_ens_deg = np.tile([np.nan], (len(d_time_min))) 
            self.gga_lon_ens_deg = np.tile([np.nan], (len(d_time_min)))
            for n in range(len(d_time_min)):
                idx = np.where(use[n,:] == True)[0]
                if len(idx) > 0:
                    idx = idx[0]
                    self.gga_lat_ens_deg[n] = GGA_lat_deg[n, idx]
                    self.gga_lon_ens_deg[n] = GGA_lon_deg[n, idx]
                    
        y_UTM, x_UTM = self.compute_UTM(self.gga_lat_ens_deg, self.gga_lon_ens_deg)
        self.utm_ens_m = (x_UTM, y_UTM)
        
        n_ensembles = self.raw_GGA_lat_deg.shape[0]
        lat = np.tile([np.nan], (n_ensembles))
        lon = np.tile([np.nan], (n_ensembles))
        self.gga_serial_time_ens = np.tile([np.nan], (n_ensembles))
        self.altitude_ens_m = np.tile([np.nan], (n_ensembles))
        self.diff_qual_ens = np.tile([np.nan], (n_ensembles))
        self.hdop_ens = np.tile([np.nan], (n_ensembles))
        self.num_sats_ens = np.tile([np.nan], (n_ensembles))
        
        #Apply method for computing velocit of ensemble
        if v_setting == 'External':
            lat = self.ext_GGA_lat_deg
            lon = self.ext_GGA_lon_deg
            self.gga_serial_time_ens = self.ext_GGA_serial_time
            self.hdop_ens = self.ext_GGA_hdop
            self.num_sats_ens = self.ext_GGA_num_sats
            self.altitude_ens_m = self.ext_GGA_altitude_m
            self.diff_qual_ens = self.ext_GGA_differential
            
        #Average all position during an ensemble
        elif v_setting == 'Average':
            lat = np.nanmean(GGA_lat_deg,1)
            lon = np.nanmean(GGA_lon_deg,1)
            self.gga_serial_time_ens = np.nanmean(GGA_serial_time,0)
            self.hdop_ens = np.nanmean(GGA_hdop,0)
            self.num_sats_ens = np.floor(np.nanmean(GGA_num_sats,0))
            self.altitude_ens_m = np.nanmean(self.raw_GGA_altitude_m,0)
            self.diff_qual_ens = np.floor(np.nanmean(self.raw_GGA_differential,0))
            
        #Use the last valid data in an ensemble
        elif v_setting == 'End':
            for n in range(n_ensembles):
                idx = np.where(np.isnan(GGA_lat_deg) == False)[0]
                if len(idx) > 0:
                    idx = idx[-1]
                    lat[n] = GGA_lat_deg[n,idx]
                    lon[n] = GGA_lon_deg[n,idx]
                    self.gga_serial_time_ens[n] = GGA_serial_time[n,idx]
                    self.altitude_ens_m[n] = GGA_altitude_m[n,idx]
                    self.diff_qual_ens[n] = GGA_differential[n,idx]
                    
                if idx < len(self.raw_GGA_hdop):
                    self.hdop_ens[n] = GGA_hdop   
                    
                if idx < len(GGA_num_sats[n]):
                    self.num_sats_ens[n] = GGA_num_sats[n, idx] 
                    
        elif v_setting == 'First':
            for n in range(n_ensembles):
                idx = 0
                
                lat[n] = GGA_lat_deg[n,idx]
                lon[n] = GGA_lon_deg[n,idx]
                self.gga_serial_time_ens[n] = GGA_serial_time[n,idx]
                self.altitude_ens_m[n] = GGA_altitude_m[n,idx]
                self.diff_qual_ens[n] = GGA_differential[n,idx]
                
                if idx < len(self.raw_GGA_hdop):
                    self.hdop_ens[n] = GGA_hdop
                    
                if idx < len(GGA_num_sats[n]):
                    self.num_sats_ens[n] = GGA_num_sats[n, idx] 
                    
        elif v_setting == 'Mindt':
            d_time = np.abs(GGA_delta_time)
            d_time_min = np.nanmin(d_time.T, 0).T
            
            use = []
            for n in range(len(d_time_min)):
                use.append(np.abs(d_time[n,:]) == d_time_min[n])
              
            use = np.array(use)  
            for n in range(len(d_time_min)):
                idx = np.where(use[n,:] == True)[0]
                if len(idx) > 0:
                    idx = idx[0]
                    lat[n] = GGA_lat_deg[n, idx]
                    lon[n] = GGA_lon_deg[n, idx]
                    self.gga_serial_time_ens[n] = GGA_serial_time[n,idx]
                    self.altitude_ens_m[n] = GGA_altitude_m[n, idx]
                    self.diff_qual_ens[n] = GGA_differential[n,idx]
                    
                if idx < len(GGA_hdop[n]):
                    self.hdop_ens[n] = GGA_hdop[n,idx]
                    
                if idx < len(GGA_num_sats[n]):
                    self.num_sats_ens = GGA_num_sats
                    
        #Identify valid values
        idx_values = np.where(np.isnan(x_UTM) == False)[0]
        if len(idx_values) > 1:
            u, v = self.gga2_vel_TRDI(lat,lon,self.gga_serial_time_ens,idx_values)
            self.gga_velocity_ens_mps = np.tile([np.nan], (2,len(lat)))
            self.gga_velocity_ens_mps[0,idx_values[1:]] = u[idx_values[1:]]
            self.gga_velocity_ens_mps[1, idx_values[1:]] = v[idx_values[1:]]
        else:
            self.gga_velocity_ens_mps = np.tile([np.nan], (2,len(lat)))
            
        
            
    def process_VTG(self, kargs=None):
        '''Processes raw vtg data to achieve a velocity for each ensemble containing data'''
        
        #Determine method used to compute ensemble velocity
        if kargs is None:
            setting = self.vtg_velocity_method
        else:
            setting = kargs[0]
            
        #Use only valid data
        VTG_speed_mps = self.raw_VTG_speed_mps
        VTG_course_deg = self.raw_VTG_course_deg
        VTG_delta_time = self.raw_VTG_delta_time
        idx = np.where(self.raw_VTG_mode_indicator == 'N')[0]
        VTG_speed_mps[idx] = np.nan
        VTG_course_deg[idx] = np.nan
        VTG_delta_time[idx] = np.nan
        
        if setting == 'Average':
            #Compute vtg elocity in x y coordinates
            direction = azdeg2rad(VTG_course_deg)
            vx, vy = pol2cart(direction,VTG_speed_mps)
            vx[vx == 0 and vy==0] = np.nan
            vy[np.isnan(vx)] = np.nan
            vx_mean = np.nanmean(vx,1)
            vy_mean = np.nanmean(vy,1)
            self.vtg_velocity_ens_mps = np.vstack([vx_mean.T, vy_mean.T])
            
        elif setting == 'End':
            n_ensembles = VTG_speed_mps.shape[0]
            vtg_vel = np.empty(n_ensembles)
            vtg_dir = np.empty(n_ensembles)
            
            for n in range(n_ensembles):
                idx = np.where(np.isnan(VTG_speed_mps))[0]  
                if len(idx) > 0:
                    idx = idx[-1]
                else:
                    idx = 0
                vtg_vel[n] = VTG_speed_mps[n,idx]      
                vtg_dir[n] = VTG_course_deg[n,idx]
                
            direction = azdeg2rad(vtg_dir)
            vx, vy = pol2cart(direction, vtg_vel)
            vx[vx == 0 and vy == 0] = np.nan
            vy[np.isnan(vx)] = np.nan
            self.vtg_velocity_ens_mps = np.vstack[vx,vy]
            
        elif setting == 'First':
            n_ensembles = VTG_speed_mps.shape[0]
            vtg_vel = np.empty(n_ensembles)
            vtg_dir = np.empty(n_ensembles)
            
            for n in range(n_ensembles):
                idx = 0
                vtg_vel[n] = VTG_speed_mps[n,idx]
                vtg_dir[n] = VTG_course_deg[n,idx]
            direction = azdeg2rad(vtg_dir)
            vx, vy = pol2cart(direction,vtg_vel)
            vx[vx == 0 and vy == 0] = np.nan
            vy[np.isnan(vx)] = np.nan
            self.vtg_velocity_ens_mps = np.vstack([vx,vy])
            
        elif setting == 'Mindt':
            d_time = np.abs(VTG_delta_time)
            d_time_min = np.nanmin(d_time.T, 0).T
            
            use = []
            vtg_speed = []
            vtg_dir = []
            
            for n in range(len(d_time_min)):
                use.append(np.abs(d_time[n,:]) == d_time_min[n])
                
            use = np.array(use)
            for n in range(len(d_time_min)):
                idx = np.where(use[n,:] == True)[0]
                if len(idx) > 0:
                    idx = idx[0]
                    vtg_speed.append(VTG_speed_mps[n,idx])
                    vtg_dir.append(VTG_course_deg[n,idx])
                else:
                    vtg_speed.append(np.nan)
                    vtg_dir.append(np.nan)
                    
        elif setting == 'External':
            direction = azdeg2rad(vtg_dir)
            vx,vy = pol2cart(direction)
            self.vtg_velocity_ens_mps = np.vstack([vx.T, vy.T])
                    
        
    def compute_UTM(self, lat_in, lon_in):
        '''Compute UTM coordinates from latitude and longitude'''
        
        lat_in[lat_in == 0] = np.nan
        lon_in[lon_in == 0] = np.nan
        
        #compute UTM coordinates
        lat2 = lat_in * (np.pi/180)
        lon2 = lon_in * (np.pi/180)
        
        y = np.tile([np.nan], lat_in.shape)
        x = np.tile([np.nan], lon_in.shape)
        idx = np.where((np.isnan(lat2) == False) & (np.isnan(lon2) == False))
        for ind in idx[0]:
            y[ind], x[ind], _, _ = utm.from_latlon(lat2[ind],lon2[ind])
        x_UTM = x.reshape(lon_in.shape)
        y_UTM = y.reshape(lat_in.shape)
        
        return (y_UTM, x_UTM)
    
    def gga2_vel_TRDI(self, lat, lon, t, idx_values):
        
        u = np.zeros(lat.shape)
        v = np.zeros(lat.shape)
        
        for n in range(1,len(idx_values)):
            lat1 = lat[idx_values[n-1]]
            lat2 = lat[idx_values[n]]
            lon1 = lon[idx_values[n-1]]
            lon2 = lon[idx_values[n]]
            t1 = t[idx_values[n-1]]
            t2 = t[idx_values[n]]
            L = ((lat1+lat2)/2) * np.pi/180
            s_L = np.sin(L)
            coefficient = 6378137*np.pi/180
            ellipticity = 1/298.257223563
            RE = coefficient * (1 + ellipticity * s_L**2)
            RN = coefficient * (1-2 * ellipticity + 3 * ellipticity * s_L**2)
            delta_x = RE * (lon2 - lon1) * np.cos(L)
            delta_y = RN * (lat2 - lat1)
            delta_time = t2 - t1
            u[idx_values[n]] = delta_x / delta_time
            v[idx_values[n]] = delta_y / delta_time
            
        return (u,v)
        