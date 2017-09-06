'''
Created on Sep 6, 2017

@author: gpetrochenkov
'''
class GPSData(object):
    '''Class containing the raw GPS data and algorithms to convert 
    that raw data to boat velocity'''
    
    def __init__(self):
        
        #raw_ properties
        raw_GGA_lat_deg = None # raw_ latitude in degress, [n,ensemble]
        raw_GGA_lon_deg = None # raw_ longitude in degrees, [n,ensemble]
        raw_GGAA_altitude_m = None # raw_ altitude in meters, [n,ensemble]
        raw_GGA_differential = None # Differential correction indicatore, [n,ensemble]
        raw_GGA_hdop = None # Horizontal dilution of precision, [n,ensemble]
        raw_GGA_utc = None # UTC time, hhmmss.ss, [n, ensemble]
        raw_GGA_serial_time = None # UTC time of gga data in seconds past midnight, [n,ensemble]
        raw_GGA_num_sats = None # Number of satellites reported in GGA sentence, [n,ensemble]
        raw_VTG_course_deg = None # Course in degress, [n, ensemble]
        raw_VTG_speed_mps = None # Speed in m/s, [n, ensemble]
        raw_VTG_delta_time = None # VTG delta time (sec)
        raw_VTG_mode_indicator = None # VTG mode indicator
        raw_GGA_delta_time = None # GGA delta time (sec)
        
        # Manufacturer assigned ensemble values
        ext_GGA_lat_deg = None # Raw latitude in degress, [1,ensemble]
        ext_GGA_lon_deg = None # Raw longitude in degrees, [1,ensemble]
        ext_GGA_altitude_m = None # Raw altitude in meters, [1,ensemble]
        ext_GGA_differential = None # Differential correction indicatore, [1,ensemble]
        ext_GGA_hdop = None # Horizontal dilution of precision, [1,ensemble]
        ext_GGA_utc = None # UTC time, hhmmss.ss, [1, ensemble]
        ext_GGA_serial_time = None # UTC time of gga data in seconds past midnight, [1,ensemble]
        ext_GGA_num_sats = None # Number of satellites reported by software [1,ensemble]
        ext_VTG_course_deg = None # Course in degress, [1, ensemble]
        ext_VTG_speed_mps = None # Speed in m/s, [1, ensemble]
       
        # User specification
        gga_position_method = None # Method used to process gga data for position ('End', 'Average' 'External')
        gga_velocity_method = None # Method used to process gga data for velocity ('End','Average' 'External')
        vtg_velocity_method = None # Method used to process vtg data for velocity ('Average' 'External)
        
        # Computed properties for ensembles
        gga_lat_ens_deg = None # Processed latitude in degrees, [1,ensemble]
        gga_lon_ens_deg = None # Processed longitude in degrees, [1,ensemble]
        utm_ens_m = None # UTM position from processed gga data, [2,ensemble]
        gga_velocity_ens_mps = None # Boat velocity computed from gga data [2,ensemble]
        gga_serial_time_ens = None # UTC time of gga data in seconds past midnight, [1,ensemble]
        vtg_velocity_ens_mps = None # Boat velocity computed from vtg data [2,ensemble]
        per_good_ens = None # Percentage of available data used to compute ensemble value
        hdop_ens = None # HDOP for each ensemble using velocity method
        num_sats_ens = None # Number of satellites for each ensemble, using velocity method
        altitude_ens_m = None # Altitude for each ensemble, using velocity method
        diff_qual_ens = None # Differential quality for each ensemble, using velocity method
        
    def populate_data(self, raw_GGA_utc,raw_GGA_lat,raw_GGA_lon,raw_GGA_alt,raw_GGA_diff,
                raw_GGA_hdop,raw_GGA_num_sats,raw_GGA_delta_time,raw_VTG_course,raw_VTG_speed,raw_VTG_delta_time,
                raw_VTG_mode_indicator,ext_GGA_utc,ext_GGA_lat,ext_GGA_lon,ext_GGA_alt,ext_GGA_diff,
                ext_GGA_hdop,ext_GGA_num_sats,ext_VTG_course,ext_VTG_speed,
                GGA_p_method,GGA_v_method,VTG_method):
        
        #raw_ properties
        raw_GGA_lat_deg = None # raw_ latitude in degress, [n,ensemble]
        raw_GGA_lon_deg = None # raw_ longitude in degrees, [n,ensemble]
        raw_GGAA_altitude_m = None # raw_ altitude in meters, [n,ensemble]
        raw_GGA_differential = None # Differential correction indicatore, [n,ensemble]
        raw_GGA_hdop = None # Horizontal dilution of precision, [n,ensemble]
        raw_GGA_utc = None # UTC time, hhmmss.ss, [n, ensemble]
        raw_GGA_serial_time = None # UTC time of gga data in seconds past midnight, [n,ensemble]
        raw_GGA_num_sats = None # Number of satellites reported in GGA sentence, [n,ensemble]
        raw_VTG_course_deg = None # Course in degress, [n, ensemble]
        raw_VTG_speed_mps = None # Speed in m/s, [n, ensemble]
        raw_VTG_delta_time = None # VTG delta time (sec)
        raw_VTG_mode_indicator = None # VTG mode indicator
        raw_GGA_delta_time = None # GGA delta time (sec)
        
        # Manufacturer assigned ensemble values
        ext_GGA_lat_deg = None # Raw latitude in degress, [1,ensemble]
        ext_GGA_lon_deg = None # Raw longitude in degrees, [1,ensemble]
        ext_GGA_altitude_m = None # Raw altitude in meters, [1,ensemble]
        ext_GGA_differential = None # Differential correction indicatore, [1,ensemble]
        ext_GGA_hdop = None # Horizontal dilution of precision, [1,ensemble]
        ext_GGA_utc = None # UTC time, hhmmss.ss, [1, ensemble]
        ext_GGA_serial_time = None # UTC time of gga data in seconds past midnight, [1,ensemble]
        ext_GGA_num_sats = None # Number of satellites reported by software [1,ensemble]
        ext_VTG_course_deg = None # Course in degress, [1, ensemble]
        ext_VTG_speed_mps = None # Speed in m/s, [1, ensemble]
       
        # User specification
        gga_position_method = None # Method used to process gga data for position ('End', 'Average' 'External')
        gga_velocity_method = None # Method used to process gga data for velocity ('End','Average' 'External')
        vtg_velocity_method = None # Method used to process vtg data for velocity ('Average' 'External)
        
        # Computed properties for ensembles
        gga_lat_ens_deg = None # Processed latitude in degrees, [1,ensemble]
        gga_lon_ens_deg = None # Processed longitude in degrees, [1,ensemble]
        utm_ens_m = None # UTM position from processed gga data, [2,ensemble]
        gga_velocity_ens_mps = None # Boat velocity computed from gga data [2,ensemble]
        gga_serial_time_ens = None # UTC time of gga data in seconds past midnight, [1,ensemble]
        vtg_velocity_ens_mps = None # Boat velocity computed from vtg data [2,ensemble]
        per_good_ens = None # Percentage of available data used to compute ensemble value
        hdop_ens = None # HDOP for each ensemble using velocity method
        num_sats_ens = None # Number of satellites for each ensemble, using velocity method
        altitude_ens_m = None # Altitude for each ensemble, using velocity method
        diff_qual_ens = None # Differential quality for each ensemble, using velocity method
        
        