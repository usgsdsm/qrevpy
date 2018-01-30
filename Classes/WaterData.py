"""
Created on Jul 31, 2017

@author: gpetrochenkov

Modified DSM 1/29/2018
    - Added docstrings
    - Cleaned up PEP8
    - Removed need for kargs in initial methods
    - Eliminated private variables and methods
    - Moved zero velocity check for SonTek data to TransectData class
"""
import numpy as np
from numpy.matlib import repmat
from MiscLibs.convenience import cosd, sind, cart2pol, pol2cart, iqr
from openpyxl.chart import print_settings
from MiscLibs.lowess import lowess
from Classes.BoatData import run_std_trim
from scipy.interpolate import interpolate


class WaterData(object):
    """Class to process and store water velocity data.

    Attributes
    ----------
    Original data provided to the class:
        raw_vel_mps: np.array(float)
            Contains the raw unfiltered velocity data in m/s.  1st index 1-4 are beams 1,2,3,4 if beam or
            u,v,w,d if otherwise.
        frequency_kHz: np.array(float)
            Defines ADCP frequency used for velocity measurement, in kHz.
        orig_coord_sys: str
            Defines the original raw data velocity coordinate system "Beam", "Inst", "Ship", "Earth".
        orig_nav_ref: str
            Defines the original taw data naviagation reference: "None", "BT", "GGA", "VTG".
        corr: np.array(float)
            Correlation values for WT, if available.
        rssi: np.array(float)
            Returned acoustic signal strength.
        rssi_units: str
            Units for returned acoustic signal strength: "Counts" "dB", "SNR".
        water_mode: str
            WaterMode for TRDI or 'Variable' for SonTek.
        blanking_distance_m: float
            Distance below transducer where data is marked invalid dur to potential ringing interference.
        cells_above_sl: np.array(bool)
            Logical array of depth cells above the sidelobe cutoff based on selected depth reference.
        cells_above_sl_bt: np.array(bool)
            Logical array of depth cells above the sidelobe cutoff based on BT
        sl_lag_effect_m: np.array(float)
            Side lobe distance due to lag and transmit length

    Data computed in this class:
        u_earth_no_ref_mps: np.array(float)
            Horizontal velocity in x-direction with no boat reference applied, in m/s.
        v_earth_no_ref_mps: np.array(float)
            Horizontal velocity in y-direction with no boat reference applied, in m/s.
        u_mps: np.array(float)
            Horizontal velocity in x-direction, earth coord, nav referenced, in m/s.
        v_mps: np.array(float)
            Horizontal velocity in y-direction, earth coord, nav referenced, in m/s.
        u_processed_mps: np.array(float)
            Horizontal velocity in x-direction, earth coord, nav referenced, filtered, and interpolated.
        v_processed_mps: np.array(float)
            Horizontal veloctiy in y-direction, earth coord, nav referenced, filtered, and interpolated.
        w_mps: np.array(float)
            Vertical velocity (+ up), in m/s.
        d_mps: np.array(float)
            Difference in vertical velocities compute from opposing beam pairs, in m/s.
        invalid_index: np.array(bool)
            Index of ensembles with no valid raw velocity data.
        num_invalid: float
            Estimated number of depth cells in ensembles with no valid raw velocity data.
        valid_data: np.array(float)
            3-D logical array of valid data
                Dim3 0 - composite
                Dim3 1 - original, cells above side lobe
                Dim3 2 - dfilter
                Dim3 3 - wfilter
                Dim3 4 - smoothFilter
                Dim3 5 - beamFilter
                Dim3 6 - excluded
                Dim3 7 - snrFilter
                Dim3 8 - validDepthFilter

    Processing settings:
        beam_filter: int
            Set 3 for 3-beam solutions, 4 for 4-beam solutions.
        d_filter: str
            Set difference velocity filter "On", "Off".
        d_filter_threshold: float
            Threshold for difference velocity filter.
        w_filter: str
            Set vertical velocity filter "On", "Off".
        w_filter_threshold: float
            Threshold for vertical velocity filter.
        excluded_dist_m: float
            Distance below transucer for which data are excluded or marked invalid, in m.
        smooth_filter: str
            Set filter based on smoothing function "On", "Off".
        smooth_speed: np.array(float)
            Smoothed mean water speed, in m/s.
        smooth_upper_limit: np.array(float)
            Smooth function upper limit of window, in m/s.
        smooth_lower_limit: np.array(float)
            Smooth funciton lower limit of window, in m/s.
        snr_filter: str
            Set SNR filter for SonTek data "On", "Off".
        snr_rng: np.array(float)
            Range of beam averaged SNR
        wt_depth_filter: np.array(bool)
            WT in ensembles with invalid depths are marked invalid.
        interpolate_ens: str
            Type of interpolation: "None", "TRDI", "Linear".
        interpolate_cells: str
            Type of cell interpolation: "None", "TRDI", "Linear"
        coord_sys: str
            Defines the velocity coordinate system "Beam", "Inst", "Ship", "Earth"
        nav_ref: str
            Defines the navigation reference: "None", "BT", "GGA", "VTG"
        sl_cutoff_percent: float
            Percent cutoff defined by cos(angle)
        sl_cutoff_number: float
            User specigied number of cells to cutoff from SonTek, not implemented, undefined
        sl_cutoff_type: str
            Type of cutoff method "Percent" or "Number".
        num_bins_filtered: float
            ???????????
    """
    # TODO check water mode for TRDI is it an array, should it be?
    # TODO filters switch on and off should probably be bool
    # TODO how is number of bins for sl cutoff implemented by SonTek?
    # TODO what is num_bins_filtered?

    def __init__(self):
        # Data input to this class
        self.raw_vel_mps = None         # Contains the raw unfiltered velocity data in m/s.  Rows 1-4 are
                                        # beams 1,2,3,4 if beam or u,v,w,d if otherwise
        self.frequency_kHz = None       # Defines ADCP frequency used for velocity measurement
        self.orig_coord_sys = None      # Defines the original velocity coordinate system Beam, Inst, Ship, Earth
        self.orig_nav_ref = None        # Defines the original taw data naviagation reference: None, BT, GGA, VTG
        self.corr = None                # Correlation values for WT, if available
        self.rssi = None                # Returned acoustic signal strength
        self.rssi_units = None          # Units for returned acoustic signal strength: Counts, dB, SNR
        self.water_mode = None          # WaterMode for TRDI or 'Variable' for SonTek
        self.blanking_distance_m = None # Distance below transducer where data is marked invalid due to potential
                                        # ringing interference
        self.cells_above_sl = None      # Bool array of depth cells above the sidelobe cutoff based
        self.cells_above_sl_bt = None   # Bool array of depth cells above the sidelobe cutoff based on BT
        self.sl_lag_effect_m = None     # Side lobe distance due to lag and transmit length
        
        # Data computed in this class
        self.u_earth_no_ref_mps = None  # Horizontal velocity in x-direction with no boat reference applied, in m/s
        self.v_earth_no_ref_mps = None  # Horizontal velocity in y-direction with no boat reference applied, in m/s
        self.u_mps = None               # Horizontal velocity in x-direction, earth coord, nav referenced, in m/s
        self.v_mps = None               # Horizontl velocity in y-direction, earth coord, nav referneced, in m/s
        self.u_processed_mps = None     # Horizontal velocity in x-direction, earth coord, nav referenced,
                                        # filtered, and interpolated
        self.v_processed_mps = None     # Horizontal veloctiy in y-direction, earth coord, nav referenced,
                                        # filtered, and interpolated
        self.w_mps = None               # Vertical velocity (+ up), in m/s
        self.d_mps = None               # Difference in vertical velocities compute from opposing beam pairs, in m/s
        self.invalid_index = None       # Index of ensembles with no valid raw velocity data
        self.num_invalid = []           # Estimated number of depth cells in ensembles with no valid raw velocity data
        self.valid_data = None          # 3-D logical array of valid data
                                            # Dim3 0 - composite
                                            # Dim3 1 - original, cells above side lobe
                                            # Dim3 2 - dfilter
                                            # Dim3 3 - wfilter
                                            # Dim3 4 - smoothFilter
                                            # Dim3 5 - beamFilter
                                            # Dim3 6 - excluded
                                            # Dim3 7 - snrFilter
                                            # Dim3 8 - validDepthFilter
                                
        # Settings
        self.beam_filter = None         # 3 for 3-beam solutions, 4 for 4-beam solutions
        self.d_filter = None            # Difference velocity filter "On", "Off"
        self.d_filter_threshold = None  # Threshold for difference velocity filter
        self.w_filter = None            # Vertical velocity filter "On", "Off"
        self.w_filter_threshold = None  # Threshold for vertical velocity filter
        self.excluded_dist_m = None     # Distance below transucer for which data are ecluded or marked invalid
        self.smooth_filter = None       # Filter based on smoothing function
        self.smooth_speed = None        # Smoothed boar speed
        self.smooth_upper_limit = None  # Smooth function upper limit of window
        self.smooth_lower_limit = None  # Smooth funciton lower limit of window
        self.snr_filter = None          # SNR filter for SonTek data
        self.snr_rng = None             # Range of beam averaged SNR
        self.wt_depth_filter = None     # WT in ensembles with invalid WT are marked invalid
        self.interpolate_ens = None     # Type of interpolation: "None", "TRDI", "Linear"
        self.interpolate_cells = None   # Type of cell interpolation: "None", "TRDI", "Linear"
        self.coord_sys = None           # Defines the velocity coordinate system "Beam", "Inst", "Ship", "Earth"
        self.nav_ref = None             # Defines the navigation reference: "None", "BT", "GGA", "VTG"
        self.sl_cutoff_percent = None   # Percent cutoff defined by cos(angle)
        self.sl_cutoff_number = None    # User specigied number of cells to cutoff above slcutoff
        self.sl_cutoff_type = None      # "Percent" or "Number"
        self.num_bins_filtered = None
        
    def populate_data(self, vel_in, freq_in, coord_sys_in, nav_ref_in, rssi_in, rssi_units_in,
                      excluded_dist_in, cells_above_sl_in, sl_cutoff_per_in, sl_cutoff_num_in,
                      sl_cutoff_type_in, sl_lag_effect_in, wm_in, blank_in, corr_in=None,
                      surface_vel_in=None, surface_rssi_in=None, surface_corr_in=None, surface_num_cells_in=None):
        
        """Populates the variables with input, computed, or default values.

        Parameters
        ----------
        vel_in: np.array(float)
            Contains the raw unfiltered velocity data in m/s.
            Rows 1-4 are beams 1,2,3,4 if beam or u,v,w,d if otherwise.
        freq_in: np.array(float)
            Defines ADCP frequency used for velocity measurement.
        coord_sys_in: str
            Defines the original raw data velcocity coordinate system "Beam", "Inst", "Ship", "Earth".
        nav_ref_in: str
            Defines the original raw data navigation reference: "None", "BT", "GGA", "VTG".
        rssi_in: np.array(float)
            Returned acoustic signal strength.
        rssi_units_in: str
            Units for returned acoustic signal strength: "Counts", "dB", "SNR".
        excluded_dist_in: float
            Distance below transducer for which data are excluded or marked invalid.
        cells_above_sl_in: np.array(bool)
            Bool array of depth cells above the sidelobe cutoff based on selected depth reference.
        sl_cutoff_per_in: float
            Percent cutoff defined by cos(angle).
        sl_cutoff_num_in: float
            User specified number of cells to cutoff above sl_cutoff.
        sl_cutoff_type_in: str
            Method used to compute cutoff "Percent" or "Number".
        sl_lag_effect_in: np.array(float)
            Lag effect for each ensemble, in m.
        wm_in: str
            Watermode for TRDI or 'Variable' for SonTek.
        blank_in: float
            Blanking distance, in m.
        corr_in: np.array(float)
            Correlation values for water track. Optional.
        surface_vel_in: np.array(float)
            Surface velocity data for RiverRay, RiverPro, RioPro. Optional.
        surface_rssi_in: np.array(float)
            Returned acoust signal strength for RiverRay, RiverPro, RioPro. Optional.
        surface_corr_in: np.array(float)
            Surface velocity correlations for RiverRay, RiverPro, RioPro. Optional.
        surface_num_cells_in: np.array(float)
            Number of surface cells in each ensemble for RiverRay, RiverPro, RioPro. Optional.
        """

        # Set object properties from input data standard for all ADCPs
        self.frequency_kHz = freq_in
        self.orig_coord_sys = coord_sys_in
        self.coord_sys = coord_sys_in
        self.orig_nav_ref = nav_ref_in
        self.nav_ref = nav_ref_in
        self.water_mode = wm_in
        self.excluded_dist_m = excluded_dist_in
        self.rssi_units = rssi_units_in
        # ============================
        # DSM need to check this for correctness
        # Set object properties that depend on the presence or absence of surface cells
        if surface_num_cells_in:
            surface_num_cells_in[np.isnan(surface_num_cells_in)] = 0
            max_cells = cells_above_sl_in.shape[0]
            num_ens = cells_above_sl_in.shape[1]
            num_reg_cells = vel_in.shape[1]
            max_surf_cells = max_cells - num_reg_cells

            # Combine surface velocity bins and regular velocity bins into one matrix
            self.raw_vel_mps = np.tile([np.nan], [4, max_cells, num_ens])
            self.rssi = np.tile([np.nan], [4, max_cells, num_ens])
            self.corr = np.tile([np.nan], [4, max_cells, num_ens])

            if max_surf_cells > 0:
                self.raw_vel_mps[:, :max_surf_cells, :] = surface_vel_in[:, :max_surf_cells, :]
                self.rssi[:, :max_surf_cells, :] = surface_rssi_in[:, :max_surf_cells, :]
                self.corr[:, :max_surf_cells, :] = surface_corr_in[:, :max_surf_cells, :]

            for i_ens in range(num_ens):
                self.raw_vel_mps[:, int(surface_num_cells_in[i_ens]):int(surface_num_cells_in[i_ens])
                    + num_reg_cells, i_ens] = vel_in[:, :num_reg_cells, i_ens]
                self.rssi[:, int(surface_num_cells_in[i_ens]):int(surface_num_cells_in[i_ens])
                    + num_reg_cells, i_ens] = rssi_in[:, :num_reg_cells, i_ens]
                self.corr[:, int(surface_num_cells_in[i_ens]):int(surface_num_cells_in[i_ens])
                    + num_reg_cells, i_ens] = corr_in[:, :num_reg_cells, i_ens]
        else:
            # No surface cells
            self.raw_vel_mps = vel_in
            self.rssi = rssi_in
            if corr_in.any():
                self.corr = corr_in
            else:
                # No correlations input
                self.corr = repmat([np.nan], rssi_in.shape[0], rssi_in.shape[1])
        # DSM changed to avoid kargs 1/26/2018
        # #No correlation data
        # if kargs is None:
        #     self.raw_vel_mps = vel_in
        #     self.rssi = rssi_in
        #     self.rssi_units = rssi_units_in
        #     self.corr = repmat([np.nan], rssi_in.shape[0], rssi_in.shape[1])
        #
        # elif len(kargs) == 1:
        #     self.raw_vel_mps = vel_in
        #     self.rssi = rssi_in
        #     self.rssi_units = rssi_units_in
        #     self.corr = kargs[0]
        #
        # elif len(kargs) > 1:
        #     self.rssi_units = rssi_units_in
        #     corr_in = kargs[0]
        #     surface_vel = kargs[1]
        #     surface_rssi = kargs[2]
        #     surface_corr = kargs[3]
        #     no_surf_cells = kargs[4]
        #     no_surf_cells[np.isnan(no_surf_cells)] = 0
        #     max_cells = cells_above_sl_in.shape[0]
        #     num_ens = cells_above_sl_in.shape[1]
        #     num_reg_cells = vel_in.shape[1]
        #     max_surf_cells = max_cells - num_reg_cells
        #
        #     #Combine surface velocity bins and regular velocity bins into one matrix
        #     self.raw_vel_mps = np.tile([np.nan],[4,max_cells, num_ens])
        #     self.rssi = np.tile([np.nan],[4,max_cells, num_ens])
        #     self.corr = np.tile([np.nan],[4,max_cells, num_ens])
        #
        #     if max_surf_cells > 0:
        #         self.raw_vel_mps[:,:max_surf_cells,:] = surface_vel[:,:max_surf_cells,:]
        #         self.rssi[:,:max_surf_cells] = surface_rssi[:,:max_surf_cells,:]
        #         self.corr[:,:max_surf_cells] = surface_corr[:,:max_surf_cells,:]
        #
        #     for i_cell in range(num_ens):
        #         self.raw_vel_mps[:,int(no_surf_cells[i_cell]):int(no_surf_cells[i_cell]) + num_reg_cells,i_cell] = vel_in[:,:num_reg_cells,i_cell]
        #         self.rssi[:,int(no_surf_cells[i_cell]):int(no_surf_cells[i_cell]) + num_reg_cells,i_cell] = rssi_in[:,:num_reg_cells,i_cell]
        #         self.corr[:,int(no_surf_cells[i_cell]):int(no_surf_cells[i_cell]) + num_reg_cells,i_cell] = corr_in[:,:num_reg_cells,i_cell]
        #
        #
        # #Set object properties from input data
        # self.frequency_kHz = freq_in
        # self.orig_coord_sys = coord_sys_in
        # self.coord_sys = coord_sys_in
        # self.orig_nav_ref = nav_ref_in
        # self.nav_ref = nav_ref_in
        self.u_mps = np.squeeze(np.copy(self.raw_vel_mps)[0, :, :])
        self.v_mps = np.squeeze(np.copy(self.raw_vel_mps)[1, :, :])
        self.w_mps = np.squeeze(np.copy(self.raw_vel_mps)[2, :, :])
        self.d_mps = np.squeeze(np.copy(self.raw_vel_mps)[3, :, :])
        
        # DSM moved to TransectData 1/26/2018
        # Because Matlab pads arrays with zeros and RR data has variable
        # number of bins, the raw data may be padded with zeros.  The next
        # four statements changes those to nan. DEBUG!---------------------
        # if len(kargs) > 1:
        #     self.u_mps[self.u_mps == 0] = np.nan
        #     self.v_mps[self.v_mps == 0] = np.nan
        #     self.w_mps[self.w_mps == 0] = np.nan
        #     self.d_mps[self.d_mps == 0] = np.nan
            
        self.water_mode = wm_in
        self.excluded_dist_m = excluded_dist_in

        # In some rare situations the blank is empty so it is set to the excluded_dist_in
        try:
            blank_in = float(blank_in)
            self.blanking_distance_m = blank_in
        except ValueError:
            self.blanking_distance_m = excluded_dist_in
            
        self.cells_above_sl = cells_above_sl_in
        self.cells_above_sl_bt = cells_above_sl_in
        self.sl_cutoff_percent = sl_cutoff_per_in
        self.sl_cutoff_number = sl_cutoff_num_in
        self.sl_cutoff_type = sl_cutoff_type_in
        self.sl_lag_effect_m = sl_lag_effect_in
        
        # Set filter defaults to no filtering and no interruption
        self.beam_filter = 3
        self.d_filter = 'Off'
        self.d_filter_threshold = 99
        self.w_filter = 'Off'
        self.w_filter_threshold = 99
        self.smooth_filter = 'Off'
        self.interpolate_ens = 'None'
        self.interpolate_cells = 'None'
        
        # Determine original valid

        # Initialize valid data property
        self.valid_data = np.tile(self.cells_above_sl, [9, 1, 1])
        
        # Find invalid raw data
        valid_vel = np.tile(self.cells_above_sl, [4, 1, 1])
        valid_vel[np.isnan(self.raw_vel_mps)] = False
        # DSM dealing with 0 velocity is dealt with in TransectData 1/26/2018
        # if len(kargs) > 1:
        #     valid_vel[self.raw_vel_mps == 0] = False
            
        # Identify invalid velocity data (less than 3 valid beams)
        valid_vel_sum = np.sum(valid_vel, axis=0)
        valid_data2 = np.copy(self.cells_above_sl)
        valid_data2[valid_vel_sum < 3] = False
        
        # Set valid_data property for original data
        self.valid_data[1, :, :] = valid_data2
        
        # Combine all filter data to composite valid data
        self.all_valid_data()
        
        # Estimate the number of cells in invalid ensembles using
        # Adjacent valid ensembles
        valid_data_2_sum = np.nansum(self.valid_data[1], 0)
        self.invalid_index = np.where(valid_data_2_sum == 0)[0]
        n_invalid = len(self.invalid_index)
        for n in range(n_invalid):
            # Find first valid ensemble
            idx1 = np.where(valid_data_2_sum[:self.invalid_index[n]] > 0)[0]
            if len(idx1) > 0:
                idx1 = idx1[0]
            else:
                idx1 = self.invalid_index[n]
                
            # Find next valid ensemble
            idx2 = np.where(valid_data_2_sum[:self.invalid_index[n]] > 0)[0]
            if len(idx2) > 0:
                idx2 = idx2[-1]
            else:
                idx2 = self.invalid_index[n]
                
            # Estimate number of cells in invalid ensemble
            self.num_invalid.append(np.floor((valid_data_2_sum[idx1]+valid_data_2_sum[idx2]) / 2))
            
        # Set processed data to non-interpolated valid data
        self.u_processed_mps = np.copy(self.u_mps)
        self.v_processed_mps = np.copy(self.v_mps)
        self.u_processed_mps[self.valid_data[0] == False] = np.nan
        self.v_processed_mps[self.valid_data[0] == False] = np.nan
        
        # Compute SNR range if SNR data is provided
        if rssi_units_in == 'SNR':
            self.compute_snr_rng()
            
    def change_coord_sys(self, new_coord_sys, sensors, adcp):
        """This function allows the coordinate system to be changed.
        Current implementation is only to allow a change to a higher order
        coordinate system Beam - Inst - Ship - Earth
        
        new_coord_sys - new coordinate system (Beam, Inst, Ship, Earth)
        sensors - obj of Sensors
        adcp - object of instrument data
        """
        
        o_coord_sys = self.orig_coord_sys[0].strip()
        
        if o_coord_sys != new_coord_sys:
            
            #Assign the transformation matrix and retrieve the sensor data
            t_matrix = adcp.t_matrix.matrix
            t_matrix_freq = adcp.frequency_hz
            
            pitch_select = getattr(sensors.pitch_deg, sensors.pitch_deg.selected)
            p =  getattr(pitch_select, '_SensorData__data')
            roll_select = getattr(sensors.roll_deg, sensors.roll_deg.selected)
            r = getattr(roll_select, '_SensorData__data')
            heading_select = getattr(sensors.heading_deg, sensors.heading_deg.selected)
            h = getattr(heading_select, 'data')
            
            #Modify the transformation matrix and heading, pitch
            #and roll values based on the original coordinate
            #system so that only the needed values ar used in
            #computing the new coordinate system.
            if o_coord_sys.strip() == 'Beam':
                orig_sys = 1
            elif o_coord_sys.strip() == 'Inst':
                orig_sys = 2
            elif o_coord_sys.strip() == 'Ship':
                orig_sys = 3
                p = np.zeros(h.shape)
                r = np.zeros(h.shape)
                t_matrix = np.eye(len(t_matrix))
            elif o_coord_sys.strip() == 'Earth':
                new_sys = 4
                
            if new_coord_sys.strip() == 'Beam':
                new_sys = 1
            elif new_coord_sys.strip() == 'Inst':
                new_sys = 2
            elif new_coord_sys.strip() == 'Ship':
                new_sys = 3
            elif new_coord_sys.strip() == 'Earth':
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
                
                vel_changed = np.tile([np.nan], self.raw_vel_mps.shape)
                n_ens = self.raw_vel_mps.shape[2]
                
                for ii in range(n_ens):
                    
                    #Compute matrix for heading, pitch, and roll
                    hpr_matrix = np.array([[((CH[ii] * CR[ii]) + (SH[ii]*SP[ii]*SR[ii])),
                                (SH[ii] * CP[ii]),
                                ((CH[ii] * SR[ii]) - SH[ii]*SP[ii]*CR[ii])],
                                [(-1 * SH[ii] * CR[ii])+(CH[ii] * SP[ii] * SR[ii]),
                                CH[ii] * CP[ii], 
                                (-1 * SH[ii] * SR[ii])-(CH[ii] * SP[ii] * CR[ii])],
                                [(-1.*CP[ii] * SR[ii]),
                                SP[ii],
                                CP[ii] * CR[ii]]])
                    
                    #Transofm beam coordinates
                    if o_coord_sys == 'Beam':
                        
                        #Determine frequency index for transformation
                        if len(t_matrix.shape) > 2 and t_matrix.shape[2] > 1:
                            idx_freq = np.where(t_matrix_freq == self.frequency_kHz[ii])
                            t_mult = np.copy(t_matrix[idx_freq])
                        else:
                            t_mult = np.copy(t_matrix)
                            
                        #Get velocity data
                        vel_beams = np.squeeze(self.raw_vel_mps[:,:,ii])
                        
                        #Apply transformation matrix for 4 beam solutions
                        temp_t = t_mult.dot(vel_beams)
                        
                        #Apply hpr_matrix
                        temp_THPR = hpr_matrix.dot(temp_t[:3])
                        temp_THPR = np.vstack([temp_THPR, temp_t[3]])
                        
                        #Check for invalid beams
                        invalid_idx = np.isnan(vel_beams)
                        
                        #Identify rows requiring 3 beam solutions
                        n_invalid_col = np.sum(invalid_idx, axis=0)
                        col_idx = np.where(n_invalid_col==1)[0]
                        
                        #Compute 3 beam solution, if necessary
                        if len(col_idx) > 0:
                            for i3 in range(len(col_idx)):
                                temp_t = repmat([np.nan], 4, 1)
                                t_mult_3_beam = t_mult
                                
                                #Id invalid beam
                                vel_3_beam = vel_beams[:,col_idx[i3]]
                                idx_3_beam = np.where(np.isnan(vel_3_beam))[0]
                        
                                #3 beam solution for non-RiverRay
                                vel_3_beam_zero = vel_3_beam
                                vel_3_beam_zero[np.isnan(vel_3_beam)] = 0
                                vel_error = t_mult[3,:].dot(vel_3_beam_zero)
                                vel_3_beam[idx_3_beam] = -1 * vel_error / t_mult[3,idx_3_beam]
                                temp_t = t_mult.dot(vel_3_beam)
                                
                                #Apply transformation matrix for 3
                                #beam solutions
                                temp_THPR[0:3,col_idx[i3]] = hpr_matrix.dot(temp_t[:3])
                                temp_THPR[3,col_idx[i3]] = np.nan
                            
                    else:
                        #Get velocity data
                        vel_raw = np.squeeze(self.raw_vel_mps[:,:,ii])
                        temp_THPR = np.array(hpr_matrix).dot(vel_raw[:3,:])
                        temp_THPR = np.vstack([temp_THPR, vel_raw[3,:]])
                        
                    #Update object
                    temp_THPR = temp_THPR.T
                    self.u_mps[:,ii] = temp_THPR[:,0]
                    self.v_mps[:,ii] = temp_THPR[:,1]
                    self.w_mps[:,ii] = temp_THPR[:,2]
                    self.d_mps[:,ii] = temp_THPR[:,3]
                    
                    #end n_ens --------------------------------------------------------
                    
                # Because of padded arrays with zeros and RR has a variable number of bins,
                # the raw data may be padded with zeros.  The next 4 statements changes
                # those to nan
                self.u_mps[self.u_mps == 0] = np.nan
                self.v_mps[self.v_mps == 0] = np.nan
                self.w_mps[self.w_mps == 0] = np.nan
                self.d_mps[self.d_mps == 0] = np.nan
                
                #Assign processed object properties
                self.u_processed_mps = np.copy(self.u_mps)
                self.v_processed_mps = np.copy(self.v_mps)
                
                #Assign coordinate system and reference properties
                self.coord_sys = new_coord_sys
                self.nav_ref = self.orig_nav_ref
                    
            else:
                
                #Reset velocity properties to raw values
                self.u_mps = self.raw_vel_mps[0]
                self.v_mps = self.raw_vel_mps[1]
                self.w_mps = self.raw_vel_mps[2]
                self.d_mps = self.raw_vel_mps[3]
                
                if adcp.manufacturer == 'TRDI':
                    self.u_mps[self.u_mps == 0] = np.nan
                    self.v_mps[self.v_mps == 0] = np.nan
                    self.w_mps[self.w_mps == 0] = np.nan
                    self.d_mps[self.d_mps == 0] = np.nan
                    
                #Assign processed properties
                self.u_processed_mps = np.copy(self.u_mps)
                self.v_processed_mps = np.copy(self.v_mps)
                
        else:
            
            #Reset velocity properties to raw values
            self.u_mps = self.raw_vel_mps[0]
            self.v_mps = self.raw_vel_mps[1]
            self.w_mps = self.raw_vel_mps[2]
            self.d_mps = self.raw_vel_mps[3]
            
            if adcp.manufacturer == 'TRDI':
                self.u_mps[self.u_mps == 0] = np.nan
                self.v_mps[self.v_mps == 0] = np.nan
                self.w_mps[self.w_mps == 0] = np.nan
                self.d_mps[self.d_mps == 0] = np.nan
                
            #Assign processed properties
            self.u_processed_mps = np.copy(self.u_mps)
            self.v_processed_mps = np.copy(self.v_mps)
            
        if new_coord_sys == 'Earth':
            self.u_earth_no_ref_mps = np.copy(self.u_mps)
            self.v_earth_no_ref_mps = np.copy(self.v_mps)
                
    def set_nav_reference(self, boat_vel):           
        """This function sets the navigation reference.  The current
        reference is first removed from the velocity and then the
        selected reference is applied"""
        
        #apply selected navigation reference
        boat_select = getattr(boat_vel, boat_vel.selected)
        if boat_select is not None:
            self.u_mps = np.add(self.u_earth_no_ref_mps, boat_select._BoatData__u_processed_mps)
            self.v_mps = np.add(self.v_earth_no_ref_mps, boat_select._BoatData__v_processed_mps)
            self.nav_ref = boat_select._BoatData__nav_ref
        else:
            self.u_mps = repmat([np.nan], self.u_earth_no_ref_mps.shape[0], self.u_earth_no_ref_mps.shape[1])
            self.v_mps = repmat([np.nan], self.v_earth_no_ref_mps.shape[0], self.v_earth_no_ref_mps.shape[1])
            if boat_vel.selected == 'bt_vel':
                self.nav_ref = 'BT'
            elif boat_vel.selected == 'gga_vel':
                self.nav_ref = 'GGA'
            elif boat_vel.selected == 'vtg_vel':
                self.nav_ref = 'VTG'
        
        valid_data2 = np.copy(self.cells_above_sl)
        valid_data2[np.isnan(self.u_mps)] = False
        self.valid_data[1] = valid_data2
        
        #Duplicate original to other filters that have yet to be applied
        self.valid_data[2:] = np.tile(self.valid_data[1], [7,1,1])
        
        #Combine all filter data and update processed properties
        self.all_valid_data()
        
    def change_mag_var(self, boat_vel, mag_var_chng):
        u_NR = self.u_earth_no_ref_mps
        v_NR = self.v_earth_no_ref_mps
        dir, mag = cart2pol(u_NR, v_NR)
        u_NR_rotated, v_NR_rotated = pol2cart(dir-np.deg2rad(mag_var_chng), mag)
        self.u_earth_no_ref_mps = u_NR_rotated
        self.v_earth_no_ref_mps = v_NR_rotated
        self.set_nav_reference(boat_vel)
        
    def change_heading_source(self, boat_vel, heading):
        u_NR = self.u_earth_no_ref_mps
        v_NR = self.v_earth_no_ref_mps
        dir, mag = cart2pol(u_NR, v_NR)
        u_NR_rotated, v_NR_rotated = pol2cart(dir-np.deg2rad(repmat(heading,len(mag),1)), mag)
        self.u_earth_no_ref_mps = u_NR_rotated
        self.v_earth_no_ref_mps = v_NR_rotated
        self.set_nav_reference(boat_vel)
            
    def apply_interpolation(self, transect, kargs=None):           
        self.u_processed_mps = np.tile([np.nan], self.u_mps.shape)
        self.v_processed_mps = np.tile([np.nan], self.v_mps.shape)
        self.u_processed_mps[self.valid_data[0]] = self.u_mps[self.valid_data[0]]
        self.v_processed_mps[self.valid_data[0]] = self.v_mps[self.valid_data[0]]
        
        #Determine interpolation methods to apply
        ens_interp = self.interpolate_ens
        cells_interp = self.interpolate_cells
        
        if kargs is not None:
            n_args = len(kargs)
            for n in np.arange(1,n_args,2):
                if kargs[n] == 'Ensembles':
                    ens_interp = kargs[n+1]
                elif kargs[n] == 'Cells':
                    cells_interp = kargs[n+1]
                    
        if ens_interp == 'None':
            #Sets invalid data to nan with no interpolation
            self.interpolate_ens_none()
        elif ens_interp == 'ExpandedT': #Expanded Ensemble Time
            #Sets interpolate to None aas the interpolation is done
            #in the clsQComp
            self.interpolate_ens_next()
        elif ens_interp == 'Hold9': #SonTek Method
            #Interpolates using SonTeks method of holding last valid for up to 9 samples
            self.interpolate_ens_hold_last_9()
        elif ens_interp == 'Hold': #Hold last valid
            #Interpolates by holding last valid indefinitely
            self.interpolate_ens_hold_last()
        elif ens_interp == 'Linear':
            #Interpolates using linear interpolation
            self.interpolate_ens_linear(transect)
        elif ens_interp == 'TRDI':
            #TRDI is applied in discharge
            self.interpolate_ens_none()
            self.interpolate_ens = ens_interp
            
            
    
            
        #Apply specified cell interpolation method
        if cells_interp == 'None':
            #Sets invalid data to nan with no interpolation
            self.interpolate_cells_none()
        elif cells_interp == 'TRDI':
            #Use TRDI method to interpolate invalid interior cells
            self.interpolate_cells_TRDI()
        elif cells_interp == 'Linear':
            #Uses linear interpolation to interpolate velocity for all
            #invalid bins including those in incvalid ensembles
            #up to 9 samples
            self.interpolate_cell_linear()    
        
    def apply_filter(self, transect, kargs = None):
        #Determine filters to apply
        
        if kargs is not None:
            n_args = len(kargs)
            n = 0
            while n < n_args:
                if kargs[n] == 'Beam':
                    n += 1
                    beam_filter_setting = kargs[n]
                    self.filter_beam(beam_filter_setting)
                elif kargs[n] == 'Difference':
                    n += 1
                    d_filter_setting = kargs[n]
                    if d_filter_setting == 'Manual':
                        n+=1
                        self.filter_diff_vel(d_filter_setting, [kargs[n]])
                    else:
                        self.filter_diff_vel(d_filter_setting)
                elif kargs[n] == 'Vertical':
                    n += 1
                    w_filter_setting = kargs[n]
                    if w_filter_setting == 'Manual':
                        n += 1
                        setting = kargs[n]
                        if np.isnan(setting):
                            setting = self.w_filter_threshold
                        self.filter_vert_vel(w_filter_setting, [setting])
                    else:
                        self.filter_vert_vel(w_filter_setting, [setting])
                elif kargs[n] == 'Other':
                    n += 1
                    self.filter_smooth(transect, [kargs[n]])
                elif kargs[n] == 'Excluded':
                    n += 1
                    self.filter_excluded(transect, [kargs[n]])
                elif kargs[n] == 'SNR':
                    n += 1
                    self.filter_snr(kargs[n])
                elif kargs[n] == 'wtDepth':
                    n += 1
                    self.filter_WT_depth(transect, [kargs[n]])
                    
                n += 1
        else:
            self.filter_beam(self.beam_filter)
            self.filter_diff_vel(self.d_filter, kargs=[self.d_filter_threshold])
            self.filter_vert_vel(self.w_filter, kargs=[self.w_filter_threshold])
            self.filter_smooth(transect, self.smooth_filter)
            self.filter_excluded(transect, self.excluded_dist_m)
            self.filter_snr(self.snr_filter)
        
        self.apply_interpolation(transect)   
        
    def sos_correction(self, transect, ratio):
        self.u_mps = np.prod(self.u_mps, float(ratio))
        self.v_mps = np.prod(self.v_mps, float(ratio))
        self.u_earth_no_ref_mps = np.prod(self.u_earth_no_ref_mps,float(ratio))
        self.v_earth_no_ref_mps = np.prod(self.v_earth_no_ref_mps,float(ratio))
        self.apply_filter(transect)   
        
    def adjust_side_lobe(self, transect):   
        depth_selected = transect.depths.selected
        cells_above_SLBT = self.cells_above_sl_bt
        
        #Compute cutoff for vertical beam depths
        if depth_selected == 'vb_depths':
            depth_selected = getattr(transect.depths, transect.depths.selected)
            sl_cutoff_VB = depth_selected.depth_processed_m - \
                depth_selected.draft_use_m * cosd(transect.adcp.beam_angle_deg) \
                - self.sl_lag_effect_m + depth_selected.draft_use_m
            cells_above_SLVB = np.round(depth_selected.depth_cell_depth_m, 2) < np.round(sl_cutoff_VB, 2)
            idx = np.where(transect.depths.bt_depths.valid_data == False)
            cells_above_SLBT[:,idx] = cells_above_SLVB[:,idx]   
            cells_above_SL = cells_above_SLBT and cells_above_SLVB
        else:
            cells_above_SL = cells_above_SLBT
            
        
                
        #Compute cutoff from interpolated depths
        #Find ensembles with no valid beam depths
        idx = np.where(np.nansum(depth_selected.valid_beams) == 0)
        
        if len(idx) > 0:
            if len(self.sl_lag_effect_m) > 1:
                sl_lag_effect_m = self.sl_lag_effect_m[idx]
            else:
                sl_lag_effect_m = self.sl_lag_effect_m
                
            sl_cutoff_int = (depth_selected.depth_processed_m[idx] - depth_selected.draft_use_m) \
                * cosd(transect.adcp.beam_angle_deg) - sl_lag_effect_m + \
                depth_selected.draft_use_m
                
            cells_above_SL[:,idx] = depth_selected.depth_cell_depth_m < sl_cutoff_int
            
        #Find ensembles with at least 1 invalid beam depth
        idx = np.where(np.nansum(depth_selected.valid_beams) < 4)
        if len(idx) > 0:
            if len(self.sl_lag_effect_m) > 1:
                sl_lag_effect_m = self.sl_lag_effect_m[idx]
            else:
                sl_lag_effect_m = self.sl_lag_effect_m[idx]
                
            sl_cutoff_int = (depth_selected.depth_processed_m[idx] - depth_selected.draft_use_m\
                * cosd(transect.adcp.beam_angle_deg)) - sl_lag_effect_m + depth_selected.draft_use_m
            cells_above_SL_Int = np.ones(cells_above_SL.shape)
            cells_above_SL_Int[:,idx] = depth_selected.depth_cell_depth_m[:,idx] < sl_cutoff_int
            
            cells_above_SL[cells_above_SL_Int == 0] = 0
        
        self.cells_above_sl = cells_above_SL
        valid_vel = not np.isnan(self._u_mps)
        self.valid_data[1,:,:] = self.cells_above_sl * valid_vel
        self.all_valid_data()
        self.compute_SNR_Rng()
        self.apply_filter(transect)
        self.apply_interpolation(transect)

    def all_valid_data(self):
        """Combines the results of all filters to determine a final set of valid data"""
        n_filters = len(self.valid_data[1:,0,0])
        sum_filters = np.nansum(self.valid_data[1:,:,:],0) / n_filters
        valid = np.tile([True], self.cells_above_sl.shape)
        valid[sum_filters < 1] = False
        self.valid_data[0] = valid
        
    def __filter_beam(self, setting):
        """The determination of invalid data depends on whether
        3-beam or 4-beam solutions are acceptable.  This function can be applied by
        specifying 3 or 4 beam solutions and setting self.beam_filter to -1
        which will trigger an automatic mode.  The automatic mode will find all 3 beam
        solutions and them compare the velocity of the 2 beam solutions to nearest 4
        beam solution before and ager the 2 beam solution.  If the 3 beam solution is
        wighin 50% of the average of the neighboring 2 beam solutions the data are deemed 
        valid if not invalid.  Thus in automatic mode only those data from 2 beam solutions
        that appear sufficiently more than the 4 beam solutions are marked invalid.  The process
        happens for each ensemble.  If the number of beams is specified manually, it is applied
        uniformly for the whole transect."""
        
        self.beam_filter = setting
        
        #in manual mode determine number of raw invalid and number of 2 beam solutions if selected
        if self.beam_filter > 0:
            
            #Find invalid raw data
            valid_vel =  np.array([self.cells_above_sl] * 4)
            valid_vel[np.isnan(self.raw_vel_mps)] = 0
            
            #Determine how many beams or transformed coordinates are valid
            valid_vel_sum = np.sum(valid_vel, 0)
            valid = self.cells_above_sl
            
            #Compare number of valid beams or velocity coordinates to filter value
            valid[(valid_vel_sum < self.beam_filter) & (valid_vel_sum > 2)] = False
            
            #Save logical of valid data to object
            self.valid_data[5,:,:] = valid
        
        else:
            #Apply automatic filter
            #Create temporary object
            
            temp = np.copy(self)
            #Apply 3 beam filter to temporary object
            temp.filter_beam(4)
            
            #determine number of ensembles (NOT NECESSARY?)
            n_ens = len(temp.valid_data[5,:,:])
            
            #create matrix of valid data with nan below sidelobe
            valid = temp.valid_data[5,:,:] 
            valid[temp.cells_above_sl == False] = np.nan
            
            #Find cells with 3 beams solutions
            idx = np.where(valid == 0)
            if len(idx) > 0:
                #Find cells with 4 beams solutions
                valid_idx = np.where(valid == 1)
                
                #Valid water u and v for cells with 4 beam solutions
                valid_u = temp.u_mps(valid == 1)
                valid_v = temp.v_mps(valid == 1)
                
                #Use griddata to estimate water speed of cells with 3 beam solutions
                F = self.scatterd_interpolant(valid_idx, valid_u)
                est_u = F
                F = self.scatterd_interpolant(valid_idx, valid_v)
                est_v = F
                
                #compute the ration of estimated value to actual 3 beam solution
                idx = np.ravel_multi_index(temp.u_mps, dims=valid_idx, order='C')
                
                if len(est_u) == 0:
                    u_ratio = 1
                else:
                    u_ratio = (temp.u_mps[idx] / est_u) - 1
                    
                if len(temp.v_mps) == 0:
                    v_ratio = 1
                else:
                    v_ratio = (temp.v_mps[idx] / est_v) - 1
                    
                #If 3-beam differs from 4-beam by more 50% mark it invalid
                num_ratio = u_ratio.shape[0]
                valid[np.isnan(valid)] = 0
                
                for n in range(num_ratio):
                    if np.abs(u_ratio[n]) < .5 or np.abs(v_ratio[n]) < .5:
                        valid[idx[n]] = 1
                        
                self.valid_data[5,:,:] = valid
            else:
                self.valid_data[5,:,:] = temp.valid_data[5,:,:]\
                
    def __filter_diff_vel(self, setting, kargs=None):
        """Applies either manual or automatic filtering of the difference (error)
        velocity.  The automatic mode is based on the following:  This filter is
        based on the assumption that the water error velocity should follow a gaussian
        distribution.  Therefore, 5 standard deviations should encompass all of the
        valid data.  The standard deviation and limits (multiplier*std dev) are computed
        in an iterative process until filtering out additional data does not change the
        computed standard deviation."""
        
        #set difference filter properties
        self.d_filter = setting
        if kargs is not None:
            self.d_filter_threshold = kargs[0]
            
        #Set multiplier
        multiplier = 5
        #Get difference data from object
        d_vel = self.d_mps
        
        #Apply selected method
        if self.d_filter == 'Manual':
            d_vel_max_ref = np.abs(self.d_filter_threshold)
            d_vel_min_ref = -1 * d_vel_max_ref
        elif self.d_filter == 'Off':
            d_vel_max_ref = np.nanmax(np.nanmax(d_vel)) + 1
            d_vel_min_ref = np.nanmin(np.nanmin(d_vel)) - 1
        elif self.d_filter == 'Auto':
            #Initialize variables
            d_vel_filtered = d_vel
            std_diff = 1
            i = -1
            #Loop until no additional data are removed
            
            self.num_bins_filtered = []
            while std_diff != 0 and i < 1000:
                i = i+1
                
                #Compute standard deviation
                d_vel_std = iqr(d_vel_filtered)
                
                #Compute maximum and minimum thresholds
                d_vel_max_ref = np.nanmedian(d_vel_filtered) + multiplier * d_vel_std
                d_vel_min_ref = np.nanmedian(d_vel_filtered) - multiplier * d_vel_std
                
                #Identify valid and invalid data
                d_vel_bad_idx = np.where(d_vel_filtered > d_vel_max_ref | d_vel_filtered < d_vel)
                d_vel_good_idx = np.where(d_vel_filtered <= d_vel_max_ref | d_vel_filtered >= d_vel)
                
                #Update filtered data array
                d_vel_filtered = d_vel_filtered[d_vel_good_idx]
                
                #Determine differences due to last filter iteration
                d_vel_std2 = iqr(d_vel_filtered)
                std_diff = d_vel_std2 - d_vel_std
                self.num_bins_filtered.append(d_vel_bad_idx.shape[1])
                
        #Set valid data row 2 for difference velocity filter results
        d_vel_bad_idx = np.where((d_vel > d_vel_max_ref) | (d_vel < d_vel_min_ref))
        valid = self.cells_above_sl
        
        valid[d_vel_bad_idx] = False
        self.valid_data[2,:,:] = valid
        
        #Set threshold property
        self.d_filter_threshold = d_vel_max_ref
        
        #Combine all filter data and update processed properties
        self.all_valid_data()
            
    def __filter_vert_vel(self, setting, kargs = None):
        """Applies either manual or automatic filter of the difference (error) velocity.  The automatic
        mode is based on the following: This filter is based on the assumption that the water error
        velocity should follow a gaussian distribution.  Therefore, 4 standard deviations should
        encompass all of the valid data.  The standard deviation and limits (multplier * standard deviation)
        are computed in an iterative process until filtering out additional data does not change
        the computed standard deviation"""
        
        #Set vertical velocity filter properties
        self.w_filter = setting
        if kargs is not None:
            self.w_filter_threshold = kargs[0]
            
            #set multiplier
            multiplier = 5 
            
            #Get difference data from object
            w_vel = self.w_mps
            
            #Apply selected method
            if self.w_filter == 'Manual':
                w_vel_max_ref = np.abs(self.w_filter_threshold)
                w_vel_min_ref = -1 * w_vel_max_ref
            elif self.w_filter == 'Off':
                w_vel_max_ref = np.nanmax(np.nanmax(w_vel)) + 1
                w_vel_min_ref = np.nanmin(np.nanmin(w_vel)) - 1
            elif self.w_filter == 'Auto':
                #Initialize variables
                w_vel_filtered = w_vel[:]
                std_diff = 1
                i = 0
                num_bins_filtered = []
                #Loop until no additional data are removed
                while std_diff != 0 and i < 1000:
                    
                    #Computed standard deviation
                    w_vel_std = iqr(w_vel_filtered)
                    
                    #Compute maximum and minimum thresholds
                    w_vel_max_ref = np.nanmedian(w_vel_filtered) + multiplier * w_vel_std
                    w_vel_min_ref = np.nanmedian(w_vel_filtered) - multiplier * w_vel_std
                    
                    #Identify valid and invalid data
                    w_vel_bad_idx = np.where((w_vel_filtered > w_vel_max_ref) or (w_vel_filtered < w_vel_min_ref))
                    w_vel_good_idx = np.where((w_vel_filtered <= w_vel_max_ref) or (w_vel_filtered >= w_vel_min_ref))
                    
                    #Update filtered data array
                    w_vel_filtered = w_vel_filtered[w_vel_good_idx]
                    
                    #Determine differences due to last filter iteration
                    w_vel_std2 = iqr(w_vel_filtered)
                    std_diff = w_vel_std2 - w_vel_std
                    num_bins_filtered.append(len(w_vel_bad_idx))
                    
                #Set valid data row 3 for difference velocity filter results
                w_vel_bad_idx = np.where((w_vel > w_vel_max_ref) or (w_vel < w_vel_min_ref))
                valid = self.cells_above_sl
                
                valid[w_vel_bad_idx] = False
                self.valid_data[3,:,:] = valid
                
                #Set threshold property
                self.w_filter_threshold = w_vel_max_ref
                
                #Combine all filter data and update processed properties
                self.all_valid_data()
                
    def __filter_smooth(self, transect, setting):
        """Filter boat speed.  Running Standard Deviation filter for water speed
           This filter employs a running trimmed standard deviation filter to
          identify and mark spikes in the water speed. First a robust Loess 
          smooth is fitted to the water speed time series and residuals between
         the raw data and the smoothed line are computed. The trimmed standard
          deviation is computed by selecting the number of residuals specified by
          "halfwidth" before the target point and after the target point, but not
          including the target point. These values are then sorted, and the points
          with the highest and lowest values are removed from the subset, and the 
          standard deviation of the trimmed subset is computed. The filter
         criteria are determined by multiplying the standard deviation by a user
         specified multiplier. This criteria defines a maximum and minimum
          acceptable residual. Data falling outside the criteria are set to nan.
          
          Reccomended filter settings are:
          filter_width = 10
          half_width = 10
          multiplier = 9
        """
        
        self.smooth_filter = setting
        
        #Compute ens_time
        ens_time = np.nancumsum(transect.datetime.ens_duration_sec)
        
        #Determine if smooth filter should be applied
        if self.smooth_filter == 'Auto':
            
            #Boat velocity components
            w_vele = self.u_mps
            w_veln = self.v_mps
            
            #Set filter parameters
            filter_width = 10
            half_width = 10
            multiplier = 9
            cycles = 3
            
            #Initialize variables
            dir = np.tile([np.nan], w_vele.shape)
            speed = np.tile([np.nan], w_vele.shape)
            speed_smooth = np.tile([np.nan], w_vele.shape)
            speed_res = np.tile([np.nan], w_vele.shape)
            speed_filtered = np.tile([np.nan], w_vele.shape)
            w_vele_filtered = w_vele
            w_veln_filtered = w_veln
            wt_bad = np.tile([np.nan], w_vele.shape)
            
            #Compute mean speed and direction of water
            w_vele_avg = np.nanmean(w_vele)
            w_veln_avg = np.nanmean(w_veln)
            dir, speed = cart2pol(w_vele_avg, w_veln_avg)
            dir = np.rad2deg(dir)
            
            #Compute residuals from a robust Loess smooth
            
            speed_smooth = lowess(ens_time, speed, filter_width / len(speed))
            speed_res = speed - speed_smooth
            
            #Apply a trimmed standard deviation filter multiple times
            speed_filtered = speed
            for i in range(cycles):
                fill_array = run_std_trim(half_width, speed_res.T)
                
                #Compute filter bounds
                upper_limit = speed_smooth + multiplier * fill_array
                lower_limit = speed_smooth - multiplier * fill_array
                
                #Apply filter to residuals
                wt_bad_idx = np.where((speed > upper_limit) or (speed < lower_limit))
                speed_res[wt_bad_idx] = np.nan
            
            valid = self.cells_above_sl
            
            valid[:, wt_bad_idx] = False
            self.valid_data[4,:,:] = valid
            self.smooth_upper_limit = upper_limit
            self.smooth_lower_limit = lower_limit
            self.smooth_speed = speed_smooth
        
        else:
            #No filter applied
            self.valid_data[4,:,:] = self.cells_above_sl
            self.smooth_upper_limit = np.nan
            self.smooth_lower_limit = np.nan
            self.smooth_speed = np.nan
            
        self.all_valid_data()
     
    def __filter_snr(self, setting):
        self.snr_filter = setting  
        
        if setting == 'Auto':
            if self.snr_rng is not None:
                bad_snr_idx = self.sn_rng > 12
                valid = self.cells_above_sl
                
                bad_snr_array = np.tile(bad_snr_idx, (valid.shape[0], 1))
                valid[bad_snr_array] = False
                self.valid_data[7,:,:] = valid
                #Combine all filter data and update processed properties
                self.all_valid_data()
        else:
            self.valid_data[7,:,:] = self.cells_above_sl
            self.all_valid_data()
        
    def __filter_wt_depth(self, transect, setting):
        self.wt_depth_filter = setting
        valid = self.cells_above_sl
        
        if setting == 'On':
            trans_select = getattr(transect.depths, transect.depths.selected)
            valid[:, np.isnan(trans_select.depth_processed_m)] = False
        self.valid_data[8,:,:] = valid
        
        self.all_valid_data()
        
    def __filter_excluded(self, transect, setting):
        """Marks all data with the cell top freater than the setting invalid"""
        trans_select = getattr(transect.depths, transect.depths.selected)
        cell_depth = trans_select.depth_cell_depth_m
        cell_size = trans_select.depth_cell_size_m
        draft = trans_select.draft_use_m
        top_cell_depth = cell_depth - 0.5 * cell_size
        threshold = np.round((setting+draft),3)
        exclude = np.round(top_cell_depth, 3) <= threshold
        valid = self.cells_above_sl
        valid[exclude] = False
        self.valid_data[6,:,:] = valid
        
        #Set threshold property
        self.excluded_dist_m = setting
        
        self.all_valid_data()

    def __interpolate_ens_next(self):
        #Set interpolation property for ensembles
        self.interpolate_ens = 'ExpandedT'
        
        #Set processed data to nan for all invalid data  
        valid = self.valid_data[0]
        self.u_processed_mps = self.u_mps
        self.v_processed_mps = self.v_mps
        self.u_processed_mps[valid[0] == False] = np.nan
        self.v_processed_mps[valid[0] == False] = np.nan
        
        #Identifying ensembles with no valid data
        valid_ens = np.any(valid)
        n_ens = len(valid_ens)
        
        #Set the invalid ensembles to the data in the next valid ensemble
        for n in np.arange(0,n_ens-1)[::-1]:
            if valid_ens[n] == False:
                self.u_processed_mps[:,n] = self.u_processed_mps[:,n+1]
                self.v_processed_mps[:,n] = self.v_processed_mps[:,n+1]
                
    def __interpolate_ens_hold_last(self):
        """Interpolates velocity data for invalid ensembles by repeating the
        the last valid data until new valid data is found"""
        
        self.interpolate_ens = 'HoldLast'
        
        valid = self.valid_data[0]
        
        #Initialize processed velocity data variables
        self.u_processed_mps = self.u_mps
        self.v_processed_mps = self.v_mps
        
        #Set invalid data to nan in processed velocity data variables
        self.u_processed_mps[valid[0] == False] = np.nan
        self.v_processed_mps[valid[0] == False] = np.nan
        
        #Determine ensembles with valid data
        valid_ens = np.any(valid)
        
        #Process each ensemble beginning with the second ensemble
        n_ens = len(valid_ens)
        
        for n in np.arange(1,n_ens):
            #If ensemble is invalid fill in with previous ensemble
            if valid_ens[n] == False:
                self.u_processed_mps[:,n] = self.u_processed_mps[:,n-1]
                self.v_processed_mps[:,n] = self.v_processed_mps[:,n-1]

    def __interpolate_ens_hold_last_9(self):
        """Interpolates velocity data for invalid ensembles by repeating the
        last valid data for up to 9 ensembles or until new valid data is
        found. If more the 9 consectutive ensembles are invalid the
        ensembles beyond the 9th remain invalid. This is for
        compatibility with SonTek RiverSurveyor Live.
        """
        
        self.interpolate_ens = 'Hold9'
        
        valid = self.valid_data[0]
        
        #Initialize processed velocity data variables
        self.u_processed_mps = self.u_mps
        self.v_processed_mps = self.v_mps
        
        #Set invalid data to nan in processed velocity data variables
        self.u_processed_mps[valid[0] == False] = np.nan
        self.v_processed_mps[valid[0] == False] = np.nan
        
        #Determine ensembles with valid data
        valid_ens = np.any(valid)
        
        #Process each ensemble beginning with the second ensemble
        n_ens = len(valid_ens)
        n_invalid = 0
        
        for n in np.arange(1,n_ens):
            #If ensemble is invalid fill in with previous ensemble
            if valid_ens[n] == False and n_invalid < 10:
                n_invalid += 1
                self.u_processed_mps[:,n] = self.u_processed_mps[:,n-1]
                self.v_processed_mps[:,n] = self.v_processed_mps[:,n-1]
            else:
                n_invalid = 0

    def __interpolate_ens_none(self):
        """Applies no interpolation for invalid ensembles"""
        
        self.interpolate_ens = 'None'
        
        valid = self.valid_data[0]
        
        #Initialize processed velocity data variables
        self.u_processed_mps = self.u_mps
        self.v_processed_mps = self.v_mps
        
        #Set invalid data to nan in processed velocity data variables
        self.u_processed_mps[valid == False] = np.nan
        self.v_processed_mps[valid == False] = np.nan

    def __interpolate_cells_none(self):
        self.interpolate_cells = 'None'
        
        valid = self.valid_data[0]
        
        self.u_processed_mps = self.u_mps
        self.v_processed_mps = self.v_mps
        
        #Set invalid data to nan in processed velocity data variables
        self.u_processed_mps[valid == False] = np.nan
        self.v_processed_mps[valid == False] = np.nan
        
    def __interpolate_ens_linear(self, transect):
        """Use linear interpolation as computed by Matlab'e scatter interpolant
        function to interpolated velocit data for ensembles with no valid velocities
        """
        
        self.interpolate_ens = 'Linear'
         
        valid = self.valid_data[0]
        
        #Initialize processed velocity data variables
        self.u_processed_mps = self.u_mps
        self.v_processed_mps = self.v_mps
        
        #Set invalid data to nan in processed velocity data variables
        self.u_processed_mps[valid[0] == False] = np.nan
        self.v_processed_mps[valid[0] == False] = np.nan
                
        #Determine ensembles with valid data
        valid_ens = np.any(valid)
        
        if np.sum(valid_ens) > 1:
            #Determine the numbe of ensembles
            n_ens = len(valid_ens)
            
            trans_select = getattr(transect.depths, transect.depths.selected)
            #compute z
            z = np.divide(np.subtract(trans_select.depth_processed_m, trans_select.depth_cell_depth_m),
                          trans_select.depth_processed_m)
            
            #Create position array
            boat_select = getattr(transect.boat_vel, transect.boat_vel.selected)
            if boat_select is not None:
                if np.nansum(boat_select.valid_data[0]) > 0:
                    boat_vel_x = boat_select.u_processed_mps
                    boat_vel_y = boat_select.v_processed_mps
                    track_x = boat_vel_x * transect.datetime.ens_duration_sec
                    track_y = boat_vel_y * transect.datetime.ens_duration_sec
                    track = np.nancumsum(np.sqrt(track_x**2 + track_y**2))
                    track_array = np.tile(track, (self.u_processed_mps.shape[0], 1))
                    
                    #Determine index of all valid data
                    valid_z = np.isnan(z) == False
                    valid_combined = valid & valid_z
                    
                    x_index = len(z)
                    xyi = np.linspace(np.nanmin(z), np.nanmax(z), x_index)
                    y_index = len(track_array)
                    yyi = np.linspace(np.nanmin(track_array), np.nanmax(track_array), y_index)
                    XI, YI = np.meshgrid(xyi,yyi)
                    
                    Fu = interpolate(np.vstack([z[valid_combined],track_array[valid_combined]], 
                                               self.u_processed_mps[valid_combined],
                                               (XI, YI)))
                    
                    Fv = interpolate(np.vstack([z[valid_combined],track_array[valid_combined]], 
                                               self.v_processed_mps[valid_combined],
                                               (XI, YI)))

    def compute_snr_rng(self):
        """Computes the range between the average snr for all beams.

        The average is computed using only data above the side lobe cutoff."""
        cells_above_sl = np.copy(float(self.cells_above_sl))
        cells_above_sl[cells_above_sl < 0.5] = np.nan
        snr_adjusted = self.rssi * cells_above_sl
        snr_average = np.squeeze(np.nanmean(snr_adjusted,0))
        self.snr_rng = np.nanmax(snr_average) - np.nanmin(snr_average)