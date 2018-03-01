import numpy as np
from numpy.matlib import repmat
import time
import os
from Classes.Pd0TRDI import Pd0TRDI
from Classes.DepthStructure import DepthStructure
from MiscLibs.convenience import cosd, arctand, tand
from Classes.WaterData import WaterData
from Classes.BoatStructure import BoatStructure
from Classes.GPSData import GPSData
from Classes.Edges import Edges
from Classes.ExtrapData import ExtrapData
from Classes.Sensors import Sensors
from Classes.SensorData import SensorData
from Classes.HeadingData import HeadingData
from Classes.DateTime import DateTime
from Classes.InstrumentData import InstrumentData
from Classes.MultiThread import MultiThread
import matplotlib.dates as mdates
from datetime import datetime
from MiscLibs.convenience import nandiff


class TransectData(object):
    """Class to hold Transect properties.

    Attributes
    ----------
    adcp: object
        Object of clsInstrument
    file_name: str
        Filename of transect data file
    w_vel: object
        Object of clsWaterData
    boat_vel: object
        Object of BoatStructure containing objects of BoatData for BT, GGA, and VTG
    gps: object
        Object of clsGPSData
    sensors: object
        Object of clsSensorData
    depths: object
        Object of clsDepthStructure containing objects of Depth data for btDepths, vbDepths, dsDepths)
    edges: object
        Object of clsEdges (left and right object of clsEdgeData)
    extrap: object
        Object of clsExtrapData
    start_edge: str
        Starting edge of transect looking downstream (Left or Right)
    date_time: object
        Object of DateTime
    checked: bool
        Setting for if transect was checked for use in mmt file assumed checked for SonTek
    in_transect_idx: np.array(int)
        Index of ensemble data associated with the moving-boat portion of the transect
    """

    def __init__(self):
        self.adcp = None  # object of clsInstrument
        self.file_name = None  # filename of transect data file
        self.w_vel = None  # object of clsWaterData
        self.boat_vel = None  # class for various boat velocity references (btVel, ggaVel, vtgVel)
        self.gps = None  # object of clsGPSData
        self.sensors = None  # object of clsSensorData
        self.depths = None  # object of clsDepthStructure for depth data including cell depths & ref depths
        self.edges = None  # object of clsEdges(left and right object of clsEdgeData)
        self.extrap = None # object of clsExtrapData
        self.start_edge = None # starting edge of transect looking downstream (Left or Right)
        self.date_time = None # object of DateTime
        self.checked = None  # transect was checked for use in mmt file assumed checked for SonTek
        self.in_transect_idx = None  # index of ensemble data associated with the moving-boat portion of the transect

    def get_data(self, source, in_file, pd0_data, mmt, mbt_idx):
        # TODO Eliminate this method
        # TODO Ensure adjusteSideLobe is applied as in the Matlab code
        if source == 'TRDI':
            self.trdi(in_file, pd0_data, mmt, mbt_idx)
            
                                #files2load idx

    def trdi(self, mmt_transect, pd0_data, mmt, mbt_idx, kargs=None):

        # TODO why is this needed
        # self.mbt = mbt_idx
        pd0 = pd0_data
        
        # Get the configuration property of the mmt_transect
        mmt_config = getattr(mmt_transect, self.active_config)
        if pd0_data.Wt is not None:
            # Get and compute ensemble beam depths
            temp_depth = np.array(pd0.Bt.depth_m)
            # Screen out invalid depths
            temp_depth[temp_depth < 0.01] = np.nan
            # Add draft
            temp_depth += mmt_config['Offsets_Transducer_Depth']
            
            # Get instrument cell data
            cell_size_all_m, cell_depth_m, sl_cutoff_per, sl_lag_effect_m = \
                TransectData.compute_instrument_cell_data(pd0)
            
            # Adjust cell depth of draft
            cell_depth_m = np.add(mmt_config['Offsets_Transducer_Depth'], cell_depth_m)
            
            # Create depth data object for BT
            self.depths = DepthStructure()
            self.depths.add_depth_object(depth_in=temp_depth,
                                         source_in='BT',
                                         freq_in=pd0.Inst.freq,
                                         draft_in=mmt_config['Offsets_Transducer_Depth'],
                                         cell_depth_in=cell_depth_m,
                                         cell_size_in=cell_size_all_m)
            
            # Compute cells above side lobe
            cells_above_sl = TransectData.side_lobe_cutoff(depths=self.depths.bt_depths.depth_orig_m,
                                                           draft=self.depths.bt_depths.draft_orig_m,
                                                           cell_depth=self.depths.bt_depths.depth_cell_depth_m,
                                                           sl_lag_effect=sl_lag_effect_m,
                                                           type='Percent',
                                                           value=1-sl_cutoff_per / 100)
            
            # Check for the presence of vertical beam data
            if np.nanmax(np.nanmax(pd0.Sensor.vert_beam_status)) > 0:
                temp_depth = pd0.Sensor.vert_beam_range_m
                
                # Screen out invalid depths
                temp_depth[temp_depth < 0.01] = np.nan
                
                # Add draft
                temp_depth = temp_depth + mmt_config['Offsets_Transducer_Depth']
                
                # Create depth data object for vertical beam
                self.depths.add_depth_object(depth_in=temp_depth,
                                             source_in='VB',
                                             freq_in=pd0.Inst.freq,
                                             draft_in=mmt_config['Offsets_Transducer_Depth'],
                                             cell_depth_in=cell_depth_m,
                                             cell_size_in=cell_size_all_m)
                                   
            # Check for the presence of depth sounder
            if np.nansum(np.nansum(pd0.Gps2.depth_m)) > 1e-5:
                temp_depth = pd0.Gps2.depth_m
                
                # Screen out invalid data
                temp_depth[temp_depth < 0.01] = np.nan
                
                # Use the last valid depth for each ensemble
                last_depth_col_idx = np.sum(np.isnan(temp_depth) == False, axis=1)-1
                last_depth_col_idx[last_depth_col_idx == -1] = 0               
                row_index = np.arange(len(temp_depth))
                last_depth = np.empty(row_index.size)
                for row in row_index:
                    last_depth[row] = temp_depth[row, last_depth_col_idx[row]]

                # Determine if mmt file has a scale factor and offset for the depth sounder
                if mmt_config['DS_Cor_Spd_Sound'] == 0:
                    scale_factor = mmt_config['DS_Scale_Factor']
                else:
                    scale_factor = pd0.sensorsdata.sos_mps / 1500.
                    
                # Apply scale factor, offset, and draft
                # Note: Only the ADCP draft is stored.  The transducer
                # draft or scaling for depth sounder data cannot be changed in QRev
                ds_depth = (last_depth * scale_factor) \
                    + mmt_config['DS_Transducer_Depth']\
                    + mmt_config['DS_Transducer_Offset']
                
                self.depths.add_depth_object(depth_in=ds_depth,
                                             source_in='DS',
                                             freq_in=pd0.Inst.freq,
                                             draft_in=mmt_config['Offsets_Transducer_Depth'],
                                             cell_depth_in=cell_depth_m,
                                             cell_size_in=cell_size_all_m)
                
            # Set depth reference to value from mmt file
            if 'Proc_River_Depth_Source' in mmt_config:
                if mmt_config['Proc_River_Depth_Source'] == 0:
                    self.depths.set_depth_reference('BT')
                    self.depths.composite_depths(self, setting=False)

                elif mmt_config['Proc_River_Depth_Source'] == 1:
                    if self.depths.ds_depths is not None:
                        self.depths.set_depth_reference('DS')
                    else:
                        self.depths.set_depth_reference('BT')
                    self.depths.composite_depths(self, setting=False)

                elif mmt_config['Proc_River_Depth_Source'] == 2:
                    if self.depths.vb_depths is not None:
                        self.depths.set_depth_reference('VB')
                    else:
                        self.depths.set_depth_reference('BT')
                    self.depths.composite_depths(self, setting=False)

                elif mmt_config['Proc_River_Depth_Source'] == 3:
                    if self.depths.vb_depths is None:
                        self.depths.set_depth_reference('BT')
                        self.depths.composite_depths(self, setting=False)
                    else:
                        self.depths.set_depth_reference('VB')
                        self.depths.composite_depths(self, setting=True)

                elif mmt_config['Proc_River_Depth_Source'] == 4:
                    if self.depths.bt_depths is not None:
                        self.depths.set_depth_reference('BT')
                        self.depths.composite_depths(self, setting=True)
                    elif self.depths.vb_depths is not None:
                        self.depths.set_depth_reference('VB')
                        self.depths.composite_depths(self, setting=True)
                    elif self.depths.ds_depths is not None:
                        self.depths.set_depth_reference('DS')
                        self.depths.composite_depths(self, setting=True)
                else:
                    self.depths.set_depth_reference('BT')
                    self.depths.composite_depths(self, setting=False)
            else:
                if mmt_config['DS_Use_Process'] > 0:
                    if self.depths.ds_depths is not None:
                        self.depths.set_depth_reference('DS')
                    else:
                        self.depths.set_depth_reference('BT')
                else:
                    self.depths.set_depth_reference('BT')
                self.depths.composite_depths(self, setting=False)
                
            # Create water_data object
            # ------------------------
            
            # Check for RiverRay and RiverPro data
            firmware = str(pd0.Inst.firm_ver[0])
            excluded_dist = 0
            if (firmware[:2] == '56') and (np.nanmax(np.isnan(pd0.Sensor.vert_beam_status))):
                excluded_dist = 0.25
                
            if (firmware[:2] == '44') or (firmware[:2] == '56'):
                # Process water velocities for RiverRay and RiverPro
                self.w_vel = WaterData()
                self.w_vel.populate_data(vel_in=pd0.Wt.vel_mps,
                                         freq_in=pd0.Inst.freq.T,
                                         coord_sys_in=pd0.Cfg.coord_sys,
                                         nav_ref_in='None',
                                         rssi_in=pd0.Wt.rssi,
                                         rssi_units_in='Counts',
                                         excluded_dist_in=excluded_dist,
                                         cells_above_sl_in=cells_above_sl,
                                         sl_cutoff_per_in=sl_cutoff_per,
                                         sl_cutoff_num_in=0,
                                         sl_cutoff_type_in='Percent',
                                         sl_lag_effect_in=sl_lag_effect_m,
                                         wm_in=pd0.Cfg.wm[0],
                                         blank_in=pd0.Cfg.wf_cm[0] / 100,
                                         corr_in=pd0.Wt.corr,
                                         surface_vel_in=pd0.Surface.vel_mps,
                                         surface_rssi_in=pd0.Surface.rssi,
                                         surface_corr_in=pd0.Surface.corr,
                                         surface_num_cells_in=pd0.Surface.no_cells)
                
            else:
                # Process water velocities for non-RiverRay ADCPs
                self.w_vel = WaterData()
                self.w_vel.populate_data(vel_in=pd0.Wt.vel_mps,
                                         freq_in=pd0.Inst.freq.T,
                                         coord_sys_in=pd0.Cfg.coord_sys,
                                         nav_ref_in='None',
                                         rssi_in=pd0.Wt.rssi,
                                         rssi_units_in='Counts',
                                         excluded_dist_in=excluded_dist,
                                         cells_above_sl_in=cells_above_sl,
                                         sl_cutoff_per_in=sl_cutoff_per,
                                         sl_cutoff_num_in=0,
                                         sl_cutoff_type_in='Percent',
                                         sl_lag_effect_in=sl_lag_effect_m,
                                         wm_in=pd0.Cfg.wm[0],
                                         blank_in=pd0.Cfg.wf_cm[0] / 100,
                                         corr_in=pd0.Wt.corr)
                
            # Initialize boat vel
            self.boat_vel = BoatStructure()
            # Apply 3-beam setting from mmt file
            if mmt_config['Proc_Use_3_Beam_Solution_For_BT'] < 0.5:
                min_beams = 4
            else:
                min_beams = 3
            self.boat_vel.add_boat_object(source='TRDI',
                                          vel_in=pd0.Bt.vel_mps,
                                          freq_in=pd0.Inst.freq.T,
                                          coord_sys_in=pd0.Cfg.coord_sys[0],
                                          nav_ref_in='BT',
                                          min_beams=min_beams,
                                          bottom_mode=pd0.Cfg.bm[0])
            
            self.boat_vel.set_nav_reference('BT')
            
            # Compute velocities from GPS Data
            # ------------------------------------
            # Stopped here 2/27/2018
            # Raw Data
            raw_gga_utc = pd0.Gps2.utc
            raw_gga_lat = pd0.Gps2.lat_deg
            raw_gga_lon = pd0.Gps2.lon_deg

            # Determine correct sign for latitude
            idx = np.where(pd0.Gps2.lat_ref == 'S')[0]
            if len(idx) > 0:
                raw_gga_lat[idx] = raw_gga_lat[idx] * -1

            # Determine correct sign for longitude
            idx = np.where(pd0.Gps2.lon_ref == 'W')
            if len(idx) > 0:
                raw_gga_lon[idx] = raw_gga_lon[idx] * -1
            
            # Assign data to local variables
            raw_gga_alt = pd0.Gps2.alt
            raw_gga_diff = pd0.Gps2.corr_qual
            raw_gga_hdop = pd0.Gps2.hdop
            raw_gga_num_sats = pd0.Gps2.num_sats
            raw_vtg_course = pd0.Gps2.course_true
            raw_vtg_speed = pd0.Gps2.speed_k_mph * 0.2777778
            raw_vtg_delta_time = pd0.Gps2.vtg_delta_time
            raw_vtg_mode_indicator = pd0.Gps2.mode_indicator
            raw_gga_delta_time = pd0.Gps2.gga_delta_time
            
            # RSL provided ensemble values, not supported for TRDI data
            ext_gga_utc = []
            ext_gga_lat = []
            ext_gga_lon = []
            ext_gga_alt = []
            ext_gga_diff = []
            ext_gga_hdop = []
            ext_gga_num_sats = []
            ext_vtg_course = []
            ext_vtg_speed = []
             
            # QRev methods GPS processing methods
            gga_p_method = 'Mindt'
            gga_v_method = 'Mindt'
            vtg_method = 'Mindt'
            
            # If valid gps data exist, process the data
            if (np.sum(np.sum(np.isnan(raw_gga_lat) == False)) > 0) \
                    or (np.sum(np.sum(np.isnan(raw_vtg_speed) == False)) > 0):
                
                # Process raw GPS data
                self.gps = GPSData()
                self.gps.populate_data(raw_gga_utc=raw_gga_utc,
                                       raw_gga_lat=raw_gga_lat,
                                       raw_gga_lon=raw_gga_lon,
                                       raw_gga_alt=raw_gga_alt,
                                       raw_gga_diff=raw_gga_diff,
                                       raw_gga_hdop=raw_gga_hdop,
                                       raw_gga_num_sats=raw_gga_num_sats,
                                       raw_gga_delta_time=raw_gga_delta_time,
                                       raw_vtg_course=raw_vtg_course,
                                       raw_vtg_speed=raw_vtg_speed,
                                       raw_vtg_delta_time=raw_vtg_delta_time,
                                       raw_vtg_mode_indicator=raw_vtg_mode_indicator,
                                       ext_gga_utc=ext_gga_utc,
                                       ext_gga_lat=ext_gga_lat,
                                       ext_gga_lon=ext_gga_lon,
                                       ext_gga_alt=ext_gga_alt,
                                       ext_gga_diff=ext_gga_diff,
                                       ext_gga_hdop=ext_gga_hdop,
                                       ext_gga_num_sats=ext_gga_num_sats,
                                       ext_vtg_course=ext_vtg_course,
                                       ext_vtg_speed=ext_vtg_speed,
                                       gga_p_method=gga_p_method,
                                       gga_v_method=gga_v_method,
                                       vtg_method=vtg_method)
                
                # If valid gga data exists create gga boat velocity object
                if np.sum(np.sum(np.isnan(raw_gga_lat) == False)) > 0:
                    self.boat_vel.add_boat_object(source='TRDI',
                                                  vel_in=self.gps.gga_velocity_ens_mps,
                                                  coord_sys_in='Earth',
                                                  nav_ref_in='GGA')

                # If valid vtg data exist create vtg boat velocity object
                if np.sum(np.sum(np.isnan(raw_vtg_speed) == False)) > 0:
                    self.boat_vel.add_boat_object(source='TRDI',
                                                  vel_in=self.gps.vtg_velocity_ens_mps,
                                                  coord_sys_in='Earth',
                                                  nav_ref_in='VTG')

            # Create Edges Object
            self.edges = Edges()
            self.edges.populate_data(rec_edge_method='Fixed', vel_method='MeasMag')
                
            # Determine number of ensembles to average
            n_ens_left = mmt_config['Q_Shore_Pings_Avg']
            # TRDI uses same number on left and right edges
            n_ens_right = n_ens_left
            
            # Set indices for ensembles in the moving-boat portion of the transect
            self.in_transect_idx = np.arange(0, pd0.Bt.vel_mps.shape[1])
            
            # Determine left and right edge distances
            if mmt_config['Edge_Begin_Left_Bank']:
                dist_left = mmt_config['Edge_Begin_Shore_Distance']
                dist_right = mmt_config['Edge_End_Shore_Distance']
                self.start_edge = 'Left'
            else:
                dist_left = mmt_config['Edge_End_Shore_Distance']
                dist_right = mmt_config['Edge_Begin_Shore_Distance']
                self.start_edge = 'Right'
                
            # Create left edge
            if mmt_config['Q_Left_Edge_Type'] == 0:
                self.edges.left.populate_data(edge_type='Triangular',
                                              distance=dist_left,
                                              number_ensembles=n_ens_left)

            elif mmt_config['Q_Left_Edge_Type'] == 1:
                self.edges.left.populate_data(edge_type='Rectangular',
                                              distance=dist_left,
                                              number_ensembles=n_ens_left)

            elif mmt_config['Q_Left_Edge_Type'] == 2:
                self.edges.left.populate_data(edge_type='Custom',
                                              distance=dist_left,
                                              number_ensembles=n_ens_left,
                                              coefficient=mmt_config['Q_Left_Edge_Coeff'])
                
            # Create right edge
            if mmt_config['Q_Right_Edge_Type'] == 0:
                self.edges.right.populate_data(edge_type='Triangular',
                                               distance=dist_right,
                                               number_ensembles=n_ens_right)

            elif mmt_config['Q_Right_Edge_Type'] == 1:
                self.edges.right.populate_data(edge_type='Rectangular',
                                               distance=dist_right,
                                               number_ensembles=n_ens_right)

            elif mmt_config['Q_Right_Edge_Type'] == 2:
                self.edges.right.populate_data(edge_type='Custom',
                                               distance=dist_right,
                                               number_ensembles=n_ens_right,
                                               coefficient=mmt_config['Q_Right_Edge_Coeff'])
                
            # Create extrap object
            # --------------------
            # Determine top method
            top = 'Power'
            if mmt_config['Q_Top_Method'] == 1:
                top = 'Constant'
            elif mmt_config['Q_Top_Method'] == 2:
                top = '3-Point'
                
            # Determine bottom method
            bot = 'Power'
            if mmt_config['Q_Bottom_Method'] == 2:
                bot = 'No Slip'
                
            self.extrap = ExtrapData()
            self.extrap.populate_data(top=top, bot=bot, exp=mmt_config['Q_Power_Curve_Coeff'])
            
            # Sensor Data
            self.sensors = Sensors() 
            
            # Heading
            
            # Internal Heading
            self.sensors.heading_deg.internal = HeadingData()
            self.sensors.heading_deg.internal.populate_data(data_in=pd0.Sensor.heading_deg.T,
                                                            source_in=pd0.Cfg.head_src[0],
                                                            magvar=mmt_config['Offsets_Magnetic_Variation'],
                                                            align=mmt_config['Ext_Heading_Offset'])

            # External Heading
            ext_heading_check = np.where(np.isnan(pd0.Gps2.heading_deg))
            if len(ext_heading_check[0]) <= 0:
                self.sensors.heading_deg.selected = 'internal'
            else:
                # Determine external heading for each ensemble
                # Using the minimum time difference
                d_time = np.abs(pd0.Gps2.hdt_delta_time)
                d_time_min = np.nanmin(d_time.T, 0).T
                use = np.tile([np.nan], d_time.shape)
                for nd_time in range(len(d_time_min)):
                    use[nd_time, :] = np.abs(d_time[nd_time, :]) == d_time_min[nd_time]
                    
                ext_heading_deg = np.tile([np.nan], (len(d_time_min)))
                for nh in range(len(d_time_min)):
                    idx = np.where(use[nh, :])[0]
                    if len(idx) > 0:
                        idx = idx[0]
                        ext_heading_deg[nh] = pd0.Gps2.heading_deg[nh, idx]
                        
                # Create external heading sensor
                self.sensors.heading_deg.external = HeadingData()
                self.sensors.heading_deg.external.populate_data(data_in=ext_heading_deg,
                                                                source_in='GPS',
                                                                magvar=mmt_config['Offsets_Magnetic_Variation'],
                                                                align=mmt_config['Ext_Heading_Offset'])

                # Determine heading source to use from mmt setting
                source_used = mmt_config['Ext_Heading_Use']
                if source_used:
                    self.sensors.heading_deg.selected = 'external'
                else:
                    self.sensors.heading_deg.selected = 'internal'

            # Pitch
            pitch = arctand(tand(pd0.Sensor.pitch_deg) * cosd(pd0.Sensor.roll_deg))
            pitch_src = pd0.Cfg.pitch_src[0]
            
            # Create pitch sensor
            self.sensors.pitch_deg.internal = SensorData()
            self.sensors.pitch_deg.internal.populate_data(data_in=pitch, source_in=pitch_src)
            self.sensors.pitch_deg.selected = 'internal'
            
            # Roll
            roll = pd0.Sensor.roll_deg.T
            roll_src = pd0.Cfg.roll_src[0]
            
            # Create Roll sensor
            self.sensors.roll_deg.internal = SensorData()
            self.sensors.roll_deg.internal.populate_data(data_in=roll, source_in=roll_src)
            self.sensors.roll_deg.selected = 'internal'
            
            # Temperature
            temperature = pd0.Sensor.temperature_deg_c.T
            temperature_src = pd0.Cfg.temp_src[0]
            
            # Create temperature sensor
            self.sensors.temperature_deg_c.internal = SensorData()
            self.sensors.temperature_deg_c.internal.populate_data(data_in=temperature, source_in=temperature_src)
            self.sensors.temperature_deg_c.selected = 'internal'
            
            # Salinity
            pd0_salinity = pd0.Sensor.salinity_ppt.T
            pd0_salinity_src = pd0.Cfg.sal_src[0]
            
            # Create salinity sensor from pd0 data
            self.sensors.salinity_ppt.internal = SensorData()
            self.sensors.salinity_ppt.internal.populate_data(data_in=pd0_salinity, source_in=pd0_salinity_src)

            # Create salinity sensor from mmt data
            mmt_salinity = mmt_config['Proc_Salinity']
            self.sensors.salinity_ppt.user = SensorData()
            self.sensors.salinity_ppt.user.populate_data(data_in=mmt_salinity, source_in='mmt')

            # Set selected salinity
            self.sensors.salinity_ppt.selected = 'internal'
            
            # Speed of Sound
            speed_of_sound = pd0.Sensor.sos_mps.T
            speed_of_sound_src = pd0.Cfg.sos_src[0]
            self.sensors.speed_of_sound_mps.internal = SensorData()
            self.sensors.speed_of_sound_mps.internal.populate_data(data_in=speed_of_sound, source_in=speed_of_sound_src)
            
            # The raw data are referenced to the internal SOS
            self.sensors.speed_of_sound_mps.selected = 'internal'
            
            # Ensemble times
            # Compute time for each ensemble in seconds
            ens_time_sec = pd0.Sensor.time[:, 0] * 3600 \
                + pd0.Sensor.time[:, 1] * 60 \
                + pd0.Sensor.time[:, 2] \
                + pd0.Sensor.time[:, 3] / 100
            
            # Compute the duration of each ensemble in seconds adjusting for lost data
            ens_delta_time = np.tile([np.nan], ens_time_sec.shape)
            idx_time = np.where(np.isnan(ens_time_sec) == False)[0]
            ens_delta_time[idx_time[1:]] = nandiff(ens_time_sec[idx_time])
            
            # Adjust for transects tha last past midnight
            idx_24hr = np.where(ens_delta_time < 0)[0]
            ens_delta_time[idx_24hr] = 24 * 3600 + ens_delta_time[idx_24hr]
            ens_delta_time = ens_delta_time.T
            
            # Start date and time
            idx = np.where(np.isnan(pd0.Sensor.time[:, 0]) == False)[0][0]
            start_year = int(pd0.Sensor.date[idx, 0])
            
            # StreamPro doesn't include y2k dates
            if start_year < 100:
                start_year = 2000 + int(pd0.Sensor.date_not_y2k[idx, 0])
                
            start_month = int(pd0.Sensor.date[idx, 1])
            start_day = int(pd0.Sensor.date[idx, 2])
            start_hour = int(pd0.Sensor.time[idx, 0])
            start_min = int(pd0.Sensor.time[idx, 1])
            start_sec = int(pd0.Sensor.time[idx, 2] + pd0.Sensor.time[idx, 3] / 100)
            
            start_dt = datetime(start_year, start_month, start_day, start_hour, start_min, start_sec)
            start_serial_time = time.mktime(start_dt.timetuple())
            start_date = time.strftime('%m/%d/%Y', time.gmtime(start_serial_time))
            
            # End data and time
            idx = np.where(np.isnan(pd0.Sensor.time[:, 0]) == False)[0][-1]
            end_year = int(pd0.Sensor.date[idx, 0])
            # StreamPro does not include Y@K dates
            if end_year < 100:
                end_year = 2000 + int(pd0.Sensor.date_not_y2k[idx, 0])
                
            end_month = int(pd0.Sensor.date[idx, 1])
            end_day = int(pd0.Sensor.date[idx, 2])
            end_hour = int(pd0.Sensor.time[idx, 0])
            end_min = int(pd0.Sensor.time[idx, 1])
            end_sec = int(pd0.Sensor.time[idx, 2] + pd0.Sensor.time[idx, 3] / 100)
            
            end_dt = datetime(end_year, end_month, end_day, end_hour, end_min, end_sec)
            end_serial_time = time.mktime(end_dt.timetuple())
            
            # Create date/time object
            self.date_time = DateTime()
            self.date_time.populate_data(date_in=start_date,
                                         start_in=start_serial_time,
                                         end_in=end_serial_time,
                                         ens_dur_in=ens_delta_time)
            
            # Transect checked for use in discharge computation
            self.checked = mmt_transect.Checked

            # STOPPED HERE This has to do with moving-boat transects

            if kargs is None:
                self.adcp = InstrumentData()
                self.adcp.populate_data('TRDI', kargs=[mmt_transect, pd0, mmt])
            else:
                self.adcp = InstrumentData()
                self.adcp.populate_data('TRDI', kargs= np.hstack([[mmt_transect, pd0, mmt], kargs]))
            
    def SonTek(self, rsdata, file_name):
        """Reads Matlab file produced by RiverSurveyor Live and populates the transect instance variables.

        Parameters
        ----------
        rsdata: object
            Object of Matlab data from SonTek Matlab files
        file_name: str
            Name of SonTek Matlab file not including path.
        """

        self.file_name = file_name

        # ADCP instrument information
        # ---------------------------
        self.adcp = InstrumentData()
        self.adcp.populate_data('SonTek', rsdata)

        # Depth
        # -----

        # Initialize depth data structure
        self.depths = DepthStructure()

        # Determine array rows and cols
        max_cells = rsdata.WaterTrack.Velocity.shape[0]
        num_ens = rsdata.WaterTrack.Velocity.shape[2]

        # Compute cell sizes and depths
        cell_size = rsdata.System.Cell_Size.reshape(1, num_ens)
        cell_size_all = np.tile(cell_size, (max_cells, 1))
        top_of_cells = rsdata.System.Cell_Start.reshape(1, num_ens)
        cell_depth = (np.tile(np.arange(1, max_cells+1, 1).reshape(max_cells, 1), (1, num_ens)) * 0.5 * cell_size_all) + np.tile(top_of_cells, (max_cells,1))

        # Prepare bottom track depth variable
        depth = rsdata.BottomTrack.BT_Beam_Depth.reshape(4, num_ens)
        depth[depth == 0] = np.nan

        # Create depth object for bottom track beams
        self.depths.add_depth_object(depth_in=depth,
                                     source_in='BT',
                                     freq_in=rsdata.BottomTrack.BT_Frequency,
                                     draft_in=rsdata.Setup.sensorDepth,
                                     cell_depth_in=cell_depth,
                                     cell_size_in=cell_size_all)
        # Prepare vertical beam depth variable
        depth = rsdata.BottomTrack.VB_Depth
        depth[depth == 0] = np.nan

        # Create depth object for vertical beam
        self.depths.add_depth_object(depth_in=depth,
                                     source_in='VB',
                                     freq_in=np.array([rsdata.Transformation_Matrices.Frequency[1]] * depth.shape[-1]),
                                     draft_in=rsdata.Setup.sensorDepth,
                                     cell_depth_in=cell_depth,
                                     cell_size_in=cell_size_all)

        # Set depth reference
        if rsdata.Setup.depthReference < 0.5:
            self.depths.set_depth_reference('VB')
        else:
            self.depths.set_depth_reference('BT')

        # Water Velocity
        # --------------

        # Rearrange arrays for consistency with WaterData class
        vel = np.swapaxes(rsdata.WaterTrack.Velocity, 1, 0)
        snr = np.swapaxes(rsdata.System.SNR, 1, 0)
        corr = np.swapaxes(rsdata.WaterTrack.Correlation, 1, 0)


        # Correct SonTek difference velocity for error in earlier transformation matrices.
        if abs(rsdata.Transformation_Matrices.Matrix[3, 0, 0]) < 0.5:
            vel[3, :, :] = vel[3, :, :] * 2

        # Apply TRDI scaling to SonTek difference velocity to convert to a TRDI compatible error velocity
        vel[3, :, :] = vel[3, :, :] / ((2**0.5) * np.tan(np.deg2rad(25)))

        # Convert velocity reference from what was used in RiverSurveyor Live to None by adding the boat velocity
        # to the reported water velocity
        boat_vel = np.swapaxes(rsdata.Summary.Boat_Vel, 1, 0)
        vel[0, :, :] = vel[0, :, :] + boat_vel[0, :]
        # Because Matlab pads arrays with zeros and RR data has variable
        # number of bins, the raw data may be padded with zeros.  The next
        # four statements changes those to nan.
        vel[vel==0] = np.nan
        ref_water = 'None'

        # The initial coordinate system must be set to earth for early versions of RiverSurveyor firmware.
        # This implementation forces all versions to use the earth coordinate system.
        if rsdata.Setup.coordinateSystem == 0:
            ref_coord = 'Beam'
            raise ValueError('Beam Coordinates are not supported for all RiverSuveyor firmware releases, use Earth coordinates.')
        elif rsdata.Setup.coordinateSystem == 1:
            ref_coord = 'Inst'
            raise ValueError('Instrument Coordinates are not supported for all RiverSuveyor firmware releases, use Earth coordinates.')
        elif rsdata.Setup.coordinateSystem == 2:
            ref_coord = 'Earth'

        # Compute side lobe cutoff using Transmit Length information if availalbe, if not it is assumed to be equal
        # to 1/2 depth_cell_size_m. The percent method is use for the side lobe cutoff computation.
        sl_cutoff_percent = rsdata.Setup.extrapolation_dDiscardPercent
        sl_cutoff_number = rsdata.Setup.extrapolation_nDiscardCells
        if 'Transmit_Length' in set(rsdata.Summary._fieldnames):
            sl_lag_effect_m = (rsdata.Summary.Transmit_Length
                               + self.depths.bt_depths.depth_cell_size_m[0, :]) / 2.0
        else:
            sl_lag_effect_m = np.copy(self.depths.bt_depths.depth_cell_depth_m[0, :])
        sl_cutoff_type = 'Percent'
        cells_above_sl = TransectData.side_lobe_cutoff(depths=self.depths.bt_depths.depth_orig_m,
                              draft=self.depths.bt_depths.draft_orig_m,
                              cell_depth=self.depths.bt_depths.depth_cell_depth_m,
                              sl_lag_effect=sl_lag_effect_m,
                              type=sl_cutoff_type,
                              value=1 - sl_cutoff_percent / 100)
        # Determine water mode
        corr_nan = np.isnan(corr)
        number_of_nan = np.count_nonzero(corr_nan)
        if number_of_nan == 0:
            wm = 'HD'
        elif corr_nan.size == number_of_nan:
            wm = 'IC'
        else:
            wm = 'Variable'

        # Determine excluded distance (Similar to SonTek's screening distance)
        excluded_distance = rsdata.Setup.screeningDistance - rsdata.Setup.sensorDepth
        if excluded_distance < 0:
            excluded_distance = 0

        # Create water velocity object
        self.w_vel = WaterData()
        self.w_vel.populate_data(vel_in=vel,
                                freq_in=rsdata.WaterTrack.WT_Frequency,
                                coord_sys_in=ref_coord,
                                nav_ref_in=ref_water,
                                rssi_in=snr,
                                rssi_units_in=rsdata.System.Units.SNR,
                                excluded_dist_in=excluded_distance,
                                cells_above_sl_in=cells_above_sl,
                                sl_cutoff_per_in=sl_cutoff_percent,
                                sl_cutoff_num_in=sl_cutoff_number,
                                sl_cutoff_type_in=sl_cutoff_type,
                                sl_lag_effect_in=sl_lag_effect_m,
                                wm_in=wm,
                                blank_in=excluded_distance,
                                corr_in=corr)

        # Bottom Track
        # ------------
        self.boat_vel = BoatStructure()
        self.boat_vel.add_boat_object(source='SonTek',
                                      vel_in=np.swapaxes(rsdata.BottomTrack.BT_Vel, 1, 0),
                                      freq_in=rsdata.BottomTrack.BT_Frequency,
                                      coord_sys_in=ref_coord,
                                      nav_ref_in='BT')

        # GPS Data
        # --------
        self.gps = GPSData()
        if np.nansum(rsdata.GPS.GPS_Quality) > 0:
            self.gps.populate_data(raw_gga_utc=rsdata.RawGPSData.GgaUTC,
                                   raw_gga_lat=rsdata.RawGPSData.GgaLatitude,
                                   raw_gga_lon=rsdata.RawGPSData.GgaLongitude,
                                   raw_gga_alt=rsdata.RawGPSData.GgaAltitude,
                                   raw_gga_diff=rsdata.RawGPSData.GgaQuality,
                                   raw_gga_hdop=np.swapaxes(np.tile(rsdata.GPS.HDOP,(rsdata.RawGPSData.GgaLatitude.shape[1], 1)), 1, 0),
                                   raw_gga_num_sats=np.swapaxes(np.tile(rsdata.GPS.Satellites,(rsdata.RawGPSData.GgaLatitude.shape[1], 1)), 1, 0),
                                   raw_gga_delta_time=None,
                                   raw_vtg_course=rsdata.RawGPSData.VtgTmgTrue,
                                   raw_vtg_speed=rsdata.RawGPSData.VtgSogMPS,
                                   raw_vtg_delta_time=None,
                                   raw_vtg_mode_indicator=rsdata.RawGPSData.VtgMode,
                                   ext_gga_utc=rsdata.GPS.Utc,
                                   ext_gga_lat=rsdata.GPS.Latitude,
                                   ext_gga_lon=rsdata.GPS.Longitude,
                                   ext_gga_alt=rsdata.GPS.Altitude,
                                   ext_gga_diff=rsdata.GPS.GPS_Quality,
                                   ext_gga_hdop=rsdata.GPS.HDOP,
                                   ext_gga_num_sats=rsdata.GPS.Satellites,
                                   ext_vtg_course=np.tile(np.nan, rsdata.GPS.Latitude.shape),
                                   ext_vtg_speed=np.tile(np.nan, rsdata.GPS.Latitude.shape),
                                   gga_p_method='End',
                                   gga_v_method='End',
                                   vtg_method='Average')

            self.boat_vel.add_boat_object(source='SonTek',
                                          vel_in=self.gps.gga_velocity_ens_mps,
                                          freq_in=None,
                                          coord_sys_in='Earth',
                                          nav_ref_in='GGA')

            self.boat_vel.add_boat_object(source='SonTek',
                                          vel_in=self.gps.vtg_velocity_ens_mps,
                                          freq_in=None,
                                          coord_sys_in='Earth',
                                          nav_ref_in='VTG')
        ref = None
        if rsdata.Setup.trackReference == 1:
            ref = 'BT'
        elif rsdata.Setup.trackReference == 2:
            ref = 'GGA'
        elif rsdata.Setup.trackReference == 3:
            ref = 'VTG'
        self.boat_vel.set_nav_reference(ref)

        # Edges
        # -----
        # Create edge object
        self.edges = Edges()
        self.edges.populate_data(rec_edge_method='Variable',
                                 vel_method='VectorProf')

        # Determine number of ensembles for each edge
        if rsdata.Setup.startEdge > 0.1:
            ensembles_right = np.nansum(rsdata.System.Step == 2)
            ensembles_left = np.nansum(rsdata.System.Step == 4)
            self.in_transect_idx = np.arange(ensembles_right + 1, num_ens - ensembles_left, 1)
            self.start_edge = 'Right'
        else:
            ensembles_right = np.nansum(rsdata.System.Step == 4)
            ensembles_left = np.nansum(rsdata.System.Step == 1)
            self.in_transect_idx = np.arange(ensembles_left + 1, num_ens - ensembles_right, 1)
            self.start_edge = 'Left'

        # Create left edge object
        edge_type = None
        if rsdata.Setup.Edges_0__Method == 2:
            edge_type = 'Triangular'
        elif rsdata.Setup.Edges_0__Method == 1:
            edge_type = 'Rectangular'
        elif rsdata.Setup.Edges_0__Method == 0:
            edge_type = 'UserQ'
        self.edges.left.populate_data(edge_type=edge_type,
                                      distance=rsdata.Setup.Edges_0__DistanceToBank,
                                      number_ensembles=ensembles_left,
                                      coefficient=None,
                                      user_discharge=rsdata.Setup.Edges_0__EstimatedQ)
        # Create right edge object
        if rsdata.Setup.Edges_1__Method == 2:
            edge_type = 'Triangular'
        elif rsdata.Setup.Edges_1__Method == 1:
            edge_type = 'Rectangular'
        elif rsdata.Setup.Edges_1__Method == 0:
            edge_type = 'UserQ'
        self.edges.right.populate_data(edge_type=edge_type,
                                      distance=rsdata.Setup.Edges_1__DistanceToBank,
                                      number_ensembles=ensembles_left,
                                      coefficient=None,
                                      user_discharge=rsdata.Setup.Edges_1__EstimatedQ)

        # Extrapolation
        # -------------

        # Top extrapolation
        if rsdata.Setup.extrapolation_Top_nFitType == 0:
            top = 'Constant'
        elif rsdata.Setup.extrapolation_Top_nFitType == 1:
            top = 'Power'
        elif rsdata.Setup.extrapolation_Top_nFitType == 2:
            top = '3-Point'
        # Bottom extrapolation
        if rsdata.Setup.extrapolation_Bottom_nFitType == 0:
            bottom = 'Constant'
        elif rsdata.Setup.extrapolation_Bottom_nFitType == 1:
            if rsdata.Setup.extrapolation_Bottom_nEntirePro > 1.1:
                bottom = 'No Slip'
            else:
                bottom = 'Power'

        # Create extrapolation object
        self.extrap = ExtrapData()
        self.extrap.populate_data(top=top,
                                  bot=bottom,
                                  exp=rsdata.Setup.extrapolation_Bottom_dExponent)

        # Sensor data
        # -----------
        self.sensors = Sensors()

        # Internal heading
        self.sensors.heading_deg.internal = HeadingData()
        # Check for firmware supporting G3 compass and associated data
        if hasattr(rsdata,'Compass'):
            # TODO need to find older file that had 3 columns in Magnetic error to test and modify code
            mag_error = rsdata.Compass.Magnetic_error
            pitch_limit = (rsdata.Compass.Maximum_Pitch, rsdata.Compass.Minimum_Pitch)
            roll_limit = (rsdata.Compass.Maximum_Roll, rsdata.Compass.Minimum_Roll)
        else:
            mag_error = None
            pitch_limit = None
            roll_limit = None
        self.sensors.heading_deg.internal.populate_data(data_in=rsdata.System.Heading,
                                                        source_in='Internal',
                                                        magvar=rsdata.Setup.magneticDeclination,
                                                        mag_error=mag_error,
                                                        pitch_limit=pitch_limit,
                                                        roll_limit=roll_limit)

        # External heading
        ext_heading = rsdata.System.GPS_Compass_Heading
        if np.nansum(np.abs(np.diff(ext_heading))) > 0:
            self.sensors.heading_deg.external.populate_data(data_in=ext_heading,
                                                            source_in='GPS',
                                                            magvar=rsdata.Setup.magneticDeclination,
                                                            align=rsdata.Setup.hdtHeadingCorrection)

        # Set selected reference
        if rsdata.Setup.headingSource > 1.1:
            self.sensors.heading_deg.selected = 'external'
        else:
            self.sensors.heading_deg.selected = 'internal'

        # Pitch and roll
        if hasattr(rsdata, 'Compass'):
            pitch = rsdata.Compass.Pitch
            roll = rsdata.Compass.Roll
        elif hasattr(rsdata.System, 'Pitch'):
            pitch = rsdata.System.Pitch
            roll = rsdata.system.Roll
        self.sensors.pitch_deg.internal = SensorData()
        self.sensors.pitch_deg.internal.populate_data(data_in=pitch, source_in='internal')
        self.sensors.pitch_deg.selected = 'internal'
        self.sensors.roll_deg.internal = SensorData()
        self.sensors.roll_deg.internal.populate_data(data_in=roll, source_in='internal')
        self.sensors.roll_deg.selected = 'internal'

        # Temperature
        if rsdata.System.Units.Temperature == 'degC':
            temperature = rsdata.System.Temperature
        else:
            temperature = (5. / 9.) * (rsdata.Temperature - 32)
        self.sensors.temperature_deg_c.internal = SensorData()
        self.sensors.temperature_deg_c.internal.populate_data(data_in=temperature, source_in='internal')
        self.sensors.temperature_deg_c.selected = 'internal'

        # Salinity
        self.sensors.salinity_ppt.user = SensorData()
        self.sensors.salinity_ppt.user.populate_data(data_in=rsdata.Setup.userSalinity, source_in='Manual')
        self.sensors.salinity_ppt.selected = 'user'
        # Matlab notes indicated that an internal sensor needed to be created for compatibility with
        # future computations
        self.sensors.salinity_ppt.internal = SensorData()
        self.sensors.salinity_ppt.internal.populate_data(data_in=rsdata.Setup.userSalinity, source_in='Manual')

        # Speed of sound
        # Not provided in SonTek data but is computed from equation used in TRDI BBSS.
        speed_of_sound = Sensors.speed_of_sound(temperature=temperature, salinity=rsdata.Setup.userSalinity)
        self.sensors.speed_of_sound_mps.internal = SensorData()
        self.sensors.speed_of_sound_mps.internal.populate_data(data_in=speed_of_sound, source_in='QRev')

        # Ensemble times
        ensemble_delta_time = np.append([0], np.diff(rsdata.System.Time))
        idx_missing = np.where(ensemble_delta_time > 1.5)
        if idx_missing[0]:
            number_missing = np.sum(ensemble_delta_time[idx_missing]) - len(idx_missing)
            error_str = self.file_name + ' is missing ' + str(number_missing) + ' samples'
            raise ValueError(error_str)

        # Date, start, end, and duration
        start_serial_time = DateTime.time_2_serial_time(time_in=rsdata.System.Time[0], source='SonTek')
        end_serial_time = DateTime.time_2_serial_time(time_in=rsdata.System.Time[-1], source='SonTek')
        meas_date = time.strftime('%m/%d/%Y', time.gmtime(start_serial_time))
        self.date_time = DateTime()
        self.date_time.populate_data(date_in=meas_date,
                                     start_in=start_serial_time,
                                     end_in=end_serial_time,
                                     ens_dur_in=ensemble_delta_time)

        # Transect checked for use in discharge computations
        self.checked = True
        # TODO refactor composite depths and set to on for SonTek data
        # Set composite depths as this is the only option in RiverSurveyor Live
        # self.depths.composite_depths('On')

    @staticmethod
    def compute_instrument_cell_data(pd0):
        
        # Number of ensembles
        num_ens = np.array(pd0.Wt.vel_mps).shape[-1]
        # Retrieve and compute cell information
        reg_cell_size = pd0.Cfg.ws_cm / 100
        reg_cell_size[reg_cell_size == 0] = np.nan
        dist_cell_m = pd0.Cfg.dist_bin1_cm / 100
        num_reg_cells = pd0.Wt.vel_mps.shape[1]
        
        #surf data are to accommodate RiverRay and RiverPro.  pd0_read sets these
        #values to nan when reading Rio Grande or StreamPro data
        
        no_surf_cells = pd0.Surface.no_cells
        no_surf_cells[np.isnan(no_surf_cells)] = 0
        max_surf_cells = np.nanmax(no_surf_cells)
        surf_cell_size = pd0.Surface.cell_size_cm / 100
        surf_cell_dist = pd0.Surface.dist_bin1_cm / 100
        
        # Compute maximum number of cells
        max_cells = int(max_surf_cells+num_reg_cells)
        
        # Combine cell size and cell range from transducer for both
        # surface and regular cells
        
        cell_depth = repmat([np.nan], max_cells, num_ens)
        cell_size_all = repmat([np.nan], max_cells, num_ens)
        
        for i in range(num_ens):
            
            #Determine number of cells to be treated as regular cells
            if np.nanmax(no_surf_cells) > 0:
                
                num_reg_cells = max_cells - no_surf_cells[i]
            else:
                num_reg_cells = max_cells
                
                
            #Surface cell are present
            if no_surf_cells[i] > 1e-5:
                cell_depth[:int(no_surf_cells[i]),i] = surf_cell_dist[i] + np.arange(0,(no_surf_cells[i] - 1)*surf_cell_size[i]+.001,surf_cell_size[i])
                cell_depth[int(no_surf_cells[i]):,i] = cell_depth[int(no_surf_cells[i]),i] \
                + (.5*surf_cell_size[i]+0.5*reg_cell_size[i]) \
                + np.arange(0, (num_reg_cells-1)*reg_cell_size[i]+0.001, reg_cell_size[i])
                cell_size_all[0:int(no_surf_cells[i]),i] = np.repeat(surf_cell_size[i],int(no_surf_cells[i]))
                cell_size_all[int(no_surf_cells[i]):,i] = np.repeat(reg_cell_size[i],int(num_reg_cells))
            else:
                
                cell_depth[:int(num_reg_cells),i] = dist_cell_m[i] + np.arange(0,(num_reg_cells)*reg_cell_size[i],reg_cell_size[i])
                cell_size_all[:,i] = np.repeat(reg_cell_size[i], num_reg_cells)
                
                
            #Firmware is used to ID RiverRay data with variable modes and lags
            firmware = str(pd0.Inst.firm_ver[0])
            
        #Compute sl_lag_effect
        lag = pd0.Cfg.lag_cm / 100
        if firmware[0:2] == '44' or firmware[0:2] == '56':
            lag_near_bottom = np.array(pd0.Cfg.lag_near_bottom)
            lag_near_bottom[lag_near_bottom == np.nan] = 0
            lag[lag_near_bottom != 0] = 0
            
        pulse_len = pd0.Cfg.xmit_pulse_cm / 100
        sl_lag_effect_m = (lag + pulse_len + reg_cell_size) / 2
        sl_cutoff_per = (1-(cosd(pd0.Inst.beam_ang[0]))) * 100
            
        return (cell_size_all, cell_depth, sl_cutoff_per, sl_lag_effect_m)
        
    def adjust_side_lobe(self):
        """Coordinates side lobe cutoff calls"""
        
        self.w_vel.adjust_side_lobe(self)
        
    def change_selection(self):
        """Changes whether the transect is is include in the final discharge"""
        
        self.checked = not self.checked
        
    def set_extrapolation(self,top,bot,exp):
        
        self.extrap.set_extrap_data(top, bot, exp)
        
    def change_q_ensembles(self, proc_method):
        """Sets in_transect_idx to all ensembles, except in the case of SonTek data
        where RSL processing is applied.
        
        Parameters
        ----------
        proc_method: str
            Processing method (WR2, RSL, QRev)
        """
        
        if proc_method == 'RSL':
            num_ens = self.boat_vel.bt_vel.u_processed_mps.shape[1]
            # Determine number of ensembles for each edge
            if self.start_edge == 'Right':
                self.in_transect_idx = np.arange(self.edges.right.num_ens_2_avg, num_ens-self.edges.left.num_ens_2_avg)
            else:
                self.in_transect_idx = np.arange(self.edges.left.num_ens_2_avg, num_ens-self.edges.right.num_ens_2_avg)
        else:
            self.in_transect_idx = np.arange(0, self.boat_vel.bt_vel.u_processed_mps.shape[0])
            
    def change_start_edge(self, setting):
        """start edge"""
        self.start_edge = setting
        
    def change_edge(self, edge, property, setting):
        """Change an edge property (velocity or coefficient)
        
        Input:
        edge: edge object to change (left or right)
        property: property to change (rec_edge_method or vel_method)
        setting: new property setting
        """
        
        self.edges.change_property(property, setting, edge)    
        
    def change_coord_sys(self, new_coord_sys):
        """Changes the coordinate system of the water and boat data.

        Current implementation only allows changes for original to higher order coordinate
        systems: Beam - Inst - Ship - Earth.
        
        Parameters
        ----------
        new_coord_sys: str
            Name of new coordinate system (Beam, Int, Ship, Earth)
        """
        self.w_vel.change_coord_sys(new_coord_sys, self.sensors, self.adcp)
        self.boat_vel.change_coord_sys(new_coord_sys, self.sensors, self.adcp)
        
    def change_nav_reference(self, update, new_nav_ref):
        """Method to set the navigation reference for the water data.
        
        Parameters
        ----------
        update: bool
            Setting to determine if water data should be updated.
        new_nav_ref: str
            New navigation reference (bt_vel, gga_vel, vtg_vel)
        """
        
        self.boat_vel.change_nav_reference(reference=new_nav_ref, transect=self)
        
        if update:
            self.update_water()
            
    def change_mag_var(self, mag_var):
        """Change magnetic variation
        
        Input:
        mag_var: magnetic variation in degrees
        """
        
        #update object
        if (self.sensorsdata.heading_deg.external is not None):
            self.sensorsdata.heading_deg.set_mag_var(mag_var, 'external')
        
        if self.sensor.heading_deg.selected == 'internal':
            old = getattr(self.sensorsdata.heading_deg, self.sensorsdata.heading_deg.selected)
            old_mag_var = old.mag_var
            mag_var_change = mag_var - old_mag_var
            self.sensorsdata.heading_deg.set_mag_var(mag_var, 'internal')
            self.boat_vel.chang_mag_var(mag_var_change)
            self.w_vel.change_mag_var(self.boat_vel,mag_var_change)
        
        self.sensorsdata.heading_deg.set_mag_var(mag_var, 'internal')
        
        self.update_water()
        
    def change_offset(self, h_offset):
        """Change the heading offset (aignment correction)
        
        Input:
        h_offset: heading offset in degrees
        """
        self.sensorsdata.heading_deg.set_align_correction(h_offset, 'internal')
        
        if self.sensorsdata.heading_deg.selected == 'external':
            old = getattr(self.sensorsdata.heading_deg, self.sensorsdata.heading_deg.selected)
            old_offset = old.align_correction
            offset_change = h_offset - old_offset
            self.boat_vel.bt_vel.change_offset(offset_change)
            self.w_vel.change_offset(self.boat_vel, offset_change)
        
        self.sensorsdata.heading_deg.set_align_correction(h_offset, 'external')
        
        self.update_water()
        
    def change_heading_source(self, h_source):
        source = getattr(self.sensorsdata.heading_deg, h_source)
        if source is not None:
            old = getattr(self.sensorsdata.heading_deg, self.sensorsdata.heading_deg.selected)
            old_heading = old.data
            new_heading = source.data
            heading_change = new_heading - old_heading
            self.sensorsdata.heading_deg.set_selected(h_source)
            self.boat_vel.bt_vel.change_heading_source(heading_change)
            self.w_vel.change_heading_source(self.boat_vel, heading_change)
            
        self.update_water()
            
    def update_water(self):
        """Method called from set_nav_reference, boat_interpolation and boat filters
        to ensure that changes in boatvel are reflected in the water data"""

        self.w_vel.set_nav_reference(self.boat_vel)
        
        # Reapply water filters and interpolations
        # Note wt_filters calls apply_filter which automatically calls
        # apply_interpolation so both filters and interpolations
        # are applied with this one call
        
        self.wt_filters()
         
    # def wt_filters(self, type=None, setting=None, threshold=None):
    def wt_filters(self, **kwargs):

        """Coordinate application of water velocity filters.
        
        Parameters
        ----------
        kwargs
            beam:
                Setting for beam filter (Auto, Off, threshold value)
            difference:
                Setting for difference filter (Auto, Off, threshold value)
            vertical:
                Setting for vertical filter (Auto, Off, threshold value)
            other:
                Setting for other filters (Off, On)
            excluded:
                Excluded distance below the transducer, in m
        """
        
        self.w_vel.apply_filter(self, **kwargs)

    def wt_interpolations(self, **kwargs):
        """Coordinate water velocity interpolation
        
        Input:
        kargs: specifies whether ensembles or cells are to be interpolated and the method
        used for interpolation.  Both setting can be included in the same call.
        
        Example:
        kargs[0]: 'Ensembles'
        kargs[1]: 'Linear'
        kargs[2]: 'Cells'
        kargs[3]: 'TRDI'
        """
        
        #Process transect
        self.w_vel.apply_interpolation(self, **kwargs)
                
    # DSM changed 1/25/2018 def side_lobe_cutoff(self, depths, draft, cell_depth, sl_lag_effect, kargs):
    @staticmethod
    def side_lobe_cutoff(depths, draft, cell_depth, sl_lag_effect, type='Percent', value=None):
        """Computes side lobe cutoff.

        The side lobe cutoff is based on the beam angle and is computed to
        ensure that the bin and any lag beyond the actual bin cutoff is
        above the side lobe cutoff.

        Parameters
        ----------
        depths: np.array
            Bottom track (all 4 beams) and vertical beam depths for each ensemble, in m.
        draft: float
            Draft of transducers, in m.
        cell_depth: np.array
            Depth to the centerline of each depth cell, in m.
        sl_lag_effect: np.array
            The extra depth below the last depth cell that must be above the side lobe cutoff, in m.
        type: str
            Method used for side lobe cutoff computation.
        value: float
            Value used in specified method to use for side lobe cutoff computation.
        """
        
     
        #Compute minimum depths for each ensemble
        min_depths = np.nanmin(depths,0)
        
        #Compute range from transducer
        range_from_xducer = min_depths - draft
        
        #Adjust fro transducter angle
        if type == 'Percent':
            coeff = value
        elif type == 'Angle':
            coeff = np.cos(np.deg2rad(value))
        
        #Compute sidelobe cutoff to centerline    
        cutoff = np.array(range_from_xducer * coeff - sl_lag_effect+draft)
        
        #Compute boolean side lobe cutoff matrix
        cells_above_sl = (cell_depth - cutoff) < 0
        return cells_above_sl
        
    def boat_interpolations(self, update, target, method=None):
        """Coordinates boat velocity interpolations.
        
        Parameters
        ----------
        target: str
            Boat velocity reference (BT or GPS)
        method: str
            Type of interpolation
        """

        # Interpolate bottom track data
        if target == 'BT':
            self.boat_vel.bt_vel.apply_interpolation(transect=self, interpolation_method=method)
            
        if target == 'GPS':
            # Interpolate GGA data
            vel = getattr(self.boat_vel, 'gga_vel')
            if vel is not None:
                self.boat_vel.gga_vel.apply_interpolation(transect=self, interpolation_method=method)
            # Interpolate VTG data
            vel = getattr(self.boat_vel, 'vtg_vel')
            if vel is not None:
                self.boat_vel.vtg_vel.apply_interpolation(transect=self, interpolation_method=method)

        # Apply composite tracks setting
        self.composite_tracks(update=False)
        
        # Update water to reflect changes in boat_vel
        if update:
            self.update_water()
            
    def composite_tracks(self, update, setting=None):
        """Coordinate application of composite tracks.
        
        Parameters
        ----------
        update: bool
            Setting to control if water data are updated.
        setting: bool
            Sets composite tracks (on or off).
        """
        
        # Determine if setting is specified
        if setting is None:
            # Process transect using saved setting
            self.boat_vel.composite_tracks(transect=self)
        else:
            # Process transect usin new setting
            self.boat_vel.composite_tracks(transect=self, setting=setting)
            
        # Update water data to reflect changes in boatvel
        if update:
            self.update_water()
            
    def boat_filters(self, update, **kargs):
        """Coordinates application of boat filters to bottom track data
        
        Input:
        kargs: specified in pairs or triplets, can be multiple groups
        kargs[n]: filter type (Beam, Difference, Vertical, Other
        kargs[n+1]: filter setting (Auto, Manual, Off)
        kargs[n+2]: Threshold if manual
        """
        
        #apply filter to transect
        self.boat_vel.bt_vel.apply_filter(self, **kargs)
        
        if self.boat_vel.selected == 'bt_vel' and update:
            self.update_water()
            
    def gps_filters(self, update, **kargs):
        """Coordinate filters for GPS based boat velocities
        
        Input:
        kargs[n]: Filter type (Differential, Altitude, HDOP, Other)
        kargs[n+1]: Filter setting (Auto, Manual, Off)
        kargs[n+2]: Threshold if manual
        """
        
        if self.boat_vel.gga_vel is not None:
            self.boat_vel.gga_vel.apply_gps_filter(self, **kargs)
        if self.boat_vel.vtg_vel is not None:
            self.boat_vel.vtg_vel.apply_gps_filter(self, **kargs)
            
        if (self.boat_vel.selected == 'VTG' or self.boat_vel.selected == 'GGA') and update == True:
            self.update_water()
            
    def set_depth_reference(self, update, setting):
        """Coordinates setting the depth reference.
        
        Parameters
        ----------
        update: bool
            Determines if associated data should be updated
        setting: str
            Depth reference (BT, VB, DS)
        """
        
        self.depths.set_depth_reference(setting)
        
        if update:
            self.process_depths(update)
            self.adjust_side_lobe()
            
    def apply_averaging_method(self, setting):
        """Method to apply the selected averaging method to the BT team depths to achieve a single
        average depth.  It is only applicable to the multiple beams used for BT, not VB or DS.
        
        Input:
        setting: averaging method (IDW, Simple)
        """
        
        self.depths.bt_depths.compute_avg_BT_depth(setting)
        
        self.process_depths(False)
            
    def process_depths(self, update=False, filter_method=None, interpolation_method=None, composite_setting=None,
                       avg_method=None, valid_method=None):
        """Method applies filter, composite, and interpolation settings to  depth objects
        so that all are updated using the same filter and interpolation settings.
        
        Parameters
        ----------
        update: bool
            Determines if water data should be updated.
        filter_method: str
            Filter method to be used (None, Smooth, TRDI).
        interpolation_method: str
            Interpolation method to be used (None, HoldLast, Smooth, Linear).
        composite_setting: bool
            Specifies use of composite depths.
        avg_method: str
            Defines averaging method: "Simple", "IDW", only applicable to bottom track.
        valid_method:
            Defines method to determine if depth is valid (QRev or TRDI).
        """
        
        # Get current settings
        depth_data = getattr(self.depths, self.depths.selected)
        if filter_method is None:
            filter_method = depth_data.filter_type

        if interpolation_method is None:
            interpolation_method = depth_data.interp_type

        if composite_setting is None:
            composite_setting = self.depths.composite

        if avg_method is None:
            avg_method = self.depths.bt_depths.avg_method

        if valid_method is None:
            valid_method = self.depths.bt_depths.valid_data_method

        self.depths.set_valid_data_method(valid_method)
        self.depths.bt_depths.set_avg_method(avg_method)
        self.depths.depth_filter(self, filter_method)
        self.depths.depth_interpolation(self, interpolation_method)
        self.depths.composite_depths(self, composite_setting)
        self.w_vel.adjust_side_lobe(self)
        
        if update:
            self.update_water()

    def change_draft(self, input):
        """Changes the draft for the specified transects and selected depth"""
        
        if self.depths.vb_depths is not None:
            self.depths.vb_depths.change_draft(input)
        if self.depths.bt_depths is not None:
            self.depths.bt_depths.change_draft(input)
            
    def set_property(self, property, setting):
        setattr(property, setting)
        
    def update_data_structure(self):
        #Salinity
        pd0_salinity = self.sensorsdata.salinity_ppt.user.data_orig
        pd0_salinity_src = self.sensorsdata.salinity_ppt.user.source
        
        #Create slainity sensor
        self.sensorsdata.add_sensor_data('salinity_ppt','internal',pd0_salinity,pd0_salinity_src)
        self.sensorsdata.set_selected('salinity_ppt', 'internal')

    def sos_user(self, kargs = None):
        """Compute new speed of sound from temperature and salinity

        Output:
        new_sos: newly computed speed of sound
        old_sos: previously used speed of sound
        """

        #Assign selected temperature data to local variable
        temp = getattr(self.sensorsdata.temperature_deg_c, self.sensorsdata.temperature_deg_c.selected)
        temperature = temp.data
        #Assign selected salinity to local variable
        sal = getattr(self.sensorsdata.salinity_ppt, self.sensorsdata.salinity_ppt.selected)
        salinity = sal.data
        old = getattr(self.sensorsdata.speed_of_sound_mps, self.sensorsdata.speed_of_sound_mps.selected)
        old_sos = old.data

        if self.sensorsdata.temperature_deg_c.selected == 'internal':
            new_sos = self.sensorsdata.speed_of_sound_mps.user.data_orig
            self.sensorsdata.speed_of_sound_mps.user.change_data(new_sos)
            self.sensorsdata.speed_of_sound_mps.user.set_source('Internal (ADCP)')
            self.sensorsdata.set_selected('speed_of_sound_mps', 'user')
        else:
            #Compute new speed of sound
            new_sos = Sensors().speed_of_sound(temperature, salinity)

            #Save new speed of sound to user sensor object with a source as computed
            if self.sensorsdata.speed_of_sound_mps.user is not None:
                self.sensorsdata.set_selected('speed_of_sound_mps', 'user')
                self.sensorsdata.speed_of_sound_mps.user.change_data(new_sos)
                self.sensorsdata.speed_of_sound_mps.user.set_source('Computed')
            else:
                self.sensorsdata.add_sensor_data('speed_of_sound_mps', 'user'. new_sos, 'Computed')
                self.sensorsdata.set_selected('speed_of_sound_mps', 'user')

        return (old_sos, new_sos)

    @staticmethod
    def raw_valid_data(transect):
        """Determines ensembles and cells with no interpolated water or boat data.

        For valid water track cells both non-interpolated valid water data and
        boat velocity data must be available. Interpolated depths are allowed.

        For valid ensembles water, boat, and depth data must all be non-interpolated.

        Parameters
        ----------
        transect: object
            Object of TransectData

        Returns
        -------
        raw_valid_ens: np.array(bool)
            Boolean array identifying raw valid ensembles.
        raw_valid_depth_cells: np.array(bool)
            Boolean array identifying raw valid depth cells.
        """

        in_transect_idx = transect.in_transect_idx

        # Determine valid water track ensembles based on water track and navigation data.
        boat_vel_select = getattr(transect.boat_vel, transect.boat_vel.selected)
        if np.nansum(boat_vel_select.u_processed_mps) > 0:
            valid_nav = boat_vel_select.valid_data[0,in_transect_idx]
        else:
            valid_nav = np.tile(False, in_transect_idx)

        valid_wt = transect.w_vel.valid_data[0, :, in_transect_idx]
        valid_wt_ens = np.any(valid_wt, 1)

        # Determine valid depths
        depths_select = getattr(transect.depths, transect.depths.selected)
        if transect.depths.composite:
            idx_na = np.where(depths_select.depth_source_ens[in_transect_idx] == 'NA')[0]
            valid_depth = not np.where(depths_select.depth_source_ens[in_transect_idx] == 'IN')[0]
            valid_depth[idx_na] = False
        else:
            valid_depth = depths_select.valid_data[in_transect_idx]
            idx = np.where(np.isnan(depths_select.depth_processed_m[in_transect_idx]))[0]
            valid_depth[idx] = False

        # Determine valid ensembles based on all data
        valid_ens = np.any(np.vstack((valid_nav, valid_wt_ens, valid_depth)))

        return valid_ens, valid_wt


# ========================================================================
# Begin multithread function included in module but not TransectData class
# Currently this is coded only for TRDI data
# ========================================================================

# DSM changed 1/23/2018 def allocate_transects(source, mmt, kargs)
def allocate_transects(source, mmt, type='Q', checked=False):
    #DEBUG, set threaded to false to get manual serial commands
    multi_threaded = True
    
    #Refactored from TransectData to iteratively create TransectData objects
        #----------------------------------------------------------------

    # Setup processing for discharge or moving-bed transects
    # DSM changed 1/23/2018 if kargs[0] == 'Q':
    if type == 'Q':
        # Identify discharge transect files to load
        transects = 'transects'
        active_config = 'active_config' 
        # DSM changed 1/23/2018 if kargs[1] == True:
        if checked == True:
            # DSM changed 1/23/2018 files_to_load = np.array([attribute.Checked for attribute in mmt.transects], dtype=bool)
            file_names = [attribute.Files[0].File for attribute in mmt.transects]
        else:
            # DSM changed 1/23/2018 files_to_load = np.array(np.ones(len(mmt.transects)), dtype=bool)
            file_names = [attribute.Files[0].File for attribute in mmt.transects]
    # DSM change 1/23/2018 elif kargs[0] == 'MB':
    elif type == 'MB':
        # Identify moving-bed transect files to load
        transects = 'mbt_transects'
        active_config = 'mbt_active_config'
        # DSM changed 1/23/2018 files_to_load =np.array([attribute.Checked for attribute in mmt.mbt_transects], dtype=bool)
        file_names = [attribute.Files[0].File for attribute in mmt.mbt_transects if attribute.Checked == 1]
        # DSM changed 1/23/2018 files_to_load_idx = np.where(files_to_load == True)[0]
    # DSM changed 1/23/2018 pathname = mmt.infile[:mmt.infile.rfind('/')]
    pathname = os.path.split(mmt.infile)[0]
    
    # Determine if any files are missing
    valid_files = []
    # DSM changed 1/23/2018
    # for x in file_names:
    #     x[0].Path = x[0].Path[x[0].Path.rfind('\\') + 1:]
    #     if os.path.exists(''.join([pathname,'/',x[0].Path])):
    #         valid_files.append((x, 1))
    #     else:
    #         valid_files.append((None, 0))
    for name in file_names:
        fullname=os.path.join(pathname, name)
        if os.path.exists(fullname):
            valid_files.append(fullname)

    # Multithread for Pd0 files
    # -------------------------
    # Seems like this section belongs in Pd0TRDI.py
    # Initialize thread variables
    pd0_data = []
    pd0_threads = []
    thread_id = 0

    # DSM 1/24/2018 could this be moved to Pd0TRDI.py as a method
    def add_pd0(file_name):
        pd0_data.append(Pd0TRDI(file_name))
        
    if multi_threaded == True:
        # DSM changed 1/24/2018 for x in valid_files:
        #     if x[1] == 1:
        #         pd0_thread = MultiThread(thread_id=thread_id, function=add_pd0, args = {'file_name': ''.join([pathname,'/',x[0][0].Path])})
        for file in valid_files:
            pd0_thread = MultiThread(thread_id=thread_id, function=add_pd0, args={'file_name': file})
            thread_id += 1
            pd0_thread.start()
            pd0_threads.append(pd0_thread)
    else:   
        # DSM changed 1/24/2018 pd0_data = [Pd0TRDI(''.join([pathname,'/',x[0][0].Path])) for x in valid_files if x[1] == 1]
        pd0_data = [Pd0TRDI(file for file in valid_files)]

    # DSM changed 1/24/2018 for x in pd0_threads:
    #     x.join()
    for thrd in pd0_threads:
        thrd.join()

    # Multithread for transect data
    # -----------------------------
    # Initialize thread variables
    processed_transects = []
    transect_threads = []
    thread_id = 0

    # DSM 1/24/2018 couldn't this be added to the TransectData class
    def add_transect(transect, source, in_file, pd0_data, mmt, mbt):
        transect.get_data(source, in_file,pd0_data, mmt, mbt)
        processed_transects.append(transect)
        
    # Process each transect
    for k in range(len(pd0_data)):
        
        transect = TransectData()
        transect.active_config = active_config

        # DSM 1/24/2018 Stopped here
        # The next line doesn't seem correct
        # the use of p_thread seems like that should have been used for PD0 and this should be t_thread
        # transect above should be empty so how do you get transect.active_config?
        transect.transects = transects
       
        mbt_idx = 0
        
        if active_config == 'mbt_active_config' or active_config == 'mbt_field_config':
            
            if multi_threaded:
                p_thread = MultiThread(thread_id = thread_id, function= add_transect, 
                                    args = {'transect': transect, 
                                            'source':'TRDI', 
                                            'in_file': mmt.mbt_transects[k], 
                                            'pd0_data': pd0_data[k], 
                                            'mmt': mmt,
                                            'mbt': mbt_idx})
                mbt_idx += 1
                p_thread.start()
                transect_threads.append(p_thread)
                
                
            else:
                add_transect(transect, 'TRDI', mmt.mbt_transects[k], pd0_data[k], mmt)
            
        else:
            if multi_threaded:
                p_thread = MultiThread(thread_id = thread_id, function= add_transect, 
                                       args = {'transect': transect, 
                                               'source':'TRDI', 
                                               'in_file': mmt.transects[k], 
                                               'pd0_data': pd0_data[k], 
                                               'mmt': mmt,
                                               'mbt': False})
                p_thread.start()
                transect_threads.append(p_thread)
                
                
            else:
                add_transect(transect, 'TRDI', mmt.transects[k], pd0_data[k], mmt)
    
    if multi_threaded:
        for x in transect_threads:
                x.join()
                
    return processed_transects   
        
def adjusted_ensemble_duration(transect, kargs=None):
        
    if transect.adcp.manufacturer == 'TRDI':
        if kargs is None:
            valid = np.isnan(transect.w_vel.u_processed_mps) == False
            valid_sum = np.sum(valid)
        else:
            valid_sum = np.isnan(transect.boat_vel.bt_vel.u_processed_mps) == False
            
        valid_ens = valid_sum > 0
        n_ens = len(valid_ens)
        ens_dur = transect.date_time.ens_duration_sec
        delta_t = np.tile([np.nan], n_ens)
        cum_dur = 0
        for j in range(n_ens):
            cum_dur = np.nansum(np.hstack([cum_dur, ens_dur[j]]))
            if valid_ens[j]:
                delta_t[j] = cum_dur
                cum_dur = 0
    else:
        delta_t = transect.date_time.ens_duration_sec
        
    return delta_t
    
                
    
            
            
          
           
                
            
        
            
            
        
        
            
        
                
        


