import numpy as np
from numpy.matlib import repmat
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
from Classes.DateTime import DateTime
from Classes.InstrumentData import InstrumentData
from Classes.MultiThread import MultiThread
import matplotlib.dates as mdates
from datetime import datetime
from MiscLibs.convenience import nandiff


class TransectData(object):
    """Class to hold Transect properties (may be removed on a refactor)
        Analogous to matlab class: clsTransectData
    """

    def __init__(self):
        self.adcp = None            # object of clsInstrument
        self.file_name = None      # filename of transect data file
        self.w_vel = None           # object of clsWaterData
        self.boat_vel = None        # class for various boat velocity references
                                    # (btVel, ggaVel, vtgVel)
        self.gps = None             # object of clsGPSData
        self.sensors = None         # object of clsSensorData
        self.depths = None      # object of clsDepthStructure for depth data including cell depths & ref depths
                                    # (btDepths, vbDepths, dsDepths)
        self.edges = None           # object of clsEdges
                                    # (left and right object of clsEdgeData)
        self.extrap = None          # object of clsExtrapData
        self.start_edge = None      # starting edge of transect looking downstream (Left or Right)
        self.datetime = None
        self.checked = None         #transect was checked for use in mmt file assumed checked for SonTek
        self.in_transect_idx = None # index of ensemble data associated with the moving-boat portion of the transect
        self.active_config = None
        self.transects = None
        self.cells_above_sl = None

    def get_data(self, source, in_file, pd0_data, mmt):
        
        if source == 'TRDI':
            self.TRDI(in_file, pd0_data, mmt)
            
                                #files2load idx
    def TRDI(self, mmt_transect, pd0_data, mmt, kargs=None):
       
       
        pd0 = pd0_data
        
        #get the configuration prooperty of the mmt_transect
        mmt_config = getattr(mmt_transect, self.active_config)
        if pd0_data.Wt is not None:
            
            #Get and compute ensemble beam depths
            temp_depth = np.array(pd0.Bt.depth_m)
            #Screen out invalid depths
            temp_depth[temp_depth<0.01] = np.nan
            #Add draft
            temp_depth += mmt_config['Offsets_Transducer_Depth']
            
            #Get instrument cell data
            cell_size_all_m, cell_depth_m, sl_cutoff_per, sl_lag_effect = self.compute_instrument_cell_data(pd0)
            
            #Adjust cell depth of draft 
            cell_depth_m = np.add(mmt_config['Offsets_Transducer_Depth'], cell_depth_m)
            
            #create depth data object for BT
            self.depths = DepthStructure()
            self.depths.add_depth_object(temp_depth, 'BT', pd0.Inst.freq, 
                                   mmt_config['Offsets_Transducer_Depth'], 
                                   kargs = [cell_depth_m, cell_size_all_m])
            
            #compute cells above side lobe
            self.side_lobe_cutoff(self.depths.bt_depths.depth_orig_m, 
                                                     self.depths.bt_depths.draft_orig_m, 
                                                     self.depths.bt_depths.depth_cell_depth_m,
                                                     sl_lag_effect,
                                                    kargs = ['Percent', 
                                                    1-sl_cutoff_per / 100])
            
            #Check for the presence of vertical beam data
            if np.nanmax(np.nanmax(pd0.Sensor.vert_beam_status)) > 0:
                
                temp_depth = pd0.Sensor.vert_beam_range_m
                
                #Screen out invalid depths
                temp_depth[temp_depth<0.01] = np.nan
                
                #Add draft
                temp_depth = temp_depth + mmt_config['Offsets_Transducer_Depth']
                
                #Create depth data object for vertical beam
                self.depths.add_depth_object(temp_depth, 'VB', pd0.Inst.freq, 
                                   mmt_config['Offsets_Transducer_Depth'],
                                   kargs = [cell_depth_m, cell_size_all_m])
                                   
            #check for the presence of depth sounder
            if np.nansum(np.nansum(pd0.Gps2.depth_m)) > 1e-5:
                temp_depth = pd0.Gps2.depth_m
                
                #Screen out invalid data
                temp_depth[temp_depth < 0.01] = np.nan
                
                #Use the last valid depth for each ensemble
                last_depth_col_idx = np.sum(np.isnan(temp_depth) == False, axis=1)-1
                last_depth_col_idx[last_depth_col_idx == -1] = 0               
                row_index = np.arange(len(temp_depth))
                last_depth=np.empty(row_index.size)
                #? DEBUG to check correctness
                for row in row_index:
                    last_depth[row]=temp_depth[row,last_depth_col_idx[row]]               
#                    idx = np.ravel_multi_index(temp_depth, dims=(temp_depth.shape[0], temp_depth.shape[1], last_depth_col_idx), order='C')
#                last_depth = temp_depth[idx]
                
                #Determine if mmt file has a scale factor and offset for the depth sounder
                if mmt_config['DS_Cor_Spd_Sound'] == 0:
                    scale_factor = mmt_config['DS_Scale_Factor']
                else:
                    scale_factor = pd0.sensors.sos_mps /  1500
                    
                #Apply scvale factor, offset. and draft
                #Note: Only the ADCP draft is stored.  The transducer
                # draft or scaling for depth sounder data cannot be changed in QRev
                ds_depth = (last_depth * scale_factor) 
                + mmt_config['DS_Transducer_Depth']
                + mmt_config['DS_Transducer_Offset']
                
                self.depths.add_depth_object(ds_depth, 'DS', pd0.Inst.freq, 
                                   mmt_config['Offsets_Transducer_Depth'],
                                   kargs = [cell_depth_m, cell_size_all_m])
                
            #Set depth reference to value from mmt file
            if 'Proc_River_Depth_Source' in mmt_config:
                if mmt_config['Proc_River_Depth_Source'] == 0:
                    self.depths.set_depth_reference('BT')
                elif mmt_config['Proc_River_Depth_Source'] == 1:
                    if self.depths.ds_depths is not None:
                        self.depths.set_depth_reference('DS')
                    else:
                        self.depths.set_depth_reference('BT')
                elif mmt_config['Proc_River_Depth_Source'] == 2:
                    self.depths.set_depth_reference('VB')
                elif mmt_config['Proc_River_Depth_Source'] == 3:
                    if self.depths.vb_depths is None:
                        self.depths.set_depth_reference('BT')
                    else:
                        self.depths.set_depth_reference('VB')
                    self.depths.composite_depths(self, kargs=['On'])
                elif mmt_config['Proc_River_Depth_Source'] == 4:
                    if self.depth.bt_depths is None:
                        self.depths.set_depth_reference('VB')
                    else:
                        self.depths.set_depth_reference('BT')
                    self.depths.composite_depths(self, kargs=['On'])
                else:
                    self.depths.set_depth_reference('BT')
                    self.depths.composite_depths(self, kargs=['On'])  
            else:
                if mmt_config['DS_Use_Process'] > 0:
                    if self.depth.ds_depths  is not None:
                        self.depths.set_depth_reference('DS')
                    else:
                        self.depths.set_depth_reference('BT')
                else:
                    self.depths.set_depth_reference('BT')
                self.depths.composite_depths(self, kargs = ['Off'])
                
            #Create water_data object -----------------------------
            
            # Check for RiverRay and RiverPro data
            firmware = str(pd0.Inst.firm_ver[0])
            excluded_dist = 0
            if firmware[:2] == '56' and np.nanmean(pd0.Inst.beams) < 5:
                excluded_dist = 0.25
                
            if firmware[:2] == '44' or firmware[:2] == '56':
                
                #Process water velocities for RiverRay and RiverPro
                self.w_vel = WaterData()
                self.w_vel.populate_data(pd0.Wt.vel_mps, 
                                         pd0.Inst.freq.T, 
                                         pd0.Cfg.coord_sys, 
                                         'None', 
                                         pd0.Wt.rssi, 
                                         'Counts', 
                                         excluded_dist, 
                                         self.cells_above_sl, 
                                         sl_cutoff_per, 
                                         0,
                                         'Percent', 
                                         sl_lag_effect, 
                                         pd0.Cfg.wm[0], 
                                         pd0.Cfg.wf_cm[0] / 100, 
                                         kargs=[pd0.Wt.corr,
                                                pd0.Surface.vel_mps,
                                                pd0.Surface.rssi,
                                                pd0.Surface.corr,
                                                pd0.Surface.no_cells])
                
            else:
                #Process water velocities for non-RiverRay ADCPs
                self.w_vel = WaterData()
                self.w_vel.populate_data(pd0.Wt.vel_mps, 
                                         pd0.Inst.freq.T, 
                                         pd0.Cfg.coord_sys, 
                                         'None', 
                                         pd0.Wt.rssi, 
                                         'Counts', 
                                         excluded_dist, 
                                         self.cells_above_sl, 
                                         sl_cutoff_per, 
                                         0,
                                         'Percent', 
                                         sl_lag_effect, 
                                         pd0.Cfg.wm[0], 
                                         pd0.Cfg.wf_cm[0] / 100, 
                                         kargs=[pd0.Wt.corr])
                
            #Initialize boat vel
            self.boat_vel = BoatStructure()
            
            self.boat_vel.add_boat_object('TRDI', pd0.Bt.vel_mps, pd0.Inst.freq.T, 
                                            pd0.Cfg.coord_sys[0], 'BT', kargs= [mmt_config['Proc_Use_3_Beam_Solution_For_BT'],
                                            pd0.Cfg.bm[0]])
            
            self.boat_vel.set_nav_reference('BT')
            
            #Compute velocities from GPS Data
            #------------------------------------
            
            #Raw Data
            raw_GGA_utc = pd0.Gps2.utc
            raw_GGA_lat = pd0.Gps2.lat_deg
            
            #Determine correct sign for latitude
            idx = np.where(pd0.Gps2.lat_ref == 'S')
            if len(idx) > 0:
                raw_GGA_lat[idx] = raw_GGA_lat[idx]  * -1
            raw_GGA_lon = pd0.Gps2.lon_deg
            
            #Determine correct sign for longitude
            idx = np.where(pd0.Gps2.lon_ref=='W')
            if len(idx) > 0:
                raw_GGA_lon[idx] = raw_GGA_lon[idx] * -1
            
            #Assign data to local variables
            raw_GGA_alt = pd0.Gps2.alt
            raw_GGA_diff = pd0.Gps2.corr_qual
            raw_GGA_hdop = pd0.Gps2.hdop
            raw_GGA_num_sats = pd0.Gps2.num_sats
            raw_VTG_course = pd0.Gps2.course_true
            raw_VTG_speed = pd0.Gps2.speed_k_mph * 0.2777778
            raw_VTG_delta_time = pd0.Gps2.vtg_delta_time
            raw_VTG_mode_indicator = pd0.Gps2.mode_indicator
            raw_GGA_delta_time = pd0.Gps2.gga_delta_time
            
            # RSL provided ensemble values, not supported for TRDI data
            ext_GGA_utc=[];
            ext_GGA_lat=[];
            ext_GGA_lon=[];
            ext_GGA_alt=[];
            ext_GGA_diff=[];
            ext_GGA_hdop=[];
            ext_GGA_num_sats=[];
            ext_VTG_course=[];
            ext_VTG_speed=[];
             
            #QRev methods GPS processing methods
            GGA_p_method = 'Mindt'
            GGA_v_method = 'Mindt'
            VTG_method = 'Mindt'
            
            #If valid gps data exist, process the data
            if np.sum(np.sum(np.isnan(raw_GGA_lat)==False)) > 0 or np.sum(np.sum(np.isnan(raw_VTG_speed) == False)) > 0:
                
                #Process raw GPS data
                self.gps = GPSData()
                self.gps.populate_data(raw_GGA_utc, raw_GGA_lat, raw_GGA_lon, raw_GGA_alt, raw_GGA_diff, raw_GGA_hdop, 
                                       raw_GGA_num_sats, raw_GGA_delta_time, raw_VTG_course, raw_VTG_speed, 
                                       raw_VTG_delta_time, raw_VTG_mode_indicator, ext_GGA_utc, ext_GGA_lat, ext_GGA_lon, 
                                       ext_GGA_alt, ext_GGA_diff, ext_GGA_hdop, ext_GGA_num_sats, ext_VTG_course, 
                                       ext_VTG_speed, GGA_p_method, GGA_v_method, VTG_method)
                
                #If valid GGA data exists create GGA boat velocity object
                if np.sum(np.sum(np.isnan(raw_GGA_lat) == False)) > 0:
                    self.boat_vel.add_boat_object('TRDI', self.gps.gga_velocity_ens_mps, np.nan, 'Earth', 'GGA')
                    
                    
                #If valid vtg data exist create VTG boat velocity object
                if np.sum(np.sum(np.isnan(raw_VTG_speed) == False)) > 0:
                    self.boat_vel.add_boat_object('TRDI', self.gps.vtg_velocity_ens_mps, np.nan, 'Earth', 'GGA' 
                                                  ,[np.nan,np.nan])
                    
                    
            #Create Edges Object
            self.edges = Edges()
            self.edges.populate_data('Fixed','MeasMag')
                
            #Determine number of ensembles to average
            nens_L = mmt_config['Q_Shore_Pings_Avg']
            #TRDI uses same number on left and right edges
            nens_R = nens_L
            
            #Set indices for ensembles in the moving-boat
            #portion of the transect
            self.in_transect_idx = np.arange(0,pd0.Bt.vel_mps.shape[1])
            
            #Determine left and right edge distances
            if mmt_config['Edge_Begin_Left_Bank'] == True:
                dist_L = mmt_config['Edge_Begin_Shore_Distance']
                dist_R = mmt_config['Edge_End_Shore_Distance']
                self.start_edge = 'Left'
            else:
                dist_L = mmt_config['Edge_End_Shore_Distance']
                dist_R = mmt_config['Edge_Begin_Shore_Distance']
                self.start_edge = 'Right'
                
            #Create edge
            if mmt_config['Q_Left_Edge_Type'] == 0:
                self.edges.create_edge('left', 'Triangular', dist_L, kargs=[nens_L])
            elif mmt_config['Q_Left_Edge_Type'] == 1:
                self.edges.create_edge('left', 'Rectangular', dist_L, kargs=[nens_L])
            elif mmt_config['Q_Left_Edge_Type'] == 2:
                coefficient = mmt_config['Q_Left_Edge_Coeff']
                self.edges.create_edge('left', 'custom', dist_L, kargs=[coefficient, nens_L])
                
            #Create edge
            if mmt_config['Q_Right_Edge_Type'] == 0:
                self.edges.create_edge('right', 'Triangular', dist_R, kargs=[nens_R])
            elif mmt_config['Q_Right_Edge_Type'] == 1:
                self.edges.create_edge('right', 'Rectangular', dist_L, kargs=[nens_L])
            elif mmt_config['Q_Right_Edge_Type'] == 2:
                coefficient = mmt_config['Q_Right_Edge_Coeff']
                self.edges.create_edge('right', 'custom', dist_L, kargs=[coefficient, nens_L])
                
            #Create extrap object
            if mmt_config['Q_Top_Method'] == 0:
                top = 'Power'
            elif mmt_config['Q_Top_Method'] == 1:
                top = 'Constant'
            elif mmt_config['Q_Top_Method'] == 2:
                top = '3-Point'
                
            #Determine bottom method
            if mmt_config['Q_Bottom_Method'] == 0:
                bot = 'Power'
            elif mmt_config['Q_Bottom_Method'] == 2:
                bot = 'No Slip'
                
            self.extrap = ExtrapData()
            self.extrap.populate_data(top, bot, mmt_config['Q_Power_Curve_Coeff']) 
            
            #Sensor Data
            
            self.sensors = Sensors() 
            
            #Heading
            
            #Internal Heading
            heading = pd0.Sensor.heading_deg.T
            heading_src = pd0.Cfg.head_src[0]
            
            #WR2 only has one set of magvar and heading offset
            magvar = mmt_config['Offsets_Magnetic_Variation']
            heading_offset = mmt_config['Ext_Heading_Offset']
            
            #Create internal heading sensor
            self.sensors.add_sensor_data('heading_deg', 'internal', heading, heading_src, kargs=[magvar, heading_offset])
            
            #External Heading
            ext_heading_check = np.where(np.isnan(pd0.Gps2.heading_deg))
            if len(ext_heading_check[0]) <= 0:
                
                self.sensors.set_selected('heading_deg', 'internal')
            else:
                #Determine external heading for each ensemble
                #Using the minimum time difference
                d_time = np.abs(pd0.Gps2.hdt_delta_time)
                d_time_min = np.nanmin(d_time.T, 0).T
                use = np.tile([np.nan], d_time.shape)
                for nd_time in range(len(d_time_min)):
                    use[nd_time,:] = np.abs(d_time[nd_time,:]) == d_time_min[nd_time]
                    
                ext_heading_deg = np.tile([np.nan], (len(d_time_min)))
                for nh in range(len(d_time_min)):
                    idx = np.where(use[nh,:]) 
                    if len(idx[0]) > 0:
                        idx = idx[0][0]
                        ext_heading_deg[nh] = pd0.Gps2.heading_deg[nh,idx]
                        
                #Create external heading sensor
                self.sensors.add_sensor_data('heading_deg', 'external', ext_heading_deg, 'GPS', kargs=[magvar,heading_offset])
            
                #Determine heading source to use from mmt setting
                source_used = mmt_config['Ext_Heading_Use']
                if source_used == True:
                    self.sensors.set_selected('heading_deg', 'external')
                else:
                    self.sensors.set_selected('heading_deg', 'internal')
                
            
            #Pitch
            pitch = arctand(tand(pd0.Sensor.pitch_deg) * cosd(pd0.Sensor.roll_deg))
            pitch_src = pd0.Cfg.pitch_src[0]
            
            #Create pitch sensor
            self.sensors.add_sensor_data('pitch_deg', 'internal', pitch, pitch_src)
            self.sensors.set_selected('pitch_deg', 'internal')
            
            #Roll
            roll = pd0.Sensor.roll_deg.T
            roll_src = pd0.Cfg.roll_src[0]
            
            #Create Roll sensor
            self.sensors.add_sensor_data('roll_deg', 'internal', roll, roll_src)
            self.sensors.set_selected('roll_deg', 'internal')
            
            #Temperature
            temperature = pd0.Sensor.temperature_deg_c.T
            temperature_src = pd0.Cfg.temp_src[0]
            
            #Create temperature sensor
            self.sensors.add_sensor_data('temperature_deg_c', 'internal', temperature, temperature_src)
            self.sensors.set_selected('temperature_deg_c', 'internal')
            
            #Salinity
            pd0_salinity = pd0.Sensor.salinity_ppt.T
            pd0_salinity_src = pd0.Cfg.sal_src[0]
            
            #Create salinity sensor
            self.sensors.add_sensor_data('salinity_ppt', 'internal', pd0_salinity, pd0_salinity_src)
            mmt_salinity = mmt_config['Proc_Salinity']
            self.sensors.add_sensor_data('salinity_ppt', 'user', mmt_salinity, 'mmt')
            self.sensors.set_selected('salinity_ppt', 'internal')
            
            #Speed of Sound
            speed_of_sound = pd0.Sensor.sos_mps.T
            speed_of_sound_src = pd0.Cfg.sos_src[0]
            self.sensors.add_sensor_data('speed_of_sound_mps','internal', speed_of_sound, speed_of_sound_src)
            
            #The raw data are referenced to the internal SOS
            self.sensors.set_selected('speed_of_sound_mps', 'internal')
            
            #Ensemble times
            #Compute time for each ensemble in seconds
            ens_time_sec = pd0.Sensor.time[:,0] * 3600 + pd0.Sensor.time[:,1] * 60 + pd0.Sensor.time[:,2] + pd0.Sensor.time[:,3] / 100
            
            #compute the duration of each ensemble in seconds
            #adjusting for lost data
            ens_delta_time = np.tile([np.nan], ens_time_sec.shape)
            idx_time = np.where(np.isnan(ens_time_sec) == False)[0]
            ens_delta_time[idx_time[1:]] = nandiff(ens_time_sec[idx_time])
            
            #Adjust for transects tha last past midnight
            idx_24hr = np.where(ens_delta_time < 0)
            ens_delta_time[idx_24hr] = 24 * 3600 + ens_delta_time[idx_24hr]
            ens_delta_time = ens_delta_time.T
            
            #Start date and time
            idx = np.where(np.isnan(pd0.Sensor.time[:,0]) == False)[0][0]
            start_year = int(pd0.Sensor.date[idx,0])
            
            #StreamPro doesn't include y2k dates
            if start_year < 100:
                start_year = 2000 + int(pd0.Sensor.date_not_y2k[idx,0])
                
            start_month = int(pd0.Sensor.date[idx,1])
            start_day = int(pd0.Sensor.date[idx,2])
            start_hour = int(pd0.Sensor.time[idx, 0])
            start_min = int(pd0.Sensor.time[idx, 1])
            start_sec = int(pd0.Sensor.time[idx,2] + pd0.Sensor.time[idx,3] / 100)
            
            start_dt = datetime(start_year, start_month, start_day, start_hour, start_min, start_sec)
            start_serial_time = mdates.date2num(start_dt)
            start_date = datetime.strftime(start_dt, '%m/%d/%Y')
            
            #End data and time
            idx = np.where(np.isnan(pd0.Sensor.time[:,0])==False)[0][-1]
            end_year = int(pd0.Sensor.date[idx,0])
            #StreamPro does not include Y@K dates
            if end_year < 100:
                end_year = 2000 + int(pd0.Sensor.date_not_y2k[idx,0])
                
            end_month = int(pd0.Sensor.date[idx,1])
            end_day = int(pd0.Sensor.date[idx,2])
            end_hour = int(pd0.Sensor.time[idx, 0])
            end_min = int(pd0.Sensor.time[idx, 1])
            end_sec = int(pd0.Sensor.time[idx, 2] + pd0.Sensor.time[idx,3] / 100)
            
            end_dt = datetime(end_year, end_month, end_day, end_hour, end_min, end_sec)
            end_serial_time = mdates.date2num(end_dt)
            end_date = datetime.strftime(start_dt, '%m/%d/%Y')
            
            #Create date/time object
            self.datetime = DateTime()
            self.datetime.populate_data(start_date, start_serial_time, end_serial_time, ens_delta_time)
            
            #Transect checked for use in discharge computation
            self.checked = mmt_transect.Checked
            
            if kargs is None:
                self.adcp = InstrumentData()
                self.adcp.populate_data('TRDI', kargs=[mmt_transect, pd0, mmt])
            else:
                self.adcp = InstrumentData()
                self.adcp.populate_data('TRDI', kargs= np.hstack([[mmt_transect, pd0, mmt], kargs]))
            
            
            
    def compute_instrument_cell_data(self, pd0):
        
        #Number of ensembles
        num_ens = np.array(pd0.Wt.vel_mps).shape[-1]
        #Retrieve and compute cell information
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
        '''Coordinates side lobe cutoff calls'''
        
        self.w_vel.adjust_side_lobe(self)
        
    def change_selection(self):
        '''Changes whether the transect is is include in the final discharge'''
        
        self.checked = not self.checked
        
        
    def set_extrapolation(self,top,bot,exp):
        
        self.extrap.set_extrap_data(top, bot, exp)
        
        
    def change_q_ensembles(self, proc_method):
        '''Sets in_transect_idx to all ensembles, except in the case of Sontek data
        where RSL processing is applied
        
        Input:
        proc_method: processing method WR2, RSL, Qrev
        '''
        
        if proc_method == 'RSL':
            num_ens = self.boat_vel.bt_vel.__u_processed_mps.shape[1]
            #Determine number of ensembles for each edge
            if self.start_edge == 'Right':
                self.in_transect_idx = np.arange(self.edges.__right.__num_ens_2_avg,num_ens-self.edges.__left.__num_ens_2_avg)
            else:
                self.in_transect_idx = np.arange(self.edges.__left.__num_ens_2_avg,num_ens-self.edges.__right.__num_ens_2_avg)
        else:
            self.in_transect_idx = np.arange(0,self.boat_vel.bt_vel.__u_processed_mps.shape[1])
            
    def change_start_edge(self, setting):
        '''start edge'''
        self.start_edge = setting
        
    def change_edge(self, edge, property, setting):
        '''Change an edge property (velocity or coefficient)
        
        Input:
        edge: edge object to change (left or right)
        property: property to change (rec_edge_method or vel_method)
        setting: new property setting
        '''
        
        self.edges.change_property(property, setting, edge)    
        
    def change_coord_sys(self, new_coord_sys):
        '''Coordinates changing the coordinate system of the water and boat data
        current implementation only allows changes for original to higher order coordinate
        systems: Beam - Inst - Ship - Earth
        
        Input:
        new_coord_sys: name of new coordinate system (Beam, Int, Ship, Earth)
        '''
        self.w_vel.change_coord_sys(new_coord_sys, self.sensors, self.adcp)
        self.boat_vel.change_coord_sys(new_coord_sys, self.sensors, self.adcp)
        
    def change_nav_reference(self, update, new_nav_ref):
        '''Method to set the naviagion reference for the water data
        
        Input:
        new_nav_ref: new navigation reference (bt_vel, gga_vel, vtg_vel)
        '''
        
        self.boat_vel.change_nav_reference(new_nav_ref, self)
        
        if update == True:
            self.update_water()
            
    def change_mag_var(self, mag_var):
        '''Change magnetic variation
        
        Input:
        mag_var: magnetic variation in degrees
        '''
        
        #update object
        if (self.sensors.heading_deg.external is not None):
            self.sensors.heading_deg.set_mag_var(mag_var, 'external')
        
        if self.sensor.heading_deg.selected == 'internal':
            old = getattr(self.sensors.heading_deg, self.sensors.heading_deg.selected)
            old_mag_var = old.mag_var
            mag_var_change = mag_var - old_mag_var
            self.sensors.heading_deg.set_mag_var(mag_var, 'internal')
            self.boat_vel.chang_mag_var(mag_var_change)
            self.w_vel.change_mag_var(self.boat_vel,mag_var_change)
        
        self.sensors.heading_deg.set_mag_var(mag_var, 'internal')
        
        self.update_water()
        
    def change_offset(self, h_offset):
        '''Change the heading offset (aignment correction)
        
        Input:
        h_offset: heading offset in degrees
        '''
        self.sensors.heading_deg.set_align_correction(h_offset, 'internal')
        
        if self.sensors.heading_deg.selected == 'external':
            old = getattr(self.sensors.heading_deg, self.sensors.heading_deg.selected)
            old_offset = old.align_correction
            offset_change = h_offset - old_offset
            self.boat_vel.bt_vel.change_offset(offset_change)
            self.w_vel.change_offset(self.boat_vel, offset_change)
        
        self.sensors.heading_deg.set_align_correction(h_offset, 'external')
        
        self.update_water()
        
    def change_heading_source(self, h_source):
        source = getattr(self.sensors.heading_deg, h_source)
        if source is not None:
            old = getattr(self.sensors.heading_deg, self.sensors.heading_deg.selected)
            old_heading = old.data
            new_heading = source.data
            heading_change = new_heading - old_heading
            self.sensors.heading_deg.set_selected(h_source)
            self.boat_vel.bt_vel.change_heading_source(heading_change)
            self.w_vel.change_heading_source(self.boat_vel, heading_change)
            
        self.update_water()
            
    def update_water(self):
        '''Method called from set_nav_reference, boat_interpolation and boat filters
        to ensure that changes in boatvel are reflected in the water data'''
        self.w_vel.set_nav_reference(self.boat_vel)
        
        #Reapply water filters and interpolations
        #Note wt_filters calls apply_filter which automatically calls
        #apply_interpolation so both filters and interpolations
        #are applied with this one call
        
        self.wt_filters()
        
         
    def wt_filters(self, kargs = None):
        '''Coordinate water velocity filters
        
        Input:
        kargs[0]: Filter Type (Beam, Difference, Vertical, Other, Excluded, SNR, WT_Depth)
        kargs[1]: Filter setting (Auto, Manual, Off)
        kargs[2]: Threshold if Manual
        '''
        
        self.w_vel.apply_filter(self, kargs)
        
    def side_lobe_cutoff(self, depths, draft, cell_depth, sl_lag_effect, kargs):
        '''Computes side lobe vutoff based on beam angle with no allowance for lag'''
        
     
        #Compute minimum depths for each ensemble
        min_depths = np.nanmin(depths,0)
        
        #Compute range from transducer
        range_from_xducer = min_depths - draft
        
        #Adjust fro transducter angle
        if kargs[0] == 'Percent':
            coeff = kargs[1]
        elif kargs[0] == 'Angle':
            coeff = cosd(kargs[1])
        
        #Compute sidelobe cutoff to centerline    
        cutoff = np.array(range_from_xducer * coeff - sl_lag_effect+draft)
        
        #Compute logical side lobe cutoff matrix
        self.cells_above_sl = (cell_depth - cutoff) < 0
        
    def boat_interpolations(self, update, target, kargs=None):
        '''Coordinates boat velocity interpolations
        
        Input:
        target: boat velocity reference (BT or GPS)
        kargs: type of interpolation
        '''
        
        if target == 'BT':
            self.boat_vel.bt_vel.apply_interpolation(self, kargs)
            
        if target == 'GPS':
            try:
                vel = getattr(self.boat_vel, 'gga_vel')
                if vel is not None:
                    self.boat_vel.gga_vel.apply_interpolation(self, kargs)
            except:
                pass
            
            try:
                vel = getattr(self.boat_vel, 'vtg_vel')
                if vel is not None:
                    self.boat_vel.vtg_vel.apply_interpolation(self, kargs)
            except:
                pass
            
        #Apply composite tracks setting
        self.composite_tracks(False)
        
        #Update water to reflect changes in boat_vel
        if update == True:
            self.update_water()
            
    def composite_tracks(self, update, kargs = None):
        '''Coordinate application of composite tracks
        
        Input:
        kargs: Optional settings for composite tracks (on or off).  If
        kargs not specified the setting currently saved in boat_vel object will be saved
        '''
        
        #Determine if setting is specified
        if kargs is None:
            #Process transect using saved setting
            self.boat_vel.composite_tracks(self)
        else:
            #Process transect usin new setting
            self.boat_vel.composite_tracks(self, kargs)
            
        #Update water data to reflect changes in boatvel
        if update == True:
            self.update_water()
            
    def boat_filters(self, update, kargs = None):
        '''Coordinates application of boat filters to bottom track data
        
        Input:
        kargs: specified in pairs or triplets, can be multiple groups
        kargs[n]: filter type (Beam, Difference, Vertical, Other
        kargs[n+1]: filter setting (Auto, Manual, Off)
        kargs[n+2]: Threshold if manual
        '''
        
        #apply filter to transect
        self.boat_vel.bt_vel.apply_filter(self, kargs)
        
        if self.boat_vel.selected == 'bt_vel' and update:
            self.update_water()
            
    def gps_filters(self, update, kargs=None):
        '''Coordinate filters for GPS based boat velocities
        
        Input:
        kargs[n]: Filter type (Differential, Altitude, HDOP, Other)
        kargs[n+1]: Filter setting (Auto, Manual, Off)
        kargs[n+2]: Threshold if manual
        '''
        
        if self.boat_vel.gga_vel is not None:
            self.boat_vel.gga_vel.apply_gps_filter(self, kargs)
        if self.boat_vel.vtg_vel is not None:
            self.boat_vel.vtg_vel.apply_gps_filter(self, kargs)
            
        if (self.boat_vel.selected == 'VTG' or self.boat_vel.selected == 'GGA') and update == True:
            self.update_water()
            
    def set_depth_reference(self, update, setting):
        '''Coordinates setting the depth reference
        
        Input:
        setting: depth reference (BT, VB, DS)
        '''
        
        self.depths.set_depth_reference(setting)
        
        if update == True:
            self.process_depths(update)
            self.adjust_side_lobe()
            
    def apply_averaging_method(self, setting):
        '''Method to apply the selected averaging method to the BT team depths to achieve a single
        average depth.  It is only applicable to the multiple beams used for BT, not VB or DS.
        
        Input:
        setting: averaging method (IDW, Simple)
        '''
        
        self.depths.bt_depths.compute_avg_BT_depth(setting)
        
        self.process_depths(False)
            
            
    def process_depths(self, update, kargs=None):
        '''Method applies filter, composite, and interpolation settings to  depth objects
        so that all are update using the same filter and interpolation settings
        
        Input:
        kargs[0] process to be applied (Filter, Composite, Interpolate)
        kargs[1]: setting for process
        '''
        
        #Get current settings
        filter = getattr(self.depths, self.depths.selected)
        filter_setting = filter.filter_type
        interp_setting = filter.interp_type
        composite_setting = self.depths.composite
        bt_avg_setting = self.depths.bt_depths.avg_method
        valid_method_setting = self.depths.bt_depths.valid_data_method
        
        #if the process and setting are provided apply those settings
        #if not simply reprocess the data using the data stored in the objects
        if kargs is not None:
            narg = len(kargs)
            for n in np.arange(0,narg,2):
                
                if kargs[n] == 'Filter':
                    filter_setting = kargs[n+1]
                elif kargs[n] == 'Composite':
                    composite_setting = kargs[n+1]
                elif kargs[n] == 'Interpolate':
                    interpolate_setting = kargs[n+1]
                elif kargs[n] == 'AvgMethod':
                    bt_avg_setting = kargs[n+1]
                elif kargs[n] == 'ValidMethod':
                    valid_method_setting = kargs[n+1]
                    
        self.depths.set_valid_data_method(valid_method_setting)
        self.depths.bt_depths.set_avg_method(bt_avg_setting)
        self.depths.depth_filter(self, filter_setting)
        self.depth_interpolation(self, interpolate_setting)
        self.depths.composite_depths(self, composite_setting)
        self.w_vel.adjust_side_lobe(self)
        
        if update == True:
            self.update_water()
            
            
    def change_draft(self, input):
        '''Changes the draft for the specified transects and selected depth'''
        
        if self.depths.vb_depths is not None:
            self.depths.vb_depths.change_draft(input)
        if self.depths.bt_depths is not None:
            self.depths.bt_depths.change_draft(input)
            
            
    def set_property(self, property, setting):
        setattr(property, setting)
        
    def update_data_structure(self):
        #Salinity
        pd0_salinity = self.sensors.salinity_ppt.user.__data_orig
        pd0_salinity_src = self.sensors.salinity_ppt.user.__source
        
        #Create slainity sensor
        self.sensors.add_sensor_data('salinity_ppt','internal',pd0_salinity,pd0_salinity_src)
        self.sensors.set_selected('salinity_ppt', 'internal')
        
    
def sos_user(self, kargs = None):
    '''Compute new speed of sound from temperature and salinity
    
    Output:
    new_sos: newly computed speed of sound
    old_sos: previously used speed of sound
    '''
    
    #Assign selected temperature data to local variable
    temp = getattr(self.sensors.temperature_deg_c, self.sensors.temperature_deg_c.selected)
    temperature = temp.data
    #Assign selected salinity to local variable
    sal = getattr(self.sensors.salinity_ppt, self.sensors.salinity_ppt.selected)
    salinity = sal.data
    old = getattr(self.sensors.speed_of_sound_mps, self.sensors.speed_of_sound_mps.selected)
    old_sos = old.data
    
    if self.sensors.temperature_deg_c.selected == 'internal':
        new_sos = self.sensors.speed_of_sound_mps.user.data_orig
        self.sensors.speed_of_sound_mps.user.change_data(new_sos)
        self.sensors.speed_of_sound_mps.user.set_source('Internal (ADCP)')
        self.sensors.set_selected('speed_of_sound_mps', 'user')
    else:
        #Compute new speed of sound
        new_sos = Sensors().speed_of_sound(temperature, salinity)
        
        #Save new speed of sound to user sensor object with a source as computed
        if self.sensors.speed_of_sound_mps.user is not None:
            self.sensors.set_selected('speed_of_sound_mps', 'user')
            self.sensors.speed_of_sound_mps.user.change_data(new_sos)
            self.sensors.speed_of_sound_mps.user.set_source('Computed')
        else:
            self.sensors.add_sensor_data('speed_of_sound_mps', 'user'. new_sos, 'Computed')
            self.sensors.set_selected('speed_of_sound_mps', 'user')
            
    return (old_sos, new_sos)
    
def allocate_transects(source, mmt, kargs): 
    
    #Refactored from TransectData to iteratively create TransectData objects
        #----------------------------------------------------------------
        if kargs[0] == 'Q':
            transects = 'transects'
            active_config = 'active_config' 
            
            if kargs[1] == True:
                files_to_load = np.array([x.Checked for x in mmt.transects], dtype=bool)
                file_names = [x.Files for x in mmt.transects]
            else:
                files_to_load = np.array(np.ones(len(mmt.transects)), dtype=bool)
                file_names = [x.Files for x in mmt.transects]
                    
                
        elif kargs[0] == 'MB':
            transects = 'mbt_transects'
            active_config = 'mbt_active_config'
            files_to_load =np.array([x.Checked for x in mmt.mbt_transects], dtype=bool)
            file_names = [x.Files for x in mmt.mbt_transects if x.Checked == 1]
        
        files_to_load_idx = np.where(files_to_load == True)[0]
          
        pathname = mmt.infile[:mmt.infile.rfind('/')]
        
        # Determine if any files are missing
        
        valid_files = []
        for x in file_names:
            x[0].Path = x[0].Path[x[0].Path.rfind('\\') + 1:]
            if os.path.exists(''.join([pathname,'/',x[0].Path])):
                valid_files.append((x, 1))
            else:
                valid_files.append((None, 0))
                
        pd0_data = []
        pd0_threads = []
        thread_id = 0
        
        def add_pd0(file_name):
            pd0_data.append(Pd0TRDI(file_name))
            
        for x in valid_files:
            if x[1] == 1:
#                 add_pd0(''.join([pathname,'/',x[0][0].Path]))
                pd0_thread = MultiThread(thread_id=thread_id, function=add_pd0, args = {'file_name': ''.join([pathname,'/',x[0][0].Path])})
                thread_id += 1
                pd0_thread.start()
                pd0_threads.append(pd0_thread)
                 
        for x in pd0_threads:
            x.join()
             
        pd0_data = [Pd0TRDI(''.join([pathname,'/',x[0][0].Path])) for x in valid_files if x[1] == 1]
        
        processed_transects = []
        transect_threads = []
        thread_id = 0
        
        def add_transect(transect, source, in_file, pd0_data, mmt):
            
            transect.get_data(source, in_file,pd0_data, mmt)
            processed_transects.append(transect)
            
        # Process each transect
        for k in range(len(pd0_data)):
            
            transect = TransectData()
            transect.active_config = active_config
            transect.transects = transects
           
            if active_config == 'mbt_active_config' or active_config == 'mbt_field_config':
                    
#                 add_transect(transect, 'TRDI', mmt.mbt_transects[k], pd0_data[k], mmt)
                p_thread = MultiThread(thread_id = thread_id, function= add_transect, 
                                       args = {'transect': transect, 
                                               'source':'TRDI', 
                                               'in_file': mmt.mbt_transects[k], 
                                               'pd0_data': pd0_data[k], 
                                               'mmt': mmt})
                p_thread.start()
                transect_threads.append(p_thread)
            else:
                add_transect(transect, 'TRDI', mmt.transects[k], pd0_data[k], mmt)
                p_thread = MultiThread(thread_id = thread_id, function= add_transect, 
                                       args = {'transect': transect, 
                                               'source':'TRDI', 
                                               'in_file': mmt.transects[k], 
                                               'pd0_data': pd0_data[k], 
                                               'mmt': mmt})
                p_thread.start()
                transect_threads.append(p_thread)
                 
         
        for x in transect_threads:
            x.join()
        
        return processed_transects   
        
def adjusted_ensemble_duration(transect, kargs=None):
        
    if transect.adcp.manufacturer == 'TRDI':
        if kargs is not None:
            valid = np.isnan(transect.w_vel.__u_processed_mps) == False
            valid_sum = np.sum(valid)
        else:
            valid_sum = np.isnan(transect.boat_vel.__u_processed_mps) == False
            
        valid_ens = valid_sum > 0
        n_ens = len(valid_ens)
        ens_dur = transect.datetime.ens_duration_sec
        delta_t = np.tile([np.nan], (1,n_ens))
        cum_dur = 0
        for j in n_ens:
            cum_dur = np.nansum(np.hstack([cum_dur, ens_dur[j]]))
            if valid_ens[j] == True:
                delta_t[j] = cum_dur
                cum_dur = 0
    else:
        delta_t = transect.datetime.ens_duration_sec
        
    return delta_t
    
                
    
            
            
          
           
                
            
        
            
            
        
        
            
        
                
        


