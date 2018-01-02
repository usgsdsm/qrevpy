from Classes.MMT_TRDI import MMT_TRDI
import numpy as np
import os
from Classes.TransectData import TransectData, allocate_transects
from Classes.Pd0TRDI import Pd0TRDI
from Classes.PreMeasurement import PreMeasurement
from Classes.MovingBedTests import MovingBedTests
from Classes.MultiThread import MultiThread
from Classes.QComp import QComp

class Measurement(object):
    """Class to hold measurement details for use in the GUI
        Analogous to matlab class: clsMeasurement
    """

    def __init__(self, in_file, source='TRDI', kargs=None):
        self.station_name = None
        self.station_number = None
        self.transects = []
        self.mb_tests = None
        self.sys_test = []
        self.compass_cal = []
        self.compass_eval = []
        self.ext_temp_chk = None
        self.extrap_fit = None
        self.processing = None
        self.comments = None
        self.discharge = None
        self.uncertainty = None
        self.intitial_settings = None
        self.qa = None
        self.userRating = None
        self.comments = []
        
        self.ext_temp_chk = {}
        if source == 'TRDI':
            self.load_trdi(in_file, **{'type': 'Q', 'checked': False})
            
            if kargs is not None and len(kargs) > 2:
                proc_type = kargs[2]
            else:
                proc_type= 'QRev'
                
        select = self.transects[0].boat_vel.selected
        if select == 'bt_vel':
            ref = 'BT'
        elif select == 'gga_vel':
            ref = 'GGA'
        elif select == 'vtg_vel':
            ref = 'VTG'
            
        if self.mb_tests is not None:
            for x in self.mb_tests:
                x.auto_use_2_correct()
                
        #update status "Save initial settings"
        self.intitial_settings = self.current_settings()
        
        #Set processing
        if proc_type == 'QRev':
            #UpdateStatus "Apply Qrev default settings"
            settings = self.QRevDefaultSettings()
            settings['Processing'] = 'Qrev'
            
            for x in self.transects:
                self.apply_settings(x, settings)
        
        elif proc_type == 'None':
            #UpdateStatus "Proicessing with no filters and interpolation
            settings = self.NoFilterInterpSettings()
            settings['Processing'] = 'None'
            
            for x in self.transects:
                self.apply_settings(x, settings)
                
        else:
            self.discharge = QComp()
        
        

    def load_trdi(self, mmt, **kargs):
        '''method to load trdi MMT file'''

        #read in the MMT file
        self.trdi = MMT_TRDI(mmt)
        mmt = self.trdi
    
        #Get properties if they exist, otherwise set them as blank strings
        self.station_name = str(self.trdi.site_info['Name'])
        
        self.station_number = str(self.trdi.site_info['Number'])
        
        self.processing = 'WR2'

#         if len(kargs) > 1:
#             self.transects = TransectData('TRDI', mmt,'Q',kargs[1])
#         else:

        self.transects = allocate_transects('TRDI', mmt, kargs=['Q', kargs['checked']])
            
        #-------------------------Done Refactor----------------------------------------------
        #Create object for pre-measurement tests
        if isinstance(mmt.qaqc, dict) or isinstance(mmt.mbt_transects):
            self.qaqc_TRDI(mmt)
        
        #Save comments from mmt file in comments
        try:
            self.comments.append('MMT Remarks: ' + mmt.site_info.remarks)
        except:
            pass
        
        for t in range(len(self.trdi.transects)):
            notes = getattr(self.trdi.transects[t], 'Notes')
            for x in range(len(notes)):
                note = 'Transect: ' + t + ' File: ' + x + ' ' + notes[x].NoteDate + ': ' + notes[x].NoteText
                self.comments.append(note)
                
        #Get externam temperature
        self.ext_temp_chk['User'] = mmt.site_info['Water_Temperature']
        
        
        #-------------------------------REFACTORED from thresholds_TRDI
        
        threshold_settings = {}
        #WT filter threshold setting
        threshold_settings['num_beam_WT'] = self.set_3_beam_WT_threshold_TRDI(mmt.transects[0])
        threshold_settings['diff_vel_threshold_WT'] = self.set_diff_vel_threshold_WT_TRDI(mmt.transects[0])
        threshold_settings['vert_vel_threshold_WT'] = self.set_vert_vel_threshold_WT_TRDI(mmt.transects[0])
        
        #BT filter threshold settings
        threshold_settings['num_beam_BT'] = self.set_3_beam_BT_threshold_TRDI(mmt.transects[0])
        threshold_settings['diff_vel_threshold_BT'] = self.set_diff_vel_threshold_BT_TRDI(mmt.transects[0])
        threshold_settings['vert_vel_threshold_BT'] = self.set_depth_screening_TRDI(mmt.transects[0])
        
        #Depth filter and averaging settings
        threshold_settings['depth_weighting'] = self.set_depth_weighting_TRDI(mmt.transects[0])
        threshold_settings['depth_valid_method'] = 'TRDI'
        threshold_settings['depth_screening'] = self.set_depth_screening_TRDI(mmt.transects[0])
        
        #--------------------------------------------DONE REFACTOR
        multi_threaded = True
        
        if multi_threaded == True:
            transect_threads = []
            
            def process_transect(mmt, transect, idx):
                transect.change_coord_sys('Earth')
                
                if 'Reference' in mmt.site_info.keys():
                    transect.change_nav_reference(0,mmt.site_info['Reference'])
                    
                self.thresholds_TRDI(transect, threshold_settings)
                
                transect.boat_interpolations(0, 'BT', 'None')
                
                if transect.gps is not None:
                    transect.boat_interpolations(0, 'GPS', 'HoldLast') 
                    
                transect.update_water()
                
                transect.wt_filters(['wt_depth', 'On'])
                transect.wt_interpolations(['Ensembles', 'None'])
                
                mmt_SOS_method = mmt.transects[idx].active_config['Proc_Speed_of_Sound_Correction']
                if mmt_SOS_method == 1:
                    self.change_SOS(idx, 'salinity', 'user')
                    self.change_SOS(idx, 'sosSrc',' Computed', [])
                elif mmt_SOS_method == 2:
                    self.change_SOS(idx, 'sosSrc', 'user', mmt.transects[idx].active_config['Proc_Fixed_Speed_of_Sound'])
             
             
            for x in range(len(self.transects)):       
                p_thread = MultiThread(thread_id = x, function= process_transect, 
                                        args = {'mmt': mmt,
                                                'transect': self.transects[x], 
                                                'idx': x, 
                                                })
                p_thread.start()
                transect_threads.append(p_thread)  
                
            for x in transect_threads:
                x.join()
                
            
        else:       
            #Update Status 'Converting to earth coordinates'
            for x in self.transects:
                x.change_coord_sys('Earth')
                
            #Update Status 'Setting navigation reference'
            if 'Reference' in mmt.site_info.keys():
                for x in self.transects:
                    x.change_nav_reference(0,mmt.site_info['Reference'])
                    
            #Set composite tracks
            #This feature not currently supported in WinRiver 2, and is turned off in Qrev by default
            
            #get threshold settings from mmt file
            #Update Status 'Applying processing settings'
            for x in self.transects:
                self.thresholds_TRDI(x, threshold_settings)
                    
            #apply boat interpolations
            #Update Status 'Applying filters and interpolation to BT
            for x in self.transects:
                x.boat_interpolations(0, 'BT', 'None')
                
            for x in self.transects:
                if x.gps is not None:
                    x.boat_interpolations(0, 'GPS', 'HoldLast') 
                 
            for x in self.transects:
                x.update_water()
                
            #Update status 'Applying filters and interpolation to WT')
            
            '''Note: this is where some computation efficiency could be gained by
            not recomputing each time but applying all filters, interpolating,
            and then computing Q'''
                
            idx = 0
            for x in self.transects:
                x.wt_filters(['wt_depth', 'On'])
                x.wt_interpolations(['Ensembles', 'None'])
                
                mmt_SOS_method = mmt.transects[idx].active_config['Proc_Speed_of_Sound_Correction']
                if mmt_SOS_method == 1:
                    self.change_SOS(idx, 'salinity', 'user')
                    self.change_SOS(idx, 'sosSrc',' Computed', [])
                elif mmt_SOS_method == 2:
                    self.change_SOS(idx, 'sosSrc', 'user', mmt.transects[idx].active_config['Proc_Fixed_Speed_of_Sound'])  
                    
                idx += 1
        
            
    def qaqc_TRDI(self, mmt):
        '''Processes qaqc test, calibrations, and evaluations
        
        Input:
        mmt: object of MMT_TRDI
        
        '''
        #ADCP Test
        if 'RG_Test_TimeStamp' in mmt.qaqc:
            for n in range(len(mmt.qaqc['RG_Test'])):
                p_m = PreMeasurement()
                p_m.populate_data(mmt.qaqc['RG_Test_TimeStamp'][n], mmt.qaqc['RG_Test'][n],'TST')
                self.sys_test.append(p_m)
        else:
            self.sys_test.append(PreMeasurement())
            
        #Compass calibration
        if 'Compass_Cal_Timestamp' in mmt.qaqc:
            for n in range(len(mmt.qaqc['Compass_Cal_Test'])):
                p_m =  PreMeasurement()
                p_m.populate_data(mmt.qaqc['Compass_Cal_Timestamp'], mmt.qaqc['Compass_Cal_Test'], 'TCC')
                self.compass_cal.append(p_m)
        else:
            self.compass_cal.append(PreMeasurement())
            
        #Compass evaluation
        if 'Compass_Eval_Timestamp' in mmt.qaqc:
            for n in range(len(mmt.qaqc['Compass_Eval_Test'])):
                p_m =  PreMeasurement()
                p_m.populate_data(mmt.qaqc['Compass_Eval_Timestamp'], mmt.qaqc['Compass_Eval_Test'], 'TCC')
                self.compass_eval.append(p_m)
        else:
            self.compass_eval.append(PreMeasurement())
            
       
        
        if len(mmt.mbt_transects) > 0:
            
            #Create transect object/s
            transects= allocate_transects('TRDI',mmt,kargs=['MB'])
            if len(transects) > 0:
                self.mb_tests = []
                
                for n in range(len(transects)):
                    mb_test = MovingBedTests()
                    mb_test.populate_data('TRDI', kargs=[transects[n], mmt.mbt_transects[n]])
                    
                    #Save notes from mmt files in comments
                    if 'NoteDate' in mmt.mbt_transects:
                        r, c = mmt.mbt_transects['NoteDate'].shape
                        for n in range(r):
                            for k in range(c):
                                if mmt.mbt_transects['NoteDate'][n,k] is not None:
                                    pass
                                
                    self.mb_tests.append(mb_test)
                                
        
            
    def thresholds_TRDI(self, transect, settings):
        '''Retrieve and pply manual filter settings from mmt file
        
        Input:
        mmt: object of MMT_TRDI
        transect: object of TransectData
        settings: threshold settings computed before processing
        '''
        
        #Apply WT settings
        wt_settings = ['Beam', settings['num_beam_WT'], 'Manual', settings['diff_vel_threshold_WT'],
                       'Vertical', 'Manual', settings['vert_vel_threshold_WT']]
        
        transect.w_vel.apply_filter(transect, kargs=wt_settings)
        
        #Apply BT settings
        bt_settings = ['Beam', settings['num_beam_BT'], 'Manual', settings['diff_vel_threshold_BT'],
                       'Vertical', 'Manual', settings['vert_vel_threshold_BT']]
        transect.boat_vel.bt_vel.apply_filter(transect, kargs=bt_settings)
    
        #Apply depth settings
        transect.depths.set_valid_data_method(settings['depth_valid_method'])
        transect.depths.depth_filter(transect, settings['depth_screening'])
        transect.depths.bt_depths.compute_avg_BT_depth(settings['depth_weighting'])
        transect.depths.composite_depths(transect)
            
            
    def set_3_beam_WT_threshold_TRDI(self, mmt_transect):
        '''Get 3-beam processing for WT from mmt file
        
        Input:
        mmt_transect: object of Transect
        
        Output:
        num_3_beam_WT_Out
        '''
        
        #check consistency
        use_3_beam_WT = mmt_transect.active_config['Proc_Use_3_Beam_Solution_For_WT']
        
        #Use setting from 1st transect for all transects
        if use_3_beam_WT == 0:
            num_beam_WT_out = 4
        else:
            num_beam_WT_out = 3
            
        return num_beam_WT_out
    
    def set_vert_vel_threshold_WT_TRDI(self, mmt_transect):
        '''Get the vertical celocity threshold for WT from mmt
        
        Input:
        mmt_transect: object of Transect
        
        Output:
        vert_vel_threshold_WT[0]: vertical celocity threshold (m/s)
        '''
        
        #Check consistency of vertical velocity threshold (m/s)
        vert_vel_threshold_WT = mmt_transect.active_config['Proc_WT_Up_Velocity_Threshold']
        
        #Use setting from 1st transect for all transects
        #Use setting from 1st transect for all transects
        if isinstance(vert_vel_threshold_WT, float):
            return vert_vel_threshold_WT 
        else:
            return vert_vel_threshold_WT[0]
        
    
    
    def set_3_beam_BT_threshold_TRDI(self, mmt_transect):
        '''Get 3-beam processing for WT from mmt file
        
        Input:
        mmt_transect: object of Transect
        
        Output:
        num_3_beam_WT_Out
        '''
        
        #check consistency
        use_3_beam_BT = mmt_transect.active_config['Proc_Use_3_Beam_Solution_For_BT']
        
        #Use setting from 1st transect for all transects
        if use_3_beam_BT == 0:
            num_beam_BT_out = 4
        else:
            num_beam_BT_out = 3
            
        return num_beam_BT_out
    
    def set_diff_vel_threshold_WT_TRDI(self, mmt_transect):
        '''Get difference (error) velocity threshold for WT from mmt
        
        Input:
        mmt_transect: object of Transect
        
        Output:
        diff_vel_threshold_WT[0]: difference velocity threshold (m/s)
        '''
        
        #Check consistency of difference (error) velocity for WT
        diff_vel_threshold_WT = mmt_transect.active_config['Proc_WT_Error_Velocity_Threshold']
        
        #Use setting from 1st transect for all transects
        if isinstance(diff_vel_threshold_WT, float):
            return diff_vel_threshold_WT 
        else:
            return diff_vel_threshold_WT[0]
    
    def set_diff_vel_threshold_BT_TRDI(self, mmt_transect):
        '''Get difference (error) velocity threshold for BT from mmt
        
        Input:
        mmt_transect: object of Transect
        
        Output:
        diff_vel_threshold_BT[0]: difference velocity threshold (m/s)
        '''
        
        #Check consistency of difference (error) velocity for BT
        diff_vel_threshold_BT = mmt_transect.active_config['Proc_BT_Error_Velocity_Threshold']
        
        #Use setting from 1st transect for all transects
        if isinstance(diff_vel_threshold_BT, float):
            return diff_vel_threshold_BT
        else:
            return diff_vel_threshold_BT[0]
    
    def set_depth_weighting_TRDI(self, mmt_transect):
        '''Get the average depth method from mmt
        
        Input:
        mmt_transect: Object of Transect
        
        Output:
        depth_weighting_setting: method for computing average depth
        '''
        
        #Check consistency of depth averaging method
        depth_weighting = mmt_transect.active_config['Proc_Use_Weighted_Mean_Depth']
        
        #Warn if setting is not consistent for all transects
        
        if isinstance(depth_weighting, int):
            if depth_weighting == 0:
                depth_weighting_setting = 'Simple'
            else:
                depth_weighting_setting = 'IDW'
        else:
            if np.abs(np.sum(np.diff(depth_weighting))) > 0.1:
                depth_weighting_setting = 'IDW'
            else:
                if depth_weighting[0] == 0:
                    depth_weighting_setting = 'Simple'
                else:
                    depth_weighting_setting = 'IDW'
                
        return depth_weighting_setting
    
    def set_depth_screening_TRDI(self, mmt_transect):
        '''Get the depth screening setting from mmt
        
        Input:
        mmt_transect: object of Transect
        
        Output:
        depth_screening_setting: depth screening setting
        '''
        
        #Check consistency of depth screening
        depth_screen = mmt_transect.active_config['Proc_Screen_Depth']
        
        #Warn if setting is not consistent for all transects
        if isinstance(depth_screen, int):
            if depth_screen == 0:
                depth_screening_setting = 'None'
            else:
                depth_screening_setting = 'TRDI'
        else:
            if np.abs(np.sum(np.diff(depth_screen))) > 0.1:
                depth_screening_setting = 'TRDI'
            else:
                #Assign setting
                if depth_screen[0] == True:
                    depth_screening_setting = 'TRDI'
                else:
                    depth_screening_setting = 'None'
                    
        
        return depth_screening_setting
        
    def change_SOS(self, transect, kargs=None):
        '''Applies a change in speed of sound to one or all transects
        and update the discharge and uncertainty computations
        
        Input:
        kargs[0]: optional, transect number to apply change to, if empty
        change is applied to all transects
        kargs[0/1]: sensor property to change
        kargs[1/2]: settinf for sensor property
        '''
        
        s = self.current_settings()
        transect.change_SOS(kargs = kargs[1:])    
        #Check to see if transect number specified
        self.apply_settings(transect, s)
        
    def apply_settings(self, transect, s, kargs = None):
        '''Applies reference, filter, and interpolation settings
        
        Input:
        s: data structure of reference, filter, and interpolation settings
        kargs: transect number to apply settings to only 1 transect
        '''
        
        #Moving-boat ensembles
        if 'processing' in s.keys():
            transect.change_q_ensembles(s['Processing'])
            self.processing = s['Processing']
        
        #Navigation reference
        if transect.boat_vel.selected != s['NavRef']:
            transect.change_nav_reference(0, s['NavRef'])
            if transect.mbt != False:
                self.mb_tests[transect.mbt_idx].auto_use_2_correct(s['NavRef'])
                
        #Composite tracks
        #Changing the nav reference applies the current setting for
        #composite tracks, check to see if a change is needed
        if transect.boat_vel.composite == s['CompTracks']:
            transect.composite_tracks(0, s['CompTracks'])
            
        #Set difference velocity BT filter
        if s.BTdFilter == 'Manual':
            BTdFilter = [s['BTdFilter'], s['BTdFilterThreshold']]
        else:
            BTdFilter = [s['BTdFilter']]
            
        #Set vertical velocity BTfilter
        if s.BTwFilter == 'Manual':
            BTwFilter = [s['BTwFilter'], s['BTwFilterThreshold']]
        else:
            BTwFilter = [s['BTwFilter']]
            
        #Apply BT settings
        bt_settings = ['Beam',s['BTbeamFilter'],'Difference',BTdFilter[:],'Vertical',BTwFilter[:],'Other',s['BTsmoothFilter']];
        transect.boat_filters(0,bt_settings)
        
        #BT Interpolation
        transect.boat_interpolations(0, 'BT', s['BTInterpolation'])
        
        #GPS filter settings
        if transect.gps is not None:
            
            if transect.boat_vel.gga_vel is not None:
                #GGA
                if s['ggaAltitudeFilter'] == 'Manual':
                    ggaAltitudeFilter = [s['ggaAltitudeFilter'], s['ggaAltitudeFilterChange']]
                else:
                    ggaAltitudeFilter = [s['ggaAltitudeFilter']]
                    
                #Set GGA HDOP Filter
                if s['GPSHDOPFilter'] == 'Manual':
                    ggaHDOPFilter = [s['GPSHDOPFilter'], s['GPSHDOPFilterMax'], s['GPSHDOPFilterChange']]
                else:
                    ggaHDOPFilter = [s['GPSHDOPFilter']]
                    
                #Apply GGA filters
                gga_settings = ['Differential', s['ggaDiffQualFilter'], 'Altitude', ggaAltitudeFilter, 'HDOP', ggaHDOPFilter, 'Other', s['GPSSmoothFilter']]
                transect.gps_filters(0, gga_settings)
            
            if transect.boat_vel.vtg_vel is not None:
                
                if s['GPSHDOPFilter'] == 'Manual':
                    vtgSettings = ['HDOP', s['GPSHDOPFilter'], s['GPSHDOPFilterMax'],s['GPSHDOPFilterChange'],'Other',s['GPSSmoothFilter']]
                else:
                    vtgSettings = ['HDOP', s['GPSHDOPFilter'], 'Other', s['GPSSmoothFilter']]
                    
                #Apply VTG filters
                transect.gps_filters(0, vtgSettings)
                
            transect.boat_interpolations(0, 'GPS', s['GPSInterpolation'])
            
        #Set depth reference
        transect.set_depth_reference(0, s['depthReference'])
        
        transect.process_depths(1, kargs = ['Filter', s['depthFilterType'],
                                            'Interpolate', s['depthInterpolation'],
                                            'Composite', s['depthComposite'],
                                            'AvgMethod', s['depthAvgMethod'],
                                            'ValidMethod', s['depthValidMethod']])
        
        #Set WT difference velocity filter
        if s['WTdFilter'] == 'Manual':
            WTdFilter = [s['WTdFilter'], s['WTdFilterThreshold']]
        else:
            WTdFilter = [s['WTdFilter']]
            
        #Set WT vertical velocity filter
        if s['WTwFilter'] == 'Manual':
            WTwFilter = [s['WTwFilter'], s['WTwFilterThreshold']]
        else:
            WTwFilter = [s['WTwFilter']]  
            
        wtSettings = ['Beam', s['WTbeamFilter'], 'Difference', WTdFilter, 
                      'Vertical', WTwFilter, 'Other', s['WTsmoothFilter'],
                      'SNR', s['WTsnrFilter'], 'wtDepth', s['WTwDepthFilter'],
                      'Excluded', s['WTExcludedDistance']]
        
        transect.wt_filters(wtSettings) 
        
        # Recompute extrapolations
        # NOTE: Extrapolations should be determined prior to WT
        # interpolations because the TRDI approach for power/power
        # using the power curve and exponent to estimate invalid cells.
        if self.extrap_fit is None or self.extrap_fit.__fit_method == 'Automatic':
            pass
        
    
            
            
        
    def current_settings(self):
        '''Saves the current settings for a measurement. Since all settings
        in QRev are consistent among all transects in a measurement onlyt the
        settings from the first transect are saved
        '''
        settings = {}
        checked = np.array([x.checked for x in self.transects])
        first_idx = np.where(checked == 1)[0][0]
        transect = self.transects[first_idx]
        
        #Navigation regerence
        settings['NavRef'] = transect.boat_vel.selected
        
        #Composite tracks
        settings['CompTracks'] = transect.boat_vel.composite
        
        #Water track settings
        settings['WTbeamFilter'] = transect.w_vel.beam_filter
        settings['WTdFilter'] = transect.w_vel.d_filter
        settings['WTdFilterThreshold'] = transect.w_vel.d_filter_threshold
        settings['WTwFilter'] = transect.w_vel.w_filter
        settings['WTwFilterThreshold'] = transect.w_vel.w_filter_threshold
        settings['WTsmoothFilter'] = transect.w_vel.smooth_filter
        settings['WTsnrFilter'] = transect.w_vel.snr_filter
        settings['WTwtDepthFilter'] = transect.w_vel.wt_depth_filter
        settings['WTensInterpolation'] = transect.w_vel.interpolate_ens
        settings['WTCellInterpolation'] = transect.w_vel.interpolate_cells
        settings['WTExcludedDistance'] = transect.w_vel.excluded_dist
        
        #Bottom track settings
        settings['BTbeamFilter'] = self.transects[first_idx].boat_vel.bt_vel._BoatData__beam_filter;
        settings['BTdFilter'] = self.transects[first_idx].boat_vel.bt_vel._BoatData__d_filter;
        settings['BTdFilterThreshold'] = self.transects[first_idx].boat_vel.bt_vel._BoatData__d_filter_threshold;
        settings['BTwFilter'] =self.transects[first_idx].boat_vel.bt_vel._BoatData__w_filter;
        settings['BTwFilterThreshold'] = self.transects[first_idx].boat_vel.bt_vel._BoatData__w_filter_threshold;
        settings['BTsmoothFilter'] = self.transects[first_idx].boat_vel.bt_vel._BoatData__smooth_filter;
        settings['BTInterpolation'] = self.transects[first_idx].boat_vel.bt_vel._BoatData__interpolate
        
        #Gps Settings
        if transect.gps is not None:
            #GGA settings
            if transect.boat_vel.gga_vel is not None:
                settings['ggaDiffQualFilter'] = transect.boat_vel.gga_vel._BoatData__gps_diff_qual_filter
                settings['ggaAltitudeFilter'] = transect.boat_vel.gga_vel._BoatData__gps_altitude_filter
                settings['ggaAltitudeFilterChange'] = transect.boat_vel.gga_vel._BoatData__gps_altitude_filter_change
                settings['GPSHDOPFilter'] = transect.boat_vel.gga_vel._BoatData__gps_HDOP_filter
                settings['GPSHDOPFilterMax'] = transect.boat_vel.gga_vel._BoatData__gps_HDOP_filter_max
                settings['GPSHDOPFilterChange'] = transect.boat_vel.gga_vel._BoatData__gps_HDOP_filter_change
                settings['GPSSmoothFilter'] = transect.boat_vel.gga_vel._BoatData__smooth_filter
                settings['GPSInterpolation'] = transect.boat_vel.gga_vel._BoatData__interpolate
            else:
                settings['ggaDiffQualFilter'] = 1
                settings['ggaAltitudeFilter'] = 'Off'
                settings['ggaAltitudeFilterChange'] = []
                
                settings['ggaSmoothFilter'] = []
                if 'GPSInterpolation' not in settings.keys():
                    settings['GPSInterpolation'] = 'None'
                if 'GPSHDOPFilter' not in settings.keys():
                    settings['GPSHDOPFilter'] = 'Off'
                    settings['GPSHDOPFilterMax'] = []
                    settings['GPSHDOPFilterChange'] = []
                if 'GPSSmoothFilter' not in settings.keys():
                    settings['GPSSmoothFilter'] = 'Off'
                    
        if transect.boat_vel.vtg_vel is not None:
            if len(transect.boat_vel.vtg_vel) > 0:
                settings['GPSHDOPFilter'] = transect.boat_vel.vtg_vel._BoatData__gps_HDOP_filter
                settings['GPSHDOPFilterMax'] = transect.boat_vel.vtg_vel._BoatData__gps_HDOP_filter_max
                settings['GPSHDOPFilterChange'] = transect.vtg_vel._BoatData__gps_HDOP_filter_change
                settings['GPSSmoothFilter'] = transect.boat_vel.vtg_vel._BoatData__smooth_filter
                settings['GPSInterpolation'] = transect.boat_vel.vtg_vel._BoatData__interpolate
            else:
                settings['vtgSmoothFilter'] = 'Off'
                if 'GPSInterpolation' not in settings.keys():
                    settings['GPSInterpolation'] = 'None'
                if 'GPSHDOPFilter' not in settings.keys():
                    settings['GPSHDOPFilter'] = 'Off'
                    settings['GPSHDOPFilterMax'] = []
                    settings['GPSHDOPFilterChange'] = []
                if 'GPSSmoothFilter' not in settings.keys():
                    settings['GPSSmoothFilter'] = 'Off'
                    
        #Depth Settings
        settings['depthAvgMethod'] = transect.depths.bt_depths.avg_method
        settings['depthValidMethod'] = transect.depths.bt_depths.valid_data_method
        
        '''Depth settings are always applied to all available depth sources.
        Only those saved in the btDepths are used here but are applied to all sources'''
        settings['depthFilterType'] = transect.depths.bt_depths.filter_type
        settings['depthReference'] = transect.depths.selected
        settings['depthComposite'] = transect.depths.composite
        select = getattr(transect.depths, transect.depths.selected)
        settings['depthInterpolation'] = select.interp_type
        
        #Extrap Settings
        settings['extrapTop'] = transect.extrap._ExtrapData__top_method
        settings['extrapBot'] = transect.extrap._ExtrapData__bot_method
        settings['extrapExp'] = transect.extrap._ExtrapData__exponent
        
        #Edge Settings
        settings['edgeVelMethod'] = transect.edges._Edges__vel_method
        settings['edgeRecEdgeMethod'] = transect.edges._Edges__rec_edge_method
        
        return settings
        
    