from Classes.MMT_TRDI import MMT_TRDI
import numpy as np
import os
from Classes.TransectData import TransectData, allocate_transects
from Classes.Pd0TRDI import Pd0TRDI
from Classes.PreMeasurement import PreMeasurement
from Classes.MovingBedTests import MovingBedTests

class Measurement(object):
    """Class to hold measurement details for use in the GUI
        Analogous to matlab class: clsMeasurement
    """

    def __init__(self, in_file, source='TRDI'):
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
        for x in self.transects:
            x.wt_filters('wt_depth', 'On')
            x.wt_interpolations('Ensembles', 'None')    
        
            
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
        transect.boat_vel.bt_vel.__apply_filter(transect, kargs=bt_settings)
    
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
        
    
    
    
    
    
    
        
        
    
            







