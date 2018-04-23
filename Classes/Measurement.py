from Classes.MMT_TRDI import MMT_TRDI
import numpy as np
import os
import scipy.io as sio
from Classes.TransectData import TransectData, allocate_transects
from Classes.Pd0TRDI import Pd0TRDI
from Classes.PreMeasurement import PreMeasurement
from Classes.MovingBedTests import MovingBedTests
from Classes.MultiThread import MultiThread
from Classes.QComp import QComp
from Classes.MatSonTek import MatSonTek
from Classes.CompassCal import CompassCal
from Classes.SystemTest import SystemTest
# from Classes.NormData import NormData
from Classes.ComputeExtrap import ComputeExtrap
from Classes.ExtrapQSensitivity import ExtrapQSensitivity

class Measurement(object):
    """Class to hold measurement details for use in the GUI
        Analogous to matlab class: clsMeasurement
    """

    # def __init__(self, in_file, source='TRDI', **kargs):
    def __init__(self, in_file, source, proc_type='QRev', checked=False):

        self.station_name = None
        self.station_number = None
        self.transects = []
        self.mb_tests = []
        self.system_test = []
        self.compass_cal = []
        self.compass_eval = []
        self.ext_temp_chk = None
        self.extrap_fit = None
        self.processing = None
        self.comments = None
        self.discharge = []
        self.uncertainty = None
        self.initial_settings = None
        self.qa = None
        self.userRating = None
        self.comments = []
        self.ext_temp_chk = {}

        if source == 'QRev':
            self.load_qrev_mat(fullname=in_file)
        else:
            if source == 'TRDI':
                self.load_trdi(in_file, checked=checked)

            elif source == 'SonTek':
                self.load_sontek(in_file)

            select = self.transects[0].boat_vel.selected
            if select == 'bt_vel':
                ref = 'BT'
            elif select == 'gga_vel':
                ref = 'GGA'
            elif select == 'vtg_vel':
                ref = 'VTG'

            if self.mb_tests is not None:
                self.mb_tests = MovingBedTests.auto_use_2_correct(moving_bed_tests=self.mb_tests)

            # Save initial settings
            self.initial_settings = self.current_settings()

            # Set processing type
            if proc_type == 'QRev':
                # Apply QRev default settings
                settings = self.QRevDefaultSettings()
                settings['Processing'] = 'QRev'
                self.apply_settings(settings)

            elif proc_type == 'None':
                #UpdateStatus "Processing with no filters and interpolation
                settings = self.NoFilterInterpSettings()
                settings['Processing'] = 'None'

                self.apply_settings(settings)

            else:
                self.discharge = QComp()
        
        

    # DSM change 1/23/2018 def load_trdi(self, mmt, **kargs):
    def load_trdi(self, mmt_file, type='Q', checked=False):
        """Method to load TRDI data.

        Parameters
        ----------
        mmt_file: str
            Full pathname to mmt file.
        type: str
            Type of data (Q: discharge, MB: moving-bed test
        checked: bool
            Determines if all files are loaded (False) or only checked (True)
        """

        mmt = MMT_TRDI(mmt_file)

        # Get properties if they exist, otherwise set them as blank strings
        self.station_name = str(mmt.site_info['Name'])
        self.station_number = str(mmt.site_info['Number'])

        # DSM Not sure this is needed 1/19/2018
        self.processing = 'WR2'

        # Create transect objects for  TRDI data
        self.transects = allocate_transects(mmt=mmt, type=type, checked=checked)

        #-------------------------Done Refactor----------------------------------------------
        #Create object for pre-measurement tests
        if isinstance(mmt.qaqc, dict) or isinstance(mmt.mbt_transects):
            self.qaqc_trdi(mmt)
        
        #Save comments from mmt file in comments
        try:
            self.comments.append('MMT Remarks: ' + mmt.site_info.remarks)
        except:
            pass


        for t in range(len(self.transects)):
            notes = getattr(mmt.transects[t], 'Notes')
            for x in range(len(notes)):
                note = 'Transect: ' + t + ' File: ' + x + ' ' + notes[x].NoteDate + ': ' + notes[x].NoteText
                self.comments.append(note)
                
        #Get external temperature
        self.ext_temp_chk['User'] = mmt.site_info['Water_Temperature']
        
        
        #-------------------------------REFACTORED from thresholds_TRDI
        
        threshold_settings = {}
        threshold_settings['wt_settings'] = {}
        threshold_settings['bt_settings'] = {}
        threshold_settings['depth_settings'] = {}
        #WT filter threshold setting
        threshold_settings['wt_settings']['beam'] = self.set_3_beam_WT_threshold_TRDI(mmt.transects[0])
        threshold_settings['wt_settings']['difference'] = 'Manual'
        threshold_settings['wt_settings']['difference_threshold'] = self.set_diff_vel_threshold_WT_TRDI(mmt.transects[0])
        threshold_settings['wt_settings']['vertical'] = 'Manual'
        threshold_settings['wt_settings']['vertical_threshold'] = self.set_vert_vel_threshold_WT_TRDI(mmt.transects[0])

        #BT filter threshold settings
        threshold_settings['bt_settings']['beam'] = self.set_3_beam_BT_threshold_TRDI(mmt.transects[0])
        threshold_settings['bt_settings']['difference'] = 'Manual'
        threshold_settings['bt_settings']['difference_threshold'] = self.set_diff_vel_threshold_BT_TRDI(mmt.transects[0])
        threshold_settings['bt_settings']['vertical'] = 'Manual'
        threshold_settings['bt_settings']['vertical_threshold'] =  self.set_vert_vel_threshold_BT_TRDI(mmt.transects[0])

        #Depth filter and averaging settings
        threshold_settings['depth_settings']['depth_weighting'] = self.set_depth_weighting_TRDI(mmt.transects[0])
        threshold_settings['depth_settings']['depth_valid_method'] = 'TRDI'
        threshold_settings['depth_settings']['depth_screening'] = self.set_depth_screening_TRDI(mmt.transects[0])
        
        #--------------------------------------------DONE REFACTOR
        multi_threaded = False
        
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

            # Determine reference used in WR2 if available
            if 'Reference' in mmt.site_info.keys():
                reference = mmt.site_info['Reference']
            else:
                reference = "BT"

            transect_idx = 0
            # Convert to earth coordinates
            for transect in self.transects:
                # Convert to earth coordinates
                transect.change_coord_sys(new_coord_sys='Earth')
                # Set navigation reference
                transect.change_nav_reference(update=False, new_nav_ref=reference)
                # Apply WR2 thresholds
                self.thresholds_TRDI(transect, threshold_settings)
                # Apply boat interpolations
                transect.boat_interpolations(update=False, target='BT', method='None')
                if transect.gps is not None:
                    transect.boat_interpolations(update=False, target='GPS', method='HoldLast')
                # Update water data for changes in boat velocity
                transect.update_water()
                # Filter water data
                transect.wt_filters(wt_depth=True)
                # Interpolate water data
                transect.wt_interpolations(target='Ensembles', interp_type='None')
                # Apply speed of sound computations as required
                mmt_sos_method = mmt.transects[transect_idx].active_config['Proc_Speed_of_Sound_Correction']
                # Speed of sound computed based on user supplied values
                if mmt_sos_method == 1:
                    transect.change_sos(parameter='salinity')
                # Speed of sound set by user
                elif mmt_sos_method == 2:
                    speed = mmt.transects[transect_idx].active_config['Proc_Fixed_Speed_of_Sound']
                    transect.change_sos(parameter='sosSrc', selected='user', speed=speed)
                transect_idx += 1

    def qaqc_trdi(self, mmt):
        """Processes qaqc test, calibrations, and evaluations
        
        Input:
        mmt: object of MMT_TRDI
        
        """
        #ADCP Test
        if 'RG_Test' in mmt.qaqc:
            for n in range(len(mmt.qaqc['RG_Test']['TestResult'])):
                p_m = PreMeasurement()
                p_m.populate_data(mmt.qaqc['RG_Test']['TestResult'][n]['TimeStamp'],
                                  mmt.qaqc['RG_Test']['TestResult'][n]['Text'],'TST')
                self.system_test.append(p_m)
        else:
            self.system_test.append(PreMeasurement())

        #Compass calibration
        if 'Compass_Calibration' in mmt.qaqc:
            for n in range(len(mmt.qaqc['Compass_Calibration']['TestResult'])):
                p_m = PreMeasurement()
                p_m.populate_data(mmt.qaqc['Compass_Calibration']['TestResult'][n]['TimeStamp'],
                                  mmt.qaqc['Compass_Calibration']['TestResult'][n]['Text'], 'TCC')
                self.compass_cal.append(p_m)
        else:
            self.compass_cal.append(PreMeasurement())
            
        #Compass evaluation
        if 'Compass_Evaluation' in mmt.qaqc:
            for n in range(len(mmt.qaqc['Compass_Evaluation']['TestResult'])):
                p_m =  PreMeasurement()
                p_m.populate_data(mmt.qaqc['Compass_Evaluation']['TestResult'][n]['TimeStamp'],
                                  mmt.qaqc['Compass_Evaluation']['TestResult'][n]['Text'], 'TCC')
                self.compass_eval.append(p_m)
        else:
            self.compass_eval.append(PreMeasurement())
            
       
        
        if len(mmt.mbt_transects) > 0:
            
            #Create transect object/s
            transects= allocate_transects('TRDI', mmt, type='MB')
            if len(transects) > 0:
                self.mb_tests = []
                
                for n in range(len(transects)):
                    mb_test = MovingBedTests()
                    # TODO need to check type for compatibility with international data
                    mb_test.populate_data('TRDI', transects[n], mmt.mbt_transects[n].moving_bed_type)
                    
                    #Save notes from mmt files in comments
                    if 'NoteDate' in mmt.mbt_transects:
                        r, c = mmt.mbt_transects['NoteDate'].shape
                        for n in range(r):
                            for k in range(c):
                                if mmt.mbt_transects['NoteDate'][n,k] is not None:
                                    pass
                                
                    self.mb_tests.append(mb_test)

    def thresholds_TRDI(self, transect, settings):
        """Retrieve and apply manual filter settings from mmt file

        Input:
        mmt: object of MMT_TRDI
        transect: object of TransectData
        settings: threshold settings computed before processing
        """

        # Apply WT settings
        transect.w_vel.apply_filter(transect, **settings['wt_settings'])

        # Apply BT settings
        transect.boat_vel.bt_vel.apply_filter(transect, **settings['bt_settings'])

        # Apply depth settings
        transect.depths.set_valid_data_method(setting=settings['depth_settings']['depth_valid_method'])
        transect.depths.depth_filter(transect=transect, filter_method=settings['depth_settings']['depth_screening'])
        transect.depths.bt_depths.compute_avg_BT_depth(method=settings['depth_settings']['depth_weighting'])
        # Apply composite depths as per setting stored in transect from TransectData.trdi
        transect.depths.composite_depths(transect)

    def load_sontek(self, fullnames):
        """Coordinates reading of all SonTek data files.

        Parameters
        ----------
        fullnames: list
            File names including path for all discharge transects converted to Matlab files.
        """

        for file in fullnames:
            # Read data file
            rsdata = MatSonTek(file)
            pathname, file_name = os.path.split(file)

            # Create transect objects for each discharge transect
            self.transects.append(TransectData())
            self.transects[-1].SonTek(rsdata ,file_name)

        # Site information pulled from last file
        if hasattr(rsdata, 'SiteInfo'):
            if hasattr(rsdata.SiteInfo, 'Site_Name'):
                self.station_name = rsdata.SiteInfo.Site_Name
            if hasattr(rsdata.SiteInfo, 'Station_Number'):
                self.station_number = rsdata.SiteInfo.Station_Number

        self.qaqc_sontek(pathname)

        for transect in self.transects:
            transect.change_coord_sys(new_coord_sys='Earth')
            transect.change_nav_reference(update=False, new_nav_ref=self.transects[0].boat_vel.selected)
            transect.boat_interpolations(update=False, target='BT', method='Hold9')
            transect.boat_interpolations(update=False, target='GPS', method='None')
            transect.apply_averaging_method(setting='Simple')
            transect.process_depths(update=False, interpolation_method='HoldLast')
            transect.update_water()
            transect.wt_filters(wt_depth=True)
            transect.wt_interpolations(target='Ensembles', interp_type='None')
            transect.wt_interpolations(target='Cells', interp_type='TRDI')

    def qaqc_sontek(self, pathname):
        """Reads and stores system tests, compass calibrations, and moving-bed tests.

        Parameters
        ----------
        pathname: str
            Path to discharge transect files.
        """

        # Compass Calibration
        compass_cal_folder = os.path.join(pathname, 'CompassCal')
        if os.path.isdir(compass_cal_folder):
            compass_cal_files=[]
            for file in os.listdir(compass_cal_folder):
                # G3 compasses
                if file.endswith('.ccal'):
                    compass_cal_files.append(file)
                # G2 compasses
                elif file.endswith('.txt'):
                    compass_cal_files.append(file)
            for file in compass_cal_files:
                time_stamp = file.split('_')[1].split('.')
                time_stamp = time_stamp[0] + ':' + time_stamp[1] + ':' + time_stamp[2]
                with open(os.path.join(compass_cal_folder, file)) as f:
                    cal_data = f.read()
                    cal = CompassCal()
                    cal.populate_data(time_stamp, cal_data, 'SSC')
                    self.compass_cal.append(cal)

        # System Test
        system_test_folder = os.path.join(pathname, 'SystemTest')
        if os.path.isdir(system_test_folder):
            for file in os.listdir(system_test_folder):
                # Find system test files.
                if file.startswith('SystemTest'):
                    with open(os.path.join(system_test_folder, file)) as f:
                        test_data = f.read()
                        test_data = test_data.replace('\x00','')
                    time_stamp = file[18:20] + ':' + file[20:22] + ':' + file[22:24]
                    sys_test = SystemTest()
                    sys_test.populate_data(time_stamp=time_stamp, data_in=test_data, data_type='SST')
                    self.system_test.append(sys_test)

        # Moving-bed tests
        self.sontek_moving_bed_tests(pathname)

    def sontek_moving_bed_tests(self, pathname):
        """Locates and processes SonTek moving-bed tests.

        Searches the pathname for Matlab files that start with Loop or SMBA.
        Processes these files as moving bed tests.

        Parameters
        ----------
        pathname: str
            Path to discharge transect files.
        """
        for file in os.listdir(pathname):
            # Find moving-bed test files.
            if file.endswith('.mat'):
                # Process Loop test
                if file.lower().startswith('loop'):
                    self.mb_tests.append(MovingBedTests())
                    self.mb_tests[-1].populate_data(source='SonTek',
                                                    file=os.path.join(pathname, file),
                                                    type='Loop')
                # Process Stationary test
                elif file.lower().startswith('smba'):
                    self.mb_tests.append(MovingBedTests())
                    self.mb_tests[-1].populate_data(source='SonTek',
                                                    file=os.path.join(pathname, file),
                                                    type='Stationary')

    def load_qrev_mat(self, fullname):
        """Loads and coordinates the mapping of existing QRev Matlab files into Python instance variables.

        Parameters
        ----------
        fullname: str
            Fullname including path to *_QRev.mat files.
        """

        # Read Matlab file and extract meas_struct
        mat_data = sio.loadmat(fullname, struct_as_record=False, squeeze_me=True)
        meas_struct=mat_data['meas_struct']

        # Assign data from meas_struct to associated instance variables in Measurement and associated objects.
        self.station_name = meas_struct.stationName
        self.station_number = meas_struct.stationNumber
        self.processing = meas_struct.processing
        self.comments = meas_struct.comments
        self.userRating = meas_struct.userRating
        self.initial_settings = meas_struct.initialSettings

        self.transects = []
        self.mb_tests = []
        self.system_test = []
        self.compass_cal = []
        self.compass_eval = []
        self.ext_temp_chk = None
        self.extrap_fit = None

        # Discharge
        if hasattr(meas_struct.discharge, 'bottom'):
            # Measurement has discharge data from only one transect
            bottom = meas_struct.discharge.bottom
            q=QComp()
            q.populate_from_qrev_mat(meas_struct.discharge)
            self.discharge.append(q)
        else:
            # Measurement has discharge data from multiple transects
            for q_data in meas_struct.discharge:
                q=QComp()
                q.populate_from_qrev_mat(q_data)
                self.discharge.append(q)

        self.uncertainty = None

        self.qa = None

        self.comments = []
        self.ext_temp_chk = {}

    def set_3_beam_WT_threshold_TRDI(self, mmt_transect):
        """Get 3-beam processing for WT from mmt file
        
        Input:
        mmt_transect: object of Transect
        
        Output:
        num_3_beam_WT_Out
        """
        
        #check consistency
        use_3_beam_WT = mmt_transect.active_config['Proc_Use_3_Beam_Solution_For_WT']
        
        #Use setting from 1st transect for all transects
        if use_3_beam_WT == 0:
            num_beam_WT_out = 4
        else:
            num_beam_WT_out = 3
            
        return num_beam_WT_out

    def set_3_beam_BT_threshold_TRDI(self, mmt_transect):
        """Get 3-beam processing for WT from mmt file

        Input:
        mmt_transect: object of Transect

        Output:
        num_3_beam_WT_Out
        """

        #check consistency
        use_3_beam_BT = mmt_transect.active_config['Proc_Use_3_Beam_Solution_For_BT']

        #Use setting from 1st transect for all transects
        if use_3_beam_BT == 0:
            num_beam_BT_out = 4
        else:
            num_beam_BT_out = 3

        return num_beam_BT_out
    
    def set_vert_vel_threshold_WT_TRDI(self, mmt_transect):
        """Get the vertical celocity threshold for WT from mmt
        
        Input:
        mmt_transect: object of Transect
        
        Output:
        vert_vel_threshold_WT[0]: vertical celocity threshold (m/s)
        """
        
        #Check consistency of vertical velocity threshold (m/s)
        vert_vel_threshold_WT = mmt_transect.active_config['Proc_WT_Up_Velocity_Threshold']
        
        #Use setting from 1st transect for all transects
        #Use setting from 1st transect for all transects
        if isinstance(vert_vel_threshold_WT, float):
            return vert_vel_threshold_WT 
        else:
            return vert_vel_threshold_WT[0]

    def set_vert_vel_threshold_BT_TRDI(self, mmt_transect):
        """Get the vertical velocity threshold for BT from mmt

        Input:
        mmt_transect: object of Transect

        Output:
        vert_vel_threshold_WT[0]: vertical celocity threshold (m/s)
        """

        # Check consistency of vertical velocity threshold (m/s)
        vert_vel_threshold_BT = mmt_transect.active_config['Proc_BT_Up_Velocity_Threshold']

        # Use setting from 1st transect for all transects
        # Use setting from 1st transect for all transects
        if isinstance(vert_vel_threshold_BT, float):
            return vert_vel_threshold_BT
        else:
            return vert_vel_threshold_BT[0]
    
    def set_diff_vel_threshold_WT_TRDI(self, mmt_transect):
        """Get difference (error) velocity threshold for WT from mmt
        
        Input:
        mmt_transect: object of Transect
        
        Output:
        diff_vel_threshold_WT[0]: difference velocity threshold (m/s)
        """
        
        #Check consistency of difference (error) velocity for WT
        diff_vel_threshold_WT = mmt_transect.active_config['Proc_WT_Error_Velocity_Threshold']
        
        #Use setting from 1st transect for all transects
        if isinstance(diff_vel_threshold_WT, float):
            return diff_vel_threshold_WT 
        else:
            return diff_vel_threshold_WT[0]
    
    def set_diff_vel_threshold_BT_TRDI(self, mmt_transect):
        """Get difference (error) velocity threshold for BT from mmt
        
        Input:
        mmt_transect: object of Transect
        
        Output:
        diff_vel_threshold_BT[0]: difference velocity threshold (m/s)
        """
        
        #Check consistency of difference (error) velocity for BT
        diff_vel_threshold_BT = mmt_transect.active_config['Proc_BT_Error_Velocity_Threshold']
        
        #Use setting from 1st transect for all transects
        if isinstance(diff_vel_threshold_BT, float):
            return diff_vel_threshold_BT
        else:
            return diff_vel_threshold_BT[0]
    
    def set_depth_weighting_TRDI(self, mmt_transect):
        """Get the average depth method from mmt
        
        Input:
        mmt_transect: Object of Transect
        
        Output:
        depth_weighting_setting: method for computing average depth
        """
        
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
        """Get the depth screening setting from mmt
        
        Input:
        mmt_transect: object of Transect
        
        Output:
        depth_screening_setting: depth screening setting
        """
        
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
        """Applies a change in speed of sound to one or all transects
        and update the discharge and uncertainty computations
        
        Input:
        kargs[0]: optional, transect number to apply change to, if empty
        change is applied to all transects
        kargs[0/1]: sensor property to change
        kargs[1/2]: settinf for sensor property
        """
        
        s = self.current_settings()
        transect.change_SOS(kargs = kargs[1:])    
        #Check to see if transect number specified
        self.apply_settings(transect, s)
        
    def apply_settings(self, settings, transect=None):
        """Applies reference, filter, and interpolation settings.
        
        Parameters
        ----------
        transect: object
            Object of TransectData
        settings: dict
            Dictionary of reference, filter, and interpolation settings
        """

        for transect in self.transects:
            # Moving-boat ensembles
            if 'Processing' in settings.keys():
                transect.change_q_ensembles(proc_method=settings['Processing'])
                self.processing = settings['Processing']

            # Navigation reference
            if transect.boat_vel.selected != settings['NavRef']:
                transect.change_nav_reference(update=False, new_nav_ref=settings['NavRef'])
                if len(self.mb_tests) > 0:
                    self.mb_tests.auto_use_2_correct(moving_bed_tests=self.mb_tests, boat_ref=settings['NavRef'])

            # Composite tracks
            # Changing the nav reference applies the current setting for
            # Composite tracks, check to see if a change is needed
            if transect.boat_vel.composite != settings['CompTracks']:
                transect.composite_tracks(update=False, setting=settings['CompTracks'])

            # Set difference velocity BT filter
            bt_kwargs = {}
            if settings['BTdFilter'] == 'Manual':
                bt_kwargs['difference'] = settings['BTdFilter']
                bt_kwargs['difference_threshold'] = settings['BTdFilterThreshold']
            else:
                bt_kwargs['difference'] = settings['BTdFilter']

            # Set vertical velocity BTfilter
            if settings['BTwFilter'] == 'Manual':
                bt_kwargs['vertical'] = settings['BTwFilter']
                bt_kwargs['vertical_threshold'] = settings['BTwFilterThreshold']
            else:
                bt_kwargs['vertical'] = settings['BTwFilter']

            # Apply BT settings
            transect.boat_filters(update=False,**bt_kwargs)

            # BT Interpolation
            transect.boat_interpolations(update=False, target='BT', method=settings['BTInterpolation'])

            # GPS filter settings
            if transect.gps is not None:
                gga_kwargs = {}
                if transect.boat_vel.gga_vel is not None:
                    # GGA
                    if settings['ggaAltitudeFilter'] == 'Manual':
                        gga_kwargs['altitude'] = settings['ggaAltitudeFilter']
                        gga_kwargs['altitude_threshold'] = settings['ggaAltitudeFilterChange']
                    else:
                        gga_kwargs['altitude'] = settings['ggaAltitudeFilter']

                    # Set GGA HDOP Filter
                    if settings['GPSHDOPFilter'] == 'Manual':
                        gga_kwargs['hdop'] = settings['GPSHDOPFilter'],
                        gga_kwargs['hdop_max_threshold'] = settings['GPSHDOPFilterMax'],
                        gga_kwargs['hdop_change_threshold'] = settings['GPSHDOPFilterChange']
                    else:
                        gga_kwargs['hdop'] = settings['GPSHDOPFilter']

                    gga_kwargs['other'] = settings['GPSSmoothFilter']
                    # Apply GGA filters
                    transect.gps_filters(update=False, **gga_kwargs)

                if transect.boat_vel.vtg_vel is not None:
                    vtg_kwargs = {}
                    if settings['GPSHDOPFilter'] == 'Manual':
                        vtg_kwargs['hdop'] = settings['GPSHDOPFilter']
                        vtg_kwargs['hdop_max_threshold'] = settings['GPSHDOPFilterMax']
                        vtg_kwargs['hdop_change_threshold'] = settings['GPSHDOPFilterChange']
                        vtg_kwargs['other'] = settings['GPSSmoothFilter']
                    else:
                        vtg_kwargs['hdop'] = settings['GPSHDOPFilter']
                        vtg_kwargs['other'] = settings['GPSSmoothFilter']

                    # Apply VTG filters
                    transect.gps_filters(update=False, **vtg_kwargs)

                transect.boat_interpolations(update=False, target='GPS', method=settings['GPSInterpolation'])

            # Set depth reference
            transect.set_depth_reference(update=False, setting=settings['depthReference'])

            transect.process_depths(update=True,
                                    filter_method=settings['depthFilterType'],
                                    interpolation_method=settings['depthInterpolation'],
                                    composite_setting=settings['depthComposite'],
                                    avg_method=settings['depthAvgMethod'],
                                    valid_method=settings['depthValidMethod'])

            # Set WT difference velocity filter
            wt_kwargs = {}
            if settings['WTdFilter'] == 'Manual':
                wt_kwargs['difference'] = settings['WTdFilter']
                wt_kwargs['difference_threshold'] = settings['WTdFilterThreshold']
            else:
                wt_kwargs['difference'] = settings['WTdFilter']

            # Set WT vertical velocity filter
            if settings['WTwFilter'] == 'Manual':
                wt_kwargs['vertical'] = settings['WTwFilter']
                wt_kwargs['vertical_threshold'] = settings['WTwFilterThreshold']
            else:
                wt_kwargs['vertical'] = settings['WTwFilter']

            wt_kwargs['beam'] = settings['WTbeamFilter']
            wt_kwargs['other'] = settings['WTsmoothFilter']
            wt_kwargs['snr'] = settings['WTsnrFilter']
            wt_kwargs['wt_depth'] = settings['WTwDepthFilter']
            wt_kwargs['excluded'] = settings['WTExcludedDistance']

            transect.wt_filters(**wt_kwargs)

            # Edge methods
            transect.edges.rec_edge_method = settings['edgeRecEdgeMethod']
            transect.edges.vel_method = settings['edgeVelMethod']

        # Recompute extrapolations
        # NOTE: Extrapolations should be determined prior to WT
        # interpolations because the TRDI approach for power/power
        # using the power curve and exponent to estimate invalid cells.

        if (self.extrap_fit is None) or (self.extrap_fit.fit_method == 'Automatic'):
            self.extrap_fit = ComputeExtrap()
            self.extrap_fit.populate_data(transects=self.transects,compute_sensitivity=False)
            top = self.extrap_fit.sel_fit[-1].top_method
            bot = self.extrap_fit.sel_fit[-1].bot_method
            exp = self.extrap_fit.sel_fit[-1].exponent
            self.change_extrapolation(top=top, bot=bot, exp=exp)
        else:
            if 'extrapTop' not in settings.keys():
                top = self.extrap_fit.sel_fit[-1].top_method
                bot = self.extrap_fit.sel_fit[-1].bot_method
                exp = self.extrap_fit.sel_fit[-1].exponent

            self.change_extrapolation(top=settings['extrapTop'], bot=settings['extrapBot'], exp=settings['extrapExp'])

        for transect in self.transects:
            # Water track interpolations
            transect.w_vel.apply_interpolation(transect=transect,
                                               target='Ensembles',
                                               interp_type=settings['WTEnsInterpolation'])
            transect.w_vel.apply_interpolation(transect=transect,
                                               target='Cells',
                                               interp_type=settings['WTCellInterpolation'])


            q = QComp()
            q.populate_data(data_in=transect)
            self.discharge.append(q)

        print('Yea!')
        # TODO add uncertainty
        # TODO add quality assurance

    def apply_settings2(self, trans_data, s, kargs=None):
        """Refactored so that the previous calculations can be multithreaded
        #Extrap Q Sensitivity the way it is now will not allow this so in this
        method trans_data is all transects instead of just one"""
        
        # Recompute extrapolations
        # NOTE: Extrapolations should be determined prior to WT
        # interpolations because the TRDI approach for power/power
        # using the power curve and exponent to estimate invalid cells.
        if self.extrap_fit is None or self.extrap_fit.fit_method == 'Automatic':
            # self.extrapfit = ComputeExtrap()
            self.extrapfit.populate_data(trans_data, 0)
            top = self.extrap_fit.sel_fit.top_method
            bot = self.extrap_fit.sel_fit.bot_method
            exp = self.extrap_fit.sel_fit.exponent

            self.set_extrapolation(trans_data, top, bot, exp)

        else:
            if 'extrapTop' in s:
                s['extrapTop'] = self.extrap_fit.sel_fit.top_method
                s['extrapBot'] = self.extrap_fit.sel_fit.bot_method
                s['extrapExp'] = self.extrap_fit.sel_fit.exponent
                
            self.set_extrapolation(trans_data, s['extrapTop'], s['extrapBot'], s['extrapExp'])
        
            self.extrap_fit.change_fit_method(trans_data, 'Manual', 'All', s['extrapTop'], s['extrapBot'], s['extrapExp'])
            
        for n in trans_data:
            n.wt_interpolations('Ensembles', s['WTEnsInterpolation'])
            n.wt_interpolations('Cells', s['WTCellInterpolation'])
            
        #Edge methods
        for n in trans_data:
            n.edges._EdgeData__vel_method = s['edgeVelMethod']
            n.edges._EdgeData__rec_edge_method = s['edgeRecEdgeMethod']
            
        #Update sensitivities for water track interpolations
        self.extrap_fit.update_q_sensitivity(trans_data)

    def current_settings(self):
        """Saves the current settings for a measurement. Since all settings
        in QRev are consistent among all transects in a measurement onlyt the
        settings from the first transect are saved
        """
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
        settings['WTExcludedDistance'] = transect.w_vel.excluded_dist_m
        
        #Bottom track settings
        settings['BTbeamFilter'] = self.transects[first_idx].boat_vel.bt_vel.beam_filter;
        settings['BTdFilter'] = self.transects[first_idx].boat_vel.bt_vel.d_filter;
        settings['BTdFilterThreshold'] = self.transects[first_idx].boat_vel.bt_vel.d_filter_threshold;
        settings['BTwFilter'] =self.transects[first_idx].boat_vel.bt_vel.w_filter;
        settings['BTwFilterThreshold'] = self.transects[first_idx].boat_vel.bt_vel.w_filter_threshold;
        settings['BTsmoothFilter'] = self.transects[first_idx].boat_vel.bt_vel.smooth_filter;
        settings['BTInterpolation'] = self.transects[first_idx].boat_vel.bt_vel.interpolate
        
        #Gps Settings
        if transect.gps is not None:
            #GGA settings
            if transect.boat_vel.gga_vel is not None:
                settings['ggaDiffQualFilter'] = transect.boat_vel.gga_vel.gps_diff_qual_filter
                settings['ggaAltitudeFilter'] = transect.boat_vel.gga_vel.gps_altitude_filter
                settings['ggaAltitudeFilterChange'] = transect.boat_vel.gga_vel.gps_altitude_filter_change
                settings['GPSHDOPFilter'] = transect.boat_vel.gga_vel.gps_HDOP_filter
                settings['GPSHDOPFilterMax'] = transect.boat_vel.gga_vel.gps_HDOP_filter_max
                settings['GPSHDOPFilterChange'] = transect.boat_vel.gga_vel.gps_HDOP_filter_change
                settings['GPSSmoothFilter'] = transect.boat_vel.gga_vel.smooth_filter
                settings['GPSInterpolation'] = transect.boat_vel.gga_vel.interpolate
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
            settings['GPSHDOPFilter'] = transect.boat_vel.vtg_vel.gps_HDOP_filter
            settings['GPSHDOPFilterMax'] = transect.boat_vel.vtg_vel.gps_HDOP_filter_max
            settings['GPSHDOPFilterChange'] = transect.boat_vel.vtg_vel.gps_HDOP_filter_change
            settings['GPSSmoothFilter'] = transect.boat_vel.vtg_vel.smooth_filter
            settings['GPSInterpolation'] = transect.boat_vel.vtg_vel.interpolate
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
        
        """Depth settings are always applied to all available depth sources.
        Only those saved in the btDepths are used here but are applied to all sources"""
        settings['depthFilterType'] = transect.depths.bt_depths.filter_type
        settings['depthReference'] = transect.depths.selected
        settings['depthComposite'] = transect.depths.composite
        select = getattr(transect.depths, transect.depths.selected)
        settings['depthInterpolation'] = select.interp_type
        
        #Extrap Settings
        settings['extrapTop'] = transect.extrap.top_method
        settings['extrapBot'] = transect.extrap.bot_method
        settings['extrapExp'] = transect.extrap.exponent
        
        #Edge Settings
        settings['edgeVelMethod'] = transect.edges.vel_method
        settings['edgeRecEdgeMethod'] = transect.edges.rec_edge_method
        
        return settings

    def QRevDefaultSettings(self):
        '''QRev default and filter settings for a measurement'''

        settings = {}

        # Naviagation reference (NEED LOGIC HERE)
        settings['NavRef'] = self.transects[0].boat_vel.selected

        # Composite tracks
        settings['CompTracks'] = False

        # Water track filter settings
        # TODO change to -1 after testing
        settings['WTbeamFilter'] = 3
        settings['WTdFilter'] = 'Auto'
        settings['WTdFilterThreshold'] = np.nan
        settings['WTwFilter'] = 'Auto'
        settings['WTwFilterThreshold'] = np.nan
        settings['WTsmoothFilter'] = False
        if self.transects[0].adcp.manufacturer == 'TRDI':
            settings['WTsnrFilter'] = 'Off'
        else:
            settings['WTsnrFilter'] = 'Auto'
        temp = [x.w_vel for x in self.transects]
        excluded_dist = np.nanmin([x.excluded_dist_m for x in temp])
        if excluded_dist < 0.158 and self.transects[0].adcp.model == 'M9':
            settings['WTExcludedDistance'] = 0.16
        else:
            settings['WTExcludedDistance'] = excluded_dist

        # Bottom track filter settings
        settings['BTbeamFilter'] = -1
        settings['BTdFilter'] = 'Auto'
        settings['BTdFilterThreshold'] = np.nan
        settings['BTwFilter'] = 'Auto'
        settings['BTwFilterThreshold'] = np.nan
        settings['BTsmoothFilter'] = False

        # GGA Filter settings
        settings['ggaDiffQualFilter'] = 2
        settings['ggaAltitudeFilter'] = 'Auto'
        settings['ggaAltitudeFilterChange'] = np.nan

        # VTG filter settings
        settings['vtgsmoothFilter'] = np.nan

        # GGA and VTG filter settings
        settings['GPSHDOPFilter'] = 'Auto'
        settings['GPSHDOPFilterMax'] = np.nan
        settings['GPSHDOPFilterChange'] = np.nan
        settings['GPSSmoothFilter'] = False

        # Depth Averaging
        settings['depthAvgMethod'] = 'IDW'
        settings['depthValidMethod'] = 'QRev'

        # Depth Reference

        # Default to 4 beam depth average
        settings['depthReference'] = 'btDepths'
        # Depth settings
        settings['depthFilterType'] = 'Smooth'
        settings['depthComposite'] = True

        # Interpolation settings
        settings = self.QRevDefaultInterp(settings)

        # Edge settings
        settings['edgeVelMethod'] = 'MeasMag'
        settings['edgeRecEdgeMethod'] = 'Fixed'

        return settings

    def QRevDefaultInterp(self, settings):
        '''Adds QRev default interpolation settings to existing settings data structure

        INPUT:
        settings: data structure of reference and filter settings

        OUTPUT:
        settings: data structure with reference, filter, and interpolation settings
        '''

        settings['BTInterpolation'] = 'Linear'
        settings['WTEnsInterpolation'] = 'Linear'
        settings['WTCellInterpolation'] = 'TRDI'
        settings['GPSInterpolation'] = 'Linear'
        settings['depthInterpolation'] = 'Linear'
        settings['WTwDepthFilter'] = True

        return settings

    def change_extrapolation(self, top=None, bot=None, exp=None):

        if top is not None:
            top = self.extrap_fit.sel_fit[-1].top_method
        if bot is not None:
            bot = self.extrap_fit.sel_fit[-1].bot_method
        if exp is not None:
            exp = self.extrap_fit.sel_fit[-1].exponent

        for transect in self.transects:
            transect.set_extrapolation(top=top, bot=bot, exp=exp)

        self.extrap_fit.q_sensitivity = ExtrapQSensitivity()
        self.extrap_fit.q_sensitivity.populate_data(transects=self.transects, extrap_fits=self.extrap_fit.sel_fit)