'''
Created on Sep 26, 2017

@author: gpetrochenkov
'''
from Classes.Measurement import Measurement
import numpy as np
from Classes.Uncertainty import Uncertainty
import re
from MiscLibs.lowess import idx

class QAData(object):
    '''QAData Class for storing QA data such as compass and system tests.  Also test_NA_drop_preserves_levels
    as the container for all methods of evaluating the quality of the measurement'''
    
    def __init__(self):
        self.__q_run_threshold_caution = None # Caution threshold for interpolated discharge for a run of invalid ensembles, in percent.
        self.__q_run_threshold_warning = None # Warning threshold for interpolated discharge for a run of invalid ensembles, in percent.
        self.__q_total_threshold_warning = None # Warning threshold for total interpolated discharge for invalid ensembles, in percent.
        self.__q_total_threshold_caution = None
        self.__transects = None # Data structure for quality assurance checks of transects.
        self.__system_test = None # Data structure for quality assurance checks of system tests.
        self.__compass = None # Data structure for quality assurance checks of compass tests and evaluations.
        self.__temperature = None # Data structure for quality assurance checks of temperature comparisons and change.
        self.__movingbed = None # Data structure for quality assurance checks of moving-bed tests and conditions.
        self.__user = None # Data structure for quality assurance checks of user input data.
        self.__depths = None # Data structure for quality assurance checks of depth data.
        self.__bt_vel = None # Data structure for quality assurance checks of bottom track velocities.
        self.__gga_vel = None # Data structure for quality assurance checks of GGA boat velocities.
        self.__vtg_vel = None # Data structure for quality assurance checks of VTG boat velocities.
        self.__w_vel = None # Data structure for quality assurance checks of water track velocities.
        self.__extrapolation = None # Data structure for quality assurance checks of extrapolations.
        self.__edges = None # Data structure for quality assurance checks of edge discharge estimates.
        self.__measurement = None #Object of class measurement
        
    def populate_data(self, meas):
        self.__measurement = meas
        self.__q_run_threshold_caution = 3
        self.__q_run_threshold_warning = 5
        self.__q_total_threshold_caution = 10
        self.__q_total_threshold_warning = 25
        
        self.__transects_qa()
        self.__system_test_qa()
        self.__compass  
       
    def transects_qa(self):
        '''Apply quality checks to transects
        
        Input:
        meas: object of Measurement'''
        meas = self.__measurement
        
        self.__transects = TestMessages('transects')
        
        checked = meas.transects[:].checked == 1
        
        num_checked = np.sum(checked)
        if num_checked >= 1:
            datetime = meas.transects[checked].datetime
            total_duration = np.nansum(datetime.transect_duration_sec)
        else:
            total_duration = 0
        self.__transects.duration = 0
        if total_duration < 720:
            self.__transects.status = 'caution'
            self.__transects.messages.append('Transects: Duration of selected transects is less than 720')
            self.__transects.duration = 1
        
        elif num_checked == 0:
            self.__transect.status = 'warning'
            self.__transect.messages.append('Transects: No transect data selected')
            self.__transects.number = 1
            
        else:
            self.__transects.number = num_checked
            
            if num_checked == 2:
                #Determine transects used to compute discharge
                checked = [ x.checked == 1 for x in meas.transects]
                
                #Assign automatically generated uncertainties to properties
                unc = Uncertainty()
                unc.ran_uncert_q(meas.discharge[checked], 'total')
                cov = unc._Uncertainty__cov
                
                if cov > 2:
                    self.__transects.status = 'caution'
                    self.__transects.messages.append('Transects: Uncertainty would be reduced by additional transects')
                    self.__transects.uncertainty = 1
                    
            
            #Check for consistent sign
            if np.sum(checked) >= 1:
                sign_check = np.unique(np.sign([x.total for x in meas.discharge[checked]]))
                self.__transects.sign = 0
                if sign_check > 1:
                    self.__transects.status = 'warning'
                    self.__transects.messages.append('Transects: Sign of total discharge is not consistent. ' + \
                                                     'One or more start banks may be incorrect;')
                    self.__transects.sign = 2
                    
            #Check for reciprocal transects
            start_edge = [x.start_edge for x in meas.transects[checked]]
            num_left = np.sum(start_edge == 'Left')
            num_right = np.sum(start_edge == 'Right')
            self.__transects.recip = 0
            if num_left != num_right:
                self.__transects.status = 'warning'
                self.__transects.messages.append('Transects: Transects selected are not reciprocal transects')
                self.__transects.recip = 2
            
        #Check for zero discharge transects
        zero_Q = np.any([x.total == 0 for x in meas.discharge[checked]])
        if zero_Q == True:
            self.__transects.status = 'warning'
            self.__transects.messages.append('TRANSECTS: One or more transects have zero discharge')
            self.__transects.number = -1
            
            
    def system_test_QA(self, meas):
        '''Apply quality checks to system tests
        
        Input:
        meas: object of Measurement
        '''
        self.__system_test = TestMessages('system_test')
        
        if meas.sys_test[0].__time_stamp is None:
            #No system test present
            self.__system_test.status = 'warning'
            self.__system_test.messages.append('SYSTEM TEST: No system test')
        else:
            #System test present
            n_test = len(meas.sys_test)
            failed_tests = np.array([])
            for n in range(0,n_test):
                failed_idx = [m.start() for m in re.finditer('FAIL', meas.sys_test[n].data)]
                failed_tests = np.append(failed_tests, len(failed_idx))
                
            self.__system_test.n_failed_tests = np.nansum(failed_tests)
            self.__system_test.valid_tests = failed_tests == 0
            
            #If test id present, determine if there are any fails
            if sum(self.__system_test.valid_tests) < n_test:
                if sum(self.__system_test.val) > 0:
                    self.__system_test.append("System Test: One or more system test sets have at least one failed test")
                    self.__system_test.status = 'caution'
                else:
                    self.__system_test.append("SYSTEM TEST: All system test have at least one test that failed")
                    self.__system_test.status = 'warning'
            else:
                self.system_test = 'good'
                
                
    def compass_QA(self, meas):
        '''Apply quality checks to compass
        
        Input:
        meas: object of class Measurement
        '''
        
        messages = []
        
        #Does ADCP have a compass?
        heading = np.unique(meas.transects[0].sensors.heading_deg.internal.data)
        self.__compass = TestMessages('compass')
        if len(heading == 1) and heading == 0:
            #No Compass
            self.__compass = 'inactive'
            self.__compass.status1 = 'default'
            self.__compass.status2 = 'default'
        else:
            #ADCP has a compass. Is it needed?
            
            #Check for loop test
            n_tests = len(meas.mb_tests)
            loop = False
            
            for n in range(n_tests):
                if meas.mb_tests[n].type == 'Loop':
                    loop = True
                    
            #Check for GPS data
            gps = False
            if meas.transects.gps is not None:
                gps = True
                
            #Compass calibration required if loop or GPS used
            if loop or gps:
                #Compass required
                
                #Check ADCP manufacturer
                if meas.transects[0].adcp.manufacturer == 'SonTek':
                    #Sontek ADCP
                    
                    #Check for calibration
                    n = len(meas.compass_cal)
                    if meas.compass_cal[n].__time_stamp:
                    
                        #No compass calibration
                        self.__compass.messages.append('COMPASS: No compass calibration')
                        self.__compass.status1 = 'warning'
                        
                    else:
                        
                        #Compass calibrated
                        self.__compass.status1 = 'good'
                        
                else:
                    #TRDI ADCP
                    
                    #Check for calibration
                    n = len(meas.compass_cal)
                    if meas.compass_cal[n].__time_stamp is None:
                        
                        #No calibration, check for evaluation
                        n = len(meas.compass_eval)
                        if meas.compass_eval[n].__time_stamp:
                            
                            #No calibration or evaluation
                            self.__compass.messages.append('COMPASS: No compass calibration')
                            self.compass.status1 = 'warning'
                            
                        else:
                            
                            #Evaluation but no calibration
                            self.__compass.messages.append('Compass: no Compass calibration')
                            self.__compass.status1 = 'caution'
                    
                    else:
                        #Compass calibrated
                        
                        #Check for evaluation
                        n = len(meas.compass_eval)
                        if meas.compass_eval[n].__time_stamp is None:
                            
                            #No Evaluation
                            self.__compass.append('Compass: No compass evaluation')
                            self.__compass.status1 = 'caution'
                            
                        else:
                            
                            #Compass evaluated
                            self.__compass.status1 = 'Good'
                            
            else:
                n_c = len(meas.compass_cal)
                n_e = len(meas.compass_eval)
                if n_c <1 and n_e < 1:
                    self.__compass.status1 = 'default'
                else:
                    self.__compass.status1 = 'good'
                    
                
            self.__compass.status2 = 'good'
            #Get data from checked transects
            checked = [x.checked == 1 for x in meas.transects]
            if sum(checked) >= 1:
                depths = [x.depths for x in meas.transects[checked]]
                sensors = [x.sensors for x in meas.transects[checked]]
                
                mag_var = np.array([])
                for n in range(len(depths)):
                    heading_deg = getattr(sensors[n].heading_deg, sensors[n].heading_deg.selected)
                    mag_var = np.append(mag_var, heading_deg.mag_var)
                    
                #Check for magvar consistency
                mag_var_check = np.unique(round(mag_var, 3))
                self.__compass.mag_var = 0
                if len(mag_var_check) > 1:
                    self.__compass.messages.append('Compass: Magnetic variation is not consistent among transects')
                    self.__compass.status2 = 'caution'
                    self.__compass.magvar = 1
                    
                #Check for magvar with GPS
                if meas.transects.gps is None:
                    mag_var_check2 = np.abs(mag_var_check) > 0
                    idx = mag_var == 0
                    if len(idx) > 0:
                        self.__compass.append('COMPASS: MAgnetic variation is 0 and GPS data are present')
                        self.__compass.status2 = 'warning'
                        self.__compass.magvar = 2
                        self.compass.mag_var_idx = idx
                        
                #Check pitch and roll
                n_transects = len(meas.transects)
                cal_roll_idx = []
                cal_pitch_idx = []
                mag_error_idx = []
                for n in range(n_transects):
                    selected_pitch = meas.transects[n].pitch_deg.selected
                    p_select = getattr(meas.transects[n].sensors.pitch_deg, selected_pitch)
                    mean_pitch = np.nanmean(p_select.data)
                    std_pitch = np.nanstd(p_select.data)
                    selected_roll = meas.transects[n].roll_deg.selected
                    r_select = getattr(meas.transects[n].roll_deg, selected_roll)
                    mean_roll = np.nanmean(r_select.data)
                    std_roll = np.nanstd(r_select.data)
                    
                    call_roll_idx, call_pitch_idx = np.array([]), np.array([])
                    #If pitch and roll limit for compass calibration exist determine
                    #what transects have violations
                    if meas.transects[n].sensors.heading_deg.internal.roll_limit is None:
                        idx_max = r_select.data > meas.transects[n].sensors.heading_deg.internal.roll_limit[0]
                        idx_min = r_select.data < meas.transects[n].sensors.heading_deg.internal.roll_limit[1]
                        
                        idx = np.array([])
                        
                        if np.sum(idx_max) > 0:
                            idx = np.append(idx, idx_max)
                            
                        if np.sum(idx_min) > 0:
                            idx = np.append(idx, idx_min)
                            
                        if np.sum(idx) > 0:
                            call_roll_idx = np.append(call_roll_idx, n)
                            
                        idx_max = p_select.data > meas.transects[n].sensors.heading_deg.internal.pitch_limit[0]
                        idx_min = p_select.data < meas.transects[n].sensors.heading_deg.internal.pitch_limit[1]
                        
                        idx = np.array([])
                        
                        if np.sum(idx_max) > 0:
                            idx = np.append(idx, idx_max)
                            
                        if np.sum(idx_min) > 0:
                            idx = np.append(idx, idx_min)
                            
                        if np.sum(idx) > 0:
                            call_pitch_idx = np.append(call_pitch_idx, n)
                            
                    #If internal compass is used check magnetic error if available and determine if any
                    #transects exceed the threshold
                    
                    mag_error_idx = np.array([])
                    if meas.transects[n].sensors.heading_deg.selected == 'internal':
                        idx = np.array([])
                        if meas.transects[n].sensors.heading_deg.internal.mag_error is not None:
                            idx = meas.transects[n].sensors.heading_deg.internal.mag_error > 2
                        if np.sum(idx) > 0:
                            mag_error_idx = np.append(mag_error_idx, n)
                            
                    
                self.__compass.mean_pitch_idx_w = np.where(np.abs(mean_pitch) > 8)
                self.__compass.mean_pitch_idx_c = np.where((np.abs(mean_pitch) > 4) & (np.abs(mean_pitch) < 8))
                
                if len(self.__compass.mean_pitch_idx_w) > 0:
                    self.__compass.status2 = 'warning'
                    self.__compass.messages.append('PITCH: One or more transects have a mean pitch > 8 deg')
                elif len(self.__compass.mean_pitch_idx_c) > 0:
                    if self.__compass.status2 == 'good':
                        self.__compass.status2 = 'caution'
                        
                    self.__compass.messages.append('Pitch: Pne or more transects have a mean pitch > 4 deg')
                    
                self.__compass.mean_roll_idx_w = np.where(np.abs(mean_roll) > 8)
                self.__compass.mean_roll_idx_c = np.where((np.abs(mean_roll) > 4) & (np.abs(mean_roll) < 8))
                
                if len(self.__compass.mean_roll_idx_w) > 0:
                    self.__compass.status2 = 'warning'
                    self.__compass.messages.append('ROLL: One or more transects have a mean roll > 8 deg')
                elif len(self.__compass.mean_roll_idx_c) > 0:
                    if self.__compass.status2 == 'good':
                        self.__compass.status2 = 'caution'
                        
                    self.__compass.messages.append('Roll: One or more transects have a mean roll > 4 deg')
                    
                self.__compass.std_pitch_idx = np.where(np.nanmean(std_pitch) > 5)
                if len(self.__compass.std_pitch_idx) > 0:
                    if self.__compass.status2 == 'good':
                        self.__compass.status2 = 'caution'
                    self.__compass.messages.append('Pitch: One or more transects have a pitch std_dev > 5 deg')
                    
                        
                self.__compass.std_roll_idx = np.where(np.nanmean(std_roll) > 5)
                if len(self.__compass.std_roll_idx) > 0:
                    if self.__compass.status2 == 'good':
                        self.__compass.status2 = 'caution'
                    self.__compass.messages.append('Roll: One or more transects have a roll std_dev > 5 deg')       
                    
                    
                self.__compass.cal_pitch_idx = cal_pitch_idx
                if len(cal_pitch_idx) > 0:
                    if self.__compass.status2 == 'good':
                        self.__compass.status2 = 'caution'
                        
                    self.__compass.messages.append('Compass: One or more transects have pitch exceeding calibration limits')
                    
                self.__compass.cal_roll_idx = cal_roll_idx
                if len(cal_roll_idx) > 0:
                    if self.__compass.status2 == 'good':
                        self.__compass.status2 = 'caution'
                        
                    self.__compass.messages.append('Compass: One or more transects have roll exceeding calibration limits') 
                    
                self.__compass.mag_error_idx = mag_error_idx
                if len(mag_error_idx) > 0:
                    if self.__compass.status2 == 'good':
                        self.__compass.status2 = 'caution'
                        
                    self.__compass.messages.append('Compass: One or more transects have a change in magnetic field exceeding 2%') 
                    
                    
                if self.__compass.status1 == 'warning' or self.__compass.status2 == 'warning':
                    self.__compass.status = 'warning'
                elif self.__compass.status1 == 'caution' or self.__comopass.status2 == 'caution':
                    self.__compass.status = 'caution'
                else:
                    self.__compass.status = 'good'
                    
                    
    def temperature_QA(self, meas):
        '''Apply quality checks to compass
        Input:
        meas: object of Measurement
        '''
        
        #Run quality assurance checks
        self.temp_QA_check(meas)
        
        self.__temperature = TestMessages('temperature')
        #Set button color
        
        n_transects = len(meas.transects)
        checked = [x.checked for x in meas.transects]
        temp = np.array([])
        
        for x in meas.transects[checked]:
            temp_select = getattr(x.sensors.temperature_deg_c, x.sensors.temperature_deg_c.selected)
            temp = np.append(temp, temp_select.data)
            
        #compute the range of temperatures for the measurement
        range = np.nanmax(temp) - np.nanmin(temp)
        #Assign a quality code to the temperature stability
        self.__temperature.status = 'good'
        if range > 2:
            self.__temperature.messages.append('TEMPERATURE: Temperature is {3.1f} degrees C which is greater than 2 degrees'.format(round(range, 1)))
            self.__temperature.status = 'warning'
        elif range > 1:
            self.__temperature.messgaes.append('Temperature: Temperature is {3.1f} degrees C which is greater than 1 degree'.format(round(range, 1)))
            self.__temperature.status = 'caution'
        
        #Assign a quality code based on whether an independent temperature was entered
        if meas.ext_temp_chk.user is None:
            self.__temperature.messages.append('Temperature: No independent temperature reading')
            self.__temperature.status = 'caution'
        
        #Assign a quality code based on the difference between the mean ADCP and the independent temperature readings
        elif meas.ext_temp_chk.adcp is None:
            diff = np.abs(meas.ext_temp_chk.user - np.nanmean(temp))
            if diff >= 2:
                self.__temperature.messages('TEMPERATURE: The difference between the average ADCP and independent temperature is ' \
                                            + '{} degrees C which is not less than 2 degrees'.format(round(diff, 1)))
                self.__temperature.status = 'warning'
                
        #Assign a quality code based on the difference between the entered ADCP and the independent temperature readings
        else:
            diff = np.abs(meas.ext_temp_chk.user - meas.ext_temp_chk.adcp)
            if diff >= 2:
                self.__temperature.messages('TEMPERATURE: The difference between the average ADCP and independent temperature is ' \
                                            + '{} degrees C which is not less than 2 degrees'.format(round(diff, 1)))
                self.__temperature.status = 'warning'   
            
         
                                                   
        

 
class TestMessages(object):
    
    def __init__(self, message_type):
        self.type = message_type
        self.status = 'good'
        self.messages = []
        self.sign = 0
        self.duration = 0
        self.number = 0
        self.uncertainty = 0   
        self.recip = 0
        
        if message_type == 'system_test':
            self.n_failed_tests = 0
            self.valid_tests = 0 
            
        if message_type == 'compass':
            self.status1 = 'good'
            self.status2 = 'good'
            self.status = 'good'
            self.mag_var = None
            self.mean_pitch_idx_w = None
            self.mean_pitch_idx_c = None
            self.mean_roll_idx_w = None
            self.mean_roll_idx_c = None
            self.mag_var_idx = None
            self.std_pitch_idx = None
            self.std_roll_idx = None
            self.cal_pitch_idx = None
            self.cal_roll_idx = None
            self.mag_error_idx = None
            
            
        