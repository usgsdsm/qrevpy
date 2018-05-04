"""
Created on Sep 26, 2017

@author: gpetrochenkov
<<<<<<< HEAD
"""
=======
'''
import numpy as np
from Classes.Uncertainty import Uncertainty
import re
>>>>>>> 6ca6c50c231afa610ed3a693864074d7104a5f20


class QAData(object):
    """QAData Class for storing QA data such as compass and system tests.  Also test_NA_drop_preserves_levels
    as the container for all methods of evaluating the quality of the measurement"""
    
    def __init__(self):
        self.q_run_threshold_caution = None # Caution threshold for interpolated discharge for a run of invalid ensembles, in percent.
        self.q_total_threshold_warning = None # Warning threshold for total interpolated discharge for invalid ensembles, in percent.
        self.q_total_threshold_caution = None
        self.transects = None # Data structure for quality assurance checks of transects.
        self.system_test = None # Data structure for quality assurance checks of system tests.
        self.compass = None # Data structure for quality assurance checks of compass tests and evaluations.
        self.temperature = None # Data structure for quality assurance checks of temperature comparisons and change.
        self.movingbed = None # Data structure for quality assurance checks of moving-bed tests and conditions.
        self.user = None # Data structure for quality assurance checks of user input data.
        self.depths = None # Data structure for quality assurance checks of depth data.
        self.bt_vel = None # Data structure for quality assurance checks of bottom track velocities.
        self.gga_vel = None # Data structure for quality assurance checks of GGA boat velocities.
        self.vtg_vel = None # Data structure for quality assurance checks of VTG boat velocities.
        self.w_vel = None # Data structure for quality assurance checks of water track velocities.
        self.extrapolation = None # Data structure for quality assurance checks of extrapolations.
        self.edges = None # Data structure for quality assurance checks of edge discharge estimates.
        self.measurement = None #Object of class measurement
        
    def populate_data(self, meas):
        self.measurement = meas
        self.q_run_threshold_caution = 3
        self.q_run_threshold_warning = 5
        self.q_total_threshold_caution = 10
        self.q_total_threshold_warning = 25
        
        self.transects_qa()
        self.system_test_qa(meas)
        self.compass_qa(meas)
       
    def transects_qa(self):
        '''Apply quality checks to transects
        
        Input:
        meas: object of Measurement'''
        meas = self.measurement
        
        self.transects = TestMessages('transects')
        
        checked = [x.checked for x in meas.transects]
        
        num_checked = np.sum(checked)
        if num_checked >= 1:
            datetime = [x.datetime for x in meas.transects if x.checked == True]
            total_duration = np.nansum([x.transect_duration_sec for x in datetime])
        else:
            total_duration = 0
        self.transects.duration = 0
        if total_duration < 720:
            self.transects.status = 'caution'
            self.transects.messages.append('Transects: Duration of selected transects is less than 720')
            self.transects.duration = 1
        
        elif num_checked == 0:
            self.transect.status = 'warning'
            self.transect.messages.append('Transects: No transect data selected')
            self.transects.number = 1
            
        else:
            self.transects.number = num_checked
            
            if num_checked == 2:
                #Determine transects used to compute discharge
                checked = [ x.checked == 1 for x in meas.transects]
                
                #Assign automatically generated uncertainties to properties
                unc = Uncertainty(self)
                unc.ran_uncert_q(meas.discharge, 'total')
                cov = unc._Uncertainty__cov
                
                if cov > 2:
                    self.transects.status = 'caution'
                    self.transects.messages.append('Transects: Uncertainty would be reduced by additional transects')
                    self.transects.uncertainty = 1
                    
            
            #Check for consistent sign
            if np.sum(checked) >= 1:
                sign_check = np.unique(np.sign([x.total for x in meas.discharge[checked]]))
                self.transects.sign = 0
                if sign_check > 1:
                    self.transects.status = 'warning'
                    self.transects.messages.append('Transects: Sign of total discharge is not consistent. ' + \
                                                     'One or more start banks may be incorrect;')
                    self.transects.sign = 2
                    
            #Check for reciprocal transects
            start_edge = [x.start_edge for x in meas.transects[checked]]
            num_left = np.sum(start_edge == 'Left')
            num_right = np.sum(start_edge == 'Right')
            self.__transects.recip = 0
            if num_left != num_right:
                self.transects.status = 'warning'
                self.transects.messages.append('Transects: Transects selected are not reciprocal transects')
                self.transects.recip = 2
            
        #Check for zero discharge transects
        zero_Q = np.any([x.total == 0 for x,y in zip(meas.discharge,checked) if y == 1])
        if zero_Q == True:
            self.transects.status = 'warning'
            self.transects.messages.append('TRANSECTS: One or more transects have zero discharge')
            self.transects.number = -1
            
            
    def system_test_qa(self, meas):
        '''Apply quality checks to system tests
        
        Input:
        meas: object of Measurement
        '''
        self.system_test = TestMessages('system_test')
        
        if meas.sys_test[0]._PreMeasurement__time_stamp is None:
            #No system test present
            self.system_test.status = 'warning'
            self.system_test.messages.append('SYSTEM TEST: No system test')
        else:
            #System test present
            n_test = len(meas.sys_test)
            failed_tests = np.array([])
            for n in range(0,n_test):
                failed_idx = [m.start() for m in re.finditer('FAIL', meas.sys_test[n].data)]
                failed_tests = np.append(failed_tests, len(failed_idx))
                
            self.system_test.n_failed_tests = np.nansum(failed_tests)
            self.system_test.valid_tests = failed_tests == 0
            
            #If test id present, determine if there are any fails
            if sum(self.__system_test.valid_tests) < n_test:
                if sum(self.system_test.val) > 0:
                    self.system_test.append("System Test: One or more system test sets have at least one failed test")
                    self.system_test.status = 'caution'
                else:
                    self.system_test.append("SYSTEM TEST: All system test have at least one test that failed")
                    self.system_test.status = 'warning'
            else:
                self.system_test = 'good'
                
                
    def compass_qa(self, meas):
        '''Apply quality checks to compass
        
        Input:
        meas: object of class Measurement
        '''
        
        messages = []
        
        #Does ADCP have a compass?
        heading = np.unique(meas.transects[0].sensors.heading_deg.internal.data)
        self.compass = TestMessages('compass')
        if len(heading) == 1 and (heading == 0):
            #No Compass
            self.compass = 'inactive'
            self.compass.status1 = 'default'
            self.compass.status2 = 'default'
        else:
            #ADCP has a compass. Is it needed?
            
            #Check for loop test
            n_tests = len(meas.mb_tests)
            loop = False
            
            for n in range(n_tests):
                if meas.mb_tests[n]._MovingBedTests__type == 'Loop':
                    loop = True
                    
            #Check for GPS data
            gps = False
            if meas.transects[0].gps is not None:
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
                        self.compass.messages.append('COMPASS: No compass calibration')
                        self.compass.status1 = 'warning'
                        
                    else:
                        
                        #Compass calibrated
                        self.compass.status1 = 'good'
                        
                else:
                    #TRDI ADCP
                    
                    #Check for calibration
                    n = len(meas.compass_cal)
                    if meas.compass_cal[n-1]._PreMeasurement__time_stamp is None:
                        
                        #No calibration, check for evaluation
                        n = len(meas.compass_eval)
                        if meas.compass_eval[n-1]._PreMeasurement__time_stamp:
                            
                            #No calibration or evaluation
                            self.compass.messages.append('COMPASS: No compass calibration')
                            self.compass.status1 = 'warning'
                            
                        else:
                            
                            #Evaluation but no calibration
                            self.compass.messages.append('Compass: no Compass calibration')
                            self.compass.status1 = 'caution'
                    
                    else:
                        #Compass calibrated
                        
                        #Check for evaluation
                        n = len(meas.compass_eval)
                        if meas.compass_eval[n].__time_stamp is None:
                            
                            #No Evaluation
                            self.compass.append('Compass: No compass evaluation')
                            self.compass.status1 = 'caution'
                            
                        else:
                            
                            #Compass evaluated
                            self.compass.status1 = 'Good'
                            
            else:
                n_c = len(meas.compass_cal)
                n_e = len(meas.compass_eval)
                if n_c <1 and n_e < 1:
                    self.compass.status1 = 'default'
                else:
                    self.compass.status1 = 'good'
                    
                
            self.compass.status2 = 'good'
            #Get data from checked transects
            checked = [x.checked == 1 for x in meas.transects]
            if sum(checked) >= 1:
                depths = [x.depths for x in meas.transects if x.checked == 1]
                sensors = [x.sensors for x in meas.transects if x.checked == 1]
                
                mag_var = np.array([])
                for n in range(len(depths)):
                    heading_deg = getattr(sensors[n].heading_deg, sensors[n].heading_deg.selected)
                    mag_var = np.append(mag_var, heading_deg.mag_var_deg)
                    
                #Check for magvar consistency
                mag_var_check = np.unique(np.round(mag_var, 3))
                self.compass.mag_var = 0
                if len(mag_var_check) > 1:
                    self.__compass.messages.append('Compass: Magnetic variation is not consistent among transects')
                    self.__compass.status2 = 'caution'
                    self.__compass.magvar = 1
                    
                #Check for magvar with GPS
                if meas.transects[0].gps is None:
                    mag_var_check2 = np.abs(mag_var_check) > 0
                    idx = mag_var == 0
                    if len(idx) > 0:
                        self.compass.messages.append('COMPASS: MAgnetic variation is 0 and GPS data are present')
                        self.compass.status2 = 'warning'
                        self.compass.magvar = 2
                        self.compass.mag_var_idx = idx
                        
                #Check pitch and roll
                n_transects = len(meas.transects)
                cal_roll_idx = []
                cal_pitch_idx = []
                mag_error_idx = []
                for n in range(n_transects):
                    selected_pitch = meas.transects[n].sensors.pitch_deg.selected
                    p_select = getattr(meas.transects[n].sensors.pitch_deg, selected_pitch)
                    mean_pitch = np.nanmean(p_select._SensorData__data)
                    std_pitch = np.nanstd(p_select._SensorData__data)
                    selected_roll = meas.transects[n].sensors.roll_deg.selected
                    r_select = getattr(meas.transects[n].sensors.roll_deg, selected_roll)
                    mean_roll = np.nanmean(r_select._SensorData__data)
                    std_roll = np.nanstd(r_select._SensorData__data)
                    
                    call_roll_idx, call_pitch_idx = np.array([]), np.array([])
                    #If pitch and roll limit for compass calibration exist determine
                    #what transects have violations
                    if meas.transects[n].sensors.heading_deg.internal.roll_limit is not None:
                        idx_max = r_select._SensorData__data > meas.transects[n].sensors.heading_deg.internal.roll_limit[0]
                        idx_min = r_select._SensorData__data < meas.transects[n].sensors.heading_deg.internal.roll_limit[1]
                        
                        idx = np.array([])
                        
                        if np.sum(idx_max) > 0:
                            idx = np.append(idx, idx_max)
                            
                        if np.sum(idx_min) > 0:
                            idx = np.append(idx, idx_min)
                            
                        if np.sum(idx) > 0:
                            call_roll_idx = np.append(call_roll_idx, n)
                            
                        idx_max = p_select._SensorData__data > meas.transects[n].sensors.heading_deg.internal.pitch_limit[0]
                        idx_min = p_select._SensorData__data < meas.transects[n].sensors.heading_deg.internal.pitch_limit[1]
                        
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
                            
                    
                self.compass.mean_pitch_idx_w = np.where(np.abs(mean_pitch) > 8)
                self.compass.mean_pitch_idx_c = np.where((np.abs(mean_pitch) > 4) & (np.abs(mean_pitch) < 8))
                
                if len(self.compass.mean_pitch_idx_w) > 0:
                    self.compass.status2 = 'warning'
                    self.compass.messages.append('PITCH: One or more transects have a mean pitch > 8 deg')
                elif len(self.compass.mean_pitch_idx_c) > 0:
                    if self.compass.status2 == 'good':
                        self.compass.status2 = 'caution'
                        
                    self.compass.messages.append('Pitch: Pne or more transects have a mean pitch > 4 deg')
                    
                self.compass.mean_roll_idx_w = np.where(np.abs(mean_roll) > 8)
                self.compass.mean_roll_idx_c = np.where((np.abs(mean_roll) > 4) & (np.abs(mean_roll) < 8))
                
                if len(self.compass.mean_roll_idx_w) > 0:
                    self.compass.status2 = 'warning'
                    self.compass.messages.append('ROLL: One or more transects have a mean roll > 8 deg')
                elif len(self.compass.mean_roll_idx_c) > 0:
                    if self.compass.status2 == 'good':
                        self.compass.status2 = 'caution'
                        
                    self.compass.messages.append('Roll: One or more transects have a mean roll > 4 deg')
                    
                self.compass.std_pitch_idx = np.where(np.nanmean(std_pitch) > 5)
                if len(self.compass.std_pitch_idx) > 0:
                    if self.compass.status2 == 'good':
                        self.compass.status2 = 'caution'
                    self.compass.messages.append('Pitch: One or more transects have a pitch std_dev > 5 deg')
                    
                        
                self.compass.std_roll_idx = np.where(np.nanmean(std_roll) > 5)
                if len(self.compass.std_roll_idx) > 0:
                    if self.compass.status2 == 'good':
                        self.compass.status2 = 'caution'
                    self.compass.messages.append('Roll: One or more transects have a roll std_dev > 5 deg')       
                    
                    
                self.compass.cal_pitch_idx = cal_pitch_idx
                if len(cal_pitch_idx) > 0:
                    if self.compass.status2 == 'good':
                        self.compass.status2 = 'caution'
                        
                    self.compass.messages.append('Compass: One or more transects have pitch exceeding calibration limits')
                    
                self.compass.cal_roll_idx = cal_roll_idx
                if len(cal_roll_idx) > 0:
                    if self.__compass.status2 == 'good':
                        self.__compass.status2 = 'caution'
                        
                    self.compass.messages.append('Compass: One or more transects have roll exceeding calibration limits') 
                    
                self.compass.mag_error_idx = mag_error_idx
                if len(mag_error_idx) > 0:
                    if self.compass.status2 == 'good':
                        self.compass.status2 = 'caution'
                        
                    self.compass.messages.append('Compass: One or more transects have a change in magnetic field exceeding 2%') 
                    
                    
                if self.compass.status1 == 'warning' or self.compass.status2 == 'warning':
                    self.compass.status = 'warning'
                elif self.compass.status1 == 'caution' or self.compass.status2 == 'caution':
                    self.compass.status = 'caution'
                else:
                    self.compass.status = 'good'
                    
                    
    def temperature_QA(self, meas):
        '''Apply quality checks to compass
        Input:
        meas: object of Measurement
        '''
        
        #Run quality assurance checks
        self.temp_QA_check(meas)
        
        self.temperature = TestMessages('temperature')
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
        self.temperature.status = 'good'
        if range > 2:
            self.temperature.messages.append('TEMPERATURE: Temperature is {3.1f} degrees C which is greater than 2 degrees'.format(round(range, 1)))
            self.temperature.status = 'warning'
        elif range > 1:
            self.temperature.messgaes.append('Temperature: Temperature is {3.1f} degrees C which is greater than 1 degree'.format(round(range, 1)))
            self.temperature.status = 'caution'
        
        #Assign a quality code based on whether an independent temperature was entered
        if meas.ext_temp_chk.user is None:
            self.temperature.messages.append('Temperature: No independent temperature reading')
            self.temperature.status = 'caution'
        
        #Assign a quality code based on the difference between the mean ADCP and the independent temperature readings
        elif meas.ext_temp_chk.adcp is None:
            diff = np.abs(meas.ext_temp_chk.user - np.nanmean(temp))
            if diff >= 2:
                self.temperature.messages('TEMPERATURE: The difference between the average ADCP and independent temperature is ' \
                                            + '{} degrees C which is not less than 2 degrees'.format(round(diff, 1)))
                self.temperature.status = 'warning'
                
        #Assign a quality code based on the difference between the entered ADCP and the independent temperature readings
        else:
            diff = np.abs(meas.ext_temp_chk.user - meas.ext_temp_chk.adcp)
            if diff >= 2:
                self.temperature.messages('TEMPERATURE: The difference between the average ADCP and independent temperature is ' \
                                            + '{} degrees C which is not less than 2 degrees'.format(round(diff, 1)))
                self.temperature.status = 'warning'   
            
         
                                                   
        

 
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
            
            
        