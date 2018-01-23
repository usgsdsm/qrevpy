'''
Created on Sep 26, 2017

@author: gpetrochenkov
'''
from Classes.Measurement import Measurement
import numpy as np
from Classes.Uncertainty import Uncertainty

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
        
        self.__transects = QATransect()
        
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
            
            
            
        
 
class QATransect(object):
    
    def __init__(self):
        self.status = 'good'
        self.messages = []
        self.sign = 0
        self.duration = 0
        self.number = 0
        self.uncertainty = 0       
        