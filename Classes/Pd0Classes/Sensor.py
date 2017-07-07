'''
Created on Jun 16, 2017

@author: gpetrochenkov
'''
import numpy as np
from numpy.matlib import repmat

class Sensor(object):
    
    def __init__(self, n_ensembles):
        
        self.ambient_temp = np.empty(n_ensembles)
        self.attitude_temp = np.empty(n_ensembles)
        self.attitude = np.empty(n_ensembles)
        self.bit_test = np.empty(n_ensembles)
        self.contam_sensor = np.empty(n_ensembles)
        self.date = np.empty([n_ensembles,3])
        self.date_y2k = np.empty([n_ensembles,4])
        self.date_not_y2k = np.empty([n_ensembles,3])
        self.error_status_word = repmat(['']*4,n_ensembles,1)
        self.heading_deg = np.empty(n_ensembles)
        self.heading_std_dev_deg =  np.empty(n_ensembles)
        self.mpt_msc = np.empty([n_ensembles,3])
        self.num = np.empty(n_ensembles)
        self.num_fact = np.empty(n_ensembles)
        self.num_tot = np.empty(n_ensembles)
        self.orient = ['']*n_ensembles
        self.pitch_std_dev_deg =  np.empty(n_ensembles)
        self.pitch_deg =  np.empty(n_ensembles)
        self.pressure_neg =  np.empty(n_ensembles)
        self.pressure_pos =  np.empty(n_ensembles)
        self.pressure_pascal = np.empty(n_ensembles)
        self.pressure_var_pascal = np.empty(n_ensembles)
        self.roll_std_dev_deg =  np.empty(n_ensembles)
        self.roll_deg = np.empty(n_ensembles)
        self.salinity_ppt =  np.empty(n_ensembles)
        self.sos_mps =  np.empty(n_ensembles)
        self.temperature_degc =  np.empty(n_ensembles)
        self.time =  np.empty([n_ensembles,4])
        self.time_y2k = np.empty([n_ensembles,4])
        self.xdcr_depth_dm = np.empty(n_ensembles)
        self.xmit_current = np.empty(n_ensembles)
        self.xmit_voltage= np.empty(n_ensembles)
        self.vert_beam_eval_amp = np.empty(n_ensembles)
        self.vert_beam_RSSI_amp = np.empty(n_ensembles)
        self.vert_beam_range_m = np.empty(n_ensembles)
        self.vert_beam_gain = ['']*n_ensembles
        self.vert_beam_status = np.empty(n_ensembles)
        
        