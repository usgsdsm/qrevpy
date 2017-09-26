'''
Created on Jun 16, 2017

@author: gpetrochenkov
'''
import numpy as np
from numpy.matlib import repmat
class Gps2(object):
    '''
    classdocs
    '''


    def __init__(self, n_ensembles, wr2):
        '''
        Constructor
        '''
        self.gga_delta_time =  np.empty([n_ensembles, 20])
        self.gga_header = [['']*20 for k in range(n_ensembles)]
        self.gga_sentence = [['']*20 for k in range(n_ensembles)]
        self.utc =   np.empty([n_ensembles, 20])
        self.lat_deg =   np.zeros([n_ensembles, 20])
        self.lat_ref =  [['']*20 for k in range(n_ensembles)]
        self.lon_deg =   np.zeros([n_ensembles, 20])
        self.lon_ref =  [['']*20 for k in range(n_ensembles)]
        self.corr_qual =   np.empty([n_ensembles, 20])
        self.num_sats =   np.empty([n_ensembles, 20])
        self.hdop =   np.empty([n_ensembles, 20])
        self.alt =   np.empty([n_ensembles, 20])
        self.alt_unit =  [['']*20 for k in range(n_ensembles)]
        self.geoid =   np.empty([n_ensembles, 20])
        self.geoid_unit = [['']*20 for k in range(n_ensembles)]
        self.d_gps_age =   np.empty([n_ensembles, 20])
        self.ref_stat_id =   np.empty([n_ensembles, 20])
        self.vtg_delta_time =   np.empty([n_ensembles, 20])
        self.vtg_header =  [['']*20 for k in range(n_ensembles)]
        self.vtg_sentence = [['']*20 for k in range(n_ensembles)]
        self.course_true =   np.empty([n_ensembles, 20])
        self.true_indicator = [['']*20 for k in range(n_ensembles)]
        self.course_mag =   np.empty([n_ensembles, 20])
        self.mag_indicator =  [['']*20 for k in range(n_ensembles)]
        self.speed_knots =  np.empty([n_ensembles, 20])
        self.knots_indicator =  [['']*20 for k in range(n_ensembles)]
        self.speed_k_mph =   np.zeros([n_ensembles, 20])
        self.kmph_indicator = [['']*20 for k in range(n_ensembles)]
        self.mode_indicator =  [['']*20 for k in range(n_ensembles)]
        self.dbt_delta_time =   np.empty([n_ensembles, 20])
        self.dbt_header = [['']*20 for k in range(n_ensembles)]
        self.depth_ft = np.empty([n_ensembles, 20])
        self.ft_indicator = [['']*20 for k in range(n_ensembles)]
        self.depth_m = np.zeros([n_ensembles,20])
        self.m_indicator = [['']*20 for k in range(n_ensembles)]
        self.depth_fath = np.empty([n_ensembles, 20])
        self.fath_indicator = [['']*20 for k in range(n_ensembles)]
        self.hdt_delta_time = np.empty([n_ensembles, 20])
        self.hdt_header = [['']*20 for k in range(n_ensembles)]
        self.heading_deg = np.empty([n_ensembles, 20])
        self.h_true_indicator =  [['']*20 for k in range(n_ensembles)]
        
        if wr2 == 1:
            self.gga_velE_mps = np.empty(n_ensembles)
            self.gga_velN_mps = np.empty(n_ensembles)
            self.vtg_velE_mps = np.empty(n_ensembles)
            self.vtg_velN_mps = np.empty(n_ensembles)
        