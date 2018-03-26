"""
Created on Jun 16, 2017

@author: gpetrochenkov
"""
import numpy as np
from numpy.matlib import repmat
class Gps2(object):
    """
    classdocs
    """


    def __init__(self, n_ensembles, wr2):
        """
        Constructor
        """
        self.gga_delta_time =  np.full([n_ensembles, 20], np.nan)
        self.gga_header = [['']*20 for k in range(n_ensembles)]
        self.gga_sentence = [['']*20 for k in range(n_ensembles)]
        self.utc =   np.full([n_ensembles, 20], np.nan)
        self.lat_deg =   np.zeros([n_ensembles, 20])
        self.lat_ref =  [['']*20 for k in range(n_ensembles)]
        self.lon_deg =   np.zeros([n_ensembles, 20])
        self.lon_ref =  [['']*20 for k in range(n_ensembles)]
        self.corr_qual =   np.full([n_ensembles, 20], np.nan)
        self.num_sats =   np.full([n_ensembles, 20], np.nan)
        self.hdop =   np.full([n_ensembles, 20], np.nan)
        self.alt =   np.full([n_ensembles, 20], np.nan)
        self.alt_unit =  [['']*20 for k in range(n_ensembles)]
        self.geoid =   np.full([n_ensembles, 20], np.nan)
        self.geoid_unit = [['']*20 for k in range(n_ensembles)]
        self.d_gps_age =   np.full([n_ensembles, 20], np.nan)
        self.ref_stat_id =   np.full([n_ensembles, 20], np.nan)
        self.vtg_delta_time =   np.full([n_ensembles, 20], np.nan)
        self.vtg_header =  [['']*20 for k in range(n_ensembles)]
        self.vtg_sentence = [['']*20 for k in range(n_ensembles)]
        self.course_true =   np.full([n_ensembles, 20], np.nan)
        self.true_indicator = [['']*20 for k in range(n_ensembles)]
        self.course_mag =   np.full([n_ensembles, 20], np.nan)
        self.mag_indicator =  [['']*20 for k in range(n_ensembles)]
        self.speed_knots =  np.full([n_ensembles, 20], np.nan)
        self.knots_indicator =  [['']*20 for k in range(n_ensembles)]
        self.speed_k_mph =   np.zeros([n_ensembles, 20])
        self.kmph_indicator = [['']*20 for k in range(n_ensembles)]
        self.mode_indicator =  [['']*20 for k in range(n_ensembles)]
        self.dbt_delta_time =   np.full([n_ensembles, 20], np.nan)
        self.dbt_header = [['']*20 for k in range(n_ensembles)]
        self.depth_ft = np.full([n_ensembles, 20], np.nan)
        self.ft_indicator = [['']*20 for k in range(n_ensembles)]
        self.depth_m = np.zeros([n_ensembles,20])
        self.m_indicator = [['']*20 for k in range(n_ensembles)]
        self.depth_fath = np.full([n_ensembles, 20], np.nan)
        self.fath_indicator = [['']*20 for k in range(n_ensembles)]
        self.hdt_delta_time = np.full([n_ensembles, 20], np.nan)
        self.hdt_header = [['']*20 for k in range(n_ensembles)]
        self.heading_deg = np.full([n_ensembles, 20], np.nan)
        self.h_true_indicator =  [['']*20 for k in range(n_ensembles)]
        
        if wr2 == 1:
            self.gga_velE_mps = np.empty(n_ensembles)
            self.gga_velN_mps = np.empty(n_ensembles)
            self.vtg_velE_mps = np.empty(n_ensembles)
            self.vtg_velN_mps = np.empty(n_ensembles)
        