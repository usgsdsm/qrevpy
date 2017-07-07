'''
Created on Jun 16, 2017

@author: gpetrochenkov
'''
import numpy as np

class Gps(object):
    '''
    classdocs
    '''


    def __init__(self, n_ensembles):
        '''
        Constructor
        '''
        self.alt_m = np.empty(n_ensembles)
        self.gga_diff = np.empty(n_ensembles)
        self.gga_hdop = np.empty(n_ensembles)
        self.gga_n_stats = np.empty(n_ensembles)
        self.gga_vel_e_mps = np.empty(n_ensembles)
        self.gga_vel_n_mps = np.empty(n_ensembles)
        self.gsa_p_dop = np.empty(n_ensembles)
        self.gsa_sat = np.empty([n_ensembles,6])
        self.gsa_v_dop = np.empty(n_ensembles)
        self.lat_deg = np.empty(n_ensembles)
        self.long_deg = np.empty(n_ensembles)
        self.vtg_vel_e_mps = np.empty(n_ensembles)
        self.vtg_vel_n_mps = np.empty(n_ensembles)