'''
Created on Jun 16, 2017

@author: gpetrochenkov
'''
import numpy as np

class Inst(object):
    
    def __init__(self, n_ensembles):
        
        self.beam_ang = np.empty(n_ensembles)
        self.beams = np.empty(n_ensembles)
        self.data_type = np.empty([n_ensembles,4], dtype=np.str)
        self.firm_ver = np.empty(n_ensembles)
        self.freq = np.empty(n_ensembles)
        self.pat = np.empty(n_ensembles, dtype=np.str)
        self.res_RDI = 0
        self.sensor_CFG = np.empty(n_ensembles)
        self.xducer = np.empty(n_ensembles, dtype=np.str)
        self.t_matrix = np.empty([4,4])
        self.demod = np.empty(n_ensembles)
    