"""
Created on Jun 19, 2017

@author: gpetrochenkov
"""
import numpy as np

class AutoMode(object):
    """
    classdocs
    """


    def __init__(self, n_ensembles):
        """
        Constructor
        """
        self.beam_count = np.empty(n_ensembles)
        self.Beam1 = Beam(n_ensembles)
        self.Beam2 = Beam(n_ensembles)
        self.Beam3 = Beam(n_ensembles)
        self.Beam4 = Beam(n_ensembles)
        self.Reserved = np.empty(n_ensembles)
        
        

class Beam(object):
    
    def __init__(self, n_ensembles):
        self.mode = np.empty(n_ensembles)
        self.depth_cm = np.empty(n_ensembles)
        self.ping_count = np.empty(n_ensembles)
        self.ping_type = np.empty(n_ensembles)
        self.cell_count = np.empty(n_ensembles)
        self.cell_size_cm = np.empty(n_ensembles)
        self.cell_mid_cm = np.empty(n_ensembles)
        self.code_repeat = np.empty(n_ensembles)
        self.trans_length_cm = np.empty(n_ensembles)
        self.lag_length_cm = np.empty(n_ensembles)
        self.transmit_bw = np.empty(n_ensembles)
        self.receive_bw = np.empty(n_ensembles)
        self.ping_interval_ms = np.empty(n_ensembles)