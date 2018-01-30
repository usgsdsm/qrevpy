"""
Created on Jun 16, 2017

@author: gpetrochenkov
"""
import numpy as np

class Bt(object):
    """
    classdocs
    """


    def __init__(self, n_ensembles, n_velocities):
        """
        Constructor
        """
        self.corr = np.empty([n_velocities,n_ensembles])
        self.depth_m = np.empty([n_velocities,n_ensembles])
        self.eval_amp = np.empty([n_velocities,n_ensembles])
        self.ext_depth_cm = np.empty(n_ensembles)
        self.pergd = np.empty([n_velocities,n_ensembles])
        self.rssi = np.empty([n_velocities,n_ensembles])
        self.vel_mps = np.empty([n_velocities,n_ensembles])