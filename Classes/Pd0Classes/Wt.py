"""
Created on Jun 16, 2017

@author: gpetrochenkov
"""
import numpy as np

class Wt(object):
    """
    classdocs
    """


    def __init__(self, n_bins,n_ensembles,n_velocities):
        self.corr = np.empty([n_velocities, n_bins, n_ensembles])
        self.pergd = np.empty([n_velocities, n_bins, n_ensembles])
        self.rssi = np.empty([n_velocities, n_bins, n_ensembles])
        self.vel_mps = np.empty([n_velocities, n_bins, n_ensembles])
