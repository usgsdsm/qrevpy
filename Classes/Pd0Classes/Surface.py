'''
Created on Jun 19, 2017

@author: gpetrochenkov
'''

import numpy as np

class Surface(object):
    
    def __init__(self, n_ensembles, n_velocities, max_surface_bins):
        self.no_cells = np.zeros(n_ensembles)
        self.cell_size_cm = np.empty(n_ensembles)
        self.dist_bin1_cm = np.empty(n_ensembles)
        self.vel_mps = np.tile([np.nan], [n_velocities, max_surface_bins, n_ensembles])
        self.corr = np.empty([n_velocities, max_surface_bins, n_ensembles])
        self.pergd = np.empty([n_velocities, max_surface_bins, n_ensembles])
        self.rssi = np.empty([n_velocities, max_surface_bins, n_ensembles])