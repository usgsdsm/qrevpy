'''
Created on Jun 16, 2017

@author: gpetrochenkov
'''
import numpy as np

class Hdr(object):
    
    def __init__(self, n_ensembles, n_types):
        
        self.bytes_per_ens =  np.empty(n_ensembles)
        self.data_offsets = np.empty([n_ensembles,n_types])
        self.n_data_types = np.empty(n_ensembles)
        self.data_ok = np.empty(n_ensembles)
        self.invalid = ['']*n_ensembles
        