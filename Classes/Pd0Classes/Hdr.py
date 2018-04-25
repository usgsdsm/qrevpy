"""
Created on Jun 16, 2017

@author: gpetrochenkov
"""
import numpy as np

class Hdr(object):
    """Class to hold header variables.

    Attributes
    ----------
    bytes_per_ens: int
    data_offsets: int
    n_data_types: int
    data_ok: int
    invalid: str
        Leader ID that was not recognized
    """
    
    def __init__(self, n_ensembles, n_types):
        """Initialize instance variables to empty arrays.

        Parameters
        ----------
        n_ensembles: int
            Number of ensembles
        n_types
            Number of data types
        """
        self.bytes_per_ens = np.empty(n_ensembles)
        self.data_offsets = np.empty([n_ensembles,n_types])
        self.n_data_types = np.empty(n_ensembles)
        self.data_ok = np.empty(n_ensembles)
        self.invalid = ['']*n_ensembles
        