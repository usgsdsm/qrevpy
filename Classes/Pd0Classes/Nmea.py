'''
Created on Jun 19, 2017

@author: gpetrochenkov
'''
import numpy as np

class Nmea(object):
    '''
    classdocs
    '''


    def __init__(self, n_ensembles):
        '''
        Constructor
        '''
        self.gga = ['']*n_ensembles
        self.gsa = ['']*n_ensembles
        self.vtg = ['']*n_ensembles
        self.raw = ['']*n_ensembles
        self.dbt = ['']*n_ensembles