'''
Created on Sep 12, 2017

@author: gpetrochenkov
'''
from Classes.EdgeData import EdgeData

class Edges(object):
    
    def __init__(self):
        
        self.__rec_edge_method = None
        self.__vel_method = None
        self.__left = EdgeData()
        self.__right = EdgeData()
        
    def populate_data(self, rec_edge_method, vel_method):
        
        self.__rec_edge_method = rec_edge_method
        self.__vel_method = vel_method
        
        
    def change_property(self, prop, setting, kargs=None):
        '''Change edge property
        
        Input:
        prop: name of property
        setting: property setting
        kargs: edge to change (left, right)
        '''
        
        if kargs is None:
            setattr(self, prop, setting)
        else:
            temp = getattr(self, kargs[0])
            temp.change_property(prop, setting)
            
    def create_edge(self, edge_loc, edge_type, dist, kargs=None):
        '''Create Edge property which is an object of EdgeData for each edge
        
        Input:
        edge_loc: left or right
        edge_type: type of edge (Triangular, Rectangular, Custom, User Q)
        dist: distance to shore
        kargs: value of cefficient for custom or discharge for User Q
        '''
        
        temp = getattr(self, '_Edges__'+edge_loc)
        temp.populate_data(edge_type, dist, kargs)
    
        
      
        
    