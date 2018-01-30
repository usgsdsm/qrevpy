"""
Created on Sep 12, 2017

@author: gpetrochenkov
"""
class EdgeData(object):
    
    def __init__(self):
        
        self.edge_type = None #Shape of edge: 'Triangular', 'Rectangular', 'Custom, 'User Q'
        self.dist_m = None #Distance to shore
        self.cust_coeff  = None #Custom coefficient provided by user
        self.num_ens_2_avg = None #Number of ensembles to average for depth and velocities
        self.user_Q_cms = None
        
    def populate_data(self, edge_type, dist, kargs = None):
        """Construct left or right edge object from provided inputs
        
        Inputs:
        edge_type: type of edge (Triangular, Rectangular, Custom, UserQ)
        dist_m: distance to shore
        kargs:
        (for custom): 0 edge coefficient
                      1 number of edge ensembles
        (for UserQ):  0 discharge supplied by user
                      1 optional number of edge ensembles
        (for Triangular and Rectangular):
                      0 number of edge ensembles
        """
        
        self.edge_type = edge_type
        self.dist_m = dist
        self.num_ens_2_avg = 10
        self.user_Q_cms = []
        
        #Set properties for custom coefficient
        if edge_type == 'Custom':
            self.cust_coef = kargs[0]
            if len(kargs) > 1:
                self.num_ens_2_avg = kargs[1]
                
        elif edge_type == 'User Q':
            self.user_Q_cms = kargs[0]
            if len(kargs) > 1:
                self.num_ens_2_avg = kargs[1]
                
        else:
            if kargs is not None:
                self.num_ens_2_avg = kargs[0]
                
    def change_property(self, prop, setting):
        """Change edge data property"""
        
        setattr(self, prop, setting)