"""
Created on Sep 12, 2017

@author: gpetrochenkov
Modified DSM 1/31/2018
    - Added numpy docstrings
    - Removed need for kargs
    - Cleaned up PEP8
"""


class EdgeData(object):
    """Class used to store edge settings.

    Attributes
    ----------
    type: str
        Shape of edge: 'Triangular', 'Rectangular', 'Custom, 'User Q'
    distance_m: float
        Distance to shore, in m.
    cust_coeff: float
        Custom coefficient provided by user.
    number_ensembles: int
        Number of ensembles to average for depth and velocities.
    user_discharge_cms: float
        User supplied discharge for edge, in cms.
    """
    
    def __init__(self):
        
<<<<<<< HEAD
        self.type = None       # Shape of edge: 'Triangular', 'Rectangular', 'Custom, 'User Q'
        self.distance_m = None          # Distance to shore
        self.cust_coeff = None     # Custom coefficient provided by user
        self.number_ensembles = None   # Number of ensembles to average for depth and velocities
        self.user_discharge_cms = None      # User supplied edge discharge.
        
    def populate_data(self, edge_type, distance=None, number_ensembles=10, coefficient=None, user_discharge=None):
        """Construct left or right edge object from provided inputs
        
        Parameters
        ----------
        edge_type: str
            Type of edge (Triangular, Rectangular, Custom, UserQ)
        distance: float
            Distance to shore, in m.
        number_ensembles: int
            Number of edge ensembles for all types but UserQ
        coefficient: float
            User supplied custom edge coefficient.
        user_discharge: float
            User supplied edge discharge, in cms.
        """

        # Set properties for custom coefficient
        self.type = edge_type
        self.distance_m = distance
        self.number_ensembles = number_ensembles
        self.user_discharge_cms = user_discharge
        self.cust_coeff = coefficient

=======
        
        self.type = None #Shape of edge: 'Triangular', 'Rectangular', 'Custom, 'User Q'
        self.dist_m = None #Distance to shore
        self.cust_coeff  = None #Custom coefficient provided by user
        self.num_ens_2_avg = None #Number of ensembles to average for depth and velocities
        self.user_Q_cms = None
        
    def populate_data(self, edge_type, dist, kargs = None):
        '''Construct left or right edge object from provided inputs
        
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
        '''
        
        self.edge_type = edge_type
        self.dist_m = dist
        self.num_ens_2_avg = 10
        self.user_Q_cms = []
        
        #Set properties for custom coefficient
        if edge_type == 'Custom':
            self.cust_coef = kargs[0]
            if len(kargs) > 1:
                self.__num_ens_2_avg = kargs[1]
                
        elif edge_type == 'User Q':
            self.user_Q_cms = kargs[0]
            if len(kargs) > 1:
                self.num_ens_2_avg = kargs[1]
                
        else:
            if kargs is not None:
                self.num_ens_2_avg = kargs[0]
                
>>>>>>> 6ca6c50c231afa610ed3a693864074d7104a5f20
    def change_property(self, prop, setting):
        """Change edge data property

        Parameters
        ----------
        prop: str
            Property to change.
        setting:
            New setting for property.
        """
        setattr(self, prop, setting)
