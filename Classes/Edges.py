from Classes.EdgeData import EdgeData


class Edges(object):
    """Class to store and process edge data.

    Attributes
    ----------
    rec_edge_method: str
        Method used to determine coef for rec. edge 'Fixed', 'Variable'.
    vel_method: str
        Method used to compute the velocity used 'MeasMag', 'VectorProf'.
    left: EdgeData
        Object of EdgeData for left edge.
    right: EdgeData
        Object of EdgeData for right edge.
    """
    
    def __init__(self):

        self.rec_edge_method = None
        self.vel_method = None
        self.left = EdgeData()
        self.right = EdgeData()
        
    def populate_data(self, rec_edge_method, vel_method):
        """Store the general methods used for edge data.

        Parameters
        ----------
        rec_edge_method: str
            Method used to determine coef for rec. edge 'Fixed', 'Variable'.
        vel_method: str
            Method used to compute the velocity used 'MeasMag', 'VectorProf'.
        """
        self.rec_edge_method = rec_edge_method
        self.vel_method = vel_method

    def change_property(self, prop, setting, edge=None):
        """Change edge property
        
        Parameters
        ----------
        prop: str
            Name of property.
        setting:
            New property setting.
        edge: str
            Edge to change (left, right)
        """
        
        if edge is None:
            setattr(self, prop, setting)
        else:
            temp = getattr(self, edge)
            temp.change_property(prop, setting)

    # DSM 1/31/2018 this method is not needed.
    # def create_edge(self, edge_loc, edge_type, dist, kargs=None):
    #     """Create Edge property which is an object of EdgeData for each edge
    #
    #     Input:
    #     edge_loc: left or right
    #     edge_type: type of edge (Triangular, Rectangular, Custom, User Q)
    #     dist: distance to shore
    #     kargs: value of cefficient for custom or discharge for User Q
    #     """
    #
    #     temp = getattr(self, '_Edges__'+edge_loc)
    #     temp.populate_data(edge_type, dist, kargs)
