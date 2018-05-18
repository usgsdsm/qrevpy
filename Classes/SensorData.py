class SensorData(object):
    """Class stores data for pitch, roll, temperature, salinity, and speed of sound and its source/

    Attributes
    ----------
    data: np.array(float)
        Data to be used in computations.
    data_orig: np.array(float)
        Original data loaded from raw data file.
    source: str
        Source of data, examples Int. Sensor, Ext. Sensor, User
    """
    
    def __init__(self):
        """Initializes class and variables."""

        self.data = None
        self.data_orig = None
        self.source = None
        
    def populate_data(self, data_in, source_in):
        """Store data in class.

        Parameters
        ----------
        data_in: np.array(float)
            Data to be stored.
        source_in: str
            Source of data to be stored.
        """

        self.data = data_in
        self.data_orig = data_in
        self.source = source_in
        
    def change_data(self, data_in):
        """Change data to be applied in computations.

        Parameters
        ----------
        data_in: np.array(float)
        """
        self.data = data_in
        
    def set_source(self, source_in):
        """Change source of data.

        Parameters
        ----------
        source_in: str
            Source of data.
        """
        self.source = source_in
