"""
Created on Sep 6, 2017

@author: gpetrochenkov

Modified DSM 2/2/2018
    - Added numpy docstrings
    - Changed time_2_serial_time to properly handle SonTek to Python serial time
    - Cleaned up PEP8
"""


class DateTime(object):
    """This stores the date and time data in Python compatible format.

    Attributes
    ----------
    date: str
        Measurement date as mm/dd/yyyy
    start_serial_time: float
        Python serial time for start of transect (seconds since 1/1/1970), timestamp
    end_serial_time: float
        Python serial time for end of transect (seconds since 1/1/1970), timestamp
    transect_duration_sec: float
        Duration of transect, in seconds.
    ens_duration_sec: np.array(float)
        Duration of each ensemble, in seconds.
    """
    
    def __init__(self):
        """Initialize class and instance variables."""

        self.date = None  # Measurement date mm/dd/yyyy
        self.start_serial_time = None  # Python serial time for start of transect, timestamp
        self.end_serial_time = None  # Python serial time for end of transect, timestamp
        self.transect_duration_sec = None  # Duration of transect in seconds
        self.ens_duration_sec = None  # Duration of each ensemble in seconds
        
    def populate_data(self, date_in, start_in, end_in, ens_dur_in):
        """Populate data in object.

        Parameters
        ----------
        date_in: str
            Measurement date as mm/dd/yyyy
        start_in: float
            Python serial time for start of transect.
        end_in: float
            Python serial time for end of transect.
        ens_dur_in: np.array(float)
            Duration of each ensemble, in seconds.
        """
        
        self.date = date_in
        self.start_serial_time = start_in
        self.end_serial_time = end_in
        self.transect_duration_sec = end_in - start_in
        self.ens_duration_sec = ens_dur_in
