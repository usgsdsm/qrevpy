"""
Created on Sep 14, 2017

@author: gpetrochenkov
Modified DSM 2/1/2018
    - Added numpy docstrings
    - Cleaned up PEP8
"""


class SensorStructure(object):
    """Class to store sensor data from various sources.

    Attributes
    ----------
    self.selected: str
        The selected sensor reference name ('internal', 'external', 'user').
    self.internal: SensorData
        Contains the data from the internal sensor, object of SensorData
    self.external: SensorData
        Contains the data from an external sensor, object of SensorData
    self.user: SensorData
        Contains user supplied value, object of SensorData
    """
    
    def __init__(self):
        """Initialize class and set variable to None."""

        self.selected = None  # The selected sensor reference name ('internal', 'external', 'user')
        self.internal = None  # Contains the data from the internal sensor
        self.external = None  # Contains the data from an external sensor
        self.user = None  # Contains user supplied value
        
    def set_selected(self, selected_name):
        """Set the selected source for the specified object

        Parameters
        ----------
        selected_name: str
            Type of data (internal, external, user).
        """
        self.selected = selected_name
