"""
Created on Sep 14, 2017

@author: gpetrochenkov
Modified DSM 2/1/2018
    - Added numpy docstrings
    - Commented out two methods
    - Cleaned up PEP8
    - Set method to static
"""

# from Classes.HeadingData import HeadingData
from Classes.SensorStructure import SensorStructure
# from Classes.SensorData import SensorData


class Sensors(object):
    """Class to store data from ADCP sensors.

    Attributes
    ----------
    heading_deg: object
        Object of HeadingData.
    pitch_deg: object
        Pitch data, object of SensorStructure
    roll_deg: object
        Roll data, object of SensorStructure
    temperature_deg_c: object
        Temperature data, object of SensorStructure
    salinity_ppt: object
        Salinity data, object of SensorStructure
    speed_of_sound_mps: object
        Speed of sound, object of SensorStructure
    """

    def __init__(self):
        """Initialize class and create variable objects"""

        self.heading_deg = SensorStructure()  # Object of HeadingData
        self.pitch_deg = SensorStructure()  # Pitch data, object of SensorStructure
        self.roll_deg = SensorStructure()  # Roll data, object of SensorStructure
        self.temperature_deg_c = SensorStructure()  # Temperature data, object of SensorStructure
        self.salinity_ppt = SensorStructure()  # Salinity data, object of SensorStructure
        self.speed_of_sound_mps = SensorStructure()  # Speed of sound, object of SensorStructure
        
    # DSM I don't think these are needed. 2/1/2018
    # def add_sensor_data(self, sensor_name, sensor_type, data,source, kargs = None):
    #     #this function will create the appropriate objects for the specified property
    #     if sensor_name == 'heading_deg':
    #         sensor = getattr(self, sensor_name)
    #         h_data = HeadingData()
    #         h_data.populate_data(data, source, kargs)
    #         setattr(sensor, sensor_type, h_data)
    #     else:
    #         sensor = getattr(self, sensor_name)
    #         s_data = SensorData()
    #         s_data.populate_data(data, source)
    #         setattr(sensor, sensor_type, s_data)
    #
    # def set_selected(self, sensor_name, selected_name):
    #     sensor = getattr(self, sensor_name)
    #     sensor.set_selected(selected_name)

    @staticmethod
    def speed_of_sound(temperature, salinity):
        """Computes speed of sound from temperature and salinity.

        Parameters
        ----------
        temperature: float or np.array(float)
            Water temperature at transducer face, in degrees C.
        salinity: float or np.array(float)
            Water salinity at transducer face, in ppt.
        """

        # Not provided in RS Matlab file computed from equation used in TRDI BBSS
        sos = 1449 + 4.6 * temperature - 0.055 * temperature**2 + 0.00029 * temperature**3 \
            + (1.34 - 0.01 * temperature) * (salinity - 35)
        return sos
