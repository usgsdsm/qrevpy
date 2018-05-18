from Classes.SensorStructure import SensorStructure
# from Classes.SensorData import SensorData


class Sensors(object):
    """Class to store data from ADCP sensors.

    Attributes
    ----------
    heading_deg: HeadingData
        Object of HeadingData.
    pitch_deg: SensorStructure
        Pitch data, object of SensorStructure
    roll_deg: SensorStructure
        Roll data, object of SensorStructure
    temperature_deg_c: SensorStructure
        Temperature data, object of SensorStructure
    salinity_ppt: SensorStructure
        Salinity data, object of SensorStructure
    speed_of_sound_mps: SensorStructure
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
