'''
Created on Sep 14, 2017

@author: gpetrochenkov
'''

from Classes.HeadingData import HeadingData
from Classes.SensorStructure import SensorStructure
from Classes.SensorData import SensorData

class Sensors(object):
    
    def __init__(self):
        
        self.heading_deg = SensorStructure() #Object of HeadingData
        self.pitch_deg = SensorStructure() #Pitch data, object of SensorStructure
        self.roll_deg = SensorStructure()  #Roll data, object of SensorStructure
        self.temperature_deg_c = SensorStructure()  #Temperature data, object of SensorStructure
        self.salinity_ppt = SensorStructure()  #Salinity data, object of SensorStructure
        self.speed_of_sound_mps = SensorStructure()  #Speed of sound, object of SensorStructure
        
    def add_sensor_data(self, sensor_name, sensor_type, data,source, kargs = None):
        #this function will create the appropriate objects for the specified property
        if sensor_name == 'heading_deg':
            sensor = getattr(self, sensor_name)
            h_data = HeadingData()
            h_data.populate_data(data, source, kargs)
            setattr(sensor, sensor_type, h_data)
        else:
            sensor = getattr(self, sensor_name)
            s_data = SensorData()
            s_data.populate_data(data, source)
            setattr(sensor, sensor_type, s_data)
            
    def set_selected(self, sensor_name, selected_name):
        sensor = getattr(self, sensor_name)
        sensor.set_selected(selected_name)
        
        
    def speed_of_sound(self, temperature, salinity):
        #Speed of sound
        #Not provided in RS Matlab file computed from equation used in TRDI BBSS
        return 1449 + 4.6 * temperature - 0.055 * temperature**2 + 0.00029 * temperature**3 \
        + (1.34 - 0.01 * temperature) * (salinity - 35)
        
    