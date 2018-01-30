"""
Created on Sep 14, 2017

@author: gpetrochenkov
"""
class SensorStructure(object):
    
    def __init__(self):
        
        self.selected = None  #The selected sensor reference name ('internal', 'external', 'user')
        self.internal = None #Contains the data from the internal sensor
        self.external = None #Contains the data from an external sensor
        self.user = None #Contains user supplied value
        
    def set_selected(self,selected_name):
        """Set the selected source for the specified object"""
        self.selected = selected_name
        