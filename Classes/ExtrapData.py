"""
Created on Sep 14, 2017

@author: gpetrochenkov
"""

class ExtrapData(object):
    
    def __init__(self):
        
        self.top_method_orig = None #Extrapolation method for top of profile: Power, Constant, 3-Point
        self.bot_method_orig = None #Extrapolation method for bottom of profile: Power, No Slip
        self.exponent_orig = None #Exponent for power of no slip methods
        self.top_method = None #Extrapolation method for top of profile: Power, Constant, 3-Point
        self.bot_method = None #Extrapolation method for bottom of profile: Power, No Slip
        self.exponent = None #Exponent for power of no slip methods
        
    def populate_data(self, top, bot, exp):
        
        self.top_method_orig = top
        self.bot_method_orig = bot
        self.top_method = top
        self.bot_method = bot
        self.exponent_orig = float(exp)
        self.exponent = float(exp)
        
    def set_extrap_data(self, top, bot, exp):
        self.top_method = top
        self.bot_method = bot
        self.exponent = exp
        
    def set_property(self, prop, setting):
        setattr(self, prop, setting)