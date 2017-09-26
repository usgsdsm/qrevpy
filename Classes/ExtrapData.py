'''
Created on Sep 14, 2017

@author: gpetrochenkov
'''

class ExtrapData(object):
    
    def __init__(self):
        
        self.__top_method_orig = None #Extrapolation method for top of profile: Power, Constant, 3-Point
        self.__bot_method_orig = None #Extrapolation method for bottom of profile: Power, No Slip
        self.__exponent_orig = None #Exponent for power of no slip methods
        self.__top_method = None #Extrapolation method for top of profile: Power, Constant, 3-Point
        self.__bot_method = None #Extrapolation method for bottom of profile: Power, No Slip
        self.__exponent = None #Exponent for power of no slip methods
        
    def populate_data(self, top, bot, exp):
        
        self.__top_method_orig = top
        self.__bot_method_orig = bot
        self.__top_method = top
        self.__bot_method = bot
        self.__exponent_orig = float(exp)
        self.__exponent = float(exp)
        
    def set_extrap_data(self, top, bot, exp):
        self.__top_method = top
        self.__bot_method = bot
        self.__exponent = exp
        
    def set_property(self, prop, setting):
        setattr(self, prop, setting)