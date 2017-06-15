from Classes.MMT_TRDI import MMT_TRDI
import numpy as np

class Measurement(object):
    """Class to hold measurement details for use in the GUI
        Analogous to matlab class: clsMeasurement
    """


    def __init__(self):
        self.station_name = None
        self.station_number = None
        self.transects = None
        self.mb_tests = None
        self.sys_test = None
        self.compass_cal = None
        self.compass_eval = None
        self.ext_temp_chk = None
        self.extrap_fit = None
        self.processing = None
        self.comments = None
        self.discharge = None
        self.uncertainty = None
        self.intitial_settings = None
        self.qa = None
        self.userRating = None

    def load_trdi(self, **kargs):
        '''method to load trdi MMT file'''

        #read in the MMT file
        self.trdi = MMT_TRDI(kargs[0])
    
        #Get properties if they exist, otherwise set them as blank strings
        self.station_name = self.trdi.site_info['Name']
        if np.isnan(self.station_name):
            self.station_name = ''

        self.station_number = self.trdi.site_info['Number']
        if np.isnan(self.station_number):
            self.station_number = ''

        self.processing = 'WR2'

#         if len(kargs) > 1:
#             self.transects = clsTransectData('TRDI', mmt,'Q',kargs[1])
#         else:
#             self.transects = clsTransectData('TRDI', mmt)








