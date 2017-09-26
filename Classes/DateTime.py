'''
Created on Sep 6, 2017

@author: gpetrochenkov
'''

class DateTime(object):
    '''This stores the date and time data'''
    
    def __init__(self):
        self.date = None # measurement date mm/dd/yyyy
        self.start_serial_time = None # python serial time for start of transect
        self.end_serial_time = None # python serial time for end of transect
        self.transect_duration_sec = None # duration of transect in seconds
        self.ens_duration_sec = None #Duration of each ensemble in seconds
        
    def populate_data(self, date_in, start_in, end_in, ens_dur_in):
        '''Populate data in object'''
        
        self.date = date_in
        self.start_serial_time = start_in
        self.end_serial_time = end_in
        self.transect_duration_sec = end_in - start_in
        self.ens_duration_sec = ens_dur_in
        
    def time_2_serial_time(self, time_in, source):
        if source == 'TRDI':
            return 719529 + time_in / (60 * 60 * 24)
        else:
            return 719529 + 10957 + time_in / (60 * 60 * 24)