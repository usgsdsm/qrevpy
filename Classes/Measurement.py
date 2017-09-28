from Classes.MMT_TRDI import MMT_TRDI
import numpy as np
import os
from Classes.TransectData import TransectData, allocate_transects
from Classes.Pd0TRDI import Pd0TRDI
from Classes.PreMeasurement import PreMeasurement
from Classes.MovingBedTests import MovingBedTests

class Measurement(object):
    """Class to hold measurement details for use in the GUI
        Analogous to matlab class: clsMeasurement
    """

    def __init__(self, in_file, source='TRDI'):
        self.station_name = None
        self.station_number = None
        self.transects = []
        self.mb_tests = None
        self.sys_test = []
        self.compass_cal = []
        self.compass_eval = []
        self.ext_temp_chk = None
        self.extrap_fit = None
        self.processing = None
        self.comments = None
        self.discharge = None
        self.uncertainty = None
        self.intitial_settings = None
        self.qa = None
        self.userRating = None
        
        if source == 'TRDI':
            self.load_trdi(in_file, **{'type': 'Q', 'checked': False})

    def load_trdi(self, mmt, **kargs):
        '''method to load trdi MMT file'''

        #read in the MMT file
        self.trdi = MMT_TRDI(mmt)
        mmt = self.trdi
    
        #Get properties if they exist, otherwise set them as blank strings
        self.station_name = str(self.trdi.site_info['Name'])
        
        self.station_number = str(self.trdi.site_info['Number'])
        
        self.processing = 'WR2'

#         if len(kargs) > 1:
#             self.transects = TransectData('TRDI', mmt,'Q',kargs[1])
#         else:

        self.transects = allocate_transects('TRDI', mmt, kargs=['Q', kargs['checked']])
            
        #-------------------------Done Refactor----------------------------------------------
        
        if isinstance(mmt.qaqc, dict) or isinstance(mmt.mbt_transects):
            self.qaqc_TRDI(mmt)
            
            
    def qaqc_TRDI(self, mmt):
        '''Processes qaqc test, calibrations, and evaluations
        
        Input:
        mmt: object of MMT_TRDI
        
        '''
        #ADCP Test
        if 'RG_Test_TimeStamp' in mmt.qaqc:
            for n in range(len(mmt.qaqc['RG_Test'])):
                p_m = PreMeasurement()
                p_m.populate_data(mmt.qaqc['RG_Test_TimeStamp'][n], mmt.qaqc['RG_Test'][n],'TST')
                self.sys_test.append(p_m)
        else:
            self.sys_test.append(PreMeasurement())
            
        #Compass calibration
        if 'Compass_Cal_Timestamp' in mmt.qaqc:
            for n in range(len(mmt.qaqc['Compass_Cal_Test'])):
                p_m =  PreMeasurement()
                p_m.populate_data(mmt.qaqc['Compass_Cal_Timestamp'], mmt.qaqc['Compass_Cal_Test'], 'TCC')
                self.compass_cal.append(p_m)
        else:
            self.compass_cal.append(PreMeasurement())
            
        #Compass evaluation
        if 'Compass_Eval_Timestamp' in mmt.qaqc:
            for n in range(len(mmt.qaqc['Compass_Eval_Test'])):
                p_m =  PreMeasurement()
                p_m.populate_data(mmt.qaqc['Compass_Eval_Timestamp'], mmt.qaqc['Compass_Eval_Test'], 'TCC')
                self.compass_eval.append(p_m)
        else:
            self.compass_eval.append(PreMeasurement())
            
       
        
        if len(mmt.mbt_transects) > 0:
            
            #Create transect object/s
            transects= allocate_transects('TRDI',mmt,kargs=['MB'])
            if len(transects) > 0:
                for n in range(len(transects)):
                    mb_test = MovingBedTests()
                    mb_test.populate_data('TRDI', kargs=[mmt, mmt.mbt_transects[n]])
                    
                    #Save notes from mmt files in comments
                    if 'NoteDate' in mmt.mbt_transects:
                        r, c = mmt.mbt_transects['NoteDate'].shape
                        for n in range(r):
                            for k in range(c):
                                if mmt.mbt_transects['NoteDate'][n,k] is not None:
                                    pass
                                
        
            
        
            
        
            
    
    
            







