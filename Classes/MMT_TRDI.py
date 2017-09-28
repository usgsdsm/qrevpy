import xmltodict
from Classes.Transect import Transect

class MMT_TRDI(object):
    """Class to read in an MMT file"""

    def __init__(self, infile):
        '''Properties are initialized in constructor
            Analogous to matlab class: clsMmtTRDI
        '''
        self.infile = infile 
        self.project = {}
        self.site_info = {}
        self.transects = []
        self.field_config = None
        self.active_config = None
        self.summary_None = None
        self.summary_BT = None
        self.summary_GGA = None
        self.summary_VTG = None
        self.qaqc = {}
        self.mbt_transects = []
        self.mbt_field_config = None
        self.mbt_active_config = None

        self.read_mmt();

    def read_mmt(self):
        '''Method to read an MMT file and assign properties'''
        
        #open up the file and convert to an ordered dictionary tree
        with open(self.infile) as fd:
            self.dict = xmltodict.parse(fd.read())
        win_river = self.dict['WinRiver']

        self.project['Name'] = win_river['Project']['@Name']
        self.project['Version'] = win_river['Project']['@Version']

        if 'Locked' in win_river['Project'].keys():
            self.project['Locked'] = win_river['Project']['Locked']
        else:
            self.project['Locked'] = None

        
        siteinfo_keys = win_river['Project']['Site_Information'].keys()

        #Iterate through all of the keys and values of site info
        for x in siteinfo_keys:
            site_data = win_river['Project']['Site_Information'][x]
            if site_data is not None:

                #Remove @ symbol from properties
                if '@' in x:
                    x = x[1:]

                if x == 'Water_Temperature':
                    self.site_info[x] = float(site_data)

                    if self.site_info[x] < -100:
                        self.site_info[x] = [] #? should be None

                else:
                    self.site_info[x] = site_data
            else:
                self.site_info[x] = None

        
        trans = win_river['Project']['Site_Discharge']['Transect']

        #Create a Transect class for each transect found under Site_Discharge
        for i in range(len(trans)):
            
            self.transects.append(Transect(trans[i]))

        #ifDischarge Summary exists create summary fields
        if 'Discharge_Summary' in win_river['Project']['Site_Discharge'].keys():
            discharge_summary = win_river['Project']['Site_Discharge']['Discharge_Summary']

            self.summary_None = self.mmtqsum(discharge_summary['None'])
            self.summary_BT = self.mmtqsum(discharge_summary['BottomTrack'])
            self.summary_GGA = self.mmtqsum(discharge_summary['GGA'])
            self.summary_VTG  = self.mmtqsum(discharge_summary['VTG'])

        #if QA_QC exists create qaqc test fields
        if 'QA_QC' in win_river['Project'].keys():
            qaqc = win_river['Project']['QA_QC']

            for key, val in enumerate(qaqc):

                key = val
                val = qaqc[key]
                
                #extract qaqc data from dictionary if the field is a test, cal, or eval
                if key in ['RG_Test', 'Compass_Calibration', 'Compass_Evaluation']:
                    self.qaqc_test(key, val)

                if key == 'Moving_Bed_Test':
                    if 'Transect' in qaqc[key].keys():
                        transects_mb = qaqc[key]['Transect']
                        self.moving_bed_test(transects_mb, qaqc[key])

            
    def qaqc_test(self, key, val):
        '''Method to extract qaqc data from dictionary'''
        self.qaqc[key] = {}
        self.qaqc[key]['Type'] = val['@Type']
        self.qaqc[key]['Status'] = val['@Status']
        self.qaqc[key]['Checked'] =  val['@Checked']
        self.qaqc[key]['TestResult'] = []
        
                       
        if type(val['TestResult']) is not list:
            val['TestResult'] = [ val['TestResult'] ]
 
        #iterate though all of the keys and values in the test result dictionary
        for test_key, test_val in enumerate(val['TestResult']):

            self.qaqc[key]['TestResult'].append({'Text': test_val['Text'], 'TimeStamp': test_val['TimeStamp']})

    
    def moving_bed_test(self, transects, moving_bed):
        '''Method to extract data from moving bed test dictionary'''
        
        #make iterable list if not already a list of dictionaries
        if isinstance(transects,dict):
            transects = [transects]

        tr_idx = 0
        #iterate through each key and value of the list
#        for x,y in enumerate(transects):
        for y in transects:
            trans = Transect(y, True)
            self.mbt_transects.append(trans)
            typeAvailable = False

            if 'Moving_Bed_Test_Summary' in moving_bed.keys():
                if 'MB_Tests' in moving_bed['Moving_Bed_Test_Summary']:
                    typeAvailable = True

            #set moving bed test type for each transect object created previously
            if typeAvailable:
                mv_type = moving_bed['Moving_Bed_Test_Summary']['MB_Tests']['Index_%d'%tr_idx]['Type']
                if mv_type == '1':
                    trans.set_moving_bed_type('Loop')
                elif mv_type == '2':
                    trans.set_moving_bed_type('Stationary')
            else:

                file_name=trans.Files[0].Path
                fidx = file_name.rfind('.')
                if file_name[fidx+1:] == 'SBT':
                    trans.set_moving_bed_type('Stationary')
                elif file_name[fidx+1:] == 'LBT':
                    trans.set_moving_bed_type('Loop')
                else:
                    print ('question')
                    
                


    def mmtqsum(self, data):
        '''Method to get the MMT Q summary data'''
        
        sum_dict = {
            'Use': [],
            'Begin_Left': [],
            'FileName': [],
            'LeftEdgeSlopeCoeff': [],
            'RightEdgeSlopeCoeff': []
            }

        #Iterate through each key and val in the dictionary data
        for key, val in enumerate(data):

            indexField = data[val]
            
            #Iterate through each key and val in the dictionary index field
            for key2, val2 in enumerate(indexField):
               
                key2 = val2
                val2 = indexField[key2]

                #append value to appropriate list
                if key2 == 'UseInSummary':
                    sum_dict['Use'].append(float(val2))
                elif key2 == "BeginLeft":
                    sum_dict['Begin_Left'].append(float(val2))
                elif key2 == 'FileName':
                    sum_dict['FileName'].append(val2)
                elif key2 == 'LeftEdgeSlopeCoeff':
                    sum_dict['LeftEdgeSlopeCoeff'].append(float(val2))
                elif key2 == 'RightEdgeSlopeCoeff':
                    sum_dict['RightEdgeSlopeCoeff'].append(float(val2))
                else:
                    if key2 not in sum_dict:
                        sum_dict[key2] = []

                    sum_dict[key2].append(float(val2))
                    
                    
       


                
                



            



            

            















