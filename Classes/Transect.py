class Transect(object):
    """Class to hold properties of MMT transect dictionary attributes"""

    def __init__(self, trans, mbt = False):
        """Constructor immediately begins extraction of data"""
        
        self.Checked = int(trans['@Checked'])
        self.Files = []

        files = trans['File']

        #Create File classes for each file associated with transect
        if type(files) is list:
            for file in files:
                self.Files.append(File(file))
        else:
            self.Files.append(File(files))

       
        self.Notes = []
        #Create Note classes for each file associated with transect
        if 'Note' in trans.keys():
            
            note = trans['Note']
            if type(note) is list:
                for n in note:
                    self.Notes.append(Note(n))
            else:
                self.Notes.append(Note(note))
        
        self.field_config = None
        self.active_config = None
        self.mbt_active_config = None
        self.mbt_field_config = None

        if mbt == False:
            #Create configuration dictionaries for each config attribute
            if type(trans['Configuration']) is list:
                for x in trans['Configuration']:
                    if int(x['@Checked']) == 0:
                        self.field_config = self.config(x)
                    if int(x['@Checked']) == 1:
                        self.active_config = self.config(x)
            else:
                if int(trans['Configuration']['@Checked']) == 0:
                        self.field_config = self.config(trans['Configuration'])
                if int(trans['Configuration']['@Checked']) == 1:
                        self.active_config = self.config(trans['Configuration'])
    
            #assign active config to field config if there is no field config
            if self.field_config is None:
                self.field_config = self.active_config
        else:
            #Create configuration dictionaries for each config attribute
            if type(trans['Configuration']) is list:
                for x in trans['Configuration']:
                    if int(x['@Checked']) == 0:
                        self.mbt_field_config = self.config(x)
                    if int(x['@Checked']) == 1:
                        self.mbt_active_config = self.config(x)
            else:
                if int(trans['Configuration']['@Checked']) == 0:
                        self.mbt_field_config = self.config(trans['Configuration'])
                if int(trans['Configuration']['@Checked']) == 1:
                        self.mbt_active_config = self.config(trans['Configuration'])
    
            #assign active config to field config if there is no field config
            if self.mbt_field_config is None:
                self.mbt_field_config = self.mbt_active_config
            
        

    def set_moving_bed_type(self, mvb_type):
        """Setter for moving bed type in the case of MBT Transects"""
        self.moving_bed_type = mvb_type

    def config(self, config):
        """Method to parse configuration file from mmt xml"""
        
        config_dict = {}

        commands = config['Commands']
        #iterate through each command and add to dictionary
        field_commands = commands.keys()
        for x in field_commands:
            f = commands[x].keys()
            if f == "#text":
                config_dict[x] = commands[x]['#text']
            else:
                for k in f:
                    config_dict[k] = commands[x][k]


        #Depth Sounder calls--------------------------------------------------------------------
        ds = config['Depth_Sounder']
        if 'Use_Depth_Sounder_In_Processing' in ds.keys():
            if ds['Use_Depth_Sounder_In_Processing']['#text'] == "YES":
                config_dict['DS_Use_Process'] = 1
            else:
                config_dict['DS_Use_Process'] = 0
        else:
            config_dict['DS_Use_Process'] = -1

        config_dict['DS_Transducer_Depth'] = float(ds['Depth_Sounder_Transducer_Depth']['#text'])
        config_dict['DS_Transducer_Offset'] = float(ds['Depth_Sounder_Transducer_Offset']['#text'])

        if ds['Depth_Sounder_Correct_Speed_of_Sound']['#text'] == 'YES':
            config_dict['DS_Cor_Spd_Sound'] = 1
        else:
            config_dict['DS_Cor_Spd_Sound'] = 0

        config_dict['DS_Scale_Factor'] = float(ds['Depth_Sounder_Scale_Factor']['#text'])

        #Ext Heading Calls -------------------------------------------------------------------------
        eh = config['Ext_Heading']
        config_dict['Ext_Heading_Offset'] = float(eh['Offset']['#text'])

        if 'Use_Ext_Heading' in eh.keys():
            if eh['Use_Ext_Heading']['#text'] == 'NO':
                config_dict['Ext_Heading_Use'] = False
            else:
                config_dict['Ext_Heading_Use'] = True
        else:
            config_dict['Ext_Heading_Use'] = True

        #GPS Calls -----------------------------------------------------------------------------------
        if 'GPS' in config.keys():
            config_dict['GPS_Time_Delay'] = config['GPS']['Time_Delay']['#text']

        #Discharge Calls -----------------------------------------------------------------------------
        dis = config['Discharge']
        config_dict['Q_Top_Method'] = float(dis['Top_Discharge_Estimate']['#text'])
        config_dict['Q_Bottom_Method'] = float(dis['Bottom_Discharge_Estimate']['#text'])
        config_dict['Q_Power_Curve_Coeff'] = float(dis['Power_Curve_Coef']['#text'])
        config_dict['Q_Cut_Top_Bins'] = float(dis['Cut_Top_Bins']['#text'])
        config_dict['Q_Bins_Above_Sidelobe'] = float(dis['Cut_Bins_Above_Sidelobe']['#text'])
        config_dict['Q_Left_Edge_Type'] = float(dis['River_Left_Edge_Type']['#text'])
        config_dict['Q_Left_Edge_Coeff'] = float(dis['Left_Edge_Slope_Coeff']['#text'])
        config_dict['Q_Right_Edge_Type'] = float(dis['River_Right_Edge_Type']['#text'])
        config_dict['Q_Right_Edge_Coeff'] = float(dis['Right_Edge_Slope_Coeff']['#text'])
        config_dict['Q_Shore_Pings_Avg'] = float(dis['Shore_Pings_Avg']['#text'])

        #Edge Estimate calls ---------------------------------------------------------------------------
        ee = config['Edge_Estimates']
        config_dict['Edge_Begin_Shore_Distance'] = ee['Begin_Shore_Distance']['#text']
        if ee['Begin_Left_Bank']['#text'] == 'YES':
            config_dict['Edge_Begin_Left_Bank'] = 1
        else:
            config_dict['Edge_Begin_Left_Bank'] = 0

        config_dict['Edge_End_Shore_Distance'] = float(ee['End_Shore_Distance']['#text'])
        
        #Offsets calls -------------------------------------------------------------------------------------
        offset = config['Offsets']
        for x in offset.keys():
            if x == 'ADCP_Transducer_Depth':
                child = "Offsets_Transducer_Depth"
            else:
                child = "Offsets_%s" % x

            config_dict[child] = float(offset[x]['#text'])

        #Processing Calls ------------------------------------------------------------------------------------
        pr = config['Processing']
        prefix = 'Proc_'
        for x in pr.keys():
            if x == 'Use_3_Beam_Solution_For_BT':
                child = prefix + 'Use_3_Beam_BT'
            if x == 'Use_3_Beam_Solution_For_WT':
                child = prefix + 'Use_3_Beam_WT'
            if x == 'BT_Error_Velocity_Threshold':
                child = prefix + 'BT_Error_Vel_Threshold'
            if x == 'WT_Error_Velocity_Threshold':
                child = prefix + 'WT_Error_Velocity_Threshold'
            if x == 'BT_Up_Velocity_Threshold':
                child = prefix + 'BT_Up_Vel_Threshold'
            if x == 'WT_Up_Velocity_Threshold':
                child = prefix + 'WT_Up_Vel_Threshold'
            if x == 'Fixed_Speed_Of_Sound':
                child = prefix + "Fixed_Speed_Of_Sound"
            if x == 'Mark_Below_Bottom_Bad':
                child = prefix + "Mark_Below_Bottom_Bad"
            if x == 'Use_Weighted_Mean':
                child = prefix + 'Use_Weighted_Mean'
            if x == 'Absorption':
                child = prefix + 'Absorption'
            else:
                child = prefix + x

            #try to cast to float otherwise assign 1 or 0 based on string value
            try:
                config_dict[child] = float(pr[x]['#text'])
            except:
                if pr[x]['#text'] == 'YES':
                    config_dict[child] = 1
                else:
                    config_dict[child] = 0
         
               
            #Recording calls----------------------------------------------------------
            rec = config['Recording']
            config_dict['Rec_Filename_Prefix'] = rec['Filename_Prefix']['#text']
            config_dict['Rec_Output_Directory'] = rec['Output_Directory']['#text']
            
            if 'Root_Directory' in rec.keys():
                config_dict['Rec_Root_Directory'] = rec['Root_Directory']['#text']
            else:
                config_dict['Rec_Root_Directory'] = None
            
           
            if rec['MeasurmentNmb'] is None:
                config_dict['Rec_MeasNmb'] = rec['MeasurmentNmb']
            else:
                config_dict['Rec_MeasNmb'] = float(rec['MeasurmentNmb'])
            config_dict['Rec_GPS'] = rec['GPS_Recording']['#text']
            config_dict['Rec_DS'] = rec['DS_Recording']['#text']
            config_dict['Rec_EH'] = rec['EH_Recording']['#text']
            config_dict['Rec_ASCII_Output'] = rec['ASCII_Output_Recording']['#text']
            config_dict['Rec_Max_File_Size'] = float(rec['Maximum_File_Size']['#text'])
            config_dict['Rec_Next_Transect_Number'] = float(rec['Next_Transect_Number']['#text'])
            config_dict['Rec_Add_Date_Time'] = float(rec['Add_Date_Time']['#text'])
            config_dict['Rec_Use_Delimiter'] = rec['Use_Delimiter']['#text']
            config_dict['Rec_Delimiter'] = rec['Custom_Delimiter']['#text']
            config_dict['Rec_Prefix'] = rec['Use_Prefix']['#text']
            config_dict['Rec_Use_MeasNmb'] = rec['Use_MeasurementNmb']['#text']
            config_dict['Rec_Use_TransectNmb'] = rec['Use_TransectNmb']['#text']
            config_dict['Rec_Use_SequenceNmb'] = rec['Use_SequenceNmb']['#text']

            #Wizard info calls -------------------------------------------------------------------------
            wiz = config['Wizard_Info']
            config_dict['Wiz_ADCP_Type'] = float(wiz['ADCP_Type'])
            config_dict['Wiz_Firmware'] = float(wiz['ADCP_FW_Version'])
            config_dict['Wiz_Use_Ext_Heading'] = wiz['Use_Ext_Heading']
            config_dict['Wiz_Use_GPS'] = wiz['Use_GPS']
            config_dict['Wiz_Use_DS'] = wiz['Use_Depth_Sounder']
            config_dict['Wiz_Max_Water_Depth'] = float(wiz['Max_Water_Depth'])
            config_dict['Wiz_Max_Water_Speed'] = float(wiz['Max_Water_Speed'])
            config_dict['Wiz_Max_Boat_Space'] = float(wiz['Max_Boat_Speed'])
            config_dict['Wiz_Material'] = float(wiz['Material'])
            config_dict['Wiz_Water_Mode'] = float(wiz['Water_Mode'])
            config_dict['Wiz_Bottom_Mode'] = float(wiz['Bottom_Mode'])
            config_dict['Wiz_Beam_Angle'] = float(wiz['Beam_Angle'])
            config_dict['Wiz_Pressure_Sensor'] = wiz['Pressure_Sensor']
            config_dict['Wiz_Water_Mode_13'] = float(wiz['Water_Mode_13_Avail'])
            config_dict['Wiz_StreamPro_Default'] = float(wiz['Use_StreamPro_Def_Cfg'])
            config_dict['Wiz_StreamPro_Bin_Size'] = float(wiz['StreamPro_Bin_Size'])
            config_dict['Wiz_StreamPro_Bin_Number'] = float(wiz['StreamPro_Bin_Num'])

            if 'Use_GPS_Internal' in wiz.keys():
                config_dict['Wiz_Use_GPS_Internal'] = wiz['Use_GPS_Internal']
            if 'Internal_GPS_Baud_Rate_Index' in wiz.keys():
                config_dict['Wiz_Internal_GPS_Baud_Rate_Index'] = float(wiz['Internal_GPS_Baud_Rate_Index'])

        return config_dict


class File(object):
    """Class for files in a transect"""

    def __init__(self, file):
        self.Path = file['@PathName']
        self.File = file['#text']
        self.Number = file['@TransectNmb']

class Note(object):
    """Class for notes in a transect"""
    def __init__(self, note):

        if 'Number' in note.keys():
            self.NoteFileNo = note['@Number']

        self.NoteDate = note['@TimeStamp']
        self.NoteText = note['@Text']


