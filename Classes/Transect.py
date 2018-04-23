class Transect(object):
    """Class to hold properties of MMT transect dictionary attributes.

    Attributes
    ----------
    Checked: int
    Files: list
    Notes: list
    """
    # TODO consider calling this MMT_Transect
    def __init__(self, trans):
        """Constructor immediately begins extraction of data"""
        
        self.Checked = int(trans['@Checked'])
        self.Files = []
        self.Notes = []
        self.field_config = None
        self.active_config = None
        self.moving_bed_type = None

        files = trans['File']

        # Create File classes for each file associated with transect
        if type(files) is list:
            for file in files:
                self.Files.append(self.file_dict(file))
        else:
            self.Files.append(self.file_dict(files))

        # Create Note classes for each file associated with transect
        if 'Note' in trans.keys():
            note = trans['Note']
            if type(note) is list:
                for n in note:
                    if type(trans['File']) is list:
                        self.Notes.append(self.note_dict(n, trans['File'][0]['@TransectNmb']))
                    else:
                        self.Notes.append(self.note_dict(n, trans['File']['@TransectNmb']))
            else:
                if type(trans['File']) is list:
                    self.Notes.append(self.note_dict(note, trans['File'][0]['@TransectNmb']))
                else:
                    self.Notes.append(self.note_dict(note, trans['File']['@TransectNmb']))

        # Create configuration dictionaries for each config attribute
        if type(trans['Configuration']) is list:
            for config in trans['Configuration']:
                if int(config['@Checked']) == 0:
                    self.field_config = self.parse_config(config)
                if int(config['@Checked']) == 1:
                    self.active_config = self.parse_config(config)
        else:
            if int(trans['Configuration']['@Checked']) == 0:
                    self.field_config = self.parse_config(trans['Configuration'])
            if int(trans['Configuration']['@Checked']) == 1:
                    self.active_config = self.parse_config(trans['Configuration'])

        # Assign active config to field config if there is no field config
        if self.field_config is None:
            self.field_config = self.active_config

    def set_moving_bed_type(self, mvb_type):
        """Setter for moving bed type in the case of MBT Transects

        Parameters
        ----------
        mvb_type: str
            Type of moving-bed test.
        """

        self.moving_bed_type = mvb_type

    @staticmethod
    def parse_config(config):
        """Method to parse configuration file from mmt xml.

        Parameters
        ----------
        config: dict
            Dictionary of configuration settings

        Returns
        -------
        config_dict: dict
            Processed dictionary of configuration settings
        """

        # Initialize dictionary for configuration
        config_dict = {}

        # Store all instrument commands
        command_groups = config['Commands']
        for group in command_groups.keys():
            config_dict[group] = []
            for key, command in command_groups[group].items():
                if key != '@Status':
                    config_dict[group].append(command)

        # Depth sounder configuration
        if 'Use_Depth_Sounder_In_Processing' in config['Depth_Sounder'].keys():
            if config['Depth_Sounder']['Use_Depth_Sounder_In_Processing']['#text'] == "YES":
                config_dict['DS_Use_Process'] = 1
            else:
                config_dict['DS_Use_Process'] = 0
        else:
            config_dict['DS_Use_Process'] = -1

        config_dict['DS_Transducer_Depth'] = float(config['Depth_Sounder']['Depth_Sounder_Transducer_Depth']['#text'])
        config_dict['DS_Transducer_Offset'] = float(config['Depth_Sounder']['Depth_Sounder_Transducer_Offset']['#text'])

        if config['Depth_Sounder']['Depth_Sounder_Correct_Speed_of_Sound']['#text'] == 'YES':
            config_dict['DS_Cor_Spd_Sound'] = 1
        else:
            config_dict['DS_Cor_Spd_Sound'] = 0

        config_dict['DS_Scale_Factor'] = float(config['Depth_Sounder']['Depth_Sounder_Scale_Factor']['#text'])

        # External heading configuration
        config_dict['Ext_Heading_Offset'] = float(config['Ext_Heading']['Offset']['#text'])

        if 'Use_Ext_Heading' in config['Ext_Heading'].keys():
            if config['Ext_Heading']['Use_Ext_Heading']['#text'] == 'NO':
                config_dict['Ext_Heading_Use'] = False
            else:
                config_dict['Ext_Heading_Use'] = True
        else:
            config_dict['Ext_Heading_Use'] = False

        # GPS configuration
        if 'GPS' in config.keys():
            config_dict['GPS_Time_Delay'] = config['GPS']['Time_Delay']['#text']

        # Discharge settings
        config_dict['Q_Top_Method'] = float(config['Discharge']['Top_Discharge_Estimate']['#text'])
        config_dict['Q_Bottom_Method'] = float(config['Discharge']['Bottom_Discharge_Estimate']['#text'])
        config_dict['Q_Power_Curve_Coeff'] = float(config['Discharge']['Power_Curve_Coef']['#text'])
        config_dict['Q_Cut_Top_Bins'] = float(config['Discharge']['Cut_Top_Bins']['#text'])
        config_dict['Q_Bins_Above_Sidelobe'] = float(config['Discharge']['Cut_Bins_Above_Sidelobe']['#text'])
        config_dict['Q_Left_Edge_Type'] = float(config['Discharge']['River_Left_Edge_Type']['#text'])
        config_dict['Q_Left_Edge_Coeff'] = float(config['Discharge']['Left_Edge_Slope_Coeff']['#text'])
        config_dict['Q_Right_Edge_Type'] = float(config['Discharge']['River_Right_Edge_Type']['#text'])
        config_dict['Q_Right_Edge_Coeff'] = float(config['Discharge']['Right_Edge_Slope_Coeff']['#text'])
        config_dict['Q_Shore_Pings_Avg'] = float(config['Discharge']['Shore_Pings_Avg']['#text'])

        # Edge estimate settings
        config_dict['Edge_Begin_Shore_Distance'] = config['Edge_Estimates']['Begin_Shore_Distance']['#text']
        if config['Edge_Estimates']['Begin_Left_Bank']['#text'] == 'YES':
            config_dict['Edge_Begin_Left_Bank'] = 1
        else:
            config_dict['Edge_Begin_Left_Bank'] = 0

        config_dict['Edge_End_Shore_Distance'] = float(config['Edge_Estimates']['End_Shore_Distance']['#text'])
        
        # Offsets
        for key in config['Offsets'].keys():
            if key == 'ADCP_Transducer_Depth':
                child = "Offsets_Transducer_Depth"
            else:
                child = "Offsets_" + key

            config_dict[child] = float(config['Offsets'][key]['#text'])

        # Processing settings
        for key in config['Processing'].keys():
            if key == 'Use_3_Beam_Solution_For_BT':
                child = 'Proc_Use_3_Beam_BT'
            elif key == 'Use_3_Beam_Solution_For_WT':
                child = 'Proc_Use_3_Beam_WT'
            elif key == 'BT_Error_Velocity_Threshold':
                child = 'Proc_BT_Error_Vel_Threshold'
            elif key == 'WT_Error_Velocity_Threshold':
                child = 'Proc_WT_Error_Velocity_Threshold'
            elif key == 'BT_Up_Velocity_Threshold':
                child = 'Proc_BT_Up_Vel_Threshold'
            elif key == 'WT_Up_Velocity_Threshold':
                child = 'Proc_WT_Up_Vel_Threshold'
            elif key == 'Fixed_Speed_Of_Sound':
                child = 'Proc_Fixed_Speed_Of_Sound'
            elif key == 'Mark_Below_Bottom_Bad':
                child = 'Proc_Mark_Below_Bottom_Bad'
            elif key == 'Use_Weighted_Mean':
                child = 'Proc_Use_Weighted_Mean'
            elif key == 'Absorption':
                child = 'Proc_Absorption'
            else:
                child = 'Proc_' + key

            # Try to cast to float otherwise assign 1 or 0 based on string value
            try:
                config_dict[child] = float(config['Processing'][key]['#text'])
            except KeyError:
                if config['Processing'][key]['#text'] == 'YES':
                    config_dict[child] = 1
                else:
                    config_dict[child] = 0

            # Recording
            config_dict['Rec_Filename_Prefix'] = config['Recording']['Filename_Prefix']['#text']
            config_dict['Rec_Output_Directory'] = config['Recording']['Output_Directory']['#text']
            
            if 'Root_Directory' in config['Recording'].keys():
                if '#text' in config['Recording']['Root_Directory']:
                    config_dict['Rec_Root_Directory'] = config['Recording']['Root_Directory']['#text']
                else:
                    config_dict['Rec_Root_Directory'] = None
            else:
                config_dict['Rec_Root_Directory'] = None

            if config['Recording']['MeasurmentNmb'] is None:
                config_dict['Rec_MeasNmb'] = config['Recording']['MeasurmentNmb']
            else:
                config_dict['Rec_MeasNmb'] = float(config['Recording']['MeasurmentNmb'])
            config_dict['Rec_GPS'] = config['Recording']['GPS_Recording']['#text']
            config_dict['Rec_DS'] = config['Recording']['DS_Recording']['#text']
            config_dict['Rec_EH'] = config['Recording']['EH_Recording']['#text']
            config_dict['Rec_ASCII_Output'] = config['Recording']['ASCII_Output_Recording']['#text']
            config_dict['Rec_Max_File_Size'] = float(config['Recording']['Maximum_File_Size']['#text'])
            config_dict['Rec_Next_Transect_Number'] = float(config['Recording']['Next_Transect_Number']['#text'])
            config_dict['Rec_Add_Date_Time'] = float(config['Recording']['Add_Date_Time']['#text'])
            config_dict['Rec_Use_Delimiter'] = config['Recording']['Use_Delimiter']['#text']
            config_dict['Rec_Delimiter'] = config['Recording']['Custom_Delimiter']['#text']
            config_dict['Rec_Prefix'] = config['Recording']['Use_Prefix']['#text']
            config_dict['Rec_Use_MeasNmb'] = config['Recording']['Use_MeasurementNmb']['#text']
            config_dict['Rec_Use_TransectNmb'] = config['Recording']['Use_TransectNmb']['#text']
            config_dict['Rec_Use_SequenceNmb'] = config['Recording']['Use_SequenceNmb']['#text']

            # Wizard settings
            config_dict['Wiz_ADCP_Type'] = float(config['Wizard_Info']['ADCP_Type'])
            config_dict['Wiz_Firmware'] = float(config['Wizard_Info']['ADCP_FW_Version'])
            config_dict['Wiz_Use_Ext_Heading'] = config['Wizard_Info']['Use_Ext_Heading']
            config_dict['Wiz_Use_GPS'] = config['Wizard_Info']['Use_GPS']
            config_dict['Wiz_Use_DS'] = config['Wizard_Info']['Use_Depth_Sounder']
            config_dict['Wiz_Max_Water_Depth'] = float(config['Wizard_Info']['Max_Water_Depth'])
            config_dict['Wiz_Max_Water_Speed'] = float(config['Wizard_Info']['Max_Water_Speed'])
            config_dict['Wiz_Max_Boat_Space'] = float(config['Wizard_Info']['Max_Boat_Speed'])
            config_dict['Wiz_Material'] = float(config['Wizard_Info']['Material'])
            config_dict['Wiz_Water_Mode'] = float(config['Wizard_Info']['Water_Mode'])
            config_dict['Wiz_Bottom_Mode'] = float(config['Wizard_Info']['Bottom_Mode'])
            config_dict['Wiz_Beam_Angle'] = float(config['Wizard_Info']['Beam_Angle'])
            config_dict['Wiz_Pressure_Sensor'] = config['Wizard_Info']['Pressure_Sensor']
            config_dict['Wiz_Water_Mode_13'] = float(config['Wizard_Info']['Water_Mode_13_Avail'])
            config_dict['Wiz_StreamPro_Default'] = float(config['Wizard_Info']['Use_StreamPro_Def_Cfg'])
            config_dict['Wiz_StreamPro_Bin_Size'] = float(config['Wizard_Info']['StreamPro_Bin_Size'])
            config_dict['Wiz_StreamPro_Bin_Number'] = float(config['Wizard_Info']['StreamPro_Bin_Num'])

            if 'Use_GPS_Internal' in config['Wizard_Info'].keys():
                config_dict['Wiz_Use_GPS_Internal'] = config['Wizard_Info']['Use_GPS_Internal']
            if 'Internal_GPS_Baud_Rate_Index' in config['Wizard_Info'].keys():
                config_dict['Wiz_Internal_GPS_Baud_Rate_Index'] = float(config['Wizard_Info']
                                                                        ['Internal_GPS_Baud_Rate_Index'])

        return config_dict

    @staticmethod
    def file_dict(file):
        """Create dictionary for file information.

        Parameters
        ----------
        file: dict
            Dictionary for file from mmt

        Returns
        -------
        transect_file: dict
            Dictionary of transect file information
                Path: str
                    Full filename of transect including path
                File: str
                    Filename of transect
                Number: str
                    Transect number assigned in WinRiver 2
        """

        transect_file = {'Path': file['@PathName'], 'File': file['#text'], 'Number': file['@TransectNmb']}
        return transect_file

    @staticmethod
    def note_dict(note, number):
        """Create dictionary for notes.

        Parameters
        ----------
        note: dict
            Dictionary from mmt for notes
        number: str
            Transect number

        Returns
        -------
        note_dict_out: dict
            Dictionary for note information
                NoteFileNo: str
                    Transect number associated with the note
                NoteDate: str
                    Date note was entered
                NoteText: str
                    Text of note
        """

        note_dict_out = {'NoteFileNo': number, 'NoteDate': note['@TimeStamp'], 'NoteText': note['@Text']}
        return note_dict_out
