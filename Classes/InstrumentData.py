"""
Created on Aug 28, 2017

@author: gpetrochenkov

Modified 1/23/2018 DSM
    - Split TRDI and SonTek into separate methods from populate_data
    - Added code for SonTek
    - Added docstrings
    - Cleaned up PEP8
"""

import numpy as np
from Classes.TransformationMatrix import TransformationMatrix


class InstrumentData(object):
    """Container for charactersistics of the ADCP used to make the measurement

    Attributes
    ----------
    serial_num: str
        Serial number of ADCP.
    manufacturer: str
        Name of manufacturer.
    model: str
        Model name of ADCP.
    firmware: str
        Firmware version in the ADCP.
    frequency_kHz:
        Frequency or frequencies used by ADCP.
    beam_angle_deg:
        Angle of the beams from vertical in degrees.
    beam_pattern:
        Pattern of the beam angles, concave or convex.
    t_matrix: np array
        Transformation matrix or matrices for the ADCP.
    configuration_commands:
        Commands used to configure the instrument.
    """
     
    def __init__(self):
        """Constructor initializes the variables to None"""
        self.serial_num = None  # Serial number of ADCP
        self.manufacturer = None  # manufacturer of ADCP (SonTek, TRDI)
        self.model = None  # model of ADCP (Rio Grande, StreamPro, RiverRay, M9, S5)
        self.firmware = None  # firmware version
        self.frequency_kHz = None  # frquency of ADCP (could be "Multi")
        self.beam_angle_deg = None  # angle of beam from vertical
        self.beam_pattern = None  # pattern of beams
        self.t_matrix = None  # object of TransformationMatrix
        self.configuration_commands = None  # configuration commands sent to ADCP
        
    def populate_data(self, manufacturer, kargs = None):
        """Manages method calls for different manufacturers.

        Parameters
        ----------
        manufacturer: str
            Name of manufacturer.
        kargs:
            Variable list of data
            mmt_transect:
            pd0:
            mmt:
        """
        if manufacturer == 'TRDI':
            self.manufacturer = manufacturer
            self.TRDI(kargs)
        elif manufacturer == 'SonTek':
            self.manufacturer = manufacturer
            self.SonTek(kargs)

    def TRDI(self, kargs):
        """Populates the variables with data from TRDI ADCPs."""
        # Assign data passed through kargs
        mmt_transect = kargs[0]
        pd0 = kargs[1]
        mmt = kargs[2]

        # Identify proper configuration to use
        config = 'field_config'
        if len(kargs) > 2:
            if kargs[2] == 'MB':
                config = 'mbt_field_config'

        #Instrument frequency
        self.frequency_kHz = pd0.Inst.freq[0]

        #Firmware
        self.firmware = pd0.Inst.firm_ver[0]

        #Instrument beam angle and pattern
        self.beam_angle_deg = pd0.Inst.beam_ang[0]
        self.beam_pattern = pd0.Inst.pat[0]

        #Instrument characteristics
        mmt_site = getattr(mmt, 'site_info')
        mmt_config = getattr(mmt_transect, config)

        self.serial_num = mmt_site['ADCPSerialNmb']

        num = float(self.firmware)
        model_switch = np.floor(num)

        #----------------------Am having trouble finding Fixed_Commands, Wizard_Commands, and User_Commands
        if model_switch == 10:
            self.model = 'Rio Grande'

#                 self.configuration_commands = ['Fixed', mmt_config['Fixed_Commands'],
#                                                'Wizard', mmt_config['Wizard_Commands'],
#                                                'User', mmt_config['User_Commands']]

        elif model_switch == 31:
            self.model = 'StreamPro'
            self.frequency_kHz = 2000
#                 self.configuration_commands = ['Fixed', mmt_config['Fixed_Commands_StreamPro'],
#                                                'Wizard', mmt_config['Wizard_Commands'],
#                                                'User', mmt_config['User_Commands']]

        elif model_switch == 44:
            self.model = 'RiverRay'
#                 self.configuration_commands = ['Fixed', mmt_config['Fixed_Commands_StreamPro'],
#                                                'Wizard', mmt_config['Wizard_Commands'],
#                                                'User', mmt_config['User_Commands']]

        elif model_switch == 56:
            self.model = 'RiverPro'
            if pd0.Cfg.n_beams < 5:
                if 'RG_Test' in mmt.qaqc.keys():
                    idx = mmt.qaqc['RG_Test'].find('RioPro')
                    if idx != -1:
                        self.model = 'RioPro'

            if 'Fixed_Commands_RiverPro' in mmt_config.keys():
                Fixed_Commands = mmt_config['Fixed_Commands_RiverPro']
            else:
                Fixed_Commands = ' '

            if 'Wizard_Commands' in mmt_config.keys():
                Wizard_Commands = mmt_config['Wizard_Commands']
            else:
                Wizard_Commands = ' '

            if 'User_Commands' in mmt_config.keys():
                User_Commands = mmt_config['User_Commands']
            else:
                User_Commands = ' '

            self.configuration_commands = ['Fixed', Fixed_Commands,
                                           'Wizard', Wizard_Commands,
                                           'User', User_Commands]

        #Obtain transformation matrix from one of the available sources
        if np.isnan(pd0.Inst.t_matrix[0,0]) == False:
            self.t_matrix =  TransformationMatrix()
            self.t_matrix.populate_data('TRDI', kargs=['pd0', pd0])
        elif self.model == 'RiverRay':
            self.t_matrix = TransformationMatrix()
            self.t_matrix.populate_data('TRDI', kargs=[self.model, 'Nominal'])
        else:
            if isinstance(mmt.qaqc, dict):
                if 'RG_Test' in mmt.qaqc.keys():

                    self.t_matrix = TransformationMatrix()
                    self.t_matrix.populate_data('TRDI', kargs=[self.model, mmt.qaqc['RG_Test']['TestResult'][0]['Text']])

                elif 'Compass_Calibration' in mmt.qaqc.keys():

                    self.t_matrix = TransformationMatrix()
                    self.t_matrix.populate_data('TRDI', kargs=[self.model, mmt.qaqc['Compass_Calibration']['TestResult'][0]['Text']])

                elif 'Compass_Eval_Timestamp' in mmt.qaqc.keys():

                    self.t_matrix = TransformationMatrix()
                    self.t_matrix.populate_data('TRDI', kargs=[self.model, mmt.qaqc['Compass_Evaluation']['TestResult'][0]['Text']])

                else:
                    self.t_matrix.populate_data('TRDI', kargs=[self.model, 'Nominal'])
            else:
                self.t_matrix = TransformationMatrix()
                self.t_matrix.populate_data('TRDI', kargs=[self.model, 'Nominal'])

    def SonTek(self, kargs):
        rs = kargs
        self.serial_num = rs.System.SerialNumber
        self.frequency_kHz = rs.Transformation_Matrices.Frequency
        if self.frequency_kHz[2]>0:
            self.model = 'M9'
        else:
            self.model = 'S5'
        if rs.SystemHW is not None:
            revision = str(rs.SystemHW.FirmwareRevision)
            if len(revision) < 2:
                revision = '0' + revision
            self.firmware = str(rs.SystemHW.FirmwareVersion) + '.' + revision
        else:
            self.firmware = ''
        self.beam_angle_deg = 25
        self.beam_pattern = 'Convex'
        self.t_matrix = TransformationMatrix()
        self.t_matrix.populate_data('SonTek', rs.Transformation_Matrices.Matrix)
        self.configuration_commands = None
