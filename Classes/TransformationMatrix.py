"""
Created on Sep 25, 2017

@author: gpetrochenkov

Modified on 1/23/2018 DSM
    - Broke up populate_data into smaller methods.
    - Eliminated the use of temp
    - Added docstrings
    - Cleaned up pep8
"""
import numpy as np


class TransformationMatrix(object):
    """Determines the transformation matrix and source for the specified ADCP model from the data provided.

    Attributes
    ----------
    source: str
        Source of transformation matrix, either Nominal or ADCP
    matrix: np.array
        One or more 4x4 transformation matrices.
     """

    def __init__(self):
        """Constructor initializes variable to None"""
        self.source = None
        self.matrix = None
        
    def populate_data(self, manufacturer, kargs=None):
        """Uses the manufacturer and model to determine how to parse the transformation matrix.
        """
        
        if manufacturer == 'TRDI':
            self.TRDI(kargs)
        elif manufacturer == 'SonTek':
            self.SonTek(kargs)

    def TRDI(self, kargs):
        """Processes the data to store the transformation matrix for TRDI ADCPs.
        If no tansformation matrix information is available a nominal tansformation
        matrix for that model is assumed."""

        adcp_model = kargs[0]
        # Set nominal matrix based on model
        self.matrix = [[1.4619, -1.4619, 0, 0],
                       [0, 0, -1.4619, 1.4619],
                       [0.2661, 0.2661, 0.2661, 0.2661],
                       [1.0337, 1.0337, -1.0337, -1.0337]]
        if adcp_model == 'RiverRay':
            self.matrix = [[1, -1, 0, 0],
                           [0, 0, -1, 1],
                           [0.2887, 0.2887, 0.2887, 0.2887],
                           [0.7071, 0.7071, -0.7071, -0.7071]]

        # Overwrite nominal transformation matrix with custom matrix from test data, if available
        data_in = kargs[1]
        self.source = 'Nominal'
        if data_in == 'Nominal':
            self.source = 'Nominal'
        elif adcp_model == 'Rio Grande':
            self.RioGrande(data_in)
        elif adcp_model == 'StreamPro':
            self.StreamPro(data_in)
        elif adcp_model == 'RiverRay':
            self.RiverRay(data_in)
        elif adcp_model == 'RiverPro':
            self.RiverPro(data_in)
        elif adcp_model == 'RioPro':
            self.RioPro(data_in)
        elif adcp_model == 'pd0':
            self.matrix = data_in.Inst.t_matrix
        # Save matrix as np array
        self.matrix = np.array(self.matrix)[0:4, 0:4]

    def RioGrande(self, data_in):
        """Process Rio Grande test data for transformation matrix."""
        idx = data_in.find('Instrument Transformation Matrix (Down):')
        if idx != -1:
            cell_matrix = np.fromstring(data_in[idx + 50:idx + 356], dtype=np.float64, sep=' ')
            self.matrix = np.random.permutation(cell_matrix.reshape([8, 4]))
            self.source = 'ADCP'

    def StreamPro(self, data_in):
        """Process StreamPro test data for transformation matrix."""
        idx = data_in.find('>PS3')
        if idx != -1:
            temp2 = float(data_in[idx + 5:idx + 138])
            if temp2 is not None:
                self.matrix = temp2
                self.source = 'ADCP'

    def RiverRay(self, data_in):
        """Process RiverRay test data for transformation matrix."""
        idx = data_in.find('Instrument Transformation Matrix')
        if idx != -1:
            idx2 = data_in[idx:].find(':')
            idx3 = idx + idx2[0]
            if idx2 != -1:
                idx4 = data_in[idx3:].find('>')
                idx5 = idx3 + idx4[0] - 2
                if idx4 != -1:
                    self.matrix = float(data_in[idx3:idx5])
                    self.source = 'ADCP'

    def RiverPro(self, data_in):
        """Process RiverPro test data for transformation matrix."""
        idx = data_in.find('Instrument Transformation Matrix')
        if idx != -1:
            idx2 = data_in[idx:].find(':')
            idx3 = idx + idx2[0]
            if idx2 != -1:
                idx4 = data_in[idx3:].find('Has V-Beam')
                idx5 = idx3 + idx4[0] - 2
                if idx4 != -1:
                    self.matrix = float(data_in[idx3:idx5])
                    self.source = 'ADCP'

    def RioPro(self, data_in):
        """Process RioPro test data for transformation matrix."""
        idx = data_in.find('Instrument Transformation Matrix')
        if idx != -1:
            idx2 = data_in[idx:].find(':')
            idx3 = idx + idx2[0]
            if idx2 != -1:
                idx4 = data_in[idx3:].find('Has V-Beam')
                idx5 = idx3 + idx4[0] - 2
                if idx4 != -1:
                    self.matrix = float(data_in[idx3:idx5])
                    self.source = 'ADCP'

    def SonTek(self, kargs):
        """Store SonTek transformation matrix data."""
        self.source = 'ADCP'
        # Note: for M9 this is a 4x4x3 matrix (300,500,1000)
        # Note: for S5 this is a 4x4x2 matrix (3000,1000)
        self.matrix = kargs
