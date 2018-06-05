import numpy as np
from Classes.SelectFit import SelectFit
from Classes.ExtrapQSensitivity import ExtrapQSensitivity
from Classes.NormData import NormData


class ComputeExtrap(object):
    """Class to compute the optimized or manually specified extrapolation methods

    Attributes
    ----------
    threshold: float
        Threshold as a percent for determining if a median is valid
    subsection: list
        Percent of discharge, does not account for transect direction
    fit_method: str
        Method used to determine fit.  Automatic or manual
    norm_data: NormData
        Object of class NormData
    sel_fit: SelectFit
        Object of class SelectFit
    q_sensitivity: ExtrapQSensitivity
        Object of class ExtrapQSensitivity
    messages: str
        Variable for messages to UserWarning

    """
    
    def __init__(self):
        """Initialize instance variables."""

        self.threshold = None  # Threshold as a percent for determining if a median is valid
        self.subsection = None  # Percent of discharge, does not account for transect direction
        self.fit_method = None  # Method used to determine fit.  Automatic or manual
        self.norm_data = []  # Object of class norm data
        self.sel_fit = []  # Object of class SelectFit
        self.q_sensitivity = None  # Object of class ExtrapQSensitivity
        self.messages = []  # Variable for messages to UserWarning
        
    def populate_data(self, transects, compute_sensitivity=True):
        """Store data in instance variables.

        Parameters
        ----------
        transects: list
            List of transects of TransectData
        compute_sensitivity: bool
            Determines is sensitivity should be computed.
        """

        self.threshold = 20
        self.subsection = [0, 100]
        self.fit_method = 'Automatic'
        self.process_profiles(transects=transects, data_type='q')
        # Compute the sensitivity of the final discharge to changes in extrapolation methods
        if compute_sensitivity:
            self.q_sensitivity = ExtrapQSensitivity()
            self.q_sensitivity.populate_data(transects=transects, extrap_fits=self.sel_fit)
            
    def process_profiles(self, transects, data_type):
        """Function that coordinates the fitting process.

        Parameters
        ----------
        transects: TransectData
            Object of TransectData
        data_type: str
            Type of data processing (q or v)
        """

        # Compute normalized data for each transect
        for transect in transects:
            norm_data = NormData()
            norm_data.populate_data(transect=transect,
                                    data_type=data_type,
                                    threshold=self.threshold,
                                    data_extent=self.subsection)
            self.norm_data.append(norm_data)

        # Compute composite normalized data
        comp_data = NormData()
        comp_data.create_composite(transects=transects, norm_data=self.norm_data, threshold=self.threshold)
        self.norm_data.append(comp_data)
        sel_fit = None

        # Compute the fit for the selected  method
        for n in range(len(self.norm_data)):
            if self.fit_method == 'Manual':
                sel_fit = SelectFit()
                sel_fit.populate_data(normalized=self.norm_data[n], fit_method=self.fit_method, transect=transects[n])
                self.sel_fit.append(sel_fit)
            else:
                sel_fit = SelectFit()
                sel_fit.populate_data(self.norm_data[n], self.fit_method)
                self.sel_fit.append(sel_fit)

        if sel_fit.top_fit_r2 is not None:
            # Evaluate if there is a potential that a 3-point top method may be appropriate
            if (sel_fit.top_fit_r2 > 0.9 or sel_fit.top_r2 > 0.9) and np.abs(sel_fit.top_max_diff) > 0.2:
                self.messages.append('The measurement profile may warrant a 3-point fit at the top')
                
    def update_q_sensitivity(self, transects):
        self.q_sensitivity = ExtrapQSensitivity()
        self.q_sensitivity.populate_data(transects, self.sel_fit)
        
    def change_fit_method(self, transect, new_fit_method, n, kargs = None):
        # TODO this function needs to be thought through only appears to be needed for setting individual transects to view in extrap window
        """Function to change the extrapolation methods associated with single transect"""
        self.fit_method = new_fit_method
        self.sel_fit = SelectFit()
        self.sel_fit.populate_data(self.norm_data, new_fit_method, kargs)
        # self.q_sensitivity = ExtrapQSensitivity()
        # self.q_sensitivity.populate_data(trans_data, self.sel_fit)
        
    def change_threshold(self, trans_data, data_type, threshold):
        """Function to change the threshold for accepting the increment median as valid.  The threshold
        is in percent of the median number of points in all increments"""
        
        self.threshold = threshold
        self.process_profiles(trans_data, data_type)
        self.q_sensitivity = ExtrapQSensitivity()
        self.q_sensitivity.populate_data(trans_data, self.sel_fit)
        
    def change_extents(self, trans_data, data_type, extents):
        """Function allows the data to be subsection by specifying the percent cumulative discharge
        for the start and end points.  Currently this function does not consider transect direction"""
        
        self.subsection = extents
        self.process_profiles(trans_data, data_type)
        self.q_sensitivity = ExtrapQSensitivity()
        self.q_sensitivity.populate_data(trans_data, self.sel_fit)
        
    def change_data_type(self, trans_data, data_type):
        self.process_profiles(trans_data, data_type)
        self.q_sensitivity = ExtrapQSensitivity(trans_data, self.selfit)
