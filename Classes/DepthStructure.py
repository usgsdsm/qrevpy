"""
Created on Jul 19, 2017

@author: gpetrochenkov

Modified DSM 1/24/2018
    - added cellDepth_in, and cellSize_in as required parameters
    - modified code to support above change
    - added docstrings
    - cleaned up PEP8
    - only modified __init__ adn add_depth_object
"""

from Classes.DepthData import DepthData
import numpy as np
from numpy.matlib import repmat

class DepthStructure(object):
    """This class creates the data structure used store depths from different sources

    Attributes
    ----------
    selected: str
        Name of object DepthData that contains depth data.
    bt_depths: DepthData
        Object of DepthData for bottom track based depths.
    vb_depths: DepthData
        Object of DepthData for vertical beam based depths.
    ds_depths: DepthData
        Object of DepthData for depth sounder based depths.
    composite: bool
        Indicates use of composite depths.
    """
    
    def __init__(self):
        """Creates object and initializes variables to None"""
        self.selected = None  # name of object DepthData that contains the depth data for q computation
        self.bt_depths = None # object of DepthData for by depth data
        self.vb_depths = None # object of DepthData for vertical beam depth data
        self.ds_depths = None # object of DepthData for depth sounder depth data
        self.composite = None # Turn composite depths "on" or "off"
        
        
    # DSM changed 1/24/2018 def add_depth_object(self, depth_in, source_in, freq_in, draft_in, kargs = None):
    def add_depth_object(self, depth_in, source_in, freq_in, draft_in, cell_depth_in, cell_size_in):
        """Adds a DepthData object to the depth structure for the specified type of depths.

        Parameters
        ----------
        depth_in: np.array
            Depth data in meters.
        source_in: str
            Specifies source of depth data: bottom track (BT), vertical beam (VB), or depth sounder (DS)
        freq_in: np.array
            Acoustic frequency in kHz of beams used to determine depth.
        draft_in:
            Draft of transducer (in meters) used to measure depths.
        cell_depth_in
            Depth of each cell in the profile. If the referenced depth does not have depth cells the depth cell
            values from the bottom track (BT) depths should be used.
        cell_size_in
            Size of each depth cell. If the referenced depth does not have depth cells the cell size from
            the bottom track (BT) depths should be used.
        """
        
        # if kargs is not None:
        #     cell_depth = kargs[0]
        #     cell_size = kargs[1]
            
        if source_in == 'BT':
            self.bt_depths = DepthData()
            # DSM changed 1/24/2018 self.bt_depths.populate_data(depth_in, source_in, freq_in, draft_in, kargs=[cell_depth, cell_size])
            self.bt_depths.populate_data(depth_in, source_in, freq_in, draft_in, cell_depth_in, cell_size_in)
        elif source_in == 'VB':
            self.vb_depths = DepthData()
            self.vb_depths.populate_data(depth_in, source_in, freq_in, draft_in, cell_depth_in, cell_size_in)
            # self.vb_depths.add_cell_data(self.bt_depths)
        elif source_in == 'DS':
            self.ds_depths = DepthData()
            self.ds_depths.populate_data(depth_in, source_in, freq_in, draft_in, cell_depth_in, cell_size_in)
            # self.ds_depths.add_cell_data(self.bt_depths)

    def set_depth_reference(self, reference):
        """Set the depth reference to used the specified depth source.

        Parameters
        ----------
        reference: str
            Depth reference identifier.
        """
        # TODO Clean up calling methods to use a consistent reference format.
        if reference == 'BT':
            self.selected = 'bt_depths'
        if reference == 'btDepths':
            self.selected = 'bt_depths'
        if reference == 'VB':
            self.selected = 'vb_depths'
        if reference == 'vbDepths':
            self.selected = 'vb_depths'
        if reference == 'DS':
            self.selected = 'ds_depths'
        if reference == 'dsDepths':
            self.selected = 'ds_depths'

    # =========================
    # DSM stopped here 1/25/2018
    # =========================
    def set_valid_data_method(self, setting):
        self.bt_depths.set_valid_data_method(setting)
        
    def composite_depths(self, transect, setting=None):
        """Depth compositing is based on the following assumptions
        
        1. If a depth sounder is available the user must have assumed the ADCP beams
        (BT or vertical) might have problems and it will be the second alternative if 
        not selected as the preferred source
        
        2. For 4-beam BT depths, if 3 beams are valid the average is considered valid.
        It may be based on interpolation of the invalid beam.  However, if only 2 beams
        are valid even though the other two beams may be interpolated and included in the average the
        average will be replaced by an alternative if available.  If no alternative is 
        available the multi-beam average based on available beams and interpolation will
        be used.

        Parameters
        ----------
        transect: object
            Transect object containing all data.
        setting: bool
            Setting to use (True) or not use (False) composite depths.
        """
        
        if setting is None:
            setting = self.composite
        else:
            self.composite = setting
            
        #The primary depth reference is the selected reference
        ref = self.selected
        
        if setting:
            #Prepare vector of valid BT averages, which are defined as
            #having at least 2 valid beams
            
            bt_valid = self.bt_depths.valid_data
            n_ensembles = bt_valid.shape[-1]
            bt_filtered = self.bt_depths.depth_processed_m
            bt_filtered[bt_filtered == 0] = np.nan
            
            #Prepare vertical beam data, using only data prior to interpolation
            if self.vb_depths is not None:
                vb_filtered = self.vb_depths.depth_processed_m
                vb_filtered[np.squeeze(self.vb_depths.valid_data) == False] = np.nan
            else:
                vb_filtered = repmat([np.nan],n_ensembles,1)

                  
            #Prepare depth sounder data, using only data prior to interpolation
            if self.ds_depths is not None:
                ds_filtered = self.ds_depths.depth_processed_m
                ds_filtered[np.squeeze(self.ds_depths.valid_data) == False] = np.nan
            else:
                ds_filtered = repmat([np.nan],n_ensembles,1)
                
            if len(bt_filtered.shape) > 1:
                comp_source = repmat([np.nan],bt_filtered.shape[0],bt_filtered.shape[1])
            else:
                comp_source = repmat([np.nan],bt_filtered.shape[0],1)
                comp_source = np.squeeze(comp_source)
            
            
            if ref == 'bt_depths':
                comp_depth = bt_filtered
                comp_source[np.isnan(comp_depth) == False] = 1
                comp_depth[np.isnan(comp_depth)] = np.squeeze(ds_filtered[np.isnan(comp_depth)])
                comp_source[(np.isnan(comp_depth) == False) & (np.isnan(comp_source) == True)] = 3
                comp_depth[np.isnan(comp_depth)] = np.squeeze(vb_filtered[np.isnan(comp_depth)])
                comp_source[(np.isnan(comp_depth == False) & (np.isnan(comp_source) == True))] = 2
                comp_depth[np.isnan(comp_depth)] = np.squeeze(self.bt_depths.depth_processed_m[np.isnan(comp_depth)])
                comp_source[(np.isnan(comp_depth) == False) & (np.isnan(comp_source) == True)] = 4
                
            elif ref == 'vb_depths':
                comp_depth = vb_filtered
                comp_source[np.isnan(comp_depth) == False] = 2
                comp_depth[np.isnan(comp_depth)] = np.squeeze(ds_filtered[np.isnan(comp_depth)])
                comp_source[(np.isnan(comp_depth) == False) & (np.isnan(comp_source) == True)] = 3
                comp_depth[np.isnan(comp_depth)] = np.squeeze(bt_filtered[np.isnan(comp_depth)])
                comp_source[(np.isnan(comp_depth) == False) & (np.isnan(comp_source) == True)] = 1
                comp_depth[np.isnan(comp_depth)] = np.squeeze(self.vb_depths.depth_processed_m[np.isnan(comp_depth)])
                comp_source[(np.isnan(comp_depth) == False) & (np.isnan(comp_source) == True)] = 4
                
            elif ref == 'ds_depths':
                comp_depth = ds_filtered
                comp_source[np.isnan(comp_depth) == False] = 3
                comp_depth[np.isnan(comp_depth)] = np.squeeze(vb_filtered[np.isnan(comp_depth)])
                comp_source[(np.isnan(comp_depth) == False) & (np.isnan(comp_source) == True)] = 2
                comp_depth[np.isnan(comp_depth)] = np.squeeze(bt_filtered[np.isnan(comp_depth)])
                comp_source[(np.isnan(comp_depth) == False) & (np.isnan(comp_source) == True)] = 1
                comp_depth[np.isnan(comp_depth)] = np.squeeze(self.ds_depths.depth_processed_m[np.isnan(comp_depth)])
                comp_source[(np.isnan(comp_depth) == False) & (np.isnan(comp_source) == True)] = 4
               
            #Save composite depth to depth_processed of selected primary reference
            
            obj = getattr(self,ref)
            obj.apply_composite(comp_depth, comp_source)
                
        else:
            obj = getattr(self,ref)
            comp_source = np.zeros(obj.depth_processed_m.shape)
            
            if ref == 'bt_depths':
                obj.valid_data[np.isnan(obj.valid_data)] = False
                comp_source[obj.valid_data] = 1
            elif ref == 'vb_depths':
                comp_source[obj.valid_data] = 2
            elif ref == 'ds_depths':
                comp_source[obj.valid_data] = 3
                
            
            obj.apply_interpolation(transect) 
            comp_depth = obj.depth_processed_m
            obj.apply_composite(comp_depth, comp_source)
            
            
    def set_draft(self, target, draft):
        """This function will change the ref_depth draft/  The associated
        depth object will also be updated because DepthData is a handle class.
        The computations are actually done in DepthData as the data are private
        to that class """
        
        if target == 'ADCP':
            self.bt_depths.change_draft(draft)
            self.vb_depths.change_draft(draft)
        else:
            self.ds_depths.change_draft(draft)    
            
    def depth_filter(self,transect,kargs):
        """Method to apply filter to all available depth sources, so that
        all sources have the same filter applied """
        
        if self.bt_depths is not None:
            self.bt_depths.apply_filter(transect, kargs)
        if self.vb_depths is not None:
            self.vb_depths.apply_filter(transect, kargs)
        if self.ds_depths is not None:
            self.ds_depths.apply_filter(transect, kargs)
            
    def depth_interpolation(self, transect, kargs):
        """Method to apply interpolation to all available depth sources, so
        that all sources have the same filter applied """
        
        if self.bt_depths is not None:
            self.bt_depths.apply_interpolation(transect, kargs)
        if self.vb_depths is not None:
            self.vb_depths.apply_interpolation(transect, kargs)
        if self.ds_depths is not None:
            self.vb_depths.apply_interpolation(transect, kargs)
            
    def sos_correction(self, ratio):
        """Correct depths for change in speed of sound"""
        
        #Bottom Track Depths
        if self.bt_depths is not None:
            self.bt_depths.sos_correction(ratio)
            
        #Vertical beam depths
        if self.vb_depths is not None:
            self.vb_depths.sos_correction(ratio)
