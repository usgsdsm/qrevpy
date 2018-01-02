'''
Created on Jul 19, 2017

@author: gpetrochenkov
'''
from Classes.DepthData import DepthData
import numpy as np
from numpy.matlib import repmat


class DepthStructure(object):
    
    def __init__(self):
        self.selected = None  # name of object DepthData that contains the depth data for q computation
        self.bt_depths = None # object of DepthData for by depth data
        self.vb_depths = None # object of DepthData for vertical beam depth data
        self.ds_depths = None # object of DepthData for depth sounder depth data
        self.composite = None # Turn composite depths "on" or "off"
        
        
    def add_depth_object(self, depth_in, source_in, freq_in, draft_in, kargs = None):
        '''Adds a DepthData object to the appropriate property
        if kargs is not empty the contents need to be split into
        variables being passed to DepthData '''
        
        if kargs is not None:
            cell_depth = kargs[0]
            cell_size = kargs[1]
            
        if source_in == 'BT':
            self.bt_depths = DepthData()
            self.bt_depths.populate_data(depth_in, source_in, freq_in, draft_in, kargs=[cell_depth, cell_size])
        elif source_in == 'VB':
            self.vb_depths = DepthData()
            self.vb_depths.populate_data(depth_in, source_in, freq_in, draft_in, None)
            self.vb_depths.add_cell_data(self.bt_depths)
        elif source_in == 'DS':
            self.ds_depths = DepthData()
            self.ds_depths.populate_data(depth_in, source_in, freq_in, draft_in, None)
            self.ds_depths.add_cell_data(self.bt_depths)
            
            
    def set_depth_reference(self, reference):
        '''This function will set the selected depth reference to the
        specified depth reference '''
        
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
            
    def set_valid_data_method(self, setting):
        self.bt_depths.set_valid_data_method(setting)
        
    def composite_depths(self, transect, kargs=None):
        '''Depth compsiting is based on the following assumptions
        
        1. If a depth sounder is available the user must have assumed the ADCP beams
        (BT or vertical) might have problems and it will be the second alternative if 
        not selected as the preferred source
        
        2. For 4-beam BT depths, if 3 beams are valid the average is considered valid.
        It may be based on interpolation of the invalid beam.  However, if only 2 beams
        are valid even though they may be interpolated and included in the average the 
        average will be replaced by an alternative if available.  If no alternative is 
        available the multi-beam average based on available beams and interpolation will
        be used '''
        
        if kargs is None:
            setting = self.composite
        else:
            setting = kargs[0]
            self.composite = setting
            
        #The primary depth reference is the selected reference
        ref = self.selected
        
        if setting == 'On':
            #Prepare vector of valid BT averages, which are defined as
            #having at least 2 valid beams
            
            bt_valid = self.bt_depths.valid_data;
            n_ensembles = bt_valid.shape[-1]
            bt_filtered = self.bt_depths.depth_processed_m
            bt_filtered[bt_filtered == 0] = np.nan
            
            #Prepare vertical beam data, using only data prior to interpolation
            if self.vb_depths is not None:
                vb_filtered = self.vb_depths.depth_processed_m
                vb_filtered[np.squeeze(self.vb_depths.valid_data) == False] = np.nan
            else:
                vb_filtered = repmat([np.nan],n_ensembles,1)

                  
            #Prepare depth sounder data, using only data prior interpolation
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
        '''This function will change the ref_depth draft/  The associated
        depth object will also be updated because DepthData is a handle class.
        The computations are actually done in DepthData as the data are private
        to that class '''
        
        if target == 'ADCP':
            self.bt_depths.change_draft(draft)
            self.vb_depths.change_draft(draft)
        else:
            self.ds_depths.change_draft(draft)    
            
    def depth_filter(self,transect,kargs):
        '''Method to apply filter to all available depth sources, so that
        all sources have the same filter applied '''
        
        if self.bt_depths is not None:
            self.bt_depths.apply_filter(transect, kargs)
        if self.vb_depths is not None:
            self.vb_depths.apply_filter(transect, kargs)
        if self.ds_depths is not None:
            self.ds_depths.apply_filter(transect, kargs)
            
    def depth_interpolation(self, transect, kargs):
        '''Method to apply interpolation to all available depth sources, so
        that all sources have the same filter applied '''
        
        if self.bt_depths is not None:
            self.bt_depths.apply_interpolation(transect, kargs)
        if self.vb_depths is not None:
            self.vb_depths.apply_interpolation(transect, kargs)
        if self.ds_depths is not None:
            self.vb_depths.apply_interpolation(transect, kargs)
            
    def sos_correction(self, ratio):
        '''Correct depths for change in speed of sound'''
        
        #Bottom Track Depths
        if self.bt_depths is not None:
            self.bt_depths.sos_correction(ratio)
            
        #Vertical beam depths
        if self.vb_depths is not None:
            self.vb_depths.sos_correction(ratio)
