'''
Created on Sep 26, 2017

@author: gpetrochenkov
'''
import numpy as np
from MiscLibs.convenience import cart2pol, pol2cart

class NormData(object):
    '''
    Class definition for normalized data. The input data are normalized
    for relative depth and relative unit discharge or velocity. The
    constuctor method allows an object to be formed without any data or
    creates a new object from OriginData or completes the properties of
    an object that is already of class NormData. The constuctor method
    also allows only a portion of the data to be used in the
    normalization process by specifying the data extent.
    David S. Mueller, 2/18/2011
    
    Modified 6/18/2011, dsm
    1) Added DisplayNos
    2) Cleaned up variable names and added comments
    
    Modified 10/17/2011, dsm
    3) Added validData index to identify data meeting threshold criteria.
    Changed threshold from a number to a percent of the median number
    of valid cells in a normalized cell. This is applied to each transect
    and to the composite so that the composite will need proportionaly
    more in each normalized cell since it is compiling all transects.
    
    Last modifications / validations 5/15/2012 dsm
    
    Modified 4/12/2013 dsm
    4) Added check for valid velocity when determining which data fit
    within a 5 increment. This affects the avgz value. See line 187.
     
    Modified 3/3/2014 dsm
    5) Modified for use in QRev
    '''
    
    def __init__(self):
        self.__file_name = None #name of transect file
        self.__cell_depth_normalized = None #normalized depth of cell
        self.__unit_normalized = None #normalized discharge or velocity for all depth cells
        self.__unit_normalized_med = None #median of normalized data within 5% partitions
        self.__unit_normalized_no = None #number of data points in each median
        self.__unit_normalized_z = None #relative depth for each median (5% increments)
        self.__unit_normalized_25 = None #value for which 25% of normalized values are smaller
        self.__unit_normalized_75 = None #value for which 75% or normalized values are larger
        self.__data_type = None #type of data (v, q, V, or Q)
        self.__data_extent = None 
        self.__valid_data = None #index of median values with point count greater than threshold cutoff
        
    def populate_data(self, trans_data, data_type, threshold, kargs = None):
        
        #If the data extent is not defined set data_extent to zero to trigger all data to be used
        if kargs is not None:
            data_extent = kargs[0]
        else:
            data_extent = [0, 100]
            
        #Determine number of transects to be proccessed
        n_cells = np.nan
        n_ens = np.nan
        
        filename = trans_data.file_name
        in_transect_idx = trans_data.in_transect_idx
        depth_select = getattr(trans_data.depths, trans_data.depths.selected)
        cell_depth = depth_select.__depth_cell_depth_m[:, in_transect_idx]
        cells_above_sl = trans_data.w_vel.__cells_above_sl
        cell_depth[cells_above_sl == False] = np.nan
        depth_ens = depth_select.__depth_processed_m[in_transect_idx]
        w_vel_x = trans_data.w_vel.__u_processed_mps[:, in_transect_idx]
        w_vel_y = trans_data.w_vel.__v_processed_mps[:, in_transect_idx]
        invalid_data = trans_data.w_vel.__u_processed_mps[:, in_transect_idx] == False
        w_vel_x[invalid_data] = np.nan
        w_vel_y[invalid_data] = np.nan
        n_cells = w_vel_y.shape[0]
        n_ens = w_vel_y.shape[1]
        
        boat_select = getattr(trans_data.boat_vel, trans_data.boat_vel.selected)
        if boat_select is not None:
            bt_vel_x = boat_select.__u_processed_mps[in_transect_idx]
            bt_vel_y = boat_select.__v_processed_mps[in_transect_idx]
        else:
            bt_vel_x = np.tile([np.nan], trans_data.boat_vel.bt_vel.__u_processed_mps[in_transect_idx].shape)
            bt_vel_y = np.tile([np.nan], trans_data.boat_vel.bt_vel.__u_processed_mps[in_transect_idx].shape)
            
        #Compute normalized cell depth by average depth in each ensemble
        norm_cell_depth = np.divide(cell_depth, depth_ens)
        
        #If datatype is discharge compute unit discharge for each cell
        if data_type == 'q':
            #Compute the cross product for each cell
            unit = np.multiply(w_vel_x, bt_vel_y) - np.multiply(w_vel_y, bt_vel_x)
        else:
            #Compute mean velocity components in each ensemble
            w_vel_mean_1 = np.nanmean(w_vel_x)
            w_vel_mean_2 = np.nanmean(w_vel_y)
            
            dir, _ = cart2pol(w_vel_mean_1, w_vel_mean_2)
            unit_vec_1, unit_vec_2 = pol2cart(dir,1)
            unit_vec = np.vstack([unit_vec_1, unit_vec_2])
            
            #Compute the velocity magnitude in the direction of the mean velocity of each
            #ensemble using the dot product
            unit = np.tile([np.nan], w_vel_x.shape)
            for i in range(w_vel_x.shape[0]):
                unit[i,:] = np.vstack([w_vel_x[i,:], w_vel_y[i,:]]).dot(unit_vec)
                                       
        #compute local
        meas = np.nansum(np.nansum(unit))
        
        #Adjust to positive value
        if meas < 0:
            unit *= -1
            
        #Compute normalize unit values
        unit_norm = np.divide(unit, np.abs(np.nanmean(unit)))
        
        #Apply extents if the have been specified
        if data_extent[0] != 0 or data_extent[1] != 100:
            
            #unit discharge is computed here because the unit norm could be based on velocity
            unit = np.multiply(w_vel_x, bt_vel_y) - np.multiply(w_vel_y, bt_vel_x)
            unit_ens = np.nansum(unit)
            unit_total = np.nancumsum(unit_ens)
            
            #Adjust so total discharge is positive
            if unit_total[:-1] < 0:
                unit_total *= -1
                
            #Apply extents
            unit_lower = unit_total[-1] * data_extent[0] / 100
            unit_upper = unit_total[-1] * data_extent[1] / 100
            idx_extent = np.where(unit_total > unit_lower and unit_total < unit_upper)
            unit_norm = unit_norm[:, idx_extent]
            norm_cell_depth = norm_cell_depth[:, idx_extent]
            n_cells, n_ens = unit_norm.shape
            
        #If whole profile is negative make positive
        idx_neg1 = np.tile([np.nan], [unit_norm.shape[1], 1])
        idx_neg2 = np.tile([np.nan], [unit_norm.shape[1], 1])
        for c in range(unit_norm.shape[1]):
            idx_neg1[c] = len(np.where(unit_norm[:, c] < 0))
            idx_neg2[c] = len(np.where(np.isnan(unit_norm[:,c]) == False))
        idx_neg = idx_neg1 == idx_neg2
        unit_norm[:, idx_neg] = unit_norm[:,idx_neg] * -1
        self.__file_name = filename
        self.__cell_depth_normalized = norm_cell_depth
        self.__unit_normalized = unit_norm
        
        #Create arrays of composited data
        max_cells = np.nanmax(n_cells)
        n_ens = np.hstack([0, n_ens])
        sum_ens = np.cumsum(n_ens)
        
        
        
            
        
                
        
        