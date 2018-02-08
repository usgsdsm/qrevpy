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
        self.__w_vel_x = None
        self.__w_vel_y = None
        self.__n_cells = None
        self.__n_ens = None
        self.__threshold = None
        self.__data_type = None
        self.__data_extent = None
        self.__valid_data = None
        
    def populate_data(self, trans_data, data_type, threshold, kargs = None):
        
        #If the data extent is not defined set data_extent to zero to trigger all data to be used
        if kargs is not None:
            self.__data_extent = kargs
        else:
            self.__data_extent = [0, 100]
            
            
        #Determine number of transects to be proccessed
        n_cells = np.nan
        n_ens = np.nan
        
        filename = trans_data.file_name
        in_transect_idx = trans_data.in_transect_idx
        depth_select = getattr(trans_data.depths, trans_data.depths.selected)
        cell_depth = depth_select.depth_cell_depth_m[:, in_transect_idx]
        cells_above_sl = trans_data.w_vel._WaterData__cells_above_sl
        cell_depth[cells_above_sl == False] = np.nan
        depth_ens = depth_select.depth_processed_m[in_transect_idx]
        self.__w_vel_x = trans_data.w_vel._WaterData__u_processed_mps[:, in_transect_idx]
        self.__w_vel_y = trans_data.w_vel._WaterData__v_processed_mps[:, in_transect_idx]
        invalid_data = trans_data.w_vel._WaterData__u_processed_mps[:, in_transect_idx] == False
        self.__w_vel_x[invalid_data] = np.nan
        self.__w_vel_y[invalid_data] = np.nan
        self.__n_cells = self.__w_vel_y.shape[0]
        self.__n_ens = self.__w_vel_y.shape[1]
        self.__threshold = threshold
        self.__data_type = data_type
        
        boat_select = getattr(trans_data.boat_vel, trans_data.boat_vel.selected)
        if boat_select is not None:
            bt_vel_x = boat_select._BoatData__u_processed_mps[in_transect_idx]
            bt_vel_y = boat_select._BoatData__v_processed_mps[in_transect_idx]
        else:
            bt_vel_x = np.tile([np.nan], trans_data.boat_vel.bt_vel._BoatData__u_processed_mps[in_transect_idx].shape)
            bt_vel_y = np.tile([np.nan], trans_data.boat_vel.bt_vel._BoatData__u_processed_mps[in_transect_idx].shape)
            
        #Compute normalized cell depth by average depth in each ensemble
        norm_cell_depth = np.divide(cell_depth, depth_ens)
        
        #If datatype is discharge compute unit discharge for each cell
        if data_type == 'q':
            #Compute the cross product for each cell
            unit = np.multiply(self.__w_vel_x, bt_vel_y) - np.multiply(self.__w_vel_y, bt_vel_x)
        else:
            #Compute mean velocity components in each ensemble
            w_vel_mean_1 = np.nanmean(self.__w_vel_x)
            w_vel_mean_2 = np.nanmean(self.__w_vel_y)
            
            dir, _ = cart2pol(w_vel_mean_1, w_vel_mean_2)
            unit_vec_1, unit_vec_2 = pol2cart(dir,1)
            unit_vec = np.vstack([unit_vec_1, unit_vec_2])
            
            #Compute the velocity magnitude in the direction of the mean velocity of each
            #ensemble using the dot product
            unit = np.tile([np.nan], self.__w_vel_x.shape)
            for i in range(self.__w_vel_x.shape[0]):
                unit[i,:] = np.vstack([self.__w_vel_x[i,:], self.__w_vel_y[i,:]]).dot(unit_vec)
                                       
        #compute local
        meas = np.nansum(np.nansum(unit))
        
        #Adjust to positive value
        if meas < 0:
            unit *= -1
            
        #Compute normalize unit values
        unit_norm = np.divide(unit, np.abs(np.nanmean(unit)))
        
        #Apply extents if the have been specified
        if self.__data_extent[0] != 0 or self.__data_extent[1] != 100:
            
            #unit discharge is computed here because the unit norm could be based on velocity
            unit = np.multiply(self.__w_vel_x, bt_vel_y) - np.multiply(self.__w_vel_y, bt_vel_x)
            unit_ens = np.nansum(unit)
            unit_total = np.nancumsum(unit_ens)
            
            #Adjust so total discharge is positive
            if unit_total[:-1] < 0:
                unit_total *= -1
                
            #Apply extents
            unit_lower = unit_total[-1] * self.__data_extent[0] / 100
            unit_upper = unit_total[-1] * self.__data_extent[1] / 100
            idx_extent = np.where(unit_total > unit_lower and unit_total < unit_upper)
            unit_norm = unit_norm[:, idx_extent]
            norm_cell_depth = norm_cell_depth[:, idx_extent]
            self.__n_cells, self.__n_ens = unit_norm.shape[0], unit_norm.shape[1]
            
        #If whole profile is negative make positive
        idx_neg1 = np.tile([np.nan], unit_norm.shape[1])
        idx_neg2 = np.tile([np.nan], unit_norm.shape[1])
        for c in range(unit_norm.shape[1]):
            idx_neg1[c] = len(np.where(unit_norm[:, c] < 0))
            idx_neg2[c] = len(np.where(np.isnan(unit_norm[:,c]) == False))
        idx_neg = idx_neg1 == idx_neg2
        unit_norm[:, idx_neg] = unit_norm[:,idx_neg] * -1
        self.__file_name = filename
        self.__cell_depth_normalized = norm_cell_depth
        self.__unit_normalized = unit_norm
        
    def get_composite_data(self, transects, norm_data):
        
        ##Create object for measurement composite
        #Create arrays of composited data
        if self.__threshold is None:
            self.__threshold = norm_data[0].__threshold
        
        maxcells= np.nanmax([x.__n_cells for x in norm_data])
        n_ens= np.hstack([[0],[x.__n_ens for x in norm_data]])
        n_cells = [x.__n_cells for x in norm_data]
        sum_ens= np.cumsum(n_ens);
        self.__unit_normalized = np.tile([np.nan], (maxcells,sum_ens[-1]))
        self.__cell_depth_normalized = np.tile([np.nan], (maxcells,sum_ens[-1]))
        
        for n in range(len(transects)):
            if transects[n].checked == 1:
                self.__unit_normalized[:n_cells[n],np.arange(sum_ens[n],sum_ens[n+1])] = norm_data[n].__unit_normalized
                self.__cell_depth_normalized[:n_cells[n],np.arange(sum_ens[n],sum_ens[n+1])] = norm_data[n].__cell_depth_normalized
            
        self.__file_Name= 'Measurement'
        
        norm_data.append(self)
        
        for n in range(len(transects) + 1):
            
            avg_interval = np.arange(0,1,.05)
            unit_norm_med = np.tile([np.nan], len(avg_interval))
            unit_norm_med_no = np.tile([np.nan], len(avg_interval))
            unit_25 = np.tile([np.nan], len(avg_interval))
            unit_75 = np.tile([np.nan], len(avg_interval))
            avgz = np.tile([np.nan], len(avg_interval))
            
            #Process each normalized increment
            for i in range(len(avg_interval) - 1):
                idx = np.where((norm_data[n].__cell_depth_normalized > avg_interval[i]) \
                               & (norm_data[n].__cell_depth_normalized <= avg_interval[i+1]) \
                               & (np.isnan(norm_data[n].__unit_normalized) == False))
                unit_norm_med[i] = np.nanmedian(norm_data[n].__unit_normalized[idx])
                unit_norm_med_no[i] = np.sum(np.isnan(norm_data[n].__unit_normalized[idx]) == False)
                unit_25[i] = np.nanpercentile(norm_data[n].__unit_normalized, [25])
                unit_75[i] = np.nanpercentile(norm_data[n].__unit_normalized, [75])
                avgz[i] = np.nanmean(norm_data[n].__cell_depth_normalized[idx])
                
            #Mark increments invalid if they do not have sufficient data
            cutoff = np.nanmedian(unit_norm_med_no[unit_norm_med_no>0]) * (self.__threshold / 1000)
            valid = np.where(unit_norm_med_no > cutoff)
            
            norm_data[n].__unit_normalized_med = unit_norm_med
            norm_data[n].__unit_normalized_no = unit_norm_med_no
            norm_data[n].__unit_normalized_25 = unit_25
            norm_data[n].__unit_normalized_75 = unit_75
            norm_data[n].__unit_normalized_z = avgz
            
            
            if norm_data[n].__data_type is None:
                norm_data[n].__data_type = norm_data[0].__data_type
            if norm_data[n].__data_extent is None:
                norm_data[n].__data_extent = norm_data[0].__data_extent
            norm_data[n].__valid_data = valid
        
        
            
        
                
        
        