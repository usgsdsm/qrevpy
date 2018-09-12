import numpy as np
from MiscLibs.common_functions import cart2pol, pol2cart
import scipy.stats as sp


class NormData(object):
    """Class creates normalized depth and unit discharge or velocity.

    The constuctor method allows an object to be formed without any data.
    populate_data method creates normalized data for a single transect.
    create_composite method creates normalized data for all check transects.
    also allows only a portion of the data to be used in the
    normalization process by specifying the data extent.

    Attributes
    ----------
    file_name: str
        Name of transect file
    cell_depth_normalized: np.array(float)
        Normalized depth of cell
    unit_normalized: np.array(float)
        Normalized discharge or velocity for all depth cells
    unit_normalized_med: np.array(float)
        Median of normalized data within 5% partitions
    unit_normalized_no: np.array(int)
        Number of data points in each median
    unit_normalized_z: np.array(float)
        Relative depth for each median (5% increments)
    unit_normalized_25: np.array(float)
        Value for which 25% of normalized values are smaller
    unit_normalized_75: np.array(float)
        Value for which 75% or normalized values are larger
    data_type: str
        Type of data (v, q, V, or Q)
    data_extent: list
        Defines percent of data from start of transect to use, default [0, 100]
    valid_data: np.array(int)
        Index of median values with point count greater than threshold cutoff
    """
    
    def __init__(self):
        """Creates object and initializes instance variables."""
        self.file_name = None  # Name of transect file
        self.cell_depth_normalized = None  # Normalized depth of cell
        self.unit_normalized = None  # Normalized discharge or velocity for all depth cells
        self.unit_normalized_med = None  # Median of normalized data within 5% partitions
        self.unit_normalized_no = None  # Number of data points in each median
        self.unit_normalized_z = None  # Relative depth for each median (5% increments)
        self.unit_normalized_25 = None  # Value for which 25% of normalized values are smaller
        self.unit_normalized_75 = None  # Value for which 75% or normalized values are larger
        self.data_type = None  # Type of data (v, q, V, or Q)
        self.data_extent = None  # Defines percent of data from start of transect to use, default [0, 100]
        self.valid_data = np.array([])  # Index of median values with point count greater than threshold cutoff
        
    def populate_data(self, transect, data_type, threshold, data_extent=None):
        """Computes the normalized values for a single transect.

        Parameters
        ----------
        transect: TransectData
            Object of TransectData
        data_type: str
            Type of data (v, q, V, or Q)
        threshold: int
            Number of data points in an increment for the increment to be valid.
        data_extent: list
            Defines percent of data from start of transect to use, default [0, 100]
        """

        # If the data extent is not defined set data_extent to zero to trigger all data to be used
        if data_extent is None:
            data_extent = [0, 100]
            
        # Determine number of transects to be processed
        # n_cells = np.nan
        # n_ens = np.nan
        
        filename = transect.file_name
        in_transect_idx = transect.in_transect_idx

        depths_selected = getattr(transect.depths, transect.depths.selected)
        cell_depth = depths_selected.depth_cell_depth_m[:, in_transect_idx]
        cells_above_sl = transect.w_vel.cells_above_sl[:, in_transect_idx]
        cell_depth[cells_above_sl == False] = np.nan
        depth_ens = depths_selected.depth_processed_m[in_transect_idx]

        w_vel_x = transect.w_vel.u_processed_mps[:, in_transect_idx]
        w_vel_y = transect.w_vel.v_processed_mps[:, in_transect_idx]

        invalid_data = np.logical_not(transect.w_vel.valid_data[0, :, in_transect_idx]).T
        w_vel_x[invalid_data] = np.nan
        w_vel_y[invalid_data] = np.nan
        
        boat_select = getattr(transect.boat_vel, transect.boat_vel.selected)
        if boat_select is not None:
            bt_vel_x = boat_select.u_processed_mps[in_transect_idx]
            bt_vel_y = boat_select.v_processed_mps[in_transect_idx]
        else:
            bt_vel_x = np.tile([np.nan], transect.boat_vel.bt_vel.u_processed_mps[in_transect_idx].shape)
            bt_vel_y = np.tile([np.nan], transect.boat_vel.bt_vel.u_processed_mps[in_transect_idx].shape)
            
        # Compute normalized cell depth by average depth in each ensemble
        norm_cell_depth = np.divide(cell_depth, depth_ens)
        norm_cell_depth[norm_cell_depth < 0] = np.nan
        # If data type is discharge compute unit discharge for each cell
        if data_type.lower() == 'q':
            # Compute the cross product for each cell
            unit = np.multiply(w_vel_x, bt_vel_y) - np.multiply(w_vel_y, bt_vel_x)
        else:
            # Compute mean velocity components in each ensemble
            w_vel_mean_1 = np.nanmean(w_vel_x, 0)
            w_vel_mean_2 = np.nanmean(w_vel_y, 0)

            # Compute a unit vector
            direction, _ = cart2pol(w_vel_mean_1, w_vel_mean_2)
            unit_vec_1, unit_vec_2 = pol2cart(direction, 1)
            unit_vec = np.vstack([unit_vec_1, unit_vec_2])
            
            # Compute the velocity magnitude in the direction of the mean velocity of each
            # ensemble using the dot product and unit vector
            unit = np.tile([np.nan], w_vel_x.shape)
            for i in range(w_vel_x.shape[0]):
                unit[i, :] = np.vstack([w_vel_x[i, :], w_vel_y[i, :]]).dot(unit_vec)
                                       
        # Compute Total
        unit_total = np.nansum(np.nansum(unit), 0)
        
        # Adjust to positive value
        if unit_total < 0:
            unit *= -1
            
        # Compute normalize unit values
        unit_norm = np.divide(unit, np.abs(np.nanmean(unit, 0)))
        
        # Apply extents if they have been specified
        if data_extent[0] != 0 or data_extent[1] != 100:
            
            # Unit discharge is computed here because the unit norm could be based on velocity
            unit = np.multiply(w_vel_x, bt_vel_y) - np.multiply(w_vel_y, bt_vel_x)
            unit_ens = np.nansum(unit)
            unit_total = np.nancumsum(unit_ens)
            
            # Adjust so total discharge is positive
            if unit_total[-1] < 0:
                unit_total *= -1
                
            # Apply extents
            unit_lower = unit_total[-1] * data_extent[0] / 100
            unit_upper = unit_total[-1] * data_extent[1] / 100
            idx_extent = np.where(np.logical_and(np.greater(unit_total, unit_lower),
                                                 np.less(unit_total, unit_upper)))[0]
            unit_norm = unit_norm[:, idx_extent]
            norm_cell_depth = norm_cell_depth[:, idx_extent]
            
        # If whole profile is negative make positive
        idx_neg1 = np.tile([np.nan], [unit_norm.shape[1], 1])
        idx_neg2 = np.tile([np.nan], [unit_norm.shape[1], 1])
        for c in range(unit_norm.shape[1]):
            idx_neg1[c] = len(np.where(unit_norm[:, c] < 0)[0])
            idx_neg2[c] = len(np.where(np.isnan(unit_norm[:, c]) == False)[0])
        idx_neg = np.squeeze(idx_neg1) == np.squeeze(idx_neg2)
        unit_norm[:, idx_neg] = unit_norm[:, idx_neg] * -1

        # Store results
        self.file_name = filename
        self.data_extent = data_extent
        self.data_type = data_type
        self.cell_depth_normalized = norm_cell_depth
        self.unit_normalized = unit_norm
        self.compute_stats(threshold)

    def compute_stats(self, threshold):
        """Computes the statistics for the normalized data.

        Parameters
        ----------
        threshold: int
            Number of data points in an increment for the increment to be valid.
        """

        # Set averaging interval
        avg_interval = np.arange(0, 1.05, .05)

        # Intialize variables to nan
        unit_norm_med = np.tile([np.nan], len(avg_interval) - 1)
        unit_norm_med_no = np.tile([np.nan], len(avg_interval) - 1)
        unit_25 = np.tile([np.nan], len(avg_interval) - 1)
        unit_75 = np.tile([np.nan], len(avg_interval) - 1)
        avgz = np.tile([np.nan], len(avg_interval) - 1)

        # Process each normalized increment
        for i in range(len(avg_interval) - 1):
            condition_1 = self.cell_depth_normalized > avg_interval[i]
            condition_2 = self.cell_depth_normalized <= avg_interval[i + 1]
            condition_3 = np.isnan(self.unit_normalized) == False
            condition_all = np.logical_and(np.logical_and(condition_1, condition_2), condition_3)

            unit_25[i], unit_norm_med[i], unit_75[i] = sp.mstats.mquantiles(self.unit_normalized[condition_all],
                                                                            alphap=0.5, betap=0.5)
            unit_norm_med_no[i] = np.sum(np.isnan(self.unit_normalized[condition_all]) == False)
            avgz[i] = 1 - np.nanmean(self.cell_depth_normalized[condition_all])

        # Mark increments invalid if they do not have sufficient data
        cutoff = np.nanmedian(unit_norm_med_no[unit_norm_med_no > 0]) * (threshold / 100)
        self.valid_data = np.where(unit_norm_med_no > cutoff)[0]

        self.unit_normalized_med = unit_norm_med
        self.unit_normalized_no = unit_norm_med_no
        self.unit_normalized_25 = unit_25
        self.unit_normalized_75 = unit_75
        self.unit_normalized_z = avgz

    def create_composite(self, transects, norm_data, threshold):
        """Compute normalized data for measurement composite.

        Parameters
        ----------
        transects: list
            List of objects of TransectData
        norm_data: list
            List of objects of NormData
        threshold: int
            Number of data points in an increment for the increment to be valid.
        """

        # Initialize lists
        n_cells = []
        n_ens = [0]

        # Determine number of cells and ensembles for each transect
        for data in norm_data:
            n_cells.append(data.unit_normalized.shape[0])
            n_ens.append(data.unit_normalized.shape[1])
        max_cells = max(n_cells)
        sum_ens = np.cumsum(n_ens)

        # Initialize normalized variables
        self.unit_normalized = np.tile([np.nan], (max_cells, sum_ens[-1]))
        self.cell_depth_normalized = np.tile([np.nan], (max_cells, sum_ens[-1]))

        # Process each transect using data from only the checked transects
        for n in range(len(transects)):
            if transects[n].checked:
                self.unit_normalized[:n_cells[n], np.arange(sum_ens[n], sum_ens[n + 1])] \
                    = norm_data[n].unit_normalized
                self.cell_depth_normalized[:n_cells[n], np.arange(sum_ens[n], sum_ens[n + 1])] \
                    = norm_data[n].cell_depth_normalized

        # Store data
        self.file_name = 'Measurement'
        self.compute_stats(threshold)
