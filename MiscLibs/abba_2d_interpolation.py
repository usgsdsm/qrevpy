"""abba_2d_interpolation

This module performs 2-D interpolation on data that is assumed to be arranged in row-column format rather
than in a random pattern. The rows represent vertical location or y-coordinate of each cell
in the data array. The columns represent a horizontal location or x-coordinate of the data.
The cell size and thus the y-coordinate of a cell can change for cell to cell or ensemble to ensemble.
The interpolation algorithm searches for the nearest valid cell above, below, before, and after
the cell to be interpolated. All cells before and after that lie completely within the target cell are
determined to be neighbors. If the target cell lies completely with in a larger cell before or after the
target that larger cell is consider a neighbor.

The methods provide the flexibility to determine neighbors based on either a raw vertical location
or a normalized location. To use a normalized location the y_normalize variable must contain the
divisor for each x location of the raw vertical locations (each column of y-coordinate).

For efficiency the data_list can contain multiple types of data that lie on the same x-y locations.
This allows multiple interpolations without having to recompute neighbors and distances.

Example
-------

For interpolating water velocities collected by an ADCP

interpolated_data = abba_idw_interpolation(data_list = [u_processed_mps, v_processed_mps]
                                           valid_data = valid_data
                                           cells_above_sl = cells_above_sl
                                           y_centers = depth_cell_depth_m
                                           y_cell_size = depth_cell_size_m
                                           y_normalize = depth_processed_m
                                           x_shiptrack = distance_along_shiptrack)

interpolated_u_values = interpolated_data[0]
interpolated_v_values = interpolated_data[1]

"""
import numpy as np


def find_neighbors(valid_data, cells_above_sl, y_cell_centers, y_cell_size, y_normalize=[]):
    """ Finds the nearest valid cells above, below, before, and after each invalid cell. The before and after
    Cells must have data in the same y range as the invalid cell.

    Parameters
    ----------
    valid_data: np.array(logical)
        Logical array indicating whether each cell is valid (true) or invalid (false)
    cells_above_sl: np.array(logical)
        Logical array indicating which cells are above the side lobe cutoff (true)
    y_cell_centers: np.array(float)
        Y coordinate corresponding to the center of the data cells
    y_cell_size: np.array(float)
        Size of each cell in the y-direction
    y_normalize: np.array(float)
        1-D array containing values that will be used to normalize the data for identifying neighbors

    Returns
    -------
    neighbors: list
        List of dictionaries providing the indices of the above, below, before, and after valid cells.
    """

    # Compute cell extents
    y_top = y_cell_centers - 0.5 * y_cell_size
    y_bottom = y_cell_centers + 0.5 * y_cell_size
    y_centers = y_cell_centers
    if len(y_normalize) > 0:
        y_top = np.round(y_top / y_normalize, 3)
        y_bottom = np.round(y_bottom / y_normalize, 3)
        y_centers = np.round(y_cell_centers / y_normalize, 3)

    # ID cells above side lobe with invalid data
    valid_data_float = valid_data.astype(float)
    valid_data_float[np.logical_not(cells_above_sl)] = np.nan
    invalid_cell_index = np.where(valid_data_float == 0)

    # Initialize output list
    neighbors = []

    # Process each index
    for cell, ens in zip(invalid_cell_index[0], invalid_cell_index[1]):
        points = []
        target = (cell, ens)

        # Identify indices of cells above and below target
        above = find_above(target, valid_data)
        if above is not None:
            points.append(above)

        below = find_below(target, valid_data)
        if below is not None:
            points.append(below)

        # Find valid cell indices for all ensembles with valid data within the target cell
        # Find cell centers that are withing the range of the target cell
        y_match1 = np.array((np.logical_and(y_centers >= y_top[target], y_centers <= y_bottom[target])))
        # Find all cells that the target falls within their range
        y_match2 = np.array((np.logical_and(y_top[target] >= y_top, y_bottom[target] <= y_bottom)))
        # Combine matches
        y_match = np.logical_or(y_match1, y_match2)
        y_match = np.logical_and(y_match, valid_data)

        # Identify indices of cells before and after target
        before = find_before(target, y_match)
        if before:
            points = points + before
        after = find_after(target, y_match)
        if after:
            points = points + after

        neighbors.append({'target': target, 'neighbors': points})

    return neighbors


def find_above(target, valid_data):
    """ Finds the nearest valid cell above the target.

    Parameters
    ----------
    target: tuple
        Indices of target cell
    valid_data: np.array(logical)

    Returns
    -------
    above_idx: tuple
        Indices of valid cell immediately above target
    """

    # Initialize cell counter
    above_idx = target[0] - 1

    # Find nearest valid cell above target
    while above_idx >= 0 and not valid_data[above_idx, target[1]]:
        above_idx = above_idx - 1
    if above_idx >= 0:
        above_idx = (above_idx, target[1])
    else:
        above_idx = None

    return above_idx


def find_below(target, valid_data):
    """ Finds the nearest valid cell below the target.

    Parameters
    ----------
    target: tuple
        Indices of target cell
    valid_data: np.array(logical)

    Returns
    -------
    below_idx: tuple
        Indices of valid cell immediately below target
    """

    # Initialize cell counter
    below_idx = target[0] + 1

    # Determine cell row index limit
    n_cells = len(valid_data[:, target[1]])-1

    # Find nearest valid cell below target
    while below_idx <= n_cells and not valid_data[below_idx, target[1]]:
        below_idx = below_idx + 1
    if below_idx <= n_cells:
        below_idx = (below_idx, target[1])
    else:
        below_idx = None

    return below_idx


def find_before(target, y_match):
    """ Finds the nearest ensemble before the target that has valid cells within the vertical range of the target

    Parameters
    ----------
    target: tuple
        Indices of target cell
    y_match: np.array(logical)
        2-D array of all cells that are within the vertical range of the target cell

    Returns
    -------
    before_idx: list
        List of tuples of indices of all cells in the nearest ensemble before that target that are within
        the vertical range of the target cell
    """

    # Initialize ensemble counter
    before_ens = target[1] - 1

    # Loop until an ensemble is found that has valid data within the vertical range of the target
    while before_ens > 0 and not np.any(y_match[:, before_ens]):
        before_ens = before_ens - 1

    # Find and store the indices all cells from the identified ensemble
    # that are within the vertical range of the target
    if before_ens >= 0:
        rows = np.where(y_match[:, before_ens])[0]
        before_idx = []
        for row in rows:
            before_idx.append((row, before_ens))
    else:
        before_idx = []

    return before_idx


def find_after(target, y_match):
    """ Finds the nearest ensemble after the target that has valid cells within the vertical range of the target

    Parameters
    ----------
    target: tuple
        Indices of target cell
    y_match: np.array(logical)
        2-D array of all cells that are within the vertical range of the target cell

    Returns
    -------
    after_idx: list
        List of tuples of indices of all cells in the nearest ensemble after that target that are within
        the vertical range of the target cell
    """

    # Initialize ensemble counter
    after_ens = target[1] + 1

    # Loop until an ensemble is found that has valid data within the vertical range of the target
    while after_ens <= y_match.shape[1]-1 and not np.any(y_match[:, after_ens]):
        after_ens = after_ens + 1

    # Find and store the indices all cells from the identified ensemble
    # that are within the vertical range of the target
    if after_ens <= y_match.shape[1]-1:
        rows = np.where(y_match[:, after_ens])[0]
        after_idx = []
        for row in rows:
            after_idx.append((row, after_ens))
    else:
        after_idx = []

    return after_idx


def compute_distances(target, neighbors, x, y):
    """ Computes distances between the target and neighbors.

    Parameters
    ----------
    target: tuple
        Indices of target cell
    neighbors: list
        List of indices of target's neighboring cells
    x: np.array(float)
        1-D array of distances between ensembles
    y: np.array(float)
        2-D array of vertical distances of cells for each ensemble

    Returns
    -------
    distances: list
        List of distances from target to each neighbor
    """

    # Intialize target location
    target_y = y[target]
    target_x = x[target[1]]

    # Compute distance from target cell to each neighbor
    distances = []
    for neighbor in neighbors:
        distances.append(np.sqrt((y[neighbor] - target_y) ** 2 + (x[neighbor[1]] - target_x) ** 2))

    return distances


def idw_interpolation(data, neighbor_indices, distances):
    """ Interpolate value using neighbors and inverse distance weighting.

    Parameters
    ----------
    data: np.array(float)
        2-D array containing data to interpolate
    neighbor_indices: list
        List of tuples defining the indices of the target's neighbors
    distances: list
        List of distances from target to each neighbor

    Returns
    -------
    interpolated_value: float
        Value of target cell interpolated from neighbors
    """

    # Compute sum of weights
    sum_of_weights = sum(distances)

    # Compute weighted sum or neighbor values
    weighted_sum = 0
    for n, index in enumerate(neighbor_indices):
        weighted_sum = weighted_sum + data[index] * distances[n]

    # Compute interpolated value
    interpolated_value = weighted_sum / sum_of_weights

    return interpolated_value


def abba_idw_interpolation(data_list, valid_data, cells_above_sl, y_centers, y_cell_size, y_normalize, x_shiptrack):
    """ Interpolates values for invalid cells using the neighboring cells above, below, before, and after and
    and inverse distance averaging.

    Parameters
    ----------
    data_list: list
        List of np.array(float) data to used for interpolation
    valid_data: np.array(logical)
        Logical array of valid data
    cells_above_sl: np.array(logical)
        Logical array of all valid cells above the side lobe cutoff
    y_centers: np.array(float)
        Y coordinate corresponding to the center of the data cells
    y_cell_size: np.array(float)
        Size of each cell in the y-direction
    y_normalize: np.array(float)
        1-D array containing values that will be used to normalize the data for identifying neighbors
    x_shiptrack: np.array(float)
        X coordinate of cumulative shiptrack

    Returns
    -------
    interp_data: list
        Indices and interpolation values for invalid cells corresponding to data list.
    """

    # Initialize output list
    interpolated_data = [[] for _ in range(len(data_list))]

    # Find neighbors associated with each target
    interpolation_points = find_neighbors(valid_data=valid_data,
                                          cells_above_sl=cells_above_sl,
                                          y_cell_centers=y_centers,
                                          y_cell_size=y_cell_size,
                                          y_normalize=y_normalize)

    # Process each target
    for point in interpolation_points:
        # Compute distance from target to neighbors
        distances = compute_distances(target=point['target'],
                                      neighbors=point['neighbors'],
                                      x=x_shiptrack,
                                      y=y_centers)

        # Interpolate target for each data set in data_list
        for n, data in enumerate(data_list):
            interpolated_value = idw_interpolation(data=data,
                                                   neighbor_indices=point['neighbors'],
                                                   distances=distances)
            interpolated_data[n].append([point['target'], interpolated_value])

    return interpolated_data