import numpy as np
from MiscLibs.abba_2d_interpolation import find_neighbors
from MiscLibs.abba_2d_interpolation import compute_distances
from MiscLibs.abba_2d_interpolation import idw_interpolation
from MiscLibs.abba_2d_interpolation import abba_idw_interpolation


def create_data_uniform_cell_size():
    """ Creates all data necessary to test interpolation for uniform cell sizes."""

    valid_data = np.array([[True, True, False, True, False],
                          [False, True, False, False, True],
                          [True, True, True, True, False],
                          [True, False, True, False, False],
                          [False, True, True, False, False],
                          [False, False, False, False, False]])

    cells_above_sl = np.array([[True, True, True, True, True],
                               [True, True, True, True, True],
                               [True, True, True, True, True],
                               [True, True, True, True, True],
                               [False, True, False, True, False],
                               [False, False, False, False, False]])

    y_cell_centers = np.array([[0.1, 0.1, 0.1, 0.1, 0.1],
                               [0.2, 0.2, 0.2, 0.2, 0.2],
                               [0.3, 0.3, 0.3, 0.3, 0.3],
                               [0.4, 0.4, 0.4, 0.4, 0.4],
                               [0.5, 0.5, 0.5, 0.5, 0.5],
                               [0.6, 0.6, 0.6, 0.6, 0.6]])

    y_cell_size = np.array([[0.1, 0.1, 0.1, 0.1, 0.1],
                            [0.1, 0.1, 0.1, 0.1, 0.1],
                            [0.1, 0.1, 0.1, 0.1, 0.1],
                            [0.1, 0.1, 0.1, 0.1, 0.1],
                            [0.1, 0.1, 0.1, 0.1, 0.1],
                            [0.1, 0.1, 0.1, 0.1, 0.1]])

    y_normalize = np.array([0.8, 0.8, 0.8, 0.8, 0.8])

    x = [1, 2, 4, 6, 7]

    data = np.array([[0.1, 0.15, 0.18, 0.12, 0.1],
                     [0.18, 0.22, 0.2, 0.12, 0.12],
                     [0.23, 0.28, 0.3, 0.23, 0.18],
                     [0.34, 0.35, 0.41, 0.34, 0.24],
                     [0.45, 0.4, 0.52, 0.45, 0.35],
                     [0.6, 0.6, 0.6, 0.6, 0.6]])

    return valid_data, cells_above_sl, y_cell_centers, y_cell_size, y_normalize, x, data


def test_neighbors_uniform():
    """ Test finding neighbors for various scenarios using uniform data."""

    # Initialize data
    valid_data, cells_above_sl, y_cell_centers, y_cell_size, y_normalize, x, data = create_data_uniform_cell_size()

    # Test uniform cell size
    result = find_neighbors(valid_data=valid_data,
                            cells_above_sl=cells_above_sl,
                            y_cell_centers=y_cell_centers,
                            y_cell_size=y_cell_size,
                            y_normalize=y_normalize)
    assert result[0]['target'] == (0, 2) and result[0]['neighbors'] == [(2, 2), (0, 1), (0, 3)]
    assert result[1]['target'] == (0, 4) and result[1]['neighbors'] == [(1, 4), (0, 3)]
    assert result[2]['target'] == (1, 0) and result[2]['neighbors'] == [(0, 0), (2, 0), (1, 1)]
    assert result[3]['target'] == (1, 2) and result[3]['neighbors'] == [(2, 2), (1, 1), (1, 4)]
    assert result[4]['target'] == (1, 3) and result[4]['neighbors'] == [(0, 3), (2, 3), (1, 1), (1, 4)]
    assert result[5]['target'] == (2, 4) and result[5]['neighbors'] == [(1, 4), (2, 3)]
    assert result[6]['target'] == (3, 1) and result[6]['neighbors'] == [(2, 1), (4, 1), (3, 0), (3, 2)]
    assert result[7]['target'] == (3, 3) and result[7]['neighbors'] == [(2, 3), (3, 2)]
    assert result[8]['target'] == (3, 4) and result[8]['neighbors'] == [(1, 4), (3, 2)]
    assert result[9]['target'] == (4, 3) and result[9]['neighbors'] == [(2, 3), (4, 2)]


def create_nonuniform_data():
    """ Created all data necessary to test interpolation for data with nonuniform cell sizes."""

    # Initialize data
    y_cell_size = np.array([[0.1, 0.05, 0.1],
                            [0.1, 0.05, 0.1],
                            [0.1, 0.05, 0.1],
                            [0.1, 0.05, 0.1],
                            [0.1, 0.05, 0.1],
                            [0.1, 0.05, 0.1],
                            [0.1, 0.05, 0.1]])

    y_cell_centers = np.array([[0.1, 0.075, 0.1],
                               [0.2, 0.125, 0.2],
                               [0.3, 0.175, 0.3],
                               [0.4, 0.225, 0.4],
                               [0.5, 0.275, 0.5],
                               [0.6, 0.325, 0.6],
                               [0.7, 0.375, 0.7]])

    y_normalize = np.array([0.8, 0.8, 0.8])

    valid_data = np.array([[True, True, True],
                           [True, True, True],
                           [True, True, True],
                           [True, True, True],
                           [True, True, True],
                           [True, True, True],
                           [True, True, True]])

    cells_above_sl = np.array([[True, True, True],
                               [True, True, True],
                               [True, True, True],
                               [True, True, True],
                               [True, True, True],
                               [False, True, False],
                               [False, True, False]])

    x = [1, 2, 5]

    data = np.array([[0.1, 0.15, 0.18],
                     [0.18, 0.22, 0.2],
                     [0.23, 0.28, 0.3],
                     [0.34, 0.35, 0.41],
                     [0.45, 0.4, 0.52],
                     [0.6, 0.6, 0.6],
                     [0.3, 0.3, 0.3]])

    return y_cell_size, y_cell_centers, y_normalize, valid_data, cells_above_sl, x, data


def test_neighbors_small_cell_invalid():
    """ Test finding neighbors for small cell size that is bordered with large cell sizes."""
    y_cell_size, y_cell_centers, y_normalize, valid_data, cells_above_sl, x, data = create_nonuniform_data()
    valid_data[2, 1] = False
    result = find_neighbors(valid_data, cells_above_sl, y_cell_centers, y_cell_size, y_normalize)
    assert result[0]['target'] == (2, 1) and result[0]['neighbors'] == [(1, 1), (3, 1), (1, 0), (1, 2)]


def test_neighbors_large_cell_invalid():
    """ Test finding neighbors for large cell size that is bordered with small cell sizes."""
    y_cell_size, y_cell_centers, y_normalize, valid_data, cells_above_sl, x, data = create_nonuniform_data()
    valid_data[1, 0] = False
    result = find_neighbors(valid_data, cells_above_sl, y_cell_centers, y_cell_size, y_normalize)
    assert result[0]['target'] == (1, 0) and result[0]['neighbors'] == [(0, 0), (2, 0), (2, 1), (3, 1)]


def test_compute_distance():
    """ Test computing distances with varying x and y components."""
    y_cell_size, y_cell_centers, y_normalize, valid_data, cells_above_sl, x, data = create_nonuniform_data()
    target = (2, 1)
    neighbors = [(1, 1), (3, 1), (1, 0), (1, 2)]
    distances = compute_distances(target=target, neighbors=neighbors, x=x, y=y_cell_centers)
    assert ['%.4f' % dist for dist in distances] == ['0.0500', '0.0500', '1.0003', '3.0001']


def test_interpolation():
    """ Test idw interpolation."""
    y_cell_size, y_cell_centers, y_normalize, valid_data, cells_above_sl, x, data = create_nonuniform_data()
    target = (2, 1)
    neighbors = [(1, 1), (3, 1), (1, 0), (1, 2)]
    distances = compute_distances(target=target, neighbors=neighbors, x=x, y=y_cell_centers)
    interpolated_value = idw_interpolation(data=data, neighbor_indices=neighbors, distances=distances)
    assert '%.4f' % interpolated_value == '0.1972'


def test_interpolate_list_of_data():
    y_cell_size, y_cell_centers, y_normalize, valid_data, cells_above_sl, x, data1 = create_nonuniform_data()
    valid_data[2, 1] = False
    valid_data[0, 2] = False
    data2 = np.array([[1.1, 1.15, 1.18],
                     [1.18, 1.22, 1.2],
                     [1.23, 1.28, 1.3],
                     [1.34, 1.35, 1.41],
                     [1.45, 1.4, 1.52],
                     [1.6, 1.6, 1.6],
                     [1.3, 1.3, 1.3]])
    interpolated_values = abba_idw_interpolation(data_list=[data1, data2],
                                                 valid_data=valid_data,
                                                 cells_above_sl=cells_above_sl,
                                                 y_centers=y_cell_centers,
                                                 y_cell_size=y_cell_size,
                                                 y_normalize=y_normalize,
                                                 x_shiptrack=x)

    # First data set
    assert interpolated_values[0][0][0] == (0, 2)
    assert '%.4f' % interpolated_values[0][0][1] == '0.1852'
    assert interpolated_values[0][1][0] == (2, 1)
    assert '%.4f' % interpolated_values[0][1][1] == '0.1972'

    # Second data set
    assert interpolated_values[1][0][0] == (0, 2)
    assert '%.4f' % interpolated_values[1][0][1] == '1.1852'
    assert interpolated_values[1][1][0] == (2, 1)
    assert '%.4f' % interpolated_values[1][1][1] == '1.1972'

