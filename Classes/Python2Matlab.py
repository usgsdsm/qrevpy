import numpy as np


class Python2Matlab(object):
    """Converts python meas class to Matlab structure.

    Attributes
    ----------
    matlab_dict: dict
        Dictionary of Matlab structures
    """

    def __init__(self, meas):
        py_2_mat_dict = self.create_py_2_mat_dict()
        self.matlab_dict = dict()

        self.matlab_dict['stationName'] = meas.station_name
        self.matlab_dict['stationNumber'] = meas.station_number
        self.matlab_dict['processing'] = meas.processing
        self.matlab_dict['extTempChk'] = meas.ext_temp_chk
        self.matlab_dict['initialSettings'] = meas.initial_settings
        # Comments really needs to convert to cell
        self.matlab_dict['comments'] = self.comment2struct(meas.comments)
        self.matlab_dict['compassCal'] = self.listobj2struct(meas.compass_cal, py_2_mat_dict)
        self.matlab_dict['compassEval'] = self.listobj2struct(meas.compass_eval, py_2_mat_dict)
        self.matlab_dict['discharge'] = self.listobj2struct(meas.discharge, py_2_mat_dict)
        self.matlab_dict['sysTest'] = self.listobj2struct(meas.system_test, py_2_mat_dict)
        self.matlab_dict['transects'] = self.listobj2struct(meas.transects, py_2_mat_dict)
        self.matlab_dict['extrapFit'] = self.listobj2struct([meas.extrap_fit], py_2_mat_dict)
        self.matlab_dict['mbTests'] = self.listobj2struct(meas.mb_tests, py_2_mat_dict)

    @staticmethod
    def listobj2struct(list_in, new_key_dict=None):
        """Converts a list of objects to a structured array.

        Parameters
        ----------
        list_in: list
            List of objects
        new_key_dict: dict
            Dictionary to translate python variable names to Matlab variable names

        Returns
        -------
        struct: np.array
            Structured array
        """

        # Create data type for each variable in object
        keys = list(vars(list_in[0]).keys())
        data_type = []
        for key in keys:
            if new_key_dict is not None and key in new_key_dict:
                data_type.append((new_key_dict[key], list))
            else:
                data_type.append((key, list))

        # Create structured array based on data type and length of list
        dt = np.dtype(data_type)
        struct = np.zeros((len(list_in),), dt)

        # Populate the structure with data from the objects
        for n, item in enumerate(list_in):
            if type(item) is list:
                struct = Python2Matlab.listobj2struct(item, new_key_dict)
            else:
                new_dict = Python2Matlab.obj2dict(item, new_key_dict)
                for key in new_dict:
                    if new_key_dict is not None and key in new_key_dict:
                        struct[new_key_dict[key]][n] = new_dict[key]
                    else:
                        struct[key][n] = new_dict[key]

        return struct

    @staticmethod
    def change_dict_keys(dict_in, new_key_dict):
        """Recursively changes the name of dictionary keys.

        dict_in: dict
            Dictionary with keys that need a name change
        new_key_dict: dict
            Dictionary to cross reference existing key to new key names
        """

        dict_out = dict()

        for key in dict_in:
            if type(dict_in[key]) is dict:
                dict_in[key] = Python2Matlab.change_dict_keys(dict_in[key], new_key_dict)

            if new_key_dict is not None and key in new_key_dict:
                dict_out[new_key_dict[key]] = dict_in[key]
            else:
                dict_out[key] = dict_in[key]

        return dict_out

    @staticmethod
    def obj2dict(obj, new_key_dict=None):
        """Converts object variables to dictionaries. Works recursively to all levels of objects.

        Parameters
        ----------
        obj: object
            Object of some class
        new_key_dict: dict
            Dictionary to translate python variable names to Matlab variable names

        Returns
        -------
        obj_dict: dict
            Dictionary of all object variables
        """
        obj_dict = vars(obj)
        new_dict = dict()
        for key in obj_dict:

            # If variable is another object convert to dictionary recursively
            if str(type(obj_dict[key]))[8:13] == 'Class':
                obj_dict[key] = Python2Matlab.obj2dict(obj_dict[key], new_key_dict)

            # If variable is a list of objects convert to dictionary
            elif type(obj_dict[key]) is list and len(obj_dict[key]) > 0 \
                    and str(type(obj_dict[key][0]))[8:13] == 'Class':
                obj_dict[key] = Python2Matlab.listobj2dict(obj_dict[key], new_key_dict)

            elif type(obj_dict[key]) is dict:
                obj_dict[key] = Python2Matlab.change_dict_keys(obj_dict[key], new_key_dict)

            # If variable is None rename as necessary and convert None to empty list
            if obj_dict[key] is None:
                if new_key_dict is not None and key in new_key_dict:
                    new_dict[new_key_dict[key]] = []
                else:
                    new_dict[key] = []
            # If varialbe is not None rename as necessary
            elif new_key_dict is not None and key in new_key_dict:
                new_dict[new_key_dict[key]] = obj_dict[key]
            else:
                new_dict[key] = obj_dict[key]

        return new_dict

    @staticmethod
    def comment2struct(comments):
        data_type = [('comment', list)]
        dt = np.dtype(data_type)
        struct = np.zeros((len(comments),), dt)
        for n, line in enumerate(comments):
            struct['comment'][n] = line
        return struct

    @staticmethod
    def listobj2dict(list_in, new_key_dict=None):
        """Converts list of objects to list of dictionaries. Works recursively to all levels of objects.

        Parameters
        ----------
        list_in: list
            List of objects of some class
        new_key_dict: dict
            Dictionary to translate python variable names to Matlab variable names

        Returns
        -------
        new_list: list
            List of dictionaries
        """
        new_list = []
        for obj in list_in:
            new_list.append(Python2Matlab.obj2dict(obj, new_key_dict))
        return new_list

    @staticmethod
    def create_py_2_mat_dict():
        """Creates a dictionary to cross reference Python names with Matlab names

        Returns
        -------
        py_2_mat_dict: dict
            Dictionary of python key to Matlab variable
        """

        py_2_mat_dict = {'Python': 'Matlab',
                         'align_correction_deg': 'alignCorrection_deg',
                         'altitude_ens_m': 'altitudeEns_m',
                         'avg_method': 'avgMethod',
                         'beam_angle_deg': 'beamAngle_deg',
                         'beam_filter': 'beamFilter',
                         'beam_pattern': 'beamPattern',
                         'blanking_distance_m': 'blankingDistance_m',
                         'boat_vel': 'boatVel',
                         'bot_diff': 'botdiff',
                         'bot_method': 'botMethod',
                         'bot_method_auto': 'botMethodAuto',
                         'bot_method_orig': 'botMethodOrig',
                         'bot_r2': 'botrsqr',
                         'bottom_ens': 'bottomEns',
                         'bottom_mode': 'bottomMode',
                         'bt_depths': 'btDepths',
                         'bt_vel': 'btVel',
                         'cell_depth_normalized': 'cellDepthNormalized',
                         'cells_above_sl': 'cellsAboveSL',
                         'cells_above_sl_bt': 'cellsAboveSL_bt',
                         'compass_cal': 'compassCal',
                         'compass_diff_deg': 'compassDiff_deg',
                         'compass_eval': 'compassEval',
                         'configuration_commands': 'configurationCommands',
                         'coord_sys': 'coordSys',
                         'corr_table': 'corrTable',
                         'correction_factor': 'correctionFactor',
                         'cov_95': 'cov95',
                         'cov_user': 'cov95User',
                         'cust_coef': 'custCoef',
                         'd_filter': 'dFilter',
                         'd_filter_threshold': 'dFilterThreshold',
                         'data_extent': 'dataExtent',
                         'data_type': 'dataType',
                         'date_time': 'dateTime',
                         'depth_beams_m': 'depthBeams_m',
                         'depth_cell_depth_m': 'depthCellDepth_m',
                         'depth_cell_depth_orig_m': 'depthCellDepthOrig_m',
                         'depth_cell_size_m': 'depthCellSize_m',
                         'depth_cell_size_orig_m': 'depthCellSizeOrig_m',
                         'depth_depth_m': 'depthCellDepth_m',
                         'depth_ens': 'depthSourceEns',
                         'depth_hz': 'depthFreq_Hz',
                         'depth_index': 'depthInvalidIndex',
                         'depth_orig_m': 'depthOrig_m',
                         'depth_processed_m': 'depthProcessed_m',
                         'depth_source': 'depthSource',
                         'depths': ' depths',
                         'diff_qual_ens': 'diffQualEns',
                         'dist_us_m': 'distUS_m',
                         'distance_m': 'dist_m',
                         'draft_orig_m': 'draftOrig_m',
                         'draft_use_m': 'draftUse_m',
                         'ds_depths': 'dsDepths',
                         'edges_95': 'edges95',
                         'edges_user': 'edges95User',
                         'end_serial_time': 'endSerialTime',
                         'ens_duration_sec': 'ensDuration_sec',
                         'excluded_dist_m': 'excludedDist',
                         'exp_method': 'expMethod',
                         'exponent_95_ci': 'exponent95confint',
                         'exponent_auto': 'exponentAuto',
                         'exponent_orig': 'exponentOrig',
                         'ext_gga_altitude_m': 'extGGAAltitude_m',
                         'ext_gga_differential': 'extGGADifferential',
                         'ext_gga_hdop': 'extGGAHDOP',
                         'ext_gga_lat_deg': 'extGGALat_deg',
                         'ext_gga_lon_deg': 'extGGALon_deg',
                         'ext_gga_num_sats': 'extGGANumSats',
                         'ext_gga_serial_time': 'extGGASerialTime',
                         'ext_gga_utc': 'extGGAUTC',
                         'ext_temp_chk': 'extTempChk',
                         'ext_vtg_course_deg': 'extVTGCourse_deg',
                         'ext_vtg_speed_mps': 'extVTGSpeed_mps',
                         'extrap_fit': 'extrapFit',
                         'extrapolation_95': 'extrapolation95',
                         'extrapolation_user': 'extrapolation95User',
                         'file_name': 'fileName',
                         'filter_type': 'filterType',
                         'fit_method': 'fitMethod',
                         'fit_r2': 'fitrsqr',
                         'flow_dir_deg': 'flowDir_deg',
                         'flow_spd_mps': 'flowSpd_mps',
                         'frequency_khz': 'frequency_hz',
                         'gga_lat_ens_deg': 'ggaLatEns_deg',
                         'gga_lon_ens_deg': 'ggaLonEns_deg',
                         'gga_position_method': 'ggaPositionMethod',
                         'gga_serial_time_ens': 'ggaSerialTimeEns',
                         'gga_vel': 'ggaVel',
                         'gga_velocity_ens_mps': 'ggaVelocityEns_mps',
                         'gga_velocity_method': 'ggaVelocityMethod',
                         'gps_HDOP_filter': 'gpsHDOPFilter',
                         'gps_HDOP_filter_change ': 'gpsHDOPFilterChange',
                         'gps_HDOP_filter_max': 'gpsHDOPFilterMax',
                         'gps_altitude_filter': 'gpsAltitudeFilter',
                         'gps_altitude_filter_change': 'gpsAltitudeFilterChange',
                         'gps_diff_qual_filter': 'gpsDiffQualFilter',
                         'hard_limit': 'hardLimit',
                         'hdop_ens': 'hdopEns',
                         'high_narrow': 'hn',
                         'high_wide': 'hw',
                         'in_transect_idx': 'inTransectIdx',
                         'initial_settings': 'initialSettings',
                         'int_cells': 'intCells',
                         'int_ens': 'intEns',
                         'interpolate_cells': 'interpolateCells',
                         'interpolate_ens': 'interpolateEns',
                         'invalid_95': 'invalid95',
                         'invalid_index': 'invalidIndex',
                         'invalid_user': 'invalid95User',
                         'left_idx': 'leftidx',
                         'low_narrow': 'ln',
                         'low_wide': 'lw',
                         'mag_error': 'magError',
                         'mag_orig_deg': 'magVarOrig_deg',
                         'mag_var_deg': 'magVar_deg',
                         'man_bot': 'manBot',
                         'man_exp': 'manExp',
                         'man_top': 'manTop',
                         'mb_dir_deg': 'mbDir_deg',
                         'mb_spd_mps': 'mbSpd_mps',
                         'mb_tests': 'mbTests',
                         'meas': 'meas_struct',
                         'middle_cells': 'middleCells',
                         'middle_ens': 'middleEns',
                         'moving_bed': 'movingBed',
                         'moving_bed_95': 'movingBed95',
                         'moving_user': 'movingBed95User',
                         'n_failed': 'nFailed',
                         'n_tests': 'nTests',
                         'nav_ref': 'navRef',
                         'near_speed_mps': 'nearBedSpeed_mps',
                         'noise_floor': 'noiseFloor',
                         'norm_data': 'normData',
                         'ns_exp': 'nsExponent',
                         'ns_exponent': 'nsexponent',
                         'num_invalid': 'numInvalid',
                         'num_sats_ens': 'numSatsEns',
                         'number_ensembles': 'numEns2Avg',
                         'orig_coord_sys': 'origCoordSys',
                         'orig_ref': 'origNavRef',
                         'orig_sys': 'origCoordSys',
                         'original_data': 'originalData',
                         'per_good_ens': 'perGoodEns',
                         'percent_invalid_bt': 'percentInvalidBT',
                         'percent_mb': 'percentMB',
                         'pitch_limit': 'pitchLimit',
                         'pp_exp': 'ppExponent',
                         'pp_exponent': 'ppexponent',
                         'processed_source': 'processedSource',
                         'q_cns_mean': 'qCNSmean',
                         'q_cns_opt_mean': 'qCNSoptmean',
                         'q_cns_opt_per_diff': 'qCNSoptperdiff',
                         'q_cns_per_diff': 'qCNSperdiff',
                         'q_man_mean': 'qManmean',
                         'q_man_per_diff': 'qManperdiff',
                         'q_ns_mean': 'q3pNSmean',
                         'q_ns_opt_mean': 'q3pNSoptmean',
                         'q_ns_opt_per_diff': 'q3pNSoptperdiff',
                         'q_ns_per_diff': 'q3pNSperdiff',
                         'q_pp_mean': 'qPPmean',
                         'q_pp_opt_mean': 'qPPoptmean',
                         'q_pp_opt_per_diff': 'qPPoptperdiff',
                         'q_pp_per_diff': 'qPPperdiff',
                         'q_run_threshold_caution': 'qRunThresholdCaution',
                         'q_run_threshold_warning': ' qRunThresholdWarning',
                         'q_sensitivity': 'qSensitivity',
                         'q_total_threshold_caution': ' qTotalThresholdWarning',
                         'q_total_threshold_warning': ' qTotalThresholdCaution',
                         'raw_gga_altitude_m': 'rawGGAAltitude_m',
                         'raw_gga_delta_time': 'rawGGADeltaTime',
                         'raw_gga_differential': 'rawGGADifferential',
                         'raw_gga_hdop': 'rawGGAHDOP',
                         'raw_gga_lat_deg': 'rawGGALat_deg',
                         'raw_gga_lon_deg': 'rawGGALon_deg',
                         'raw_gga_serial_time': 'rawGGASerialTime',
                         'raw_gga_utc': 'rawGGAUTC',
                         'raw_gga_num_sats': 'rawGGANumSats',
                         'raw_vel_mps': 'rawVel_mps',
                         'raw_vtg_course_deg': 'rawVTGCourse_deg',
                         'raw_vtg_delta_time': 'rawVTGDeltatTime',
                         'raw_vtg_mode_indicator': 'rawVTGModeIndicator',
                         'raw_vtg_speed_mps': 'rawVTGSpeed_mps',
                         'rec_method': 'recEdgeMethod',
                         'right_idx': 'rightidx',
                         'roll_limit': 'rollLimit',
                         'rssi_units': 'rssiUnits',
                         'sel_fit': 'selFit',
                         'serial_num': 'serialNum',
                         'sl_lag_effect_m': 'slLagEffect_m',
                         'sl_cutoff_number': 'slCutoffNum',
                         'sl_cutoff_percent': 'slCutoffPer',
                         'sl_cutoff_type': 'slCutoffType',
                         'smooth_depth': 'smoothDepth',
                         'smooth_filter': 'smoothFilter',
                         'smooth_lower_limit': 'smoothLowerLimit',
                         'smooth_speed': 'smoothSpeed',
                         'smooth_upper_limit': 'smoothUpperLimit',
                         'snr_filter': 'snrFilter',
                         'speed_of_sound_mps': 'speedOfSound_mps',
                         'snr_rng': 'snrRng',
                         'start_edge': 'startEdge',
                         'start_serial_time': 'startSerialTime',
                         'station_name': 'stationName',
                         'station_number': 'stationNumber',
                         'stationary_cs_track': 'stationaryCSTrack',
                         'stationary_mb_vel': 'stationaryMBVel',
                         'stationary_us_track': 'stationaryUSTrack',
                         'system_test': 'sysTest',
                         'system_tst': ' systemTest',
                         'systematic_user': 'systematicUser',
                         't_matrix': 'tMatrix',
                         'temperature': ' temperature',
                         'temperature_c': 'temperature_degC',
                         'test_quality': 'testQuality',
                         'time_stamp': 'timeStamp',
                         'top_ens': 'topEns',
                         'top_fit_r2': 'topfitr2',
                         'top_max_diff': 'topmaxdiff',
                         'top_method': 'topMethod',
                         'top_method_auto': 'topMethodAuto',
                         'top_method_orig': 'topMethodOrig',
                         'top_r2': 'topr2',
                         'total_95': 'total95',
                         'total_uncorrected': 'totalUncorrected',
                         'total_user': 'total95User',
                         'transect_duration_sec': 'transectDuration_sec',
                         'u_auto': 'uAuto',
                         'u_processed_mps': 'uProcessed_mps',
                         'u_earth_no_ref_mps': 'uEarthNoRef_mps',
                         'unit_normalize_dz': 'unitNormalizedz',
                         'unit_normalized': 'unitNormalized',
                         'unit_normalized_25': 'unitNormalized25',
                         'unit_normalized_75': 'unitNormalized75',
                         'unit_normalized_med': 'unitNormalizedMed',
                         'unit_normalized_no': 'unitNormalizedNo',
                         'use_2_correct': 'use2Correct',
                         'user_discharge_cms': 'userQ_cms',
                         'user_rating': 'userRating',
                         'user_valid': 'userValid',
                         'utm_ens_m': 'UTMEns_m',
                         'v_processed_mps': 'vProcessed_mps',
                         'v_earth_no_ref_mps': 'vEarthNoRef_mps',
                         'valid_beams': 'validBeams',
                         'valid_data': 'validData',
                         'valid_method': 'validDataMethod',
                         'vb_depths': 'vbDepths',
                         'vel_method': 'velMethod',
                         'vtg_vel': 'vtgVel',
                         'vtg_velocity_ens_mps': 'vtgVelocityEns_mps',
                         'vtg_velocity_method': 'vtgVelocityMethod',
                         'w_filter': 'wFilter',
                         'w_filter_threshold': 'wFilterThreshold',
                         'w_vel': 'wVel',
                         'water_mode': 'waterMode',
                         'wt_depth_filter': 'wtDepthFilter',
                         'z_auto': 'zAuto'}
        return py_2_mat_dict
