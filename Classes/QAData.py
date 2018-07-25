import numpy as np
from Classes.Uncertainty import Uncertainty
from Classes.QComp import QComp


class QAData(object):
    """Evaluates and stores quality assurance characteristics and messages.

    Attributes
    ----------
    q_run_threshold_caution: int
        Caution threshold for interpolated discharge for a run of invalid ensembles, in percent.
    q_run_threshold_warning: int
        Warning threshold for interpolated discharge for a run of invalid ensembles, in percent.
    q_total_threshold_caution: int
        Caution threshold for total interpolated discharge for invalid ensembles, in percent.
    q_total_threshold_warning: int
        Warning threshold for total interpolated discharge for invalid ensembles, in percent.
    transects: dict
        Dictionary of quality assurance checks for transects
    system_tst: dict
        Dictionary of quality assurance checks on the system test(s)
    compass: dict
        Dictionary of quality assurance checks on compass calibration and evaluations
    temperature: dict
        Dictionary of quality assurance checks on temperature comparions and variation
    moving_bed: dict
        Dictionary of quality assurance checks on moving-bed tests
    user: dict
        Dictionary of quality assurance checks on user input data
    bt_vel: dict
        Dictionary of quality assurance checks on bottom track velocities
    gga_vel: dict
        Dictionary of quality assurance checks on gga boat velocities
    vtg_vel: dict
        Dictionary of quality assurance checks on vtg boat velocities
    w_vel: dict
        Dictionary of quality assurance checks on water track velocities
    extrapolation: dict
        Dictionary of quality assurance checks on extrapolations
    edges: dict
        Dictionary of quality assurance checks on edges
    """
    
    def __init__(self, meas):
        """Checks the measurement for all quality assurance issues.

        Parameters
        ----------
        meas: Measurement
            Object of class Measurement
        """

        # Set default thresholds
        self.q_run_threshold_caution = 3
        self.q_run_threshold_warning = 5
        self.q_total_threshold_caution = 10
        self.q_total_threshold_warning = 25

        # Initialize instance variables
        self.transects = dict()
        self.system_tst = dict()
        self.compass = dict()
        self.temperature = dict()
        self.moving_bed = dict()
        self.user = dict()
        self.depths = dict()
        self.bt_vel = dict()
        self.gga_vel = dict()
        self.vtg_vel = dict()
        self.w_vel = dict()
        self.extrapolation = dict()
        self.edges = dict()

        # Apply QA checks
        self.transects_qa(meas)
        self.system_tst_qa(meas)
        self.compass_qa(meas)
        self.temperature_qa(meas)
        self.moving_bed_qa(meas)
        self.user_qa(meas)
        self.depths_qa(meas)
        self.boat_qa(meas)
        self.water_qa(meas)
        self.extrapolation_qa(meas)
        self.edges_qa(meas)

    def transects_qa(self, meas):
        """Apply quality checks to transects

        Parameters
        ----------
        meas: Measurement
            Object of class Measurement
        """

        # Assume good results
        self.transects['status'] = 'good'

        # Initialize keys
        self.transects['messages'] = []
        self.transects['recip'] = 0
        self.transects['sign'] = 0
        self.transects['duration'] = 0
        self.transects['number'] = 0
        self.transects['uncertainty'] = 0

        checked = []
        discharges = []
        start_edge = []
        for n in range(len(meas.transects)):
            checked.append(meas.transects[n].checked)
            if meas.transects[n].checked:
                discharges.append(meas.discharge[n])
                start_edge.append(meas.transects[n].start_edge)

        num_checked = np.nansum(np.asarray(checked))

        # Check duration
        total_duration = 0
        if num_checked >= 1:
            for transect in meas.transects:
                if transect.checked:
                    total_duration += transect.date_time.transect_duration_sec

        # Check duration against USGS policy
        if total_duration < 720:
            self.transects['status'] = 'caution'
            self.transects['messages'].append(
                ['Transects: Duration of selected transects is less than 720 seconds;', 2, 0])
            self.transects['duration'] = 1

        # Check transects for missing ensembles
        for transect in meas.transects:
            if transect.checked:

                # Determine number of missing ensembles
                if transect.adcp.manufacturer == 'SonTek':
                    # Determine number of missing ensembles for SonTek data
                    idx_missing = np.where(transect.date_time.ens_duration_sec > 1.5)[0]
                    if len(idx_missing) > 0:
                        average_ensemble_duration = (np.nansum(transect.date_time.ens_duration_sec)
                                                     - np.nansum(transect.date_time.ens_duration_sec[idx_missing])
                                                     / (len(transect.date_time.ens_duration_sec) - len(idx_missing)))
                        num_missing = np.round(np.nansum(transect.date_time.ens_duration_sec[idx_missing])
                                               / average_ensemble_duration) - len(idx_missing)
                    else:
                        num_missing = 0
                else:
                    # Determine number of lost ensembles for TRDI data
                    idx_missing = np.where(np.isnan(transect.date_time.ens_duration_sec) == True)[0]
                    num_missing = len(idx_missing) - 1

                # Save caution message
                if num_missing > 0:
                    self.transects['messages'].append([['Transects: ' + str(transect.file_name) + ' is missing'
                                                       + str(int(num_missing)) + ' ensembles;'], 2, 0])
                    self.transects['status'] = 'caution'

        # Check number of transects checked
        if num_checked == 0:
            # No transects selected
            self.transects['status'] = 'warning'
            self.transects['messages'].append(['TRANSECTS: No transects selected;', 1, 0])
            self.transects['number'] = 2
        elif num_checked == 1:
            # Only one transect selected
            self.transects['status'] = 'caution'
            self.transects['messages'].append(['Transects: Only one transect selected;', 2, 0])
            self.transects['number'] = 2
        else:
            self.transects['number'] = num_checked
            if num_checked == 2:
                # Only 2 transects selected
                cov, _ = Uncertainty.uncertainty_q_random(discharges, 'total')
                # Check uncertainty
                if cov > 2:
                    self.transects['status'] = 'caution'
                    self.transects['messages'].append(
                        ['Transects: Uncertainty would be reduced by additional transects;', 2, 0])

            # Check for consistent sign
            q_positive = []
            for q in discharges:
                if q.total >= 0:
                    q_positive.append(True)
                else:
                    q_positive.append(False)
            if len(np.unique(q_positive)) > 1:
                self.transects['status'] = 'warning'
                self.transects['messages'].append(
                    ['Transects: Sign of total Q is not consistent. One or more start banks may be incorrect;', 1, 0])

            # Check for reciprocal transects
            num_left = start_edge.count('Left')
            num_right = start_edge.count('Right')

            if not num_left == num_right:
                self.transects['status'] = 'warning'
                self.transects['messages'].append(['Transects: Transects selected are not reciprocal transects;', 1, 0])

        # Check for zero discharge transects
        q_zero = False
        for q in discharges:
            if q.total == 0:
                q_zero = True
        if q_zero:
            self.transects['status'] = 'warning'
            self.transects['messages'].append(['TRANSECTS: One or more transects have zero Q;', 1, 0])

    def system_tst_qa(self, meas):
        """Apply QA checks to system test.

        Parameters
        ----------
        meas: Measurement
            Object of class Measurement
        """

        self.system_tst['messages'] = []
        self.system_tst['status'] = 'good'

        # Determine is a system test was recorded
        if not meas.system_test:
            # No system test data recorded
            self.system_tst['status'] = 'warning'
            self.system_tst['messages'].append(['SYSTEM TEST: No system test;', 1, 3])
        else:

            pt3_fail = False
            num_tests_with_failure = 0

            for test in meas.system_test:
                if hasattr(test, 'result'):
                    if 'pt3' in test.result and test.result['pt3'] is not None:

                        # Check hard_limit, high gain, wide bandwidth
                        if 'hard_limit' in test.result['pt3']:
                            if 'high_wide' in test.result['pt3']['hard_limit']:
                                corr_table = test.result['pt3']['hard_limit']['high_wide']['corr_table']
                                if len(corr_table) > 0:
                                    # All lags past lag 2 should be less than 50% of lag 0
                                    qa_threshold = corr_table[0, :] * 0.5
                                    all_lag_check = np.greater(corr_table[3::, :], qa_threshold)

                                    # Lag 7 should be less than 25% of lag 0
                                    lag_7_check = np.greater(corr_table[7, :], corr_table[0, :] * 0.25)

                                    # If either condition is met for any beam the test fails
                                    if np.sum(np.sum(all_lag_check)) + np.sum(lag_7_check) > 1:
                                        pt3_fail = True

                    if test.result['n_failed'] is not None and test.result['n_failed'] > 0:
                        num_tests_with_failure += 1

            if pt3_fail:
                self.system_tst['status'] = 'caution'
                self.system_tst['messages'].append(
                    ['System Test: One or more PT3 tests in the system test indicate potential EMI;', 2, 3])

            # Check for failed tests
            if num_tests_with_failure == len(meas.system_test):
                # All tests had a failure
                self.system_tst['status'] = 'warning'
                self.system_tst['messages'].append(
                    ['SYSTEM TEST: All system test sets have at least one test that failed;', 1, 3])
            elif num_tests_with_failure > 0:
                self.system_tst['status'] = 'caution'
                self.system_tst['messages'].append(
                    ['System Test: One or more system test sets have at least one test that failed;', 2, 3])

    def compass_qa(self, meas):
        """Apply QA checks to compass calibration and evaluation.

        Parameters
        ----------
        meas: Measurement
            Object of class Measurement
        """

        self.compass['messages'] = []

        checked = []
        for transect in meas.transects:
            checked.append(transect.checked)

        if np.any(checked):
            heading = np.unique(meas.transects[checked.index(1)].sensors.heading_deg.internal.data)
        else:
            heading = np.array([0])

        # Intialize variable as if ADCP has no compass
        self.compass['status'] = 'inactive'
        self.compass['status1'] = 'good'
        self.compass['status2'] = 'good'
        self.compass['magvar'] = 0
        self.compass['magvar_idx'] = 0

        if len(heading) > 1 and np.any(np.not_equal(heading, 0)):
            # ADCP has a compass
            # A compass calibration is required is a loop test or GPS are used

            # Check for loop test
            loop = False
            for test in meas.mb_tests:
                if test.type == 'Loop':
                    loop = True

            # Check for GPS data
            gps = False
            if meas.transects[checked.index(True)].boat_vel.gga_vel is not None or \
               meas.transects[checked.index(True)].boat_vel.vtg_vel is not None:
                gps = True

            if gps or loop:
                # Compass calibration is required

                # Determine the ADCP manufacturer
                if meas.transects[checked.index(True)].adcp.manufacturer == 'SonTek':
                    # SonTek ADCP
                    if not meas.compass_cal:
                        # No compass calibration
                        self.compass['status1'] = 'warning'
                        self.compass['messages'].append(['COMPASS: No compass calibration;', 1, 4])
                    elif meas.compass_cal[-1].result['compass']['error'] == 'N/A':
                        # If the error cannot be decoded from the calibration assume the calibration is good
                        self.compass['status1'] = 'good'
                    else:
                        if meas.compass_cal[-1].result['compass']['error'] <= 0.2:
                            self.compass['status1'] = 'good'
                        else:
                            self.compass['status1'] = 'caution'
                            self.compass['messages'].append(['COMPASS: Calibration result > 0.2 deg;', 2, 4])

                elif meas.transects[checked.index(True)].adcp.manufacturer == 'TRDI':
                    # TRDI ADCP
                    if not meas.compass_cal:
                        # No compass calibration
                        if not meas.compass_eval:
                            # No calibration or evaluation
                            self.compass['status1'] = 'warning'
                            self.compass['messages'].append(['COMPASS: No compass calibration or evaluation;', 1, 4])
                        else:
                            # No calibration but an evaluation was completed
                            self.compass['status1'] = 'caution'
                            self.compass['messages'].append(['Compass: No compass calibration;', 2, 4])
                    else:
                        # Compass was calibrated
                        if not meas.compass_eval:
                            # No compass evaluation
                            self.compass['status1'] = 'caution'
                            self.compass['messages'].append(['Compass: No compass evaluation;', 2, 4])
                        else:
                            # Check results of evaluation
                            try:
                                if float(meas.compass_cal[-1].result['compass']['error']) <= 1:
                                    self.compass['status1'] = 'good'
                                else:
                                    self.compass['status1'] = 'caution'
                                    self.compass['messages'].append(['Compass: Evaluation result > 1 deg;', 2, 4])
                            except ValueError:
                                self.compass['status1'] = 'good'
            else:
                # Compass not required
                if (not meas.compass_cal) and (not meas.compass_eval):
                    # No compass calibration or evaluation
                    self.compass['status1'] = 'default'
                else:
                    # Compass was calibrated and evaluated
                    self.compass['status1'] = 'good'

            # Check for consistent magvar
            magvar = []
            mag_error_exceeded = []
            pitch_mean = []
            pitch_std = []
            pitch_exceeded = []
            roll_mean = []
            roll_std = []
            roll_exceeded = []
            for transect in meas.transects:
                if transect.checked:
                    heading_source_selected = getattr(
                        transect.sensors.heading_deg, transect.sensors.heading_deg.selected)
                    pitch_source_selected = getattr(transect.sensors.pitch_deg, transect.sensors.pitch_deg.selected)
                    roll_source_selected = getattr(transect.sensors.roll_deg, transect.sensors.roll_deg.selected)

                    magvar.append(heading_source_selected.mag_var_deg)

                    pitch_mean.append(np.nanmean(pitch_source_selected.data))
                    pitch_std.append(np.nanstd(pitch_source_selected.data))
                    roll_mean.append(np.nanmean(roll_source_selected.data))
                    roll_std.append(np.nanstd(roll_source_selected.data))

                    # SonTek G3 compass provides pitch, roll, and magnetic error parameters that can be checked
                    if meas.transects[checked.index(True)].adcp.manufacturer == 'SonTek':
                        if heading_source_selected.pitch_limit is not None:
                            # Check for bug in SonTek data where pitch and roll was n x 3 use n x 1
                            if len(pitch_source_selected.data.shape) == 1:
                                pitch_data = pitch_source_selected.data
                            else:
                                pitch_data = pitch_source_selected.data[:, 0]
                            idx_max = np.where(pitch_data > heading_source_selected.pitch_limit[0])[0]
                            idx_min = np.where(pitch_data < heading_source_selected.pitch_limit[1])[0]
                            if len(idx_max) > 0 or len(idx_min) > 0:
                                pitch_exceeded.append(True)
                            else:
                                pitch_exceeded.append(False)

                        if heading_source_selected.roll_limit is not None:
                            if len(roll_source_selected.data.shape) == 1:
                                roll_data = roll_source_selected.data
                            else:
                                roll_data = roll_source_selected.data[:, 0]
                            idx_max = np.where(roll_data > heading_source_selected.pitch_limit[0])[0]
                            idx_min = np.where(roll_data < heading_source_selected.pitch_limit[1])[0]
                            if len(idx_max) > 0 or len(idx_min) > 0:
                                roll_exceeded.append(True)
                            else:
                                roll_exceeded.append(False)

                        if heading_source_selected.mag_error is not None:
                            idx_max = np.where(heading_source_selected.mag_error > 2)[0]
                            if len(idx_max) > 0:
                                mag_error_exceeded.append(True)
                            else:
                                mag_error_exceeded.append(False)

            if len(np.unique(magvar)) > 1:
                self.compass['status2'] = 'caution'
                self.compass['messages'].append(
                    ['Compass: Magnetic variation is not consistent among transects;', 2, 4])
                self.compass['magvar'] = 1

            # Check that magvar was set if GPS data are available
            if gps:
                if 0 in magvar:
                    self.compass['status2'] = 'warning'
                    self.compass['messages'].append(
                        ['COMPASS: Magnetic variation is 0 and GPS data are present;', 1, 4])
                    self.compass['magvar'] = 2
                    self.compass['magvar_idx'] = magvar.index(0)

            # Check pitch mean
            if np.any(np.asarray(pitch_mean) > 8):
                self.compass['status2'] = 'warning'
                self.compass['messages'].append(['PITCH: One or more transects have a mean pitch > 8 deg;', 1, 4])
            elif np.any(np.asarray(pitch_mean) > 4):
                if self.compass['status2'] == 'good':
                    self.compass['status2'] = 'caution'
                self.compass['messages'].append(['Pitch: One or more transects have a mean pitch > 4 deg;', 2, 4])

            # Check roll mean
            if np.any(np.asarray(roll_mean) > 8):
                self.compass['status2'] = 'warning'
                self.compass['messages'].append(['ROLL: One or more transects have a mean roll > 8 deg;', 1, 4])
            elif np.any(np.asarray(roll_mean) > 4):
                if self.compass['status2'] == 'good':
                    self.compass['status2'] = 'caution'
                self.compass['messages'].append(['Roll: One or more transects have a mean roll > 4 deg;', 2, 4])

            # Check pitch standard deviation
            if np.any(np.asarray(pitch_std) > 5):
                if self.compass['status2'] == 'good':
                    self.compass['status2'] = 'caution'
                self.compass['messages'].append(['Pitch: One or more transects have a pitch std dev > 5 deg;', 2, 4])

            # Check roll standard deviation
            if np.any(np.asarray(roll_std) > 5):
                if self.compass['status2'] == 'good':
                    self.compass['status2'] = 'caution'
                self.compass['messages'].append(['Roll: One or more transects have a roll std dev > 5 deg;', 2, 4])

            # Additional checks for SonTek G3 compass
            if meas.transects[checked.index(True)].adcp.manufacturer == 'SonTek':
                # Check if pitch limits were exceeded
                if any(pitch_exceeded):
                    if self.compass['status2'] == 'good':
                        self.compass['status2'] = 'caution'
                    self.compass['messages'].append(
                        ['Compass: One or more transects have pitch exceeding calibration limits;', 2, 4])

                # Check if roll limits were exceeded
                if any(roll_exceeded):
                    if self.compass['status2'] == 'good':
                        self.compass['status2'] = 'caution'
                    self.compass['messages'].append(
                        ['Compass: One or more transects have roll exceeding calibration limits;', 2, 4])

                # Check if magnetic error was exceeded
                if any(mag_error_exceeded):
                    if self.compass['status2'] == 'good':
                        self.compass['status2'] = 'caution'
                    self.compass['messages'].append(
                        ['Compass: One or more transects have a change in mag field exceeding 2%;', 2, 4])

            if self.compass['status1'] == 'warning' or self.compass['status2'] == 'warning':
                self.compass['status'] = 'warning'
            elif self.compass['status1'] == 'caution' or self.compass['status2'] == 'caution':
                self.compass['status'] = 'caution'
            else:
                self.compass['status'] = 'good'

    def temperature_qa(self, meas):
        """Apply QA checks to temperature.

        Parameters
        ----------
        meas: Measurement
            Object of class Measurement
        """

        self.temperature['messages'] = []
        check = [0, 0]

        # Create array of all temperatures
        temp = np.array([])
        checked = []
        for transect in meas.transects:
            if transect.checked:
                checked.append(transect.checked)
                temp_selected = getattr(transect.sensors.temperature_deg_c, transect.sensors.temperature_deg_c.selected)
                if len(temp) == 0:
                    temp = temp_selected.data
                else:
                    temp = np.hstack((temp, temp_selected.data))

        # Check temperature range
        if np.any(checked):
            temp_range = np.nanmax(temp) - np.nanmin(temp)
        else:
            temp_range = 0

        if temp_range > 2:
            check[0] = 3
            self.temperature['messages'].append(['TEMPERATURE: Temperature range is '
                                                + '%3.1f % temp_range'
                                                + 'degrees C which is greater than 2 degrees;', 1, 5])
        elif temp_range > 1:
            check[0] = 2
            self.temperature['messages'].append(['TEMPERATURE: Temperature range is '
                                                 + '%3.1f % temp_range'
                                                 + 'degrees C which is greater than 1 degrees;', 2, 5])
        else:
            check[0] = 1

        # Check for independent temperature reading
        if 'user' in meas.ext_temp_chk:
            try:
                user = float(meas.ext_temp_chk['user'])
            except ValueError:
                user = None
            if user is None:
                # No independent temperature reading
                check[1] = 2
                self.temperature['messages'].append(['Temperature: No independent temperature reading;', 2, 5])
            elif 'adcp' in meas.ext_temp_chk:
                # Compare user to manually entered ADCP temperature
                diff = np.abs(user - meas.ext_temp_chk['adcp'])
                if diff < 2:
                    check[1] = 1
                else:
                    check[1] = 3
                    self.temperature['messages'].append(
                        ['TEMP.: The difference between ADCP and reference is > 2:  ' + '%3.1f % diff' + ' C;', 1, 5])
            else:
                # Compare user to mean of all temperature data
                diff = np.abs(user - np.nanmean(temp))
                if diff < 2:
                    check[1] = 1
                else:
                    check[1] = 3
                    self.temperature['messages'].append(
                        ['TEMP.: The difference between ADCP and reference is > 2:  ' + '%3.1f % diff' + ' C;', 1, 5])

        # Assign temperature status
        max_check = max(check)
        if max_check == 1:
            self.temperature['status'] = 'good'
        elif max_check == 2:
            self.temperature['status'] = 'caution'
        elif max_check == 3:
            self.temperature['status'] = 'warning'

    def moving_bed_qa(self, meas):
        """Applies quality checks to moving-bed tests.

        Parameters
        ----------
        meas: Measurement
            Object of class Measurement
        """

        self.moving_bed['messages'] = []
        self.moving_bed['code'] = 0

        # Are there moving-bed tests?
        if len(meas.mb_tests) < 1:
            # No moving-bed test
            self.moving_bed['messages'].append(['MOVING-BED TEST: No moving bed test;', 1, 6])
            self.moving_bed['status'] = 'warning'
            self.moving_bed['code'] = 3

        else:
            # Moving-bed tests available
            mb_data = meas.mb_tests

            # Are tests valid according to the user
            user_valid_test = []
            file_names = []
            idx_selected = []
            test_quality = []
            mb_tests = []
            mb = []
            mb_test_type = []
            loop = []
            for n, test in enumerate(mb_data):
                if test.user_valid:
                    user_valid_test.append(True)
                    file_names.append(test.transect.file_name)
                    if test.type == 'Loop' and not test.test_quality == 'Errors':
                        loop.append(test.moving_bed)
                    # Selected test
                    if test.selected:
                        idx_selected.append(n)
                        test_quality.append(test.test_quality)
                        mb_tests.append(test)
                        mb.append(test.moving_bed)
                        mb_test_type.append(test.type)
                else:
                    user_valid_test.append(False)

            if not any(user_valid_test):
                # No valid test according to user
                self.moving_bed['messages'].append(['MOVING-BED TEST: No valid moving-bed test based on user input;',
                                                    1, 6])
                self.moving_bed['status'] = 'warning'
                self.moving_bed['code'] = 3
            else:
                # Check for duplicate valid moving-bed tests
                if len(np.unique(file_names)) < len(file_names):
                    self.moving_bed['messages'].append([
                        'MOVING-BED TEST: Duplicate moving-bed test files marked valid;', 1, 6])
                    self.moving_bed['status'] = 'warning'
                    self.moving_bed['code'] = 3

            if self.moving_bed['code'] == 0:
                # Check test quality
                if len(test_quality) > 0 and sum(np.array(test_quality) == 'Good') > 0:
                    self.moving_bed['status'] = 'good'
                    self.moving_bed['code'] = 1

                    # Check if there is a moving-bed
                    if any(mb):
                        # Moving-bed present
                        self.moving_bed['messages'].append(
                            ['Moving-Bed Test: A moving-bed is present, use GPS or moving-bed correction;', 2, 6])
                        self.moving_bed['code'] = 2
                        self.moving_bed['status'] = 'caution'

                        # Check for test type
                        if sum(np.array(mb_test_type) == 'Stationary'):
                            # Check for GPS or 3 stationary tests
                            if len(mb_tests) < 3:
                                gps = []
                                for transect in meas.transects:
                                    if transect.checked:
                                        if transect.gps is None:
                                            gps.append(False)
                                        else:
                                            gps.append(True)
                                if not all(gps):
                                    # GPS not available for all selected transects
                                    self.moving_bed['messages'].append([
                                        'Moving-Bed Test: '
                                        + 'Less than 3 stationary tests available for moving-bed correction;',
                                        2, 6])

                elif len(test_quality) > 0 and sum(np.array(test_quality) == 'Warnings') > 0:
                    # Quality check has warnings
                    self.moving_bed['messages'].append(['Moving-Bed Test: The moving-bed test(s) has warnings, '
                                                        + 'please review tests to determine validity;', 2, 6])
                    self.moving_bed['status'] = 'caution'
                    self.moving_bed['code'] = 2

                elif len(test_quality) > 0 and sum(np.array(test_quality) == 'Manual') > 0:
                    # Manual override used
                    self.moving_bed['messages'].append(['MOVING-BED TEST: '
                                                        + 'The user has manually forced the use of some tests;', 1, 6])
                    self.moving_bed['status'] = 'warning'
                    self.moving_bed['code'] = 3

                else:
                    # Test has critical errors
                    self.moving_bed['messages'].append(['MOVING-BED TEST: The moving-bed test(s) have critical errors '
                                                        + 'and will not be used;', 1, 6])
                    self.moving_bed['status'] = 'warning'
                    self.moving_bed['code'] = 3

                # Check multiple loops for consistency
                if len(np.unique(loop)) > 1:
                    self.moving_bed['messages'].append(['Moving-Bed Test: Results of valid loops are not consistent, '
                                                        + 'review moving-bed tests;', 2, 6])
                    if self.moving_bed['code'] < 3:
                        self.moving_bed['code'] = 2
                        self.moving_bed['status'] = 'caution'

    def user_qa(self, meas):
        """Apply quality checks to user input data.

        Parameters
        ----------
        meas: Measurement
            Object of class Measurement
        """

        self.user['messages'] = []
        self.user['status'] = 'good'

        # Check for Station Name
        self.user['sta_name'] = False
        if meas.station_name is None:
            self.user['messages'].append(['Site Info: Station name not entered;', 2, 2])
            self.user['status'] = 'caution'
            self.user['sta_name'] = True

        # Check for Station Number
        self.user['sta_number'] = False
        if meas.station_number is None:
            self.user['messages'].append(['Site Info: Station number not entered;', 2, 2])
            self.user['status'] = 'caution'
            self.user['sta_name'] = True

    def depths_qa(self, meas):
        """Apply quality checks to depth data.

        Parameters
        ----------
        meas: Measurement
            Object of class Measurement
        """

        # Initialize variables
        n_transects = len(meas.transects)
        self.depths['q_total'] = np.tile(np.nan, n_transects)
        self.depths['q_max_run'] = np.tile(np.nan, n_transects)
        self.depths['q_total_caution'] = np.tile(False, n_transects)
        self.depths['q_run_caution'] = np.tile(False, n_transects)
        self.depths['q_total_warning'] = np.tile(False, n_transects)
        self.depths['q_run_warning'] = np.tile(False, n_transects)
        self.depths['all_invalid'] = np.tile(False, n_transects)
        self.depths['messages'] = []
        self.depths['status'] = 'good'
        self.depths['draft'] = 0
        checked = []
        drafts = []
        for n, transect in enumerate(meas.transects):
            checked.append(transect.checked)
            if transect.checked:
                in_transect_idx = transect.in_transect_idx

                depths_selected = getattr(transect.depths, transect.depths.selected)
                drafts.append(depths_selected.draft_use_m)

                # Determine valid measured depths
                if transect.depths.composite:
                    depth_na = depths_selected.depth_source_ens[in_transect_idx] != 'NA'
                    depth_in = depths_selected.depth_source_ens[in_transect_idx] != 'IN'
                    depth_valid = np.all(np.vstack((depth_na, depth_in)), 0)
                else:
                    depth_valid_temp = depths_selected.valid_data[in_transect_idx]
                    depth_nan = depths_selected.depth_processed_m[in_transect_idx] != np.nan
                    depth_valid = np.all(np.vstack((depth_nan, depth_valid_temp)), 0)

                if not np.any(depth_valid):
                    self.depths['all_invalid'][n] = True

                # Compute QA characteristics
                q_total, q_max_run, number_invalid_ensembles = QAData.invalid_qa(depth_valid, meas.discharge[n])
                self.depths['q_total'][n] = q_total
                self.depths['q_max_run'][n] = q_max_run

                # Compute percentage compared to total
                q_total_percent = np.abs((q_total / meas.discharge[n].total) * 100)
                q_max_run_percent = np.abs((q_max_run / meas.discharge[n].total) * 100)

                # Apply total interpolated discharge threshold
                if q_total_percent > self.q_total_threshold_warning:
                    self.depths['q_total_warning'][n] = True
                elif q_total_percent > self.q_total_threshold_caution:
                    self.depths['q_total_caution'][n] = True

                # Apply interpolated discharge run thresholds
                if q_max_run_percent > self.q_run_threshold_warning:
                    self.depths['q_run_warning'][n] = True
                elif q_max_run_percent > self.q_run_threshold_caution:
                    self.depths['q_run_caution'][n] = True

        if checked:

            # Create array of all unique draft values
            draft_check = np.unique(np.round(drafts, 3))

            # Check draft consistency
            if len(draft_check) > 1:
                self.depths['status'] = 'caution'
                self.depths['draft'] = 1
                self.depths['messages'].append(['Depth: Transducer depth is not consistent among transects;', 2, 10])

            # Check for zero draft
            if np.any(np.less(draft_check, 0.01)):
                self.depths['status'] = 'warning'
                self.depths['draft'] = 2
                self.depths['messages'].append(['DEPTH: Transducer depth is too shallow, likely 0;', 1, 10])

            # Check consecutive interpolated discharge criteria
            if np.any(self.depths['q_run_warning']):
                self.depths['messages'].append(['DEPTH: Int. Q for consecutive invalid ensembles exceeds '
                                                + '%2.0f % self.q_run_threshold_warning' + '%;', 1, 10])
                self.depths['status'] = 'warning'
            elif np.any(self.depths['q_run_caution']):
                self.depths['messages'].append(['Depth: Int. Q for consecutive invalid ensembles exceeds '
                                                + '%2.0f % self.q_run_threshold_caution' + '%;', 2, 10])
            self.depths['status'] = 'caution'

            # Check total interpolated discharge criteria
            if np.any(self.depths['q_total_warning']):
                self.depths['messages'].append(['DEPTH: Int. Q for invalid ensembles in a transect exceeds '
                                                + '%2.0f % self.q_total_threshold_warning' + '%;', 1, 10])
                self.depths['status'] = 'warning'
            elif np.any(self.depths['q_total_caution']):
                self.depths['messages'].append(['Depth: Int. Q for invalid ensembles in a transect exceeds '
                                                + '%2.0f % self.q_total_threshold_caution' + '%;', 2, 10])
                self.depths['status'] = 'caution'

            # Check if all depths are invalid
            if np.any(self.depths['all_invalid']):
                self.depths['messages'].append(['DEPTH: There are no valid depths for one or more transects.', 2, 10])
                self.depths['status'] = 'warning'

        else:
            self.depths['status'] = 'inactive'

    def boat_qa(self, meas):
        """Apply quality checks to boat data.

        Parameters
        ----------
        meas: Measurement
            Object of class Measurement
        """

        # Initialize variables
        n_transects = len(meas.transects)
        data_type = {'BT': {'class': 'bt_vel', 'warning': 'BT-', 'caution': 'bt-',
                            'filter': [('All: ', 0), ('Original: ', 1), ('ErrorVel: ', 2),
                                       ('VertVel: ', 3), ('Other: ', 4), ('3Beams: ', 5)]},
                     'GGA': {'class': 'gga_vel', 'warning': 'GGA-', 'caution': 'gga-',
                             'filter': [('All: ', 0), ('Original: ', 1), ('DGPS: ', 2),
                                        ('Altitude: ', 3), ('Other: ', 4), ('HDOP: ', 5)]},
                     'VTG': {'class': 'vtg_vel', 'warning': 'VTG-', 'caution': 'vtg-',
                             'filter': [('All: ', 0), ('Original: ', 1), ('HDOP: ', 5)]}}

        for dt_key, dt_value in data_type.items():
            boat = getattr(self, dt_value['class'])
            # Initialize dictionaries for each data type
            boat['q_total_caution'] = np.tile(False, (n_transects, 6))
            boat['q_max_run_caution'] = np.tile(False, (n_transects, 6))
            boat['q_total_warning'] = np.tile(False, (n_transects, 6))
            boat['q_max_run_warning'] = np.tile(False, (n_transects, 6))
            boat['all_invalid'] = np.tile(False, n_transects)
            boat['q_total'] = np.tile(np.nan, (n_transects, 6))
            boat['q_max_run'] = np.tile(np.nan, (n_transects, 6))
            boat['messages'] = []

            status_switch = 0
            avg_speed_check = 0

            # Check the results of each filter
            for dt_filter in dt_value['filter']:
                boat['status'] = 'inactive'

                # Quality check each transect
                for n, transect in enumerate(meas.transects):

                    # Evaluate on transects used in the discharge computation
                    if transect.checked:

                        in_transect_idx = transect.in_transect_idx

                        # Check to see if data are available for the data_type
                        if getattr(transect.boat_vel, dt_value['class']) is not None:
                            boat['status'] = 'good'

                            # Compute quality characteristics
                            valid = getattr(transect.boat_vel, dt_value['class']).valid_data[dt_filter[1],
                                                                                             in_transect_idx]
                            q_total, q_max_run, number_invalid_ens = QAData.invalid_qa(valid, meas.discharge[n])
                            boat['q_total'][n, dt_filter[1]] = q_total
                            boat['q_max_run'][n, dt_filter[1]] = q_max_run

                            # Compute percentage compared to total
                            q_total_percent = np.abs((q_total / meas.discharge[n].total) * 100)
                            q_max_run_percent = np.abs((q_max_run / meas.discharge[n].total) * 100)

                            # Check if all invalid
                            if dt_filter[1] == 0 and not np.any(valid):
                                boat['all_invalid'][n] = True

                            # Apply total interpolated discharge threshold
                            if q_total_percent > self.q_total_threshold_warning:
                                boat['q_total_warning'][n, dt_filter[1]] = True
                            elif q_total_percent > self.q_total_threshold_caution:
                                boat['q_total_caution'][n, dt_filter[1]] = True

                            # Apply interpolated discharge run thresholds
                            if q_max_run_percent > self.q_run_threshold_warning:
                                boat['q_max_run_warning'][n, dt_filter[1]] = True
                            elif q_max_run_percent > self.q_run_threshold_caution:
                                boat['q_max_run_caution'][n, dt_filter[1]] = True

                            # Check boat velocity for vtg data
                            if dt_key is 'VTG' and transect.boat_vel.selected is 'vtg_vel' and avg_speed_check == 0:
                                avg_speed = np.nanmean((transect.boat_vel.vtg_vel.u_mps**2
                                                        + transect.boat_vel.vtg_vel.v_mps**2)**0.5)
                                if avg_speed < 0.24:
                                    boat['q_total_caution'][n, dt_filter[1]] = True
                                    boat['messages'].append(
                                        ['vtg-AvgSpeed: VTG data may not be accurate for average boat speed less than'
                                         + '0.24 m/s (0.8 ft/s);', 2, 8])
                                    avg_speed_check = 1

                # Create message for consecutive invalid discharge
                if boat['q_max_run_warning'][:, dt_filter[1]].any():
                    if dt_key is 'BT':
                        module_code = 7
                    else:
                        module_code = 8
                    boat['messages'].append(
                        [dt_value['warning'] + dt_filter[0] +
                            'Int. Q for consecutive invalid ensembles exceeds ' +
                            '%3.1f' % self.q_run_threshold_warning + '%;', 1, module_code])
                    status_switch = 2
                elif boat['q_max_run_caution'][:, dt_filter[1]].any():
                    if dt_key is 'BT':
                        module_code = 7
                    else:
                        module_code = 8
                    boat['messages'].append(
                        [dt_value['caution'] + dt_filter[0] +
                            'Int. Q for consecutive invalid ensembles exceeds ' +
                            '%3.1f' % self.q_run_threshold_caution + '%;', 2, module_code])
                    if status_switch < 1:
                        status_switch = 1

                # Create message for total invalid discharge
                if boat['q_total_warning'][:, dt_filter[1]].any():
                    if dt_key is 'BT':
                        module_code = 7
                    else:
                        module_code = 8
                    boat['messages'].append(
                        [dt_value['warning'] + dt_filter[0] +
                            'Int. Q for invalid ensembles in a transect exceeds ' +
                            '%3.1f' % self.q_total_threshold_warning + '%;', 1, module_code])
                    status_switch = 2
                elif boat['q_max_run_caution'][:, dt_filter[1]].any():
                    if dt_key is 'BT':
                        module_code = 7
                    else:
                        module_code = 8
                    boat['messages'].append(
                        [dt_value['caution'] + dt_filter[0] +
                            'Int. Q for invalid ensembles in a transect exceeds ' +
                            '%3.1f' % self.q_total_threshold_caution + '%;', 2, module_code])
                    if status_switch < 1:
                        status_switch = 1

            # Create message for all data invalid
            if boat['all_invalid'].any():
                boat['status'] = 'warning'
                if dt_key is 'BT':
                    module_code = 7
                else:
                    module_code = 8
                boat['messages'].append(
                    [dt_value['warning'] + dt_value['filter'][0][0] +
                        'There are no valid data for one or more transects.;', 1, module_code])

            # Set status
            if status_switch == 2:
                boat['status'] = 'warning'
            elif status_switch == 1:
                boat['status'] = 'caution'

            setattr(self, dt_value['class'], boat)

    def water_qa(self, meas):
        """Apply quality checks to water data.

        Parameters
        ----------
        meas: Measurement
            Object of class Measurement
        """

        # Initialize filter labels and indices
        prefix = ['All: ', 'Original: ', 'ErrorVel: ', 'VertVel: ', 'Other: ', '3Beams: ', 'SNR:']
        if meas.transects[0].adcp.manufacturer is 'TRDI':
            filter_index = [0, 1, 2, 3, 4, 5]
        else:
            filter_index = [0, 1, 2, 3, 4, 5, 7]

        n_transects = len(meas.transects)
        n_filters = len(filter_index) + 1
        # Initialize dictionaries for each data type
        self.w_vel['q_total_caution'] = np.tile(False, (n_transects, n_filters))
        self.w_vel['q_max_run_caution'] = np.tile(False, (n_transects, n_filters))
        self.w_vel['q_total_warning'] = np.tile(False, (n_transects, n_filters))
        self.w_vel['q_max_run_warning'] = np.tile(False, (n_transects, n_filters))
        self.w_vel['all_invalid'] = np.tile(False, n_transects)
        self.w_vel['q_total'] = np.tile(np.nan, (n_transects, n_filters))
        self.w_vel['q_max_run'] = np.tile(np.nan, (n_transects, n_filters))
        self.w_vel['messages'] = []
        status_switch = 0

        # TODO if meas had a property checked as list it would save creating that list multiple times
        checked = []
        for transect in meas.transects:
            checked.append(transect.checked)

        # At least one transect is being used to compute discharge
        if any(checked):
            # Loop through filters
            for prefix_idx, filter_idx in enumerate(filter_index):
                # Loop through transects
                for n, transect in enumerate(meas.transects):
                    if transect.checked:
                        valid_original = np.any(transect.w_vel.valid_data[1, :, transect.in_transect_idx].T, 0)

                        # Determine what data each filter have marked invalid. Original invalid data are excluded
                        valid = np.any(transect.w_vel.valid_data[filter_idx, :, transect.in_transect_idx].T, 0)
                        if filter_idx > 1:
                            valid_int = valid.astype(int) - valid_original.astype(int)
                            valid = valid_int != -1

                        # Check if all data are invalid
                        if filter_idx == 0:
                            if np.nansum(valid.astype(int)) < 1:
                                self.w_vel['all_invalid'][n] = True
                        # TODO seems like the rest of this should be under else of all invalid or multiple messages
                        # generated.

                        # Compute characteristics
                        q_total, q_max_run, number_invalid_ens = QAData.invalid_qa(valid, meas.discharge[n])
                        self.w_vel['q_total'][n, filter_idx] = q_total
                        self.w_vel['q_max_run'][n, filter_idx] = q_max_run

                        # Compute percentage compared to total
                        q_total_percent = np.abs((q_total / meas.discharge[n].total) * 100)
                        q_max_run_percent = np.abs((q_max_run / meas.discharge[n].total) * 100)

                        # Check total invalid discharge in ensembles for warning
                        if q_total_percent > self.q_total_threshold_warning:
                            self.w_vel['q_total_warning'][n, filter_idx] = True

                        # Apply run or cluster thresholds
                        if q_max_run_percent > self.q_run_threshold_warning:
                            self.w_vel['q_max_run_warning'][n, filter_idx] = True
                        elif q_max_run_percent > self.q_run_threshold_caution:
                            self.w_vel['q_max_run_caution'][n, filter_idx] = True

                        # Compute percent discharge interpolated for both cells and ensembles
                        # This approach doesn't exclude original data
                        valid_cells = transect.w_vel.valid_data[filter_idx, :, transect.in_transect_idx].T
                        q_invalid_total = np.nansum(meas.discharge[n].middle_cells[np.logical_not(valid_cells)]) \
                            + np.nansum(meas.discharge[n].top_ens[np.logical_not(valid)]) \
                            + np.nansum(meas.discharge[n].bottom_ens[np.logical_not(valid)])
                        q_invalid_total_percent = (q_invalid_total / meas.discharge[n].total) * 100

                        if q_invalid_total_percent > self.q_total_threshold_caution:
                            self.w_vel['q_total_caution'][n, filter_idx] = True

                # Generate messages for ensemble run or clusters
                if np.any(self.w_vel['q_max_run_warning'][:, filter_idx]):
                    self.w_vel['messages'].append(['WT-' + prefix[prefix_idx]
                                                   + 'Int. Q for consecutive invalid ensembles exceeds '
                                                   + '%3.0f' % self.q_run_threshold_warning
                                                   + '%;', 1, 11])
                    status_switch = 2
                elif np.any(self.w_vel['q_max_run_caution'][:, filter_idx]):
                    self.w_vel['messages'].append(['wt-' + prefix[prefix_idx]
                                                   + 'Int. Q for consecutive invalid ensembles exceeds '
                                                   + '%3.0f' % self.q_run_threshold_caution
                                                   + '%;', 2, 11])
                    if status_switch < 1:
                        status_switch = 1

                # Generate message for total_invalid Q
                if np.any(self.w_vel['q_total_warning'][:, filter_idx]):
                    self.w_vel['messages'].append(['WT-' + prefix[prefix_idx]
                                                   + 'Int. Q for invalid ensembles in a transect exceeds '
                                                   + '%3.0f' % self.q_total_threshold_warning
                                                   + '%;', 1, 11])
                    status_switch = 2
                elif np.any(self.w_vel['q_total_caution'][:, filter_idx]):
                    self.w_vel['messages'].append(['wt-' + prefix[prefix_idx]
                                                   + 'Int. Q for invalid cells and ensembles in a transect exceeds '
                                                   + '%3.0f' % self.q_total_threshold_caution
                                                   + '%;', 2, 11])
                    if status_switch < 1:
                        status_switch = 1

            # Generate message for all invalid
            if np.any(self.w_vel['all_invalid']):
                self.w_vel['messages'].append(['WT-', prefix[0], 'There are no valid data for one or more transects.',
                                               1, 11])
                status_switch = 2

            # Set status
            self.w_vel['status'] = 'good'
            if status_switch == 2:
                self.w_vel['status'] = 'warning'
            elif status_switch == 1:
                self.w_vel['status'] = 'caution'
        else:
            self.w_vel['status'] = 'inactive'

    def extrapolation_qa(self, meas):
        """Apply quality checks to extrapolation methods

        Parameters
        ----------
        meas: Measurement
            Object of class Measurement
        """

        self.extrapolation['messages'] = []

        checked = []
        discharges = []
        for n, transect in enumerate(meas.transects):
            checked.append(transect.checked)
            if transect.checked:
                discharges.append(meas.discharge[n])

        if any(checked):
            self.extrapolation['status'] = 'good'
            extrap_uncertainty = Uncertainty.uncertainty_extrapolation(meas, discharges)

            if np.abs(extrap_uncertainty) > 2:
                self.extrapolation['messages'].append(['Extrapolation: The extrapolation uncertainty is more than '
                                                       + '2 percent;', 2, 12])
                self.extrapolation['messages'].append(['    Carefully review the extrapolation;', 2, 12])
                self.extrapolation['status'] = 'caution'
        else:
            self.extrapolation['status'] = 'inactive'

    def edges_qa(self, meas):
        """Apply quality checks to edge estimates

        Parameters
        ----------
        meas: Measurement
            Object of class Measurement
        """

        # Intialize variables
        self.edges['messages'] = []
        checked = []
        left_q = []
        right_q = []
        total_q = []
        edge_dist_left = []
        edge_dist_right = []
        dist_moved_left = []
        dist_moved_right = []
        dist_made_good = []
        left_type = []
        right_type = []
        for n, transect in enumerate(meas.transects):
            checked.append(transect.checked)
            if transect.checked:
                left_q.append(meas.discharge[n].left)
                right_q.append(meas.discharge[n].right)
                total_q.append(meas.discharge[n].total)
                dmr, dml, dmg = QAData.edge_distance_moved(transect)
                dist_moved_right.append(dmr)
                dist_moved_left.append(dml)
                dist_made_good.append(dmg)
                edge_dist_left.append(transect.edges.left.distance_m)
                edge_dist_right.append(transect.edges.right.distance_m)
                left_type.append(transect.edges.left.type)
                right_type.append(transect.edges.right.type)

        if any(checked):
            # Set default status to good
            self.edges['status'] = 'good'

            # Check left edge q > 5%
            self.edges['left_q'] = 0
            left_q_percent = (np.nanmean(left_q) / np.nanmean(total_q)) * 100
            if np.abs(left_q_percent) > 5:
                self.edges['status'] = 'caution'
                self.edges['messages'].append(['Edges: Left edge Q is greater than 5%;', 2, 13])
                self.edges['left_q'] = 1

                # Check right edge q > 5%
                self.edges['right_q'] = 0
                right_q_percent = (np.nanmean(right_q) / np.nanmean(total_q)) * 100
                if np.abs(right_q_percent) > 5:
                    self.edges['status'] = 'caution'
                    self.edges['messages'].append(['Edges: Right edge Q is greater than 5%;', 2, 13])
                    self.edges['right_q'] = 1

                # Check for consistent sign
                q_positive = []
                self.edges['left_sign'] = 0
                for q in left_q:
                    if q >= 0:
                        q_positive.append(True)
                    else:
                        q_positive.append(False)
                if len(np.unique(q_positive)) > 1 and left_q_percent > 0.5:
                    self.edges['status'] = 'caution'
                    self.edges['messages'].append(['Edges: Sign of left edge Q is not consistent;', 2, 13])
                    self.edges['left_sign'] = 1

                q_positive = []
                self.edges['right_sign'] = 0
                for q in right_q:
                    if q >= 0:
                        q_positive.append(True)
                    else:
                        q_positive.append(False)
                if len(np.unique(q_positive)) > 1 and right_q_percent > 0.5:
                    self.edges['status'] = 'caution'
                    self.edges['messages'].append(['Edges: Sign of right edge Q is not consistent;', 2, 13])
                    self.edges['right_sign'] = 1

                # Check distance moved
                dmg_5_percent = 0.05 * np.nanmean(dist_made_good)
                avg_right_edge_dist = np.nanmean(edge_dist_right)
                right_threshold = np.nanmin([dmg_5_percent, avg_right_edge_dist])
                self.edges['right_dist_moved_idx'] = np.where(dist_moved_right > right_threshold)[0]
                if np.any(self.edges['right_dist_moved_idx']):
                    self.edges['status'] = 'caution'
                    self.edges['messages'].append(['Edges: Excessive boat movement in right edge ensembles;', 2, 13])

                avg_left_edge_dist = np.nanmean(edge_dist_left)
                left_threshold = np.nanmin([dmg_5_percent, avg_left_edge_dist])
                self.edges['left_dist_moved_idx'] = np.where(dist_moved_left > left_threshold)[0]
                if np.any(self.edges['left_dist_moved_idx']):
                    self.edges['status'] = 'caution'
                    self.edges['messages'].append(['Edges: Excessive boat movement in left edge ensembles;', 2, 13])

                # Check for edge ensembles marked invalid due to excluded distance
                for transect in meas.transects:
                    if transect.checked:
                        ens_sum_excluded_data = np.nansum(transect.w_vel.valid_data[6, :, :], 0)
                        cells_above_sl = np.nansum(transect.w_vel.cells_above_sl, 0)
                        ens_excluded_data = np.not_equal(ens_sum_excluded_data, cells_above_sl)
                        if any(ens_excluded_data):
                            self.edges['status'] = 'caution'
                            self.edges['messages'].append(['Edges: The excluded distance caused invalid ensembles '
                                                           + 'in an edge, check edge distance;', 2, 13])
                            break

                # Check edges for zero discharge
                self.edges['left_zero'] = 0
                left_zero_idx = np.where(left_q == 0)[0]
                if left_zero_idx:
                    self.edges['status'] = 'warning'
                    self.edges['messages'].append(['EDGES: Left edge has zero Q;', 1, 13])
                    self.edges['left_zero'] = 2

                self.edges['right_zero'] = 0
                right_zero_idx = np.where(right_q == 0)[0]
                if right_zero_idx:
                    self.edges['status'] = 'warning'
                    self.edges['messages'].append(['EDGES: Right edge has zero Q;', 1, 13])
                    self.edges['right_zero'] = 2

                # Check consistent edge type
                self.edges['left_type'] = 0
                if len(np.unique(left_type)) > 1:
                    self.edges['status'] = 'warning'
                    self.edges['messages'].append(['EDGES: Left edge type is not consistent;', 1, 13])
                    self.edges['left_type'] = 2

                self.edges['right_type'] = 0
                if len(np.unique(right_type)) > 1:
                    self.edges['status'] = 'warning'
                    self.edges['messages'].append(['EDGES: Right edge type is not consistent;', 1, 13])
                    self.edges['right_type'] = 2
        else:
            self.edges['status'] = 'inactive'

    @staticmethod
    def invalid_qa(valid, discharge):
        """Computes the total invalid discharge in ensembles that have invalid data. The function also computes
        the maximum run or cluster of ensembles with the maximum interpolated discharge.

        Parameters
        ----------
        valid: np.array(bool)
            Array identifying valid and invalid ensembles.
        discharge: QComp
            Object of class QComp

        Returns
        -------
        q_invalid_total: float
            Total interpolated discharge in invalid ensembles
        q_invalid_max_run: float
            Maximum interpolated discharge in a run or cluster of invalid ensembles
        ens_invalid: int
            Total number of invalid ensembles
        """

        # Create bool for invalid data
        invalid = np.logical_not(valid)
        q_invalid_total = np.nansum(discharge.middle_ens[invalid]) + np.nansum(discharge.top_ens[invalid]) \
            + np.nansum(discharge.bottom_ens[invalid])

        # Compute total number of invalid ensembles
        ens_invalid = np.sum(invalid)

        # Compute the indices of where changes occur

        valid_int = np.insert(valid.astype(int), 0, -1)
        valid_int = np.append(valid_int, -1)
        valid_run = np.where(np.diff(valid_int) != 0)[0]
        run_length = np.diff(valid_run)
        run_length0 = run_length[(valid[0] == 1)::2]

        n_runs = len(run_length0)

        if valid[0] is True:
            n_start = 1
        else:
            n_start = 0

        n_end = len(valid_run)-1

        if n_runs > 1:
            m = 0
            q_invalid_run = []
            for n in range(n_start, n_end, 2):
                m += 1
                idx_start = valid_run[n]
                idx_end = valid_run[n+1]
                q_invalid_run.append(np.nansum(discharge.middle_ens[idx_start:idx_end])
                                     + np.nansum(discharge.top_ens[idx_start:idx_end])
                                     + np.nansum(discharge.bottom_ens[idx_start:idx_end]))

            # Determine the maximum discharge in a single run
            q_invalid_max_run = np.nanmax(np.abs(q_invalid_run))

        else:
            q_invalid_max_run = 0

        return q_invalid_total, q_invalid_max_run, ens_invalid

    @staticmethod
    def edge_distance_moved(transect):
        """Computes the boat movement during edge ensemble collection.

        Parameters
        ----------
        transect: Transect
            Object of class Transect

        Returns
        -------
        right_dist_moved: float
            Distance in m moved during collection of right edge samples
        left_dist_moved: float
            Distance in m moved during collection of left edge samples
        dmg: float
            Distance made good for the entire transect
        """

        boat_selected = getattr(transect.boat_vel, transect.boat_vel.selected)
        ens_duration = transect.date_time.ens_duration_sec

        # Get boat velocities
        if boat_selected is not None:
            u_processed = boat_selected.u_processed_mps
            v_processed = boat_selected.v_processed_mps
        else:
            u_processed = np.tile(np.nan, transect.boat_vel.bt_vel.u_processed_mps.shape)
            v_processed = np.tile(np.nan, transect.boat_vel.bt_vel.v_processed_mps.shape)

        # Compute boat coordinates
        x_processed = np.nancumsum(u_processed * ens_duration)
        y_processed = np.nancumsum(v_processed * ens_duration)
        dmg = (x_processed[-1]**2 + y_processed[-1]**2)**0.5

        # Compute left distance moved
        # TODO should be a dist moved function
        left_edge_idx = QComp.edge_ensembles('left', transect)
        if len(left_edge_idx) > 0:
            boat_x = x_processed[left_edge_idx[-1]] - x_processed[left_edge_idx[0]]
            boat_y = y_processed[left_edge_idx[-1]] - y_processed[left_edge_idx[0]]
            left_dist_moved = (boat_x**2 + boat_y**2)**0.5
        else:
            left_dist_moved = np.nan

        # Compute right distance moved
        right_edge_idx = QComp.edge_ensembles('right', transect)
        if len(right_edge_idx) > 0:
            boat_x = x_processed[right_edge_idx[-1]] - x_processed[right_edge_idx[0]]
            boat_y = y_processed[right_edge_idx[-1]] - y_processed[right_edge_idx[0]]
            right_dist_moved = (boat_x ** 2 + boat_y ** 2) ** 0.5
        else:
            right_dist_moved = np.nan

        return right_dist_moved, left_dist_moved, dmg
