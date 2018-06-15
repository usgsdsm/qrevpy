import numpy as np
from Classes.Uncertainty import Uncertainty


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
        self.user = dict()
        self.depths = dict()
        self.boat = dict()
        self.water = dict()
        self.extrapolation = dict()
        self.edges = dict()

        # Apply QA checks
        self.transects_qa(meas)
        self.system_tst_qa(meas)
        self.compass_qa(meas)
        self.temperature_qa(meas)
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
            self.transects.duration = 1

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
            self.transects.number = 2
        elif num_checked == 1:
            # Only one transect selected
            self.transects['status'] = 'caution'
            self.transects['messages'].append(['Transects: Only one transect selected;', 2, 0])
            self.transects.number = 2
        else:
            self.transects.number = num_checked
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
            if np.all(q_positive) > 1:
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
                if test.result['pt3'] is not None:

                    # Check hard_limit, high gain, wide bandwidth
                    if 'hard_limit' in test.result['pt3']:
                        if 'high_wide' in test.result['pt3']['hard_limit']:
                            corr_table = test.result['pt3']['hard_limit']['high_wide']['corr_table']
                            # All lags past lag 2 should be less than 50% of lag 0
                            qa_threshold = corr_table[0, :] * 0.5
                            all_lag_check = np.greater(corr_table[3::, :], qa_threshold)

                            # Lag 7 should be less than 25% of lag 0
                            lag_7_check = np.greater(corr_table[7, :], corr_table[0, :] * 0.25)

                            # If either condition is met for any beam the test fails
                            if np.sum(np.sum(all_lag_check)) + np.sum(lag_7_check) > 1:
                                pt3_fail = True

                if test.result['n_failed'] > 0:
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

        heading = np.unique(meas.transects[checked.index(True)].sensors.heading_deg.internal.data)

        if len(heading) == 1 and heading == 0:
            # ADCP has no compass
            self.compass['status'] = 'inactive'
            self.compass['status1'] = 'good'
            self.compass['status2'] = 'good'
            self.compass['magvar'] = None
            self.compass['magvar_idx'] = None
        else:
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
                            if meas.compass_cal[-1].result['compass']['error'] <= 1
                                self.compass['status1'] = 'good'
                            else:
                                self.compass['status1'] = 'caution'
                                self.compass['messages'].append(['Compass: Evaluation result > 1 deg;', 2, 4])
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
                            idx_max = np.where(pitch_source_selected.data > pitch_source_selected.pitch_limit[0])[0]
                            idx_min = np.where(pitch_source_selected.data < pitch_source_selected.pitch_limit[1])[0]
                            if not idx_max or not idx_min:
                                pitch_exceeded.append(True)
                            else:
                                pitch_exceeded.append(False)

                        if heading_source_selected.roll_limit is not None:
                            idx_max = np.where(roll_source_selected.data > roll_source_selected.pitch_limit[0])[0]
                            idx_min = np.where(roll_source_selected.data < roll_source_selected.pitch_limit[1])[0]
                            if not idx_max or not idx_min:
                                roll_exceeded.append(True)
                            else:
                                roll_exceeded.append(False)

                        if heading_source_selected.mag_error is not None:
                            idx_max = np.where(heading_source_selected.mag_error > 2)[0]
                            if not idx_max:
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
            elif np.any(np.asarray(pitch_mean) > 4 ):
                if self.compass['status2'] == 'good':
                    self.compass['status2'] = 'caution'
                self.compass['messages'].append(['Pitch: One or more transects have a mean pitch > 4 deg;', 2, 4])

            # Check roll mean
            if np.any(np.asarray(roll_mean) > 8):
                self.compass['status2'] = 'warning'
                self.compass['messages'].append(['ROLL: One or more transects have a mean roll > 8 deg;', 1, 4])
            elif np.any(np.asarray(roll_mean) > 4 ):
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

        # Create array of all temperatures
        temp = []
        for transect in meas.transects:
            if transect.checked:
                temp_selected = getattr(transect.sensors.temperature_deg_c, transect.sensors.temperature_deg_c.selected)
                temp.append(temp_selected.data)
        np.asarray(temp)

        # Check temperature range
        temp_range = np.max(temp) - np.min(temp)
        if temp_range > 2:
            check = [3]
            self.temperature['messages'].append(['TEMPERATURE: Temperature range is '
                                            + '%3.1f % temp_range'
                                            +  'degrees C which is greater than 2 degrees;', 1, 5])
        elif temp_range > 1:
            check = [2]
            self.temperature['messages'].append(['TEMPERATURE: Temperature range is '
                                                 + '%3.1f % temp_range'
                                                 + 'degrees C which is greater than 1 degrees;', 2, 5])
        else:
            check = [1]

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
                if  diff < 2:
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