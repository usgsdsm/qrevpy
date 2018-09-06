"""
The module consists of multiple classes.

Classes
-------
SystemTest:
    The container for all other classes in this module.
PT3:
    Contains the results of PT3 tests for TRDI ADCPs.
CorrContainer:
    Hold the various correlation test results.
CorrData
    Contines the correlation table and noise level data for a test.

Most of this code was original programmed by gpetrochenkov in a PreMeasurement class.
DSM split that class into CompassCal and SystemTest. DSM made the code modifications, added
numpy docstrings and cleaned up PEP8 recommendations.
"""

import numpy as np
import re


class SystemTest(object):
    """Class to store all system test data.

    Attributes
    ----------
    time_stamp: str
        Time stamp of test (hh:mm:ss)
    data: str
        All output from the test.
    number_of_tests: int
        Number of individual tests checked.
    number_failed: int
        Number of individual tests that failed.
    emi_tests: PT3
        Object of PT3 tests.
    """
    def __init__(self):
        """Initialize class and instance variables."""
        self.time_stamp = None
        self.data = None
        self.number_of_tests = None
        self.number_failed = None
        self.emi_tests = None

    def populate_data(self, time_stamp, data_in):
        """Store data in SystemTest.

        Parameters
        ----------
        time_stamp: str
            Time stamp of test data.
        data_in: str
            All output from the test.
        """

        self.time_stamp = time_stamp
        self.data = data_in
        self.number_of_tests, self.number_failed = SystemTest.sys_test_read(data_in)
        # PT3 tests are only found in TRDI data so test for presence of PT3 test.
        if data_in.find('PT3') > 0:
            # Some TRDI test files have been observed to be incomplete. Thus, a
            # try / except is used so that processing can continue even if the
            # PT3 test data is corrupted.
            try:
                self.emi_tests = PT3()
                self.emi_tests.populate_data(data_in)
            except ValueError:
                pass

    @staticmethod
    def sys_test_read(data):
        """Method for getting the system test data

        Parameters
        ----------
            data: str
                Test data from system test.

        Returns
        -------
            n_test: int
                Number of tests completed.
            n_failed: int
                Number of tests that failed."""

        # match regex for number of tests and number of failures
        num_tests = re.findall('(Fail|FAIL|F A I L|Pass|PASS|NOT DETECTED|P A S S)', data)
        num_fails = re.findall('(Fail|FAIL|F A I L)', data)

        # read in the numbers as doubles
        n_tests = len(num_tests)
        n_failed = len(num_fails)

        return n_tests, n_failed


class PT3(object):
    """Class to store PT3 test data.

    Attributes
    ----------
    hard_limit: CorrContainer
        Data from hard limit tests.
    linear: CorrContainer
        Data from linear tests.
    """

    def __init__(self):
        """Initialize object and instance variables."""

        self.hard_limit = None
        self.linear = None

    def populate_data(self, data):
        """Method for getting the data for the correlation matrices.

        Parameters
        ----------
        data: str
            All test data.
        """

        # match regex for correlation tables
        match = re.findall('Lag.*?0', data)

        # count the number or correlation tables to read in
        correl_count = 0
        for x in range(len(match)):
            match2 = re.findall('Bm1', match[x])
            correl_count += len(match2)

        # Initialize correlation matricces
        corr_hlimit_hgain_wband = None
        corr_hlimit_lgain_wband = None
        corr_hlimit_hgain_nband = None
        corr_hlimit_lgain_nband = None
        corr_linear_hgain_wband = None
        corr_linear_lgain_wband = None
        corr_linear_hgain_nband = None
        corr_linear_lgain_nband = None

        # Get matching string for correlation tables
        lag_match = re.findall('Lag.*?(High|H-Gain|Sin)', data)

        # Read depending on the number of correlation tables computed
        # at the beginning of the function
        if correl_count == 1:

            # Get strings to parse numbers from
            num_match = re.findall('0.*?(High|H-Gain|Sin)', lag_match[0])

            # Get all numbers from matches
            numbers = re.findall('\d+\.*\d*', ''.join(num_match))

            # Convert all numbers to double and add them to an array
            nums = []
            for x in range(len(numbers)):
                if x % 5 != 0:
                    nums.append(float(numbers[x]))

            nums = np.array(nums)
            corr_hlimit_hgain_wband = nums.reshape([4, 8]).T

        elif correl_count == 4:

            # For all string in the correlation regular expression
            for x in range(len(lag_match)):

                # Count the Bm1 string to know how many tables to read
                # DSM added [x] 2/3/2018
                bm_count = re.findall('Bm1', lag_match[x])

                # Get string to parse numbers from
                num_match = re.findall('0.*?(High|H-Gain|Sin)', ''.join(bm_count))

                # Get all numbers from matches
                numbers = re.findall('\d+\.*\d*', ''.join(num_match))

                if len(bm_count) == 2:

                    # read in all strings as DOUBLESLASH
                    nums = []
                    for y in range(len(numbers)):

                        if y % 9 != 0:
                            nums.append(float(numbers[y]))

                    nums = np.array(nums)
                    nums = nums.reshape(nums, [8, 8]).T

                    # assign matrix slices to corresponding variables
                    if corr_hlimit_hgain_wband is None:
                        corr_hlimit_hgain_wband = nums[:, 0:4]
                        corr_hlimit_lgain_wband = nums[:, 4:8]
                    else:
                        corr_hlimit_hgain_nband = nums[:, 0:4]
                        corr_hlimit_lgain_nband = nums[:, 4:8]

                else:

                    # read in all strings as doubles
                    nums = []
                    for y in range(len(numbers)):

                        if y % 17 == 0:
                            nums.append(float(numbers[y]))

                    # reshape numbers to appropriate sized matrix
                    nums = np.array(nums)
                    nums = nums.reshape(16, 8).T

                    # assign matrix slices to corresponding variables
                    corr_hlimit_hgain_wband = nums[:, 0:4]
                    corr_hlimit_lgain_wband = nums[:, 4:8]
                    corr_hlimit_hgain_nband = nums[:, 8:12]
                    corr_hlimit_lgain_nband = nums[:, 12:16]

        else:
            # for all string in the correlation regular expresssion

            for x in range(len(lag_match)):

                # Count the Bm1 strings to know how many tables to read
                bm_count = re.findall('Bm1', lag_match[x])

                # Get string to parse numbers from
                num_match = re.findall('0.*?(High|H-Gain|Sin)', lag_match[x])

                # Get all numbers from matches
                numbers = re.findall('\d+\.*\d*', ''.join(num_match))

                if len(bm_count) == 2:

                    # Read in all strings as doubles
                    nums = []
                    for y in range(len(numbers)):
                        if y % 9 != 0:
                            nums.append(float(numbers[y]))

                    # Reshape numbers to the appropriate sized matrix
                    nums = np.array(nums)
                    nums = nums.reshape([8, 8]).T

                    # assign matrix slices to corresponding variables
                    if corr_hlimit_hgain_wband is None:
                        corr_hlimit_hgain_wband = nums[:, 0:4]
                        corr_hlimit_lgain_wband = nums[:, 4:8]
                    elif corr_hlimit_hgain_nband is None:
                        corr_hlimit_hgain_nband = nums[:, 0:4]
                        corr_hlimit_lgain_nband = nums[:, 4:8]
                    elif corr_linear_hgain_wband is None:
                        corr_linear_hgain_wband = nums[:, 0:4]
                        corr_linear_lgain_wband = nums[:, 4:8]
                    else:
                        corr_linear_hgain_nband = nums[:, 0:4]
                        corr_linear_lgain_nband = nums[:, 4:8]

                else:

                    # read in all strings as doubles
                    nums = []

                    for y in range(len(numbers)):

                        if y % 17 != 0:
                            nums.append(float(numbers[y]))

                    nums = np.array(nums)
                    nums = nums.reshape([16, 8]).T

                    # Assign matrix slices to corresponding variables
                    if corr_hlimit_hgain_wband is None:
                        corr_hlimit_hgain_wband = nums[:, 0:4]
                        corr_hlimit_lgain_wband = nums[:, 4:8]
                        corr_hlimit_hgain_nband = nums[:, 8:12]
                        corr_hlimit_lgain_nband = nums[:, 12:16]
                    else:
                        corr_linear_hgain_wband = nums[:, 0:4]
                        corr_linear_lgain_wband = nums[:, 4:8]
                        corr_linear_hgain_nband = nums[:, 8:12]
                        corr_linear_lgain_nband = nums[:, 12:16]

        # Get all string for the noise floor data
        rssi_match = re.findall('(High Gain RSSI|RSSI Noise Floor \(counts\)|RSSI \(counts\)).*?(Low|>)', data)
        rssi_nums = []

        # Read in noise floor strings as doubles
        for x in range(len(rssi_match)):
            numbers = re.findall('\d+\.*\d*', rssi_match[x])
            for y in range(len(numbers)):
                rssi_nums.append(float(numbers[y]))

        # Get all strings for the sin duty cycle data
        sin_match = re.findall('(SIN Duty Cycle|Sin Duty Cycle \(percent\)|Sin Duty\(%\)).*?(COS|Cos)', data)
        sin_nums = []

        # read in sdc strings as doubles
        for x in range(len(sin_match)):
            numbers = re.findall('\d+\.*\d*', sin_match[x])
            for y in range(len(numbers)):
                sin_nums.append(float(numbers[y]))

        # Get all strings for the cos duty cycle data
        cos_match = re.findall('(COS Duty Cycle|Cos Duty Cycle \(percent\)|Cos Duty\(%\)).*?(Receive|RSSI)', data)
        cos_nums = []

        # Read in sdc strings as doubles
        for x in range(len(cos_match)):
            numbers = re.findall('\d+\.*\d*', cos_match[x])
            for y in range(len(numbers)):
                cos_nums.append(float(numbers[y]))

        # Create CorrData objects for each type of test.
        if corr_hlimit_hgain_wband is not None:
            hlim_hw = CorrData(corr_hlimit_hgain_wband, rssi_nums[0:4], sin_nums[0:4], cos_nums[0:4])
        else:
            hlim_hw = None

        if corr_hlimit_lgain_wband is not None:
            hlim_lw = CorrData(corr_hlimit_lgain_wband, rssi_nums[4:8], sin_nums[4:8], cos_nums[4:8])
        else:
            hlim_lw = None

        if corr_hlimit_hgain_nband is not None:
            hlim_hn = CorrData(corr_hlimit_hgain_nband, rssi_nums[8:12], sin_nums[8:12], cos_nums[8:12])
        else:
            hlim_hn = None

        if corr_hlimit_lgain_nband is not None:
            hlim_ln = CorrData(corr_hlimit_lgain_nband, rssi_nums[12:16], sin_nums[12:16], cos_nums[12:16])
        else:
            hlim_ln = None

        if corr_linear_hgain_wband is not None:
            lin_hw = CorrData(corr_linear_hgain_wband, rssi_nums[16:20])
        else:
            lin_hw = None

        if corr_linear_lgain_wband is not None:
            lin_lw = CorrData(corr_linear_lgain_wband, rssi_nums[20:24])
        else:
            lin_lw = None

        if corr_linear_hgain_nband is not None:
            lin_hn = CorrData(corr_linear_hgain_nband, rssi_nums[24:28])
        else:
            lin_hn = None

        if corr_linear_lgain_nband is not None:
            lin_ln = CorrData(corr_linear_lgain_nband, rssi_nums[28:32])
        else:
            lin_ln = None

        # Create data container classes
        self.hard_limit = CorrContainer(hlim_hw, hlim_lw, hlim_hn, hlim_ln)
        self.linear = CorrContainer(lin_hw, lin_lw, lin_hn, lin_ln)


class CorrContainer(object):
    """Stores data for various types of PT3 tests."""

    def __init__(self, hw, lw=None, hn=None, ln=None):
        """Initialize class and instance variables."""

        self.hw = hw
        if lw is not None:
            self.lw = lw
        if hn is not None:
            self.hn = hn
        if ln is not None:
            self.ln = ln


class CorrData(object):
    """Store results of a PT3 test."""

    def __init__(self, corr_table, noise_floor, sdc=None, cdc=None):
        """Initalize class and assign instance variables.

        Parameters
        ----------
        corr_table: np.array(float)
            Correlation table from PT3 test.
        noise_floor: np.array(float)
            Noise floor values for each beam.
        sdc: np.array(float)
            Sine duty cycle for each beam.
        cdc: np.array(float)
            Cosine duty cycle for each beam.
        """

        self.corr_table = corr_table
        self.noise_floor = noise_floor

        if sdc is not None:
            self.sdc = sdc
            self.cdc = cdc
