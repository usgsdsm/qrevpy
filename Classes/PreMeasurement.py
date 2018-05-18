import numpy as np
import re
import copy


class PreMeasurement(object):
    
    def __init__(self):
        self.time_stamp = None
        self.data = None
        self.result = {}
        
    def populate_data(self, time_stamp, data_in, data_type):

        # Store time stamp and data
        self.time_stamp = time_stamp
        self.data = data_in

        # Process data depending on data type and store result
        if data_type[1] == 'C':
            self.compass_read()
        # elif data_type == 'TCE':
        #     self.compass_read()
        elif data_type == 'TST':
            self.sys_test_read()
            self.pt3_data()
        # elif data_type == 'SCC':
        #     self.compass_read()
        elif data_type == 'SST':
            self.sys_test_read()

    def compass_read(self):
        """Method for getting compass evaluation data"""

        # Match regex for compass evaluation error:
        splits = re.split('(Total error:|Double Cycle Errors:|Error from calibration:)', self.data)
        if len(splits) > 1:
            error = re.search('\d+\.*\d*', splits[-1])[0]
        else:
            error = 'N/A'
        self.result['compass'] = {'error': error}
            
    def sys_test_read(self):
        """Method for reading the system test data"""
        
        # Match regex for number of tests and number of failures
        num_tests = re.findall('(Fail|FAIL|F A I L|Pass|PASS|NOT DETECTED|P A S S)', self.data)
        num_fails = re.findall('(Fail|FAIL|F A I L)', self.data)
        
        # Store results
        self.result['n_tests'] = len(num_tests)
        self.result['n_failed'] = len(num_fails)
        
    def pt3_data(self):
        """Method for processing the data in the correlation matrices."""
        data_types = {'corr_table': np.array([]), 'sdc': np.array([]), 'cdc': np.array([]), 'noise_floor': np.array([])}
        test_types = {'high_wide': data_types.copy(), 'high_narrow': data_types.copy(), 'low_wide': data_types.copy(),
                      'low_narrow': data_types.copy()}
        pt3 = {'hard_limit': copy.deepcopy(test_types), 'linear': copy.deepcopy(test_types)}
        # pt3 = {'hard_limit':
        #            {'high_wide': data_types.copy(),
        #             'high_narrow': data_types.copy(),
        #             'low_wide': data_types.copy(),
        #             'low_narrow': data_types.copy()},
        #        'linear':
        #            {'high_wide': data_types.copy(),
        #             'high_narrow': data_types.copy(),
        #             'low_wide': data_types.copy(),
        #             'low_narrow': data_types.copy()}}

        # Match regex for correlation tables
        matches = re.findall('Lag.*?0', self.data, re.DOTALL)

        # Count the number or correlation tables to process
        correl_count = 0
        for match in matches:
            bm1_matches = re.findall('Bm1', match)
            correl_count += len(bm1_matches)
            
        # # Initialize correlation matrices
        # corr_hlimit_hgain_wband = None
        # corr_hlimit_lgain_wband = None
        # corr_hlimit_hgain_nband = None
        # corr_hlimit_lgain_nband = None
        # corr_linear_hgain_wband = None
        # corr_linear_lgain_wband = None
        # corr_linear_hgain_nband = None
        # corr_linear_lgain_nband = None
        
        # Correlation table match
        lag_matches = re.findall('Lag.*?^\s*$', self.data, re.MULTILINE | re.DOTALL)

        # Sin match
        sin_match = re.findall('((Sin|SIN).*?^\s*$)', self.data, re.MULTILINE | re.DOTALL)[0][0]
        sin_array = np.array(re.findall('\d+\.*\d*', sin_match))

        # Cos match
        cos_match = re.findall('((Cos|COS).*?^\s*$)', self.data, re.MULTILINE | re.DOTALL)[0][0]
        cos_array = np.array(re.findall('\d+\.*\d*', cos_match))

        # RSSI match
        rssi_array = np.array([])
        rssi_matches = re.findall('RSSI.*?^\s*$', self.data, re.MULTILINE | re.DOTALL)
        for rssi_match in rssi_matches:
            rssi_array = np.hstack((rssi_array, np.array(re.findall('\d+\.*\d*', rssi_match))))

        # Process each set of correlation tables
        for n, lag_match in enumerate(lag_matches):

            # Count the Bm1 string to know how many tables to read
            bm_count = len(re.findall('Bm1', lag_match))

            # Extract the table into list
            numbers = re.findall('\d+\.*\d*', lag_match)

            # Create array from data in table
            corr_data = np.array(numbers[(bm_count * 4):(bm_count * 44)]).reshape([8, (bm_count * 4) + 1])[:, 1::]

            # Only one pt3 test. Typical of Rio Grande and Streampro
            if bm_count == 1:

                # Assign matrix slices to corresponding variables
                # corr_hlimit_hgain_wband = corr_data
                pt3['hard_limit']['high_wide']['corr_table'] = corr_data
                pt3['hard_limit']['high_wide']['sdc'] = sin_array[0:4]
                pt3['hard_limit']['high_wide']['cdc'] = cos_array[0:4]
                pt3['hard_limit']['high_wide']['noise_floor'] = rssi_array

            # 4 tests arranged in groups of 2. All data are hard limited.
            elif bm_count == 2 and correl_count == 4:

                # Hard limited wide bandwidth (n=0)
                if n == 0:

                    pt3['hard_limit']['high_wide']['corr_table'] = corr_data[:, 0:4]
                    pt3['hard_limit']['high_wide']['sdc'] = sin_array[n * 4: (n + 1) * 4]
                    pt3['hard_limit']['high_wide']['cdc'] = cos_array[n * 4: (n + 1) * 4]
                    pt3['hard_limit']['high_wide']['noise_floor'] = rssi_array[n * 4: (n + 1) * 4]

                    pt3['hard_limit']['low_wide']['corr_table'] = corr_data[:, 4::]
                    pt3['hard_limit']['low_wide']['sdc'] = sin_array[(n + 1) * 4: (n + 2) * 4]
                    pt3['hard_limit']['low_wide']['cdc'] = cos_array[(n + 1) * 4: (n + 2) * 4]
                    pt3['hard_limit']['low_wide']['noise_floor'] = rssi_array[(n + 1) * 4: (n + 2) * 4]

                # Hard limited narrow bandwidth (n=1)
                elif n == 1:

                    pt3['hard_limit']['high_narrow']['corr_table'] = corr_data[:, 0:4]
                    pt3['hard_limit']['high_narrow']['sdc'] = sin_array[(n + 1) * 4: (n + 2) * 4]
                    pt3['hard_limit']['high_narrow']['cdc'] = cos_array[(n + 1) * 4: (n + 2) * 4]
                    pt3['hard_limit']['high_narrow']['noise_floor'] = rssi_array[(n + 1) * 4: (n + 2) * 4]

                    pt3['hard_limit']['low_narrow']['corr_table'] = corr_data[:, 4::]
                    pt3['hard_limit']['low_narrow']['sdc'] = sin_array[(n + 2) * 4: (n + 3) * 4]
                    pt3['hard_limit']['low_narrow']['cdc'] = cos_array[(n + 2) * 4: (n + 3) * 4]
                    pt3['hard_limit']['low_narrow']['noise_floor'] = rssi_array[(n + 2) * 4: (n + 3) * 4]

            # 8 tests arranged in sets of 2. The linear is 1st followed by the hard limit.
            elif bm_count == 2 and correl_count == 8:

                # Hard limit bandwidth (n=0)
                if n == 0:

                    pt3['hard_limit']['high_wide']['corr_table'] = corr_data[:, 0:4]
                    pt3['hard_limit']['high_wide']['sdc'] = sin_array[n * 4: (n + 1) * 4]
                    pt3['hard_limit']['high_wide']['cdc'] = cos_array[n * 4: (n + 1) * 4]
                    pt3['hard_limit']['high_wide']['noise_floor'] = rssi_array[n * 4: (n + 1) * 4]

                    pt3['hard_limit']['low_wide']['corr_table'] = corr_data[:, 4::]
                    pt3['hard_limit']['low_wide']['sdc'] = sin_array[(n + 1) * 4: (n + 2) * 4]
                    pt3['hard_limit']['low_wide']['cdc'] = cos_array[(n + 1) * 4: (n + 2) * 4]
                    pt3['hard_limit']['low_wide']['noise_floor'] = rssi_array[(n + 1) * 4: (n + 2) * 4]

                # Hard limit narrow bandwidth (n=1)
                elif n == 1:

                    pt3['hard_limit']['high_narrow']['corr_table'] = corr_data[:, 0:4]
                    pt3['hard_limit']['high_narrow']['sdc'] = sin_array[(n + 1) * 4: (n + 2) * 4]
                    pt3['hard_limit']['high_narrow']['cdc'] = cos_array[(n + 1) * 4: (n + 2) * 4]
                    pt3['hard_limit']['high_narrow']['noise_floor'] = rssi_array[(n + 1) * 4: (n + 2) * 4]

                    pt3['hard_limit']['low_narrow']['corr_table'] = corr_data[:, 4::]
                    pt3['hard_limit']['low_narrow']['sdc'] = sin_array[(n + 2) * 4: (n + 3) * 4]
                    pt3['hard_limit']['low_narrow']['cdc'] = cos_array[(n + 2) * 4: (n + 3) * 4]
                    pt3['hard_limit']['low_narrow']['noise_floor'] = rssi_array[(n + 2) * 4: (n + 3) * 4]

                # Linear wide bandwidth (n=2)
                elif n == 2:

                    pt3['linear']['high_wide']['corr_table'] = corr_data[:, 0:4]
                    pt3['linear']['high_wide']['noise_floor'] = rssi_array[(n + 2) * 4: (n + 3) * 4]

                    pt3['linear']['low_wide']['corr_table'] = corr_data[:, 4::]
                    pt3['linear']['low_wide']['noise_floor'] = rssi_array[(n + 3) * 4: (n + 4) * 4]

                # Linear narrow bandwidth (n=3)
                elif n == 3:

                    pt3['linear']['high_narrow']['corr_table'] = corr_data[:, 0:4]
                    pt3['linear']['high_narrow']['noise_floor'] = rssi_array[(n + 3) * 4: (n + 4) * 4]

                    pt3['linear']['low_narrow']['corr_table'] = corr_data[:, 4::]
                    pt3['linear']['low_narrow']['noise_floor'] = rssi_array[(n + 4) * 4: (n + 5) * 4]

            # 8 tests in groups of 4. Hard limit is the first group then the linear.
            elif bm_count == 4:

                # Hard limit data (n=0)
                if n == 0:

                    pt3['hard_limit']['high_wide']['corr_table'] = corr_data[:, 0:4]
                    pt3['hard_limit']['high_wide']['sdc'] = sin_array[n * 4: (n + 1) * 4]
                    pt3['hard_limit']['high_wide']['cdc'] = cos_array[n * 4: (n + 1) * 4]
                    pt3['hard_limit']['high_wide']['noise_floor'] = rssi_array[n * 4: (n + 1) * 4]

                    pt3['hard_limit']['low_wide']['corr_table'] = corr_data[:, 4:8]
                    pt3['hard_limit']['low_wide']['sdc'] = sin_array[(n + 1) * 4: (n + 2) * 4]
                    pt3['hard_limit']['low_wide']['cdc'] = cos_array[(n + 1) * 4: (n + 2) * 4]
                    pt3['hard_limit']['low_wide']['noise_floor'] = rssi_array[(n + 1) * 4: (n + 2) * 4]

                    pt3['hard_limit']['high_narrow']['corr_table'] = corr_data[:, 8:12]
                    pt3['hard_limit']['high_narrow']['sdc'] = sin_array[(n + 2) * 4: (n + 3) * 4]
                    pt3['hard_limit']['high_narrow']['cdc'] = cos_array[(n + 2) * 4: (n + 3) * 4]
                    pt3['hard_limit']['high_narrow']['noise_floor'] = rssi_array[(n + 2) * 4: (n + 3) * 4]

                    pt3['hard_limit']['low_narrow']['corr_table'] = corr_data[:, 12::]
                    pt3['hard_limit']['low_narrow']['sdc'] = sin_array[(n + 3) * 4: (n + 4) * 4]
                    pt3['hard_limit']['low_narrow']['cdc'] = cos_array[(n + 3) * 4: (n + 4) * 4]
                    pt3['hard_limit']['low_narrow']['noise_floor'] = rssi_array[(n + 3) * 4: (n + 4) * 4]

                # Linear data (n=1)
                else:
                    pt3['linear']['high_wide']['corr_table'] = corr_data[:, 0:4]
                    pt3['linear']['high_wide']['noise_floor'] = rssi_array[(n + 3) * 4: (n + 4) * 4]

                    pt3['linear']['low_wide']['corr_table'] = corr_data[:, 4:8]
                    pt3['linear']['low_wide']['noise_floor'] = rssi_array[(n + 4) * 4: (n + 5) * 4]

                    pt3['linear']['high_narrow']['corr_table'] = corr_data[:, 8:12]
                    pt3['linear']['high_narrow']['noise_floor'] = rssi_array[(n + 5) * 4: (n + 6) * 4]

                    pt3['linear']['low_narrow']['corr_table'] = corr_data[:, 12::]
                    pt3['linear']['low_narrow']['noise_floor'] = rssi_array[(n + 6) * 4: (n + 7) * 4]
        self.result['pt3'] = pt3
