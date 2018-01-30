"""
Created on Jul 12, 2017

@author: gpetrochenkov
"""

import numpy as np
import re

class PT3(object):
    
    def __init__(self, hard_limit, linear = None):
        
        self.hard_limit = hard_limit
        if linear is not None:
            self.linear = linear
            
class CorrContainer(object):
    
    def __init__(self, hw, lw = None, hn = None, ln=None):
        
        self.hw = hw
        if lw is not None:
            self.lw = lw
        if hn is not None:
            self.hn = hn
        if ln is not None:
            self.ln = ln
            
class CorrData(object):
    
    def __init__(self, corr_table, noise_floor, sdc=None, cdc=None):
        
        self.corr_table = corr_table
        self.noise_floor = noise_floor
        
        if sdc is not None:
            self.sdc = sdc
            self.cdc = cdc
            
            

class PreMeasurement(object):
    
    def __init__(self):
        self.time_stamp = None
        self.data = None
        self.result = None
        self.compass = None
        self.sys_test = None
        self.pt3 = None
        self.error = None
        self.n_tests = None
        self.n_failed = None
        
        
    def populate_data(self, time_stamp, data_in, data_type):
        
        self.time_stamp = time_stamp
        self.data = data_in
        self.result = self.read_result(data_type)
        
        
    def read_result(self, data_type):
        """Function to aprse the data from the tests"""
        
        if data_type == 'TCC':
            self.compass_read(self.data)
        elif data_type == 'TCE':
            self.compass_read(self.data)
        elif data_type == 'TST':
            self.sys_test_read(self.data)
            self.pt3_data(self.data)
        elif data_type == 'SCC':
            self.compass_read(self.data)
        elif data_type == 'SST':
            self.sys_test_read(self.data)  
            
    def compass_read(self, data):
        """Method for getting compass evaluation data"""
        
        #match regex for compass evaluation error:
        splits = re.split('(Total error:|Double Cycle Errors:|Error from calibration:)', data)
        if len(splits) > 0:
            error = re.search('\d+\.*\d*')
            
            self.error = error
        else:
            self.error = 'N/A'
            
    def sys_test_read(self, data):
        """Method for getting the system test data"""
        
        #match regex for number of tests and number of failures
        num_tests = re.findall('(Fail|FAIL|F A I L|Pass|PASS|NOT DETECTED|P A S S)', data)
        num_fails = re.findall('(Fail|FAIL|F A I L)', data)
        
        #read in the numbers as doubles
        self.n_tests = len(num_tests)
        self.n_failed = len(num_fails)
        
    def pt3_data(self, data):
        """Method for getting the data for the correlation matrices"""
        
        #match regex for correlation tables
        match = re.findall('Lag.*?0', data)
        
        #count the number or correlation tables to read in
        correl_count = 0
        for x in range(len(match)):
            match2 = re.findall('Bm1', match[x])
            correl_count += len(match2)
            
        #Initialize correlation matricces
        corr_hlimit_hgain_wband = None
        corr_hlimit_lgain_wband = None
        corr_hlimit_hgain_nband = None
        corr_hlimit_lgain_nband = None
        corr_linear_hgain_wband = None
        corr_linear_lgain_wband = None
        corr_linear_hgain_nband = None
        corr_linear_lgain_nband = None
        
        #get matching string for correlation tables
        lag_match = re.findall('Lag.*?(High|H-Gain|Sin)', data)
        
        #Read depending on the number of correlation tables computed
        #at the beginning of the function
        if correl_count == 1:
            
            #Get strings to parse numbers from
            num_match = re.findall('0.*?(High|H-Gain|Sin)', lag_match[0])
            
            #get all numbers from matches
            numbers = re.findall('\d+\.*\d*', ''.join(num_match))
            
            #convert all numbers to double and add them to an array
            nums = []
            for x in range(len(numbers)):
                if x % 5 != 0:
                    nums.append(float(numbers[x]))
                    
        
            nums = np.array(nums)
            corr_hlimit_hgain_wband = nums.reshape([4,8]).T
            
        elif correl_count == 4:    
        
            #For all string in the correlation regular expression
            for x in range(len(lag_match)):
                
                #Count the Bm1 string to know how many tables to read
                bm_count = re.findall('Bm1', lag_match)
                
                #get string to parse numbers from
                num_match = re.findall('0.*?(High|H-Gain|Sin)', ''.join(bm_count))
                
                #get all numbers from matches
                numbers = re.findall('\d+\.*\d*', ''.join(num_match))
                
                if len(bm_count) == 2:
                    
                    #read in all strings as DOUBLESLASH
                    nums = []
                    for y in range(len(numbers)):
                        
                        if y % 9 != 0:
                            nums.append(float(numbers[y]))
                        
                    nums = np.array(nums)    
                    nums = nums.reshape(nums, [8,8]).T
                    
                    #assign matrix slices to corresponding variables
                    if corr_hlimit_hgain_wband is None:
                        corr_hlimit_hgain_wband = nums[:, 0:4]
                        corr_hlimit_lgain_wband = nums[:, 4:8]
                    else:
                        corr_hlimit_hgain_nband = nums[:, 0:4]
                        corr_hlimit_lgain_nband = nums[:, 4:8]
                        
                else:
                    
                    #read in all strings as doubles
                    nums = []
                    for y in range(len(numbers)):
                        
                        if y % 17 == 0:
                            nums.append(float(numbers[y]))
                            
                    #reshape numbers to appropriate sized matrix
                    nums = np.array(nums)
                    nums = nums.reshape(16,8).T
                    
                    #assign matrix slices to corresponding variables
                    corr_hlimit_hgain_wband = nums[:, 0:4]
                    corr_hlimit_lgain_wband = nums[:, 4:8]
                    corr_hlimit_hgain_nband = nums[:, 8:12]
                    corr_hlimit_lgain_nband = nums[:, 12:16]
                    
        else:
            #for all string in the correlation regular expresssion
            
            for x in range(len(lag_match)):
                
                #count the Bm1 strings to know how many tables to read
                bm_count = re.findall('Bm1', lag_match[x])
                
                #get string to parse numbers from
                num_match = re.findall('0.*?(High|H-Gain|Sin)', lag_match[x])
                
                #get all numbers from matches
                numbers = re.findall('\d+\.*\d*', ''.join(num_match))
                
                if len(bm_count) == 2:
                    
                    #read in all strings as doubles
                    nums = []
                    for y in range(len(numbers)):
                        if y % 9 != 0:
                            nums.append(float(numbers[y]))
                            
                    #reshape numbers to the appropriate sized matrix
                    nums = np.array(nums)
                    nums = nums.reshape([8,8]).T
                    
                    #assign matrix slices to corresponding variables
                    if corr_hlimit_hgain_wband is None:
                        corr_hlimit_hgain_wband = nums[:,0:4]
                        corr_hlimit_lgain_wband = nums[:,4:8]
                    elif corr_hlimit_hgain_nband is None:
                        corr_hlimit_hgain_nband = nums[:,0:4]
                        corr_hlimit_lgain_nband = nums[:,4:8]
                    elif corr_linear_hgain_wband is None:
                        corr_linear_hgain_wband = nums[:,0:4]
                        corr_linear_lgain_wband = nums[:,4:8]
                    else:
                        corr_linear_hgain_nband = nums[:,0:4]
                        corr_linear_lgain_nband = nums[:,4:8]
                        
                else:
                    
                    #read in all strings as doubles
                    nums = []
                    
                    for y in range(len(numbers)):
                        
                        if y % 17 != 0:
                            nums.append(float(numbers[y]))
                    
                    nums = np.array(nums)
                    nums.reshape([16,8]).T
                    
                    #Assign matrix slices to corresponding variables
                    if corr_hlimit_hgain_wband is None:
                        corr_hlimit_hgain_wband = numbers[:,0:4]
                        corr_hlimit_lgain_wband = numbers[:,4:8]
                        corr_hlimit_hgain_wband = numbers[:,8:12]
                        corr_hlimit_lgain_nband = numbers[:,12:16]
                    else:
                        corr_linear_hgain_wband = numbers[:,0:4]
                        corr_linear_lgain_wband = numbers[:,4:8]
                        corr_linear_hgain_nband = numbers[:,8:12]
                        corr_linear_lgain_nband = numbers[:,12:16]
                        
                        
        #Get all string for the noise floor data
        rssi_match = re.findall('(High Gain RSSI|RSSI Noise Floor \(counts\)|RSSI \(counts\)).*?(Low|>)', data)
        rssi_nums = []
        
        #Read in noise floor strings as doubles
        for x in range(len(rssi_match)):
            numbers = re.findall('\d+\.*\d*', rssi_match[x])
            for y in range(len(numbers)):
                rssi_nums.append(float(numbers[y]))
                
        #Get all strings for the sin duty cycle data
        sin_match = re.findall('(SIN Duty Cycle|Sin Duty Cycle \(percent\)|Sin Duty\(%\)).*?(COS|Cos)', data)
        sin_nums = []
        
        #read in sdc strings as doubles
        for x in range(len(sin_match)):
            numbers = re.findall('\d+\.*\d*', sin_match[x])
            for y in range(len(numbers)):
                sin_nums.append(float(numbers[y]))
                
        #Get all strings for the cos duty cycle data
        cos_match = re.findall('(COS Duty Cycle|Cos Duty Cycle \(percent\)|Cos Duty\(%\)).*?(Receive|RSSI)', data)
        cos_nums = []
        
        #read in sdc strings as doubles
        for x in range(len(cos_match)):
            numbers = re.findall('\d+\.*\d*', cos_match[x])
            for y in range(len(numbers)):
                cos_nums.append(float(numbers[y]))
                
        
        #8 Statements assigning data to pt3 dictionary if there is data
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
           
        #Create data container classes 
        hlimit = CorrContainer(hlim_hw, hlim_lw, hlim_hn, hlim_ln)
        linear = CorrContainer(lin_hw, lin_lw, lin_hn, lin_ln)
        
        #Create PT3 class
        self.pt3 = PT3(hlimit, linear)
        
              
        
                            