'''
Created on Sep 26, 2017

@author: gpetrochenkov
'''

class SelectFit(object):
    '''Class definition for data class that contains all of the automated
    extrapolation method selection information. This inherits all the
    properties from FitData and addes the automated selection data to it.
    David S. Mueller 6/17/2011
    
    
    Last modificaitons / validation 5/15/2012
    7/24/2013
    DSM modified to catch C shaped profile lines 63-66, 123
    DSM modified to fix manual fits and how it displays line 172
    7/25/2013
    DSM fixed bug in line 172, changed 0.1667 to obj.exponent
    
    4/4/2014
    DSM modified for use in QRev'''
    
    def __init__(self):
        self.__fit_method = None #User selected method Automatic or Manual
        self.__bot_method_auto = None #Selected extrapolation for top
        self.__top_method_auto = None #Selected extrapolation for bottom
        self.__exponent_auto = None #Selected exponent
        self.__top_fit_r2 = None #Top fit custom r^2
        self.__top_max_diff = None #Maximum difference between power and 3-pt at top
        self.__bot_diff = None #Difference between power and no slop at z = 0.1
        self.__bot_r2 = None #Bottom fit r^2
        self.__fit_r2 = None #Selected fit of selected power/no slip fit
        self.__ns_exponent = None # No slip optimized exponent
        self.__pp_exponent = None #Power Power optimized exponent
        self.__top_r2 = None
        
    
