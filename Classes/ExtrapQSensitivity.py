'''
Created on Sep 26, 2017

@author: gpetrochenkov
'''


class ExtrapQSensitivity(object):
    '''Class to compute the sensitivity of the discharge to various extrapolation methods'''
    
    def __init__(self):
        
        self.__q_pp_mean = None #discharge power power 1/6
        self.__q_pp_optmean = None #discharge power power optimized
        self.__q_cns_mean = None #discharge constant no RoutingSlipDelivery
        self.__q_cns_optmean = None #discharge constant optimized no slip
        self.__q_3p_NS_mean = None #discharge 3-pt no slip
        self.__q_3p_NS_optmean = None #discharge 3-pt optimized no slip 
        self.__q_PP_per_diff = None #power power 1/6 difference from reference
        self.__q_PP_opt_per_diff = None #power power optimized percent difference from reference
        self.__q_cns_per_diff = None #constant no slip percent difference from reference
        self.__q_cns_opt_per_diff = None #constant optimized no slip percent difference from reference
        self.__q_3p_NS_per_diff = None #3-point no skip percent difference from reference
        self.__q_3p_ns_opt_per_diff = None #3-point optimized no slip percent difference from reference
        self.__pp_exponent = None #optimized power power exponent
        self.__ns_exponent = None #optimized no slip Exponent
        self.__man_top = None #manually specified top method
        self.__man_bot = None #manually specified bottom method
        self.__man_exp = None #manually specified exponent
        self.__q_nan_mean = None #mean discharge for manually specified extrapolations 
        self.__q_man_per_diff = None #manually specified extrapolations percent difference from reference
        
    def populate_data(self, trans_data, extrap_fits):
        pass
        