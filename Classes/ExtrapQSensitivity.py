<<<<<<< HEAD
from Classes.QComp import QComp
import numpy as np
from MiscLibs.common_functions import get_object_values
=======
'''
Created on Sep 26, 2017

@author: gpetrochenkov
'''

>>>>>>> 6ca6c50c231afa610ed3a693864074d7104a5f20

class ExtrapQSensitivity(object):
    """Class to compute the sensitivity of the discharge to various extrapolation methods.

    Attributes
    ----------
    q_pp_mean: float
        Discharge power power 1/6
    q_pp_opt_mean: float
        Discharge power power optimized
    q_cns_mean: float
        Discharge constant no RoutingSlipDelivery
    q_cns_opt_mean: float
        Discharge constant optimized no slip
    q_3p_ns_mean: float
        Discharge 3-pt no slip
    q_3p_ns_opt_mean: float
        Discharge 3-pt optimized no slip
    q_PP_per_diff: float
        Power power 1/6 difference from reference
    q_PP_opt_per_diff: float
        Power power optimized percent difference from reference
    q_cns_per_diff: float
        Constant no slip percent difference from reference
    q_cns_opt_per_diff: float
        Constant optimized no slip percent difference from reference
    q_3p_ns_per_diff: float
        3-point no skip percent difference from reference
    q_3p_ns_opt_per_diff: float
        3-point optimized no slip percent difference from reference
    pp_exponent: float
        Optimized power power exponent
    ns_exponent: float
        Optimized no slip Exponent
    man_top: str
        Manually specified top method
    man_bot: str
        Manually specified bottom method
    man_exp: float
        Manually specified exponent
    q_man_mean: float
        Mean discharge for manually specified extrapolations
    q_man_per_diff: float
        Manually specified extrapolations percent difference from reference
    """
    
    def __init__(self):
        """Initialize object and instance variables."""

        self.q_pp_mean = None  # Discharge power power 1/6
        self.q_pp_opt_mean = None  # discharge power power optimized
        self.q_cns_mean = None  # Discharge constant no RoutingSlipDelivery
        self.q_cns_opt_mean = None  # Discharge constant optimized no slip
        self.q_3p_ns_mean = None  # Discharge 3-pt no slip
        self.q_3p_ns_opt_mean = None  # Discharge 3-pt optimized no slip
        self.q_PP_per_diff = None  # Power power 1/6 difference from reference
        self.q_PP_opt_per_diff = None  # Power power optimized percent difference from reference
        self.q_cns_per_diff = None  # Constant no slip percent difference from reference
        self.q_cns_opt_per_diff = None  # Constant optimized no slip percent difference from reference
        self.q_3p_nS_per_diff = None  # 3-point no skip percent difference from reference
        self.q_3p_ns_opt_per_diff = None  # 3-point optimized no slip percent difference from reference
        self.pp_exponent = None  # Optimized power power exponent
        self.ns_exponent = None  # Optimized no slip Exponent
        self.man_top = None  # Manually specified top method
        self.man_bot = None  # Manually specified bottom method
        self.man_exp = None  # Manually specified exponent
        self.q_man_mean = None  # Mean discharge for manually specified extrapolations
        self.q_man_per_diff = None  # Manually specified extrapolations percent difference from reference
        
<<<<<<< HEAD
    def populate_data(self, transects, extrap_fits):
        """Compute means and percent differences.

        Parameters
        ----------
        transects: list
            List of objects of TransectData
        extrap_fits: object
            Object of SelectFit
        """
        q_pp = []
        q_pp_opt = []
        q_cns = []
        q_cns_opt = []
        q_3p_ns = []
        q_3p_ns_opt = []
        checked = []
        self.pp_exponent = extrap_fits[-1].pp_exponent
        self.ns_exponent = extrap_fits[-1].ns_exponent

        # Compute discharges for each transect for possible extrapolation combinations
        for transect in transects:
            checked.append(transect.checked)
            q = QComp()

            q.populate_data(data_in=transect, top_method='Power', bot_method='Power', exponent=0.1667)
            q_pp.append(q)

            q.populate_data(data_in=transect, top_method='Power', bot_method='Power', exponent=self.pp_exponent)
            q_pp_opt.append(q)

            q.populate_data(data_in=transect, top_method='Constant', bot_method='No Slip', exponent=0.1667)
            q_cns.append(q)

            q.populate_data(data_in=transect, top_method='Constant', bot_method='No Slip', exponent=self.ns_exponent)
            q_cns_opt.append(q)

            q.populate_data(data_in=transect, top_method='3-Point', bot_method='No Slip', exponent=0.1667)
            q_3p_ns.append(q)

            q.populate_data(data_in=transect, top_method='3-Point', bot_method='No Slip', exponent=self.ns_exponent)
            q_3p_ns_opt.append(q)

        # Compute mean discharge for each combination
        self.q_pp_mean = np.nanmean(get_object_values(q_pp, 'total', checked))
        self.q_pp_opt_mean = np.nanmean(get_object_values(q_pp_opt, 'total', checked))
        self.q_cns_mean = np.nanmean(get_object_values(q_cns, 'total', checked))
        self.q_cns_opt_mean = np.nanmean(get_object_values(q_cns_opt, 'total', checked))
        self.q_3p_ns_mean = np.nanmean(get_object_values(q_3p_ns, 'total', checked))
        self.q_3p_ns_opt_mean = np.nanmean(get_object_values(q_3p_ns_opt, 'total', checked))

        self.compute_percent_diff(extrap_fits=extrap_fits, transects=transects)



    def compute_percent_diff(self, extrap_fits, transects=None):
        """Computes the percent difference for each of the extrapolation options as compared to selected method.

        Parameters
        ----------
        extrap_fits: object
            Object of SelectFit
        """
        # Determine which mean is the reference
        if extrap_fits[-1].fit_method == 'Manual':
            self.man_top = extrap_fits[-1].top_method
            self.man_bot = extrap_fits[-1].bot_method
            self.man_exp = extrap_fits[-1].exponent

            if transects is not None:
                q_man = []
                checked = []
                for transect in transects:
                    q = QComp()
                    checked.append(transect.checked)

                    q.populate_data(data_in=transect,
                                    top_method=self.man_top,
                                    bot_method=self.man_bot,
                                    exponent=self.man_exp)
                    q_man.append(q)
                self.q_man_mean = np.nanmean(q_man[checked].total)
            reference_mean = self.q_man_mean

        else:
            if extrap_fits[-1].top_method_auto == 'Power':
                if np.abs(extrap_fits[-1].exponent_auto - 0.1667) < 0.0001:
                    reference_mean = self.q_pp_mean
                else:
                    reference_mean = self.q_pp_opt_mean
            elif extrap_fits[-1].top_method_auto == 'Constant':
                if np.abs(extrap_fits[-1].exponent_auto - 0.1667) < 0.0001:
                    reference_mean = self.q_cns_mean
                else:
                    reference_mean = self.q_cns_opt_mean
            else:
                if np.abs(extrap_fits[-1].exponent_auto - 0.1667) < 0.0001:
                    reference_mean = self.q_3p_ns_mean
                else:
                    reference_mean = self.q_3p_ns_opt_mean

        # Compute percent difference from reference
        self.q_pp_per_diff = ((self.q_pp_mean - reference_mean) / reference_mean) * 100
        self.q_pp_opt_per_diff = ((self.q_pp_opt_mean - reference_mean) / reference_mean) * 100
        self.q_cns_per_diff = ((self.q_cns_mean - reference_mean) / reference_mean) * 100
        self.q_cns_opt_per_diff = ((self.q_cns_opt_mean - reference_mean) / reference_mean) * 100
        self.q_3p_ns_per_diff = ((self.q_3p_ns_mean - reference_mean) / reference_mean) * 100
        self.q_3p_ns_opt_per_diff = ((self.q_3p_ns_opt_mean - reference_mean) / reference_mean) * 100

        if extrap_fits[-1].fit_method == 'Manual':
            self.q_man_per_diff =((self.q_man_mean - reference_mean) / reference_mean) * 100
=======
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
        
>>>>>>> 6ca6c50c231afa610ed3a693864074d7104a5f20
