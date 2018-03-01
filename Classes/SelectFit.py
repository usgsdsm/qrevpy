import numpy as np
import statsmodels.api as sm
from Classes.FitData import FitData

class SelectFit(object):
    """Class definition for data class that contains all of the automated
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
    DSM modified for use in QRev"""

    def __init__(self):
        self.fit_method = None  # User selected method Automatic or Manual
        self.top_method = None
        self.bot_method = None
        self.exponent = None
        self.exp_method = None
        self.u = None
        self.u_auto = None
        self.z = None
        self.z_auto = None
        self.residuals = None
        self.coef = None
        self.bot_method_auto = None  # Selected extrapolation for top
        self.top_method_auto = None  # Selected extrapolation for bottom
        self.exponent_auto = None  # Selected exponent
        self.top_fit_r2 = None  # Top fit custom r^2
        self.top_max_diff = None  # Maximum difference between power and 3-pt at top
        self.bot_diff = None  # Difference between power and no slop at z = 0.1
        self.bot_r2 = None  # Bottom fit r^2
        self.fit_r2 = None  # Selected fit of selected power/no slip fit
        self.ns_exponent = None  # No slip optimized exponent
        self.pp_exponent = None  # Power Power optimized exponent
        self.top_r2 = None
        self.rsqr = None
        self.exponent_95_ci = None

    def populate_data(self, normalized, fit_method, kargs=None):

        valid_data = np.squeeze(normalized.valid_data)

        # Store data in properties to object
        self.fit_method = fit_method

        # Compute power fit with optimized exponent as reference to determine
        # if constant no slip will be more appropriate
        ppobj = FitData()
        ppobj.populate_data(normalized, 'Power', 'Power', 'optimize')

        # Store results in obhect
        self.pp_exponent = ppobj.exponent
        self.residuals = ppobj.residuals
        self.rsqr = ppobj.r_squared
        self.exponent_95_ci = ppobj.exponent_95_ci

        # Begin automatic fit

        # More than 6 cells are required to compute an optimized fit.  For fewer
        # than 7 cells the default power/power fit is selected due to lack of sufficient
        # data for a good analysis
        if len(self.residuals) > 6:
            # Compute the difference between the top two cells of data and the optimized power fit
            top2 = np.nansum(normalized.unit_normalized_med[valid_data[-2:]] - ppobj.coef *
                             normalized.unit_normalized_z[valid_data[-2:]]) ** ppobj.exponent

            # Compute the difference between the bottom two cells of data and the optimized power fit
            bot2 = np.nansum(normalized.unit_normalized_med[valid_data[:2]] \
                             - ppobj.coef * normalized.unit_normalized_z[valid_data[:2]]) \
                   ** ppobj.exponent

            # Compute the difference between the middle two cells of data and the optimized power fit
            mid1 = int(np.floor(len(np.isnan(valid_data) == False) / 2))
            mid2 = np.nansum(normalized.unit_normalized_med[valid_data[mid1:mid1 + 2]] \
                             - ppobj.coef * normalized.unit_normalized_z[
                                 valid_data[mid1:mid1 + 2]]) \
                   ** ppobj._FitData__exponent

            self.top_method_auto = 'Power'
            self.bot_method_auto = 'Power'

            # Evaluate difference in data and power fit at water surface using a linear fit throught the top 4
            # median cells and save results
            x = normalized.unit_normalized_med[valid_data[:4]]
            #             x = sm.add_constant(x)
            y = normalized.unit_normalized_z[valid_data[:4]]
            lin_fit = sm.OLS(x, y)
            result = lin_fit.fit()
            dsmfitr2 = 1 - (np.sum(result.resid ** 2) / np.mean(np.abs(result.resid)))
            self.top_r2 = dsmfitr2
            self.top_r2 = result.rsquared

            # Evaluate overall fit
            # If the optimized power fit does not have an r^2 better than 0.8 or if the optimized
            # exponent if 0.1667 falls within the 95% confidence interval of the optimized fit,
            # there is insufficient justification to change the exponent from 0.1667
            if ppobj.rsqr < 0.8 or (0.1667 > self.exponent_95_ci[0] and 0.1667 < self.exponent_95_ci):
                # If an optimized exponent cannot be justified the linear fit is used to determine if a constant
                # fit at the top is a better alternative than a power fit.  If the power fit is the better
                # alternative the exponent is set to the default 0.1667 and the data is refit
                if np.abs(self.top_fit_r2 < 0.8 or self.top_r2 < 0.9):
                    fd = FitData()
                    ppobj = fd(normalized, 'Power', 'Power', 'Manual', 0.1667)

            # Evaluate fit of top and bottom portions of the profile
            # Set save selected exponent and associated fit statistics
            self.exponent_auto = ppobj.exponent
            self.fit_r2 = ppobj.r_squared

            # Compute the difference at the water surface between a linear fit of the top 4 measured cells
            # and the best selected power fit of the whole profile
            self.top_max_diff = ppobj.u[-1] - np.sum(result.params)

            # Evaluate the difference at the bottom between power using the whole profile and power using
            # only the bottom third
            ns_fd = FitData()
            ns_fd.populate_data(normalized, 'Constant', 'No Slip', 'Optimize')
            self.ns_exponent = ns_fd.exponent
            self.bot_r2 = ns_fd.r_squared
            self.bot_diff = ppobj.u[np.round(ppobj.z, 2) == 0.1] \
                              - ns_fd.u[np.round(ns_fd.z, 2) == 0.1]

            # Begin automatic selection logic
            # -----------------------------------

            # A constant no slip fit condition is selected if:
            #
            # 1)The top of the power fit doesn't fit the data well.
            # This is determined to be the situation when
            # (a) the difference at the water surface between the
            # linear fit and the power fit is greater than 10% and
            # (b) the difference is either positive or the difference
            # of the top measured cell differs from the best
            # selected power fit by more than 5%.
            # OR
            # 2) The bottom of the power fit doesn't fit the data
            # well. This is determined to be the situation when (a)
            # the difference between and optimized no slip fit
            # and the selected best power fit of the whole profile
            # is greater than 10% and (b) the optimized on slip fit has
            # and r^2 greater than 0.6.
            # OR
            # 3) Flow is bidirectional. The sign of the top of the
            # profile is different from the sign of the bottom of
            # the profile.
            # OR
            # 4) The profile is C-shaped. This is determined by
            # (a) the sign of the top and bottom difference from
            # the best selected power fit being different than the
            # sign of the middle difference from the best selected
            # power fit and (b) the combined difference of the top
            # and bottom difference from the best selected power
            # fit being greater than 10%.

            if (np.abs(self.top_max_diff > 0.1) \
                and (self.top_max_diff > 0 or np.abs(
                        normalized.unit_normalized_med[valid_data[0]] - ppobj._FitData__u[-1]) > 0.05) \
                or ((np.abs(self.bot_diff) > 0.1) and self.__bot_r2 > 0.6) \
                or (np.sign(normalized.unit_normalized_med)[valid_data[0]] != np.sign(
                        normalized.unit_normalized_med[valid_data[-1]]))) \
                    or np.sign(bot2) * np.sign(top2) == np.sign(mid2) and np.abs(bot2 + top2) > 0.1:

                # Set the bottom to no slip
                self.bot_method_auto = 'No Slip'
                # If the no slip fit with an optimized exponent does not have r^2 better than 0.8 use the default 0.1667 for the no slip exponent
                if ns_fd.r_squared > 0.8:
                    self.exponent_auto = ns_fd.exponent
                    self.fit_r2 = ns_fd.rsquared
                else:
                    self.exponent_auto = 0.1667
                    self.fit_r2 = np.nan

                # Use the no slip 95% confidence intervals if they are available
                if ns_fd.exponent_95_ci is not None and np.all(
                        np.isnan(ns_fd.exponent_95_ci) == False):
                    self.exponent_95_ci[0] = ns_fd.exponent_95_ci[0]
                    self.exponent_95_ci[1] = ns_fd.exponent_95_ci[1]
                else:
                    self.exponent_95_ci[0] = np.nan
                    self.exponent_95_ci[1] = np.nan

                # Set the top method to constant
                self.top_method_auto = 'Constant'

            else:

                # Leave the fit to power/power and set the best selected optimized exponent as the automatic fit exponent
                self.__exponent_auto = ppobj.exponent

        else:

            # if the data are insufficient for a valid analysis use the power/power fit with the default 0.1667 exponent
            self.__top_method_auto = 'Power'
            self.__bot_method_auto = 'Power'
            self.__exponent_auto = 0.1667
            self.__ns_exponent = 0.1667

        # Update the fit uysing the automatically selected methods
        update_fd = FitData()
        update_fd.populate_data(normalized, self.top_method_auto, self.bot_method_auto, 'Manual',
                                [self.exponent_auto])
        update_auto = update_fd

        if fit_method == 'Manual':

            if len(kargs) == 1:
                trans_data = kargs[0]
                update_fd = FitData()
                update_fd.populate_data(normalized, trans_data.extrap.top_method,
                                        trans_data.extrap.bot_method, 'Manual',
                                        [trans_data.extrap.exponent])
            else:
                update_fd = FitData()
                update_fd.populate_data(normalized, kargs[1], kargs[2], 'Manual', kargs[3])

        # Store fit data in object
        self.top_method = update_fd.top_method
        self.bot_method = update_fd.bot_method
        self.exponent = update_fd.exponent
        self.coef = update_fd.coef
        self.u = update_fd.u
        self.u_auto = update_auto.u_auto
        self.z_auto = update_auto.z_auto
        self.z = update_fd.z
        self.exp_method = update_fd.exp_method
        self.residuals = update_fd.residuals
        
    
