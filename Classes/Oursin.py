import numpy as np
import pandas as pd
from scipy.stats import t
import copy
#copy.copy()
#copy.deepcopy()

class Oursin(object):
    """Computes the uncertainty of a measurement.

    Attributes
    ----------
    cov: float
        Coefficient of variation for all transects used in dicharge computation
    cov_95: float
        Coefficient of variation inflated to the 95% percent coverage
    invalid_95: float
        Estimated 95% uncertainty for dicharge in invalid bins and ensembles
    edges_95: float
        Estimated 95% uncertainty for the computed edge discharges
    extrapolation_95: float
        Estimated 95% uncertainty in discharge due to top and bottom extrapolations
    moving_bed_95: float
        Estimated 95% uncertainty due to moving-bed tests and conditions
    systematic: float
        Systematic error estimated at 1.5% at 1 sigma
    total_95: float
        Estimated 95% uncertainty in discharge using automated values
    cov_95_user: float
        User provided estimate of coefficient of variation inflated to the 95% percent coverage
    invalid_95_user: float
        User provided estimate of95% uncertainty for discharge in invalid bins and ensembles
    edges_95_user: float
        User provided estimate of 95% uncertainty for the computed edge discharges
    extrapolation_95_user: float
        User provided estimate of 95% uncertainty in discharge due to top and bottom extrapolations
    moving_bed_95_user: float
        User provided estimate of 95% uncertainty due to moving-bed tests and conditions
    systematic_user: float
        User provided estimate of systematic error estimated at 1.5% at 1 sigma
    total_95_user: float
        Estimated 95% uncertainty in discharge using user provide values to override automated values
    """
    #TODO : update info above

    def __init__(self):
        """Initializes the instance variables."""
        self.mode = None
        
        # Paramaters user
        self.exp_min_user=None 
        self.exp_max_user=None
        self.edg_prct_left_user=None
        self.edg_prct_right_user=None
        self.draft_prct_user=None
        self.draft_cm_user=None
        self.ev_user=None
        self.u_syst_user=None
        
        # Paramaters defaults        
        self.exp_min=None 
        self.exp_max=None
        self.edg_prct_left=None
        self.edg_prct_right=None
        self.draft_prct=None
        self.draft_cm=None
        self.WTev=None
        self.BTev=None
        self.u_syst=None        
        
        # By transects
        self.u_combined_by_transect = []
        self.U_total_by_transect = []
        self.variance_terms_by_transect= None
        # Overall transects 
        self.u_combined_total = None
        self.U_total_95 = None
        self.u_terms = []
        
        self.contrib_terms_total = None
        self.contrib_legend = None
        self.nb_transects = None
        
        #@@@ relative uncertainty at 68% level of confidence
        # in a list : one element for each checked transect
        self.u_meas_list = [] # uncertainty of measured area (-> PlanADCP)
        self.u_syst_list = [] # systematic uncertainty (minimum)
        self.u_ens_list = [] # uncertainty due to number of ensembles
        self.u_proj_list = [] # uncertainty due to projection neglected now
        self.u_ev_list = [] # error velocity
        self.u_badcell_list = [] # bad cells
        self.u_badens_list = [] # bad ensembles
        self.u_top_list = []
        self.u_bot_list = []
        self.u_left_list = []
        self.u_right_list = []
        
        self.simu1 = []
        self.simu1_left =  []
        self.simu1_right =  []   
        self.simu1_top = []
        self.simu1_bot = []
        self.simu2 = []
        self.simu2_top = []
        self.simu2_bot = []
        self.simu3 = []
        self.simu3_top = []
        self.simu3_bot = []        
        self.simu4 = []
        self.simu4_top = []
        self.simu4_bot = []        
        self.simu5 = []
        self.simu5_top = []
        self.simu5_bot = []        
        self.simu6 = []
        self.simu6_top = []
        self.simu6_bot = []        
        self.simu7 = []
        self.simu7_top = []
        self.simu7_bot = []        
        self.simu8 = []
        self.simu8_top = []
        self.simu8_bot = []        
        self.simu9 = []       
        self.simu9_left =[]
        self.simu9_right = [] 
        self.simu10 = []   
        self.simu10_left = []
        self.simu10_right =[]
        self.simu11 = []
        self.simu11_top = []
        self.simu11_left = []
        self.simu11_right = [] 
        self.simu12 = []
        self.simu12_top = []
        self.simu12_left = []
        self.simu12_right = []         
        self.simu13 = []
        self.simu14 = []
        self.simu15 = []
        self.simu16 = []
        self.simu16_left = []
        self.simu16_right = []
        self.simu9 =  []
        self.simu10 = []
        self.simu11 = []
        self.simu12 = []
        self.simu13 = []
        self.simu14 = []
        self.simu15 = []
        self.simu16 = []
        self.simu17 = []
        self.simu18 = []
        self.simu19 = []
        self.simu20 = []
        self.simu21 = []
        self.simu22 = []
        self.simu23 = []
        self.simu24 = []
        self.simu25 = []
        
    def compute_oursin(self, meas, mode , exp_min_user=None, exp_max_user=None, 
                       edg_prct_left_user=None,edg_prct_right_user=None,
                       draft_prct_user=None,draft_cm_user=None,ev_user=None,
                       u_syst_user=None):
        """Computes the uncertainty for the components of the discharge measurement
        using measurement data or user provided values.

        Parameters
        ----------
        meas: Measurement
            Object of class Measurement
        
        TODO : ajouter les paramètres de l'utilsateur
        """
        
        # Use user paramters, or default parameters otherwise
        if mode=="Default":
            self.mode="Default"
            self.exp_min=0.25
            self.exp_max=0.125
            self.edg_prct_left=0.20
            self.edg_prct_right=0.20
            self.draft_prct=None
            self.draft_cm=0.02
            self.WTev=0.1
            self.u_syst=1.31           
                        
        elif mode=="Manual":
            self.mode="Manual"
            
            # TODO si la valeur est entrée en prct : calculer la valeur abs
            
        else:
            self.mode="Expert"
            
            if "Power" in [meas.extrap_fit.sel_fit[x].bot_method_auto for x in range(len(meas.transects)) if meas.transects[x].checked==1]:
                self.exp_min=list(min([meas.extrap_fit.sel_fit[x].exponent_95_ci[0] for x in range(len(meas.transects)) if meas.transects[x].checked==1 and meas.extrap_fit.sel_fit[x].top_method_auto=="Power"]))[0]
                self.exp_max=list(max([meas.extrap_fit.sel_fit[x].exponent_95_ci[1] for x in range(len(meas.transects)) if meas.transects[x].checked==1 and meas.extrap_fit.sel_fit[x].top_method_auto=="Power"]))[0]
            elif "No Slip" in [meas.extrap_fit.sel_fit[x].bot_method_auto for x in range(len(meas.transects)) if meas.transects[x].checked==1]:
                self.exp_min=list(min([meas.extrap_fit.sel_fit[x].exponent_95_ci[0] for x in range(len(meas.transects)) if meas.transects[x].checked==1 and meas.extrap_fit.sel_fit[x].bot_method_auto=="No Slip"]))[0]
                self.exp_max=list(max([meas.extrap_fit.sel_fit[x].exponent_95_ci[1] for x in range(len(meas.transects)) if meas.transects[x].checked==1 and meas.extrap_fit.sel_fit[x].bot_method_auto=="No Slip"]))[0]
            else:
                self.exp_min=0.25
                self.exp_max=0.125
                
            self.edg_prct_left=0.20
            self.edg_prct_right=0.20
            
            self.draft_prct=None
            self.draft_cm=0.02
            
            self.WTev=meas.current_settings()['WTdFilterThreshold']
            self.BTev=meas.current_settings()['BTdFilterThreshold']
            
            self.u_syst=1.31     
            

        # TODO : mode user
        # TODO : mode default (draft 2cm, + or - 20% edge and coeffs)
        # TODO : automatic mode based on statistical results when fitting exponents...
        # TODO : extraire le best fit de Extrap
        
        
        # Use only checked discharges
        checked = []
        discharges = []
        simu1 = []
        for n in range(len(meas.transects)):
            checked.append(meas.transects[n].checked)
            if meas.transects[n].checked:
                discharges.append(meas.discharge[n])
                simu1.append(meas.discharge[n])
                self.u_syst_list.append(1.31)
                self.u_ens_list.append(32*len(meas.discharge[n].middle_ens)**(-0.88))

        self.nb_transects = len(self.u_syst_list)
        
        #ind_checked = [x for x in range(len(meas.transects)) if meas.transects[x].checked==1]

        # Run all simulations for OURSIN
        self.run_oursin_simulations(meas) # it should compute all simulations
        
        # TODO : sortir :
        # - best method Extrap
        # - best exponent
        # - left & right edge distance
        # - draft
        # - error velocity computed by QRev
        
        # TODO :
        # - Analyser comment le fit est fait dans Extrap
        # - Analyser comment le seuil error velocity est calculé
        
        # Compute uncertainties based on simulations using staticmethods
        # TODO refaire le calcul pour u_meas à partir des courbes de Plan ADCP 
        self.u_meas_compute(meas) # = [] # uncertainty of measured area (-> PlanADCP) 
        #self.u_meas_list = self.u_ens_list
        self.u_ev_list = 100*Oursin.apply_u_rect([self.simu1,self.simu17])/self.simu1
        self.u_badens_list = 100*Oursin.apply_u_rect([self.simu13,self.simu14,self.simu15])/self.simu1 
        self.u_badcell_list = 100*Oursin.apply_u_rect([self.simu18,self.simu19,self.simu20,self.simu21,self.simu22,self.simu23,self.simu24,self.simu25])/self.simu1
        self.u_top_list  = 100*Oursin.apply_u_rect([self.simu1_top,self.simu2_top,self.simu3_top,self.simu4_top,self.simu5_top,self.simu6_top,self.simu11_top,self.simu12_top])/self.simu1        
        self.u_bot_list  = 100*Oursin.apply_u_rect([self.simu1_bot,self.simu3_bot,self.simu4_bot,self.simu7_bot,self.simu8_bot])/self.simu1
        self.u_left_list  = 100*Oursin.apply_u_rect([self.simu1_left,self.simu9_left,self.simu10_left,self.simu11_left,self.simu12_left,self.simu16_left])/self.simu1
        self.u_right_list  = 100*Oursin.apply_u_rect([self.simu1_right,self.simu9_right,self.simu10_right,self.simu11_right,self.simu12_right,self.simu16_right])/self.simu1
         
        # Recherche des paramètres ou calcul par défault
        # TODO : mode expert : results = user        

        # Estimate the total measurement uncertainty and the combined uncertainty for each transects
        self.compute_combined_uncertainty()
  
    def func_pow(x, a=None, b=None, c=None):
        return a * np.exp(-b * x) + c

    def u_meas_compute(self,meas):
        u_meas_list=[]
        ind_transect = 0
        u_prct_dzi = 0.5 # k=1
        r_corr = 1 # Correlation between errors of contiguous cells
        
        for transect in meas.transects:
            #print(ind_transect)
            # ---- Computation of v_boatm by ensemble
            v_bm_ens_abs = np.absolute(meas.transects[ind_transect].boat_vel.bt_vel.d_mps)
            v_bm_ens_abs[np.isnan(v_bm_ens_abs)] = 0.00    
            d_ens = meas.transects[ind_transect].depths.bt_depths.depth_processed_m # depth by ens 
            freq_v_boat_adcp = np.asarray(list(meas.transects[ind_transect].boat_vel.bt_vel.frequency_khz))    
            u_prct_v_boatm_ens = 0.01*( (0.03*v_bm_ens_abs)+((1+(0.3*v_bm_ens_abs))/(1+(0.001*d_ens*freq_v_boat_adcp))) ) / v_bm_ens_abs
            u_prct_v_boatm_ens[u_prct_v_boatm_ens == np.inf] = 0 # replace inf by 0
            
            #-----Cell size [m] for each ensemble ARRAY ## meas.transects[ind_transect].depths.bt_depths.depth_cell_size_m # meas.transects[ind_transect].depths.bt_depths.depth_cell_size_m.shape
            h_cell_mean_ens = np.mean(meas.transects[ind_transect].depths.bt_depths.depth_cell_size_m,axis=0)
            
            #-----Mean velocity [m/s] for each ensemble # ARRAY
            wv_mean_ens = meas.discharge[ind_transect].middle_ens/meas.transects[ind_transect].depths.bt_depths.depth_processed_m

            #-----Computation of v_water by ensembles
            if meas.transects[0].adcp.manufacturer == 'SonTek': # OTHER fit: a=7.412, b=4.178, c=4.749
                freq_adcp = None
                mode_adcp = None
                water_ping = None
                #u_prct_v_water_ens=[0.01*self.func_pow(s,7.412,4.178,4.749) for s in list(h_cell_mean_ens)]               
                u_prct_v_water_ens=np.array([0.01*7.412*np.exp(-4.178 * xi) +4.749 for xi in h_cell_mean_ens])    
                #u_prct_v_water_ens=[0.01*7.412*np.exp(-4.178 * s) +4.749 for s in list(h_cell_mean_ens)] 
            elif meas.transects[0].adcp.manufacturer =='TRDI' : 
                # meas.transects[0].adcp.frequency_khz
                # meas.transects[0].adcp.configuration_commands
                #meas.transects[0].adcp.model
                
                #meas.transects[0].w_vel.water_mode
    
                #-----Criterion 1: water ping
#                if meas.transects[0].adcp.configuration_commands == None:
#                    water_ping = 1
#                    print('Water ping unknown')
                if [True for s in list(meas.transects[0].adcp.configuration_commands) if "WP" in s]:
                    water_ping = int([s for s in list(meas.transects[0].adcp.configuration_commands) if "WP" in s][0].split('WP')[1])
                else:
                    water_ping = 1
                    
                #-----Criterion 2: water mode adcp
                if int(meas.transects[0].w_vel.water_mode)==1:
                    mode_adcp = 1
                    if water_ping ==1:
#                        u_prct_v_water_ens=[0.01*self.func_pow(s,a=53.156, 5.553, 3.098) for s in list(h_cell_mean_ens)]
#                        u_prct_v_water_ens=np.array([0.01*self.func_pow(xi,53.156, 5.553, 3.098) for xi in h_cell_mean_ens])
                        u_prct_v_water_ens=np.array([0.01*53.156*np.exp(-5.553 * xi) +3.098 for xi in h_cell_mean_ens])
                    else:
#                        u_prct_v_water_ens=[0.01*self.func_pow(s,a=23.805, 5.618, 1.445) for s in list(h_cell_mean_ens)]
#                        u_prct_v_water_ens=np.array([0.01*self.func_pow(xi,23.805, 5.618, 1.445) for xi in h_cell_mean_ens])
                        u_prct_v_water_ens=np.array([0.01*23.805*np.exp(-5.618 * xi) +1.445 for xi in h_cell_mean_ens])
                        
                elif int(meas.transects[0].w_vel.water_mode)==5:
                    mode_adcp = 5
                    if water_ping ==1:
#                        u_prct_v_water_ens=[0.01*self.func_pow(s,a=0.719, 10.534, 0.179) for s in list(h_cell_mean_ens)]
#                        u_prct_v_water_ens=np.array([0.01*self.func_pow(xi,0.719, 10.534, 0.179) for xi in h_cell_mean_ens])
                        u_prct_v_water_ens=np.array([0.01*0.719*np.exp(-10.534 * xi) +0.179 for xi in h_cell_mean_ens])
                    else:
#                        u_prct_v_water_ens=[0.01*self.func_pow(s,a=0.307, 9.538, 0.080) for s in list(h_cell_mean_ens)]
#                        u_prct_v_water_ens=np.array([0.01*self.func_pow(xi,0.307, 9.538, 0.080) for xi in h_cell_mean_ens])
                        u_prct_v_water_ens=np.array([0.01*0.307*np.exp(-9.538 * xi) +0.080 for xi in h_cell_mean_ens])
                        
                elif int(meas.transects[0].w_vel.water_mode)==11:
                    mode_adcp = 11
                    if water_ping ==1:
#                        u_prct_v_water_ens=[0.01*self.func_pow(s,1.512, 9.907, 0.419) for s in list(h_cell_mean_ens)]
#                        u_prct_v_water_ens=np.array([0.01*self.func_pow(xi,1.512, 9.907, 0.419) for xi in h_cell_mean_ens])
                        u_prct_v_water_ens=np.array([0.01*1.512*np.exp(-9.907 * xi) +0.419 for xi in h_cell_mean_ens])
                    else:
#                        u_prct_v_water_ens=[0.01*self.func_pow(s,0.719, 10.534, 0.179) for s in list(h_cell_mean_ens)]
#                        u_prct_v_water_ens=np.array([0.01*self.func_pow(xi,0.719, 10.534, 0.179) for xi in h_cell_mean_ens])
                        u_prct_v_water_ens=np.array([0.01*0.719*np.exp(-10.534 * xi) +0.179 for xi in h_cell_mean_ens])
                        
                elif int(meas.transects[0].w_vel.water_mode)==12:
                    mode_adcp = 12 
                    if water_ping ==1:
#                        u_prct_v_water_ens=[0.01*self.func_pow(s,21.735, 5.620, 1.320) for s in list(h_cell_mean_ens)]
#                        u_prct_v_water_ens=np.array([0.01*self.func_pow(xi,21.735, 5.620, 1.320) for xi in h_cell_mean_ens])
                        u_prct_v_water_ens=np.array([0.01*21.735*np.exp(-5.620 * xi) +1.320 for xi in h_cell_mean_ens])
                    else:
#                        u_prct_v_water_ens=[0.01*self.func_pow(s,9.729, 5.634, 0.594) for s in list(h_cell_mean_ens)]
#                        u_prct_v_water_ens=np.array([0.01*self.func_pow(xi,9.729, 5.634, 0.594) for xi in h_cell_mean_ens])
                        u_prct_v_water_ens=np.array([0.01*9.729*np.exp(-5.634 * xi) +0.594 for xi in h_cell_mean_ens])
                        
                else: # elif: int(meas.transects[0].w_vel.water_mode)==8 or OTHER water mode
                    mode_adcp = 8
                    if water_ping ==1:
#                        u_prct_v_water_ens=[0.01*self.func_pow(s,8.298, 10.008, 2.099) for s in list(h_cell_mean_ens)]
#                        u_prct_v_water_ens=np.array([0.01*self.func_pow(xi,8.298, 10.008, 2.099) for xi in h_cell_mean_ens])
                        u_prct_v_water_ens=np.array([0.01*8.298*np.exp(-10.008 * xi) +2.099 for xi in h_cell_mean_ens])
                    else:
#                        u_prct_v_water_ens=[0.01*self.func_pow(s,3.719, 10.043, 0.939) for s in list(h_cell_mean_ens)]
#                        u_prct_v_water_ens=np.array([0.01*self.func_pow(xi,3.719, 10.043, 0.939) for xi in h_cell_mean_ens])
                        u_prct_v_water_ens=np.array([0.01*3.719*np.exp(-10.043 * xi) +0.939 for xi in h_cell_mean_ens])
                #-----Criterion 3: frequence meas.transects[ind_transect].w_vel.frequency
                # not taken into account now due to possible variation in frequences with sontek
                # we should wait for a more explicit equation
                


                #-----Equation 
#                RG_600_1ping_Mode1 fit: a=53.572, 2.875, c=3.146
#                RG_600_1ping_Mode5 fit: a=1.387, b=12.760, c=0.255
#                RG_600_1ping_Mode8 fit: a=10.946, b=10.674, c=2.531
#                RG_600_1ping_Mode11 fit: a=2.053, b=10.063, c=0.569
#                RG_600_1ping_Mode12 fit: a=30.460, b=4.608, c=2.415
#                RG_600_5ping_Mode1 fit: a=30.461, b=4.058, c=2.306
#                RG_600_5ping_Mode5 fit: a=0.816, b=15.352, c=0.115
#                RG_600_5ping_Mode8 fit: a=4.842, b=10.598, c=1.133
#                RG_600_5ping_Mode11 fit: a=0.910, b=9.805, c=0.251
#                RG_600_5ping_Mode12 fit: a=13.615, b=4.607, c=1.081
                
#                RG_1200_1ping_Mode1 fit: a=53.156, b=5.553, c=3.098 ##
#                RG_1200_1ping_Mode5 fit: a=0.719, b=10.534, c=0.179 ##
#                RG_1200_1ping_Mode8 fit: a=8.298, b=10.008, c=2.099 ##
#                RG_1200_1ping_Mode11 fit: a=1.512, b=9.907, c=0.419 ##
#                RG_1200_1ping_Mode12 fit: a=21.735, b=5.620, c=1.320 ##
                
#                RG_1200_5ping_Mode1 fit: a=23.805, b=5.618, c=1.445 ##
#                RG_1200_5ping_Mode5 fit: a=0.307, b=9.538, c=0.080 ##
#                RG_1200_5ping_Mode8 fit: a=3.719, b=10.043, c=0.939
#                RG_1200_5ping_Mode11 fit: a=0.719, b=10.534, c=0.179 ##
#                RG_1200_5ping_Mode12 fit: a=9.729, b=5.634, c=0.594 ##
    
            else:
#                OTHER fit: a=7.412, b=4.178, c=4.749
                 u_prct_v_water_ens=[0.01*self.func_pow(s,a=7.412, b=4.178, c=4.749) for s in list(h_cell_mean_ens)]                    
                
            u_prct_v_water_ens[u_prct_v_water_ens == np.inf] = 0 # replace inf by 0
            u_prct_v_water_ens[u_prct_v_water_ens == -np.inf] = 0 # replace -inf by 0
            
            # ----- Computation of u_meas
            Q_2_tran = meas.discharge[ind_transect].total**2
            q_2_ens  = meas.discharge[ind_transect].middle_ens**2
            n_cell_ens = meas.transects[ind_transect].w_vel.cells_above_sl.sum(axis=0)
    
            u_2_meas = q_2_ens*(u_prct_v_boatm_ens**2 + (1/n_cell_ens)*(u_prct_v_water_ens**2 + u_prct_v_water_ens**2 * r_corr * (n_cell_ens-1) + u_prct_dzi**2 ))
            u_2_meas = np.nan_to_num(u_2_meas)
            u_2_prct_meas = u_2_meas.sum()/Q_2_tran
            u_prct_meas = u_2_prct_meas**0.5
            
            u_meas_list.append(u_prct_meas)
            #print('OK append')
            ind_transect += 1    
        
        self.u_meas_list= u_meas_list  
    
    def run_oursin_simulations(self,meas):
        coeff_user = 1/6
        prct_dist_edge = 0.20
        draft_diff = 0.02 #m draft error
        WTdFilterThreshold_test = 0.2 #0.1 # 1.067
        BTdFilterThreshold_test = 0.1
        
        ind_checked = [x for x in range(len(meas.transects)) if meas.transects[x].checked==1]
        
        # TODO : add options (run or not the simulation, if not the simu1 is considered instead)
        simextrap = copy.deepcopy(meas)
        simedge = copy.deepcopy(meas)
        simdraft = copy.deepcopy(meas)
        simbadcells = copy.deepcopy(meas)
        simbadens = copy.deepcopy(meas)
        simerrorvel = copy.deepcopy(meas)
        
        #-----Simulation 1 is the reference
        self.simu1 = [meas.discharge[x].total for x in range(len(meas.transects)) if meas.transects[x].checked==1]
        self.simu1_left = [meas.discharge[x].left for x in range(len(meas.transects)) if meas.transects[x].checked==1]
        self.simu1_right = [meas.discharge[x].right for x in range(len(meas.transects)) if meas.transects[x].checked==1]
        self.simu1_top = [meas.discharge[x].top for x in range(len(meas.transects)) if meas.transects[x].checked==1]
        self.simu1_bot = [meas.discharge[x].bottom for x in range(len(meas.transects)) if meas.transects[x].checked==1]
        
        #-----Simulations based on change_extrapolations
        # TODO : mettre ici les résultats de la best simulation retenue
        # TODO : ajouter en variable les paramètres de l'utilisateur (coeff 1/4, 1/8 ou pourcentage par 20% de moins que le coeff retenu)
        # TODO : différencier le cas d'un PP d'un CNS 
        # TODO : négliger la simulation 2, c'est le PP mean en principe (prendre directement le résultat fourni par QRev)

        simextrap.change_extrapolation(method='Manual', top='Power', bot='Power', exp=(1/6))
        self.simu2 = [simextrap.discharge[x].total for x in range(len(simextrap.transects)) if simextrap.transects[x].checked==1]
        self.simu2_top = [simextrap.discharge[x].top for x in range(len(simextrap.transects)) if simextrap.transects[x].checked==1]
        self.simu2_bot = [simextrap.discharge[x].bottom for x in range(len(simextrap.transects)) if simextrap.transects[x].checked==1]        
        
        simextrap.change_extrapolation(method='Manual', top='Power', bot='Power', exp=(1/4))
        self.simu3 = [simextrap.discharge[x].total for x in range(len(simextrap.transects)) if simextrap.transects[x].checked==1]
        self.simu3_top = [simextrap.discharge[x].top for x in range(len(simextrap.transects)) if simextrap.transects[x].checked==1]
        self.simu3_bot = [simextrap.discharge[x].bottom for x in range(len(simextrap.transects)) if simextrap.transects[x].checked==1]   
        
        simextrap.change_extrapolation(method='Manual', top='Power', bot='Power', exp=(1/8))
        self.simu4 = [simextrap.discharge[x].total for x in range(len(simextrap.transects)) if simextrap.transects[x].checked==1]
        self.simu4_top = [simextrap.discharge[x].top for x in range(len(simextrap.transects)) if simextrap.transects[x].checked==1]
        self.simu4_bot = [simextrap.discharge[x].bottom for x in range(len(simextrap.transects)) if simextrap.transects[x].checked==1]   
        
        simextrap.change_extrapolation(method='Manual', top='Constant', bot='No Slip', exp=coeff_user)
        self.simu5 = [simextrap.discharge[x].total for x in range(len(simextrap.transects)) if simextrap.transects[x].checked==1]
        self.simu5_top = [simextrap.discharge[x].top for x in range(len(simextrap.transects)) if simextrap.transects[x].checked==1]
        self.simu5_bot = [simextrap.discharge[x].bottom for x in range(len(simextrap.transects)) if simextrap.transects[x].checked==1]     
        
        simextrap.change_extrapolation(method='Manual', top='Constant', bot='No Slip', exp=coeff_user)
        self.simu6 = [simextrap.discharge[x].total for x in range(len(simextrap.transects)) if simextrap.transects[x].checked==1]
        self.simu6_top = [simextrap.discharge[x].top for x in range(len(simextrap.transects)) if simextrap.transects[x].checked==1]
        self.simu6_bot = [simextrap.discharge[x].bottom for x in range(len(simextrap.transects)) if simextrap.transects[x].checked==1]        
        
        simextrap.change_extrapolation(method='Manual', top='Power', bot='Power', exp=(1/6))
        self.simu7 = [simextrap.discharge[x].total for x in range(len(simextrap.transects)) if simextrap.transects[x].checked==1]
        self.simu7_top = [simextrap.discharge[x].top for x in range(len(simextrap.transects)) if simextrap.transects[x].checked==1]
        self.simu7_bot = [simextrap.discharge[x].bottom for x in range(len(simextrap.transects)) if simextrap.transects[x].checked==1] 
        
        simextrap.change_extrapolation(method='Manual', top='Constant', bot='No Slip', exp=coeff_user)
        self.simu8 = [simextrap.discharge[x].total for x in range(len(simextrap.transects)) if simextrap.transects[x].checked==1]
        self.simu8_top = [simextrap.discharge[x].top for x in range(len(simextrap.transects)) if simextrap.transects[x].checked==1]
        self.simu8_bot = [simextrap.discharge[x].bottom for x in range(len(simextrap.transects)) if simextrap.transects[x].checked==1]         
        
        del simextrap
        
        #-----Simulations 9 and 10 (edges) shloud be run transect by transect
        #-----Simulation 11 and 12 (draft) """"
        #-----Simulation 13 to 15 bad ensembles
        #-----Simulation 18 to 25 bad cells """"     
        
        for ind_transect in ind_checked:
            init_dist_right = simedge.transects[ind_transect].edges.right.distance_m
            init_dist_left = simedge.transects[ind_transect].edges.left.distance_m
            # Simulation 9 (underestimation case)
            simedge.transects[ind_transect].edges.left.distance_m = (1-prct_dist_edge)*init_dist_left
            simedge.transects[ind_transect].edges.right.distance_m = (1-prct_dist_edge)*init_dist_right
            simedge.transects[ind_transect].edges.left.type = 'Triangular'
            simedge.transects[ind_transect].edges.right.type = 'Triangular'
            
            simedge.compute_discharge()
            self.simu9.append(simedge.discharge[ind_transect].total)
            self.simu9_right.append(simedge.discharge[ind_transect].right)
            self.simu9_left.append(simedge.discharge[ind_transect].left)
            
            # Simulation 10 (overestimation case)
            simedge.transects[ind_transect].edges.left.distance_m = (1+prct_dist_edge)*init_dist_left
            simedge.transects[ind_transect].edges.right.distance_m = (1+prct_dist_edge)*init_dist_right
            simedge.transects[ind_transect].edges.left.type = 'Rectangular'
            simedge.transects[ind_transect].edges.right.type = 'Rectangular'
            
            simedge.compute_discharge()
            self.simu10.append(simedge.discharge[ind_transect].total)
            self.simu10_right.append(simedge.discharge[ind_transect].right)
            self.simu10_left.append(simedge.discharge[ind_transect].left)
            
            # Simulation 11 -0.02m
            simdraft.transects[ind_transect].change_draft(meas.transects[ind_transect].depths.bt_depths.draft_orig_m-draft_diff)
            
            simdraft.compute_discharge()
            self.simu11.append(simdraft.discharge[ind_transect].total)
            self.simu11_top.append(simedge.discharge[ind_transect].top)
            self.simu11_left.append(simedge.discharge[ind_transect].left)
            self.simu11_right.append(simedge.discharge[ind_transect].right)
           
            # Simulation 12 +0.02m
            simdraft.transects[ind_transect].change_draft(meas.transects[ind_transect].depths.bt_depths.draft_orig_m+draft_diff)   
            
            simdraft.compute_discharge()
            self.simu12.append(simdraft.discharge[ind_transect].total)
            self.simu12_top.append(simedge.discharge[ind_transect].top)
            self.simu12_left.append(simedge.discharge[ind_transect].left)
            self.simu12_right.append(simedge.discharge[ind_transect].right)
           
            
            # Simulation 13 Bad ensembles (last ens)
            simbadens.transects[ind_transect].w_vel.apply_interpolation(simbadens.transects[ind_transect],ens_interp='HoldLast', cells_interp= 'TRDI')# Enter desired interpolation method here)
            simbadens.compute_discharge()
            self.simu13.append(simbadens.discharge[ind_transect].total)
           
            # Simulation 14 Bad ensembles (next ens)
            simbadens.transects[ind_transect].w_vel.apply_interpolation(simbadens.transects[ind_transect],ens_interp='ExpandedT', cells_interp= 'TRDI')# Enter desired interpolation method here)
            simbadens.compute_discharge()
            self.simu14.append(simbadens.discharge[ind_transect].total)
            
            # Simulation 15 Bad ensembles (linear interpolation)
            simbadens.transects[ind_transect].w_vel.apply_interpolation(simbadens.transects[ind_transect],ens_interp='Hold', cells_interp= 'TRDI')# Enter desired interpolation method here)
            simbadens.compute_discharge()
            self.simu15.append(simbadens.discharge[ind_transect].total) 
            
            # simu 18
            simbadens.transects[ind_transect].w_vel.apply_interpolation(simbadens.transects[ind_transect],ens_interp='HoldLast', cells_interp= 'TRDI')# Enter desired interpolation method here)
            simbadcells.compute_discharge()
            self.simu18.append(simbadcells.discharge[ind_transect].total)
           
            # simu 19
            simbadens.transects[ind_transect].w_vel.apply_interpolation(simbadens.transects[ind_transect],ens_interp='Hold9', cells_interp= 'TRDI')# Enter desired interpolation method here)
            simbadcells.compute_discharge()
            self.simu19.append(simbadcells.discharge[ind_transect].total)
            
            # simu 20
            simbadens.transects[ind_transect].w_vel.apply_interpolation(simbadens.transects[ind_transect],ens_interp='Linear', cells_interp= 'TRDI')# Enter desired interpolation method here)
            simbadcells.compute_discharge()
            self.simu20.append(simbadcells.discharge[ind_transect].total) 
         
            # simu 21
            simbadens.transects[ind_transect].w_vel.apply_interpolation(simbadens.transects[ind_transect],ens_interp='HoldLast', cells_interp= 'TRDI')# Enter desired interpolation method here)
            simbadcells.compute_discharge()
            self.simu21.append(simbadcells.discharge[ind_transect].total) 
            
            # simu 22
            simbadens.transects[ind_transect].w_vel.apply_interpolation(simbadens.transects[ind_transect],ens_interp='HoldLast', cells_interp= 'TRDI')# Enter desired interpolation method here)
            simbadcells.compute_discharge()
            self.simu22.append(simbadcells.discharge[ind_transect].total) 
            
            # simu 23
            simbadens.transects[ind_transect].w_vel.apply_interpolation(simbadens.transects[ind_transect],ens_interp='HoldLast', cells_interp= 'TRDI')# Enter desired interpolation method here)
            simbadcells.compute_discharge()
            self.simu23.append(simbadcells.discharge[ind_transect].total) 
            
            # simu 24
            simbadens.transects[ind_transect].w_vel.apply_interpolation(simbadens.transects[ind_transect],ens_interp='HoldLast', cells_interp= 'TRDI')# Enter desired interpolation method here)
            simbadcells.compute_discharge()
            self.simu24.append(simbadcells.discharge[ind_transect].total) 
            
            # simu 25
            simbadens.transects[ind_transect].w_vel.apply_interpolation(simbadens.transects[ind_transect],ens_interp='HoldLast', cells_interp= 'Linear')# Enter desired interpolation method here)
            simbadcells.compute_discharge()
            self.simu25.append(simbadcells.discharge[ind_transect].total) 
            
        del simedge 
        del simdraft           
        del simbadens
        del simbadcells
                
        #-----Simulation 16 Neglected at the moment
        self.simu16 = self.simu1
        self.simu16_left = self.simu1_left
        self.simu16_right = self.simu1_right
        
        #-----Simulation 17 (and 26) error velocity
        # TODO : add another simulation to produce a min-max interval
        settings = simerrorvel.current_settings()
        settings['WTdFilter'] = 'Manual'
        settings['BTdFilter'] = 'Manual'
        settings['WTdFilterThreshold']=WTdFilterThreshold_test
        settings['BTdFilterThreshold']=BTdFilterThreshold_test
        simerrorvel.apply_settings(settings)
        self.simu17=[simerrorvel.discharge[x].total for x in ind_checked]
        
        del simerrorvel
        # END DEF

    def compute_combined_uncertainty(self):
        TABLE_u = pd.DataFrame()
        TABLE_u['u_meas']=self.u_meas_list
        TABLE_u['u_syst']=self.u_syst_list        
        TABLE_u['u_ens']=self.u_ens_list
        TABLE_u['u_ev']=self.u_ev_list
        TABLE_u['u_badcell']=self.u_badcell_list
        TABLE_u['u_badens']=self.u_badens_list
        TABLE_u['u_top']=self.u_top_list
        TABLE_u['u_bot']=self.u_bot_list
        TABLE_u['u_left']=self.u_left_list
        TABLE_u['u_right']=self.u_right_list
        
        TABLE_u2 = TABLE_u.applymap(lambda x:x**2)
        
        self.u_combined_by_transect = [x**0.5 for x in list(TABLE_u2.sum(axis=1))]
        self.U_total_by_transect = [2*x for x in self.u_combined_by_transect]
    
        # Uncertainty for the measurement
        u2_alea = TABLE_u2['u_meas'].mean()
        u2_syst = list(TABLE_u2.drop(['u_meas'], axis=1).mean())
        u2_c_gauging = (1/self.nb_transects)*u2_alea + sum(u2_syst)      
        
        u_terms = []
        u_terms.append(((1/self.nb_transects)*u2_alea)**0.5)
        u_terms=u_terms+[x**0.5 for x in u2_syst]
        
        contrib_terms_total = []
        contrib_terms_total.append(((1/self.nb_transects)*u2_alea)/u2_c_gauging)
        contrib_terms_total=contrib_terms_total+[x/u2_c_gauging for x in u2_syst]
        
        self.u_terms = u_terms
        self.u_combined_total = (u2_c_gauging**0.5)
        self.U_total_95 = 2*self.u_combined_total
           
        self.variance_terms_by_transect= TABLE_u2
        
        # TABLE_u2.div(TABLE_u2.sum(axis=1), axis=0)
        self.contrib_terms_total = [100*x for x in contrib_terms_total]
        self.contrib_legend = list(TABLE_u.columns)
    
    @staticmethod
    def apply_u_rect(list_simu):
        simu = pd.DataFrame()
        for i in range(len(list_simu)):
            name_col='simu' + str(i)
            simu[name_col]=list_simu[i]
        u_rect = (simu.max(axis=1)-simu.min(axis=1))/(2*(3**0.5))
        
        return(u_rect)
    
    @staticmethod
    def get_array_attr(list_in, prop):
        # Create array of specified attribute
        data = []
        for item in list_in:
            data.append(getattr(item, prop))
        np.asarray(data)
        return data
    