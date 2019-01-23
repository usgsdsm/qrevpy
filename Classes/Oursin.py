# -*- coding: utf-8 -*-
"""
Created on Wed Jan 23 09:59:49 2019
@author: aurelien.despax
Contact: aurelien.despax@irstea.fr
"""

import numpy as np
import pandas as pd
from scipy.stats import t
import copy

class Oursin(object):
    """Similar to Uncertainty object /// Computes the uncertainty of a measurement.

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
        self.u_meas_list = self.u_ens_list # = [] # uncertainty of measured area (-> PlanADCP)
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
            simbadens.transects[ind_transect].w_vel.apply_interpolation(simbadens.transects[ind_transect],target='Ensembles', interp_type= 'HoldLast')# Enter desired interpolation method here)
            simbadens.compute_discharge()
            self.simu13.append(simbadens.discharge[ind_transect].total)
           
            # Simulation 14 Bad ensembles (next ens)
            simbadens.transects[ind_transect].w_vel.apply_interpolation(simbadens.transects[ind_transect],target='Ensembles', interp_type= 'ExpandedT')# Enter desired interpolation method here)
            simbadens.compute_discharge()
            self.simu14.append(simbadens.discharge[ind_transect].total)
            
            # Simulation 15 Bad ensembles (linear interpolation)
            simbadens.transects[ind_transect].w_vel.apply_interpolation(simbadens.transects[ind_transect],target='Ensembles', interp_type= 'Linear')# Enter desired interpolation method here)
            simbadens.compute_discharge()
            self.simu15.append(simbadens.discharge[ind_transect].total) 
            
            # simu 18
            simbadcells.transects[ind_transect].w_vel.apply_interpolation(simbadcells.transects[ind_transect],target='Cells', interp_type= 'TRDI')# Enter desired interpolation method here)
            simbadcells.compute_discharge()
            self.simu18.append(simbadcells.discharge[ind_transect].total)
           
            # simu 19
            simbadcells.transects[ind_transect].w_vel.apply_interpolation(simbadcells.transects[ind_transect],target='Cells', interp_type= 'None')# Enter desired interpolation method here)
            simbadcells.compute_discharge()
            self.simu19.append(simbadcells.discharge[ind_transect].total)
            
            # simu 20
            simbadcells.transects[ind_transect].w_vel.apply_interpolation(simbadcells.transects[ind_transect],target='Cells', interp_type= 'Linear')# Enter desired interpolation method here)
            simbadcells.compute_discharge()
            self.simu20.append(simbadcells.discharge[ind_transect].total) 
         
            # simu 21
            simbadcells.transects[ind_transect].w_vel.apply_interpolation(simbadcells.transects[ind_transect],target='Cells', interp_type= 'Linear')# Enter desired interpolation method here)
            simbadcells.compute_discharge()
            self.simu21.append(simbadcells.discharge[ind_transect].total) 
            
            # simu 22
            simbadcells.transects[ind_transect].w_vel.apply_interpolation(simbadcells.transects[ind_transect],target='Cells', interp_type= 'Linear')# Enter desired interpolation method here)
            simbadcells.compute_discharge()
            self.simu22.append(simbadcells.discharge[ind_transect].total) 
            
            # simu 23
            simbadcells.transects[ind_transect].w_vel.apply_interpolation(simbadcells.transects[ind_transect],target='Cells', interp_type= 'Linear')# Enter desired interpolation method here)
            simbadcells.compute_discharge()
            self.simu23.append(simbadcells.discharge[ind_transect].total) 
            
            # simu 24
            simbadcells.transects[ind_transect].w_vel.apply_interpolation(simbadcells.transects[ind_transect],target='Cells', interp_type= 'Linear')# Enter desired interpolation method here)
            simbadcells.compute_discharge()
            self.simu24.append(simbadcells.discharge[ind_transect].total) 
            
            # simu 25
            simbadcells.transects[ind_transect].w_vel.apply_interpolation(simbadcells.transects[ind_transect],target='Cells', interp_type= 'Linear')# Enter desired interpolation method here)
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
    