"""
Created on Jun 16, 2017

@author: gpetrochenkov
"""
import numpy as np

class Cfg(object):
    
    def __init__(self, n_ensembles):
        
        self.ba = np.empty(n_ensembles)
        self.bc = np.empty(n_ensembles)
        self.be_mmps = np.empty(n_ensembles)
        self.bg = np.empty(n_ensembles)
        self.bm = np.empty(n_ensembles)
        self.bp = np.empty(n_ensembles)
        self.bx_dm = np.empty(n_ensembles)
        self.code_reps = np.empty(n_ensembles)
        self.coord_sys = [''] * n_ensembles
        self.cpu_ser_no = np.empty([n_ensembles,8])
        self.cq = np.empty(n_ensembles)
        self.cx = np.empty(n_ensembles)
        self.dist_bin1_cm = np.empty(n_ensembles)
        self.ea_deg = np.empty(n_ensembles)
        self.eb_deg = np.empty(n_ensembles)
        self.sensor_avail = [''] * n_ensembles
        self.ex = [''] * n_ensembles
        self.ez = [''] * n_ensembles
        self.head_src = [''] * n_ensembles
        self.lag_cm = np.empty(n_ensembles)
        self.map_bins = [''] * n_ensembles
        self.n_beams = np.empty(n_ensembles)
        self.pitch_src = [''] * n_ensembles
        self.ref_lay_end_cell = np.empty(n_ensembles)
        self.ref_lay_str_cell = np.empty(n_ensembles)
        self.roll_src = [''] * n_ensembles
        self.sal_src = [''] * n_ensembles
        self.wm = np.empty(n_ensembles)
        self.sos_src = [''] * n_ensembles
        self.temp_src = [''] * n_ensembles
        self.tp_sec = np.empty(n_ensembles)
        self.use_3beam = [''] * n_ensembles
        self.use_pr = [''] * n_ensembles
        self.wa = np.empty(n_ensembles)
        self.wb = np.empty(n_ensembles)
        self.wc = np.empty(n_ensembles)
        self.we_mmps = np.empty(n_ensembles)
        self.wf_cm = np.empty(n_ensembles)
        self.wg_per = np.empty(n_ensembles)
        self.wj = np.empty(n_ensembles)
        self.wn = np.empty(n_ensembles)
        self.wp = np.empty(n_ensembles)
        self.ws_cm = np.empty(n_ensembles)
        self.xdcr_dep_srs =[''] * n_ensembles
        self.xmit_pulse_cm = np.empty(n_ensembles)
        self.lag_near_bottom = np.empty(n_ensembles)
        
        
        
        
        
    