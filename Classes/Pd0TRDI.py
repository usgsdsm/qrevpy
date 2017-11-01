'''
Created on Jun 14, 2017

@author: gpetrochenkov
'''

import numpy as np
from numpy.matlib import repmat
import struct
import os
import re
from Classes.Pd0Classes.Hdr import Hdr
from Classes.Pd0Classes.Inst import Inst
from Classes.Pd0Classes.Cfg import Cfg
from Classes.Pd0Classes.Sensor import Sensor
from Classes.Pd0Classes.Wt import Wt
from Classes.Pd0Classes.Bt import Bt
from Classes.Pd0Classes.Gps import Gps
from Classes.Pd0Classes.Gps2 import Gps2
from Classes.Pd0Classes.Surface import Surface
from Classes.Pd0Classes.AutoMode import AutoMode
from Classes.Pd0Classes.Nmea import Nmea
from scipy.special._ufuncs import sindg, cosdg

class Pd0TRDI(object):
    '''Class to extract data from PD0 files
        Analogous to matlab class: clsPd0TRDI
    '''
    
    def __init__(self, file_name):
        '''Constructor initializing properties'''
        
        self.file_name = None
        self.WR2 = 0
        self.Hdr = None
        self.Inst = None
        self.Cfg = None
        self.Sensor = None
        self.Wt = None
        self.Bt = None
        self.Gps = None
        self.Gps2 = None
        self.Surface = None
        self.AutoMode = None
        self.Nmea = None
        
        self.pd_read(file_name, None)
        
    def create_classes(self, **kargs):
#         print kargs['n_ensembles']
        self.Hdr = Hdr(kargs['n_ensembles'],kargs['n_types'])
        self.Inst = Inst(kargs['n_ensembles'])
        self.Cfg = Cfg(kargs['n_ensembles'])
        self.Sensor = Sensor(kargs['n_ensembles'])
        self.Wt = Wt(kargs['n_bins'],kargs['n_ensembles'],kargs['n_velocities'])
        self.Bt = Bt(kargs['n_ensembles'],kargs['n_velocities'])
        self.Gps = Gps(kargs['n_ensembles'])
        self.Gps2 = Gps2(kargs['n_ensembles'], self.WR2)
        self.Surface = Surface(kargs['n_ensembles'],kargs['n_velocities'],kargs['max_surface_bins'])
        self.AutoMode = AutoMode(kargs['n_ensembles'])
        self.Nmea = Nmea(kargs['n_ensembles'])
        

    def pd_read(self, fullname, args=0):
        n_velocities = 4
        max_surface_bins = 5
        
        if os.path.exists(fullname):
            file_info = os.path.getsize(fullname)
            
            with open(fullname,'rb') as f:
              
                leader_id =  hex(np.fromfile(f,np.uint16, count = 1)[0])
                if leader_id != '0x7f7f':
                    while leader_id != '0x7f7f':
                        f.seek(-1,1)
                        leader_id = hex(np.fromfile(f,np.uint16, count = 1)[0])
                  
                initial_pos = f.tell()-2
                     
                bytes_per_ens = np.fromfile(f, dtype=np.uint16, count=1)[0]
                f.seek(1,1)
                n_types =  np.fromfile(f,np.uint8, count = 1)[0]
                offset = np.fromfile(f,np.uint16, count = 1)[0]
                f.seek(initial_pos+offset+8, 0)
                n_beams = np.fromfile(f,np.uint8,count=1)[0]
                n_bins = np.fromfile(f, np.uint8,count=1)[0]
             
                n_ensembles = self.number_of_ensembles(f, file_info)
#                 print 'n of ens', n_ensembles
            
               
                
                data = {'n_ensembles': n_ensembles,
                        'n_types': n_types,
                        'n_bins': n_bins,
                        'n_beams': n_beams,
                        'max_surface_bins': max_surface_bins,
                        'n_velocities': n_velocities}
                
                self.WR2 = args
                
                self.create_classes(**data)
                
                i_ens = -1
                end_file_check = 0
                end_file = file_info
                i_data_types = 0
                n_data_types = 1
                RR_BtDepth_correction = repmat([np.nan], n_beams, n_ensembles)
                
                f.seek(initial_pos, 0)
                  
                while end_file_check < end_file:
                    
                    
                    leader_id = hex(np.fromfile(f,np.uint16, count = 1)[0])     
                    
                    if i_data_types >= n_data_types and leader_id != '0x7f7f':
                        leader_id = '0x9999'
                          
                
                    if leader_id == '0x7f7f':
                          
                        
                        i2022 = 0
                        file_loc = f.tell() - 2
                          
                        if file_loc+bytes_per_ens > end_file and i_ens >= n_ensembles:
                            end_file_check = end_file+1
                        else:
                            i_data_types = 0
                            store_file_loc = f.tell()
                            bytes_per_ens = np.fromfile(f, np.uint16, count=1)[0]
                            
                            if self.check_sum(f,file_loc, bytes_per_ens):
                                f.seek(file_loc+5,0)
                                n_data_types = np.fromfile(f, np.uint8, count=1)[0]
                                data_offsets = np.fromfile(f, np.uint16, count=n_data_types)
                                  
                                while i_data_types+1 <= n_data_types and leader_id != '0x80':
                                    f.seek(data_offsets[i_data_types]+file_loc,0)
                                    leader_id = hex(np.fromfile(f, np.uint16, count=1)[0])
                                    i_data_types += 1
                                      
                                if i_ens > -1 and leader_id == '0x80':
                                    ens_num = np.fromfile(f, np.uint16, count=1)[0]
                                    ens_num_diff = ens_num - self.Sensor.num[i_ens]
                                    if ens_num_diff > 1:
                                        for nn in range(0,int(ens_num_diff-1)):
                                            if i_ens < n_ensembles:
                                                self.Sensor.num[i_ens] = self.Sensor.num[i_ens-1]+1
                                                i_ens += 1
                                    elif ens_num_diff < 1:
                                        i_ens -= 1
                                          
                            else:
                                self.bad_check_sum(f, file_loc)
                                i_data_types = -1
                                  
                            f.seek(store_file_loc, 0)
                            i_data_types = 0
#                            j100,j101,j102,j103 = 0,0,0,0
                            j100,j101,j102,j103 = -1,-1,-1,-1
                            i_ens += 1
                              
                            self.Hdr.bytes_per_ens[i_ens] = np.fromfile(f,np.uint16, count=1)[0]
                              
                            #if checksum is valid read header data
                            if self.check_sum(f,file_loc, int(self.Hdr.bytes_per_ens[i_ens])):
                                f.seek(file_loc+5,0)
                                self.Hdr.n_data_types[i_ens] = np.fromfile(f, np.uint8, count=1)[0]
                                test = np.fromfile(f, np.uint16,count=int(self.Hdr.n_data_types[i_ens]))
                                
                            
                                if test.shape[0] > self.Hdr.data_offsets.shape[1]:
                                    self.Hdr.data_offsets.resize(n_ensembles, test.shape[0])
#                                 
                                self.Hdr.data_offsets[i_ens,0:int(self.Hdr.n_data_types[i_ens])] = test[0:int(self.Hdr.n_data_types[i_ens])]
                                  
                                self.end_reading(f, file_loc, i_data_types, i_ens, bytes_per_ens)
                                    
                            else:
                                self.bad_check_sum(f, file_loc)
                                i_data_types = -1
                                
                    
                    elif leader_id == '0x0':   
                        
                        i_data_types += 1
                        self.Inst.firm_ver[i_ens] = np.fromfile(f, np.uint8, count=1)[0]
                        self.Inst.firm_ver[i_ens] = self.Inst.firm_ver[i_ens] + np.fromfile(f, np.uint8, count=1)[0] / 100
                        
                        bitls = np.fromfile(f,np.uint8, count=1)[0]
                        bitls = "{0:08b}".format(bitls)
                        
                        val = int(bitls[5:], 2)
                        if val == 0: self.Inst.freq[i_ens] = 75
                        elif val == 1: self.Inst.freq[i_ens] = 150
                        elif val == 2: self.Inst.freq[i_ens] = 300
                        elif val == 3: self.Inst.freq[i_ens] = 600
                        elif val == 4: self.Inst.freq[i_ens] = 1200
                        elif val == 5: self.Inst.freq[i_ens] = 2400
                        else: self.Inst.freq[i_ens] = np.nan
                        
                        val = int(bitls[4], 2)
                        if val == 0: self.Inst.pat[i_ens] = 'Concave'
                        elif val == 1: self.Inst.pat[i_ens] = 'Convex'
                        else: self.Inst.pat[i_ens] = 'N/a'
                        
                        self.Inst.sensor_CFG[i_ens] = int(bitls[2:3],2) + 1
                         
                        val = int(bitls[1], 2)
                        if val == 0:self.Inst.xducer[i_ens]='Not Attached'
                        elif val == 1:self.Inst.xducer[i_ens]='Attached' 
                        else:self.Inst.xducer[i_ens] = 'N/a'
                        
                        val = int(bitls[0], 2)
                        if val == 0: self.Sensor.orient[i_ens] = 'Down'
                        elif val == 1: self.Sensor.orient[i_ens] = 'Up'
                        else: self.Sensor.orient[i_ens] = 'N/a'
                        
                        bitms = np.fromfile(f, np.uint8, count=1)[0]
                        bitms= "{0:08b}".format(bitms)
                        
                        val = int(bitms[6:], 2)
                        if val == 0:self.Inst.beam_ang[i_ens] = 15
                        elif val == 1:self.Inst.beam_ang[i_ens] = 20
                        elif val == 2:self.Inst.beam_ang[i_ens] = 30
                        elif val == 3:self.Inst.beam_ang[i_ens] = np.nan
                        else:self.Inst.beam_ang[i_ens] = np.nan
                        
                        val = int(bitms[:4], 2)
                        if val == 4:self.Inst.beams[i_ens] = 4
                        elif val == 5:
                            self.Inst.beams[i_ens] = 5
                            self.Inst.demod[i_ens] = 1
                        elif val == 15:
                            self.Inst.beams[i_ens] = 5
                            self.Inst.demod[i_ens] = 2
                        else:
                            self.Inst.beams[i_ens] = np.nan
                            self.Inst.demod[i_ens] = np.nan
                            
                        
                        val = np.fromfile(f, np.uint8, count=1)[0]
                        if val == 0:self.Inst.data_type[i_ens] = 'Real'
                        else:self.Inst.data_type[i_ens] = 'Simu'
                        
                        
                        #Position file pointer and read configuration information
                        f.seek(1,1)
                        self.Cfg.n_beams[i_ens] = np.fromfile(f, np.uint8, count=1)[0]
                        self.Cfg.wn[i_ens] = np.fromfile(f, np.uint8, count=1)[0]
                        self.Cfg.wp[i_ens] = np.fromfile(f, np.uint16, count=1)[0]
                        self.Cfg.ws_cm[i_ens] = np.fromfile(f, np.uint16, count=1)[0]
                        self.Cfg.wf_cm[i_ens] = np.fromfile(f, np.uint16, count=1)[0]
                        self.Cfg.wm[i_ens] = np.fromfile(f, np.uint8, count=1)[0]
                        self.Cfg.wc[i_ens] = np.fromfile(f, np.uint8, count=1)[0]
                        self.Cfg.code_reps[i_ens] = np.fromfile(f, np.uint8, count=1)[0]
                        self.Cfg.wg_per[i_ens] = np.fromfile(f, np.uint8, count=1)[0]
                        self.Cfg.we_mmps[i_ens] = np.fromfile(f, np.uint16, count=1)[0]
                        self.Cfg.tp_sec[i_ens] = np.sum(np.fromfile(f, np.uint8, count=3) * np.array([60, 1, .01]))
                        self.Cfg.ex[i_ens] = "{0:08b}".format(ord(f.read(1)))
                        
                        
                        val = int(self.Cfg.ex[i_ens][3:5], 2)
                        if val == 0: self.Cfg.coord_sys[i_ens] = 'Beam'
                        elif val == 1: self.Cfg.coord_sys[i_ens] = 'Inst'
                        elif val == 2: self.Cfg.coord_sys[i_ens] = 'Ship'
                        elif val == 3: self.Cfg.coord_sys[i_ens] = 'Earth'
                        else: self.Cfg.coord_sys[i_ens] = "N/a"
                        
                        val = int(self.Cfg.ex[i_ens][5], 2)
                        if val == 0: self.Cfg.use_pr = 'No'
                        elif val == 1: self.Cfg.use_pr = 'Yes'
                        else: self.Cfg.use_pr = 'N/a'
                        
                        val = int(self.Cfg.ex[i_ens][6], 2)
                        if val == 0: self.Cfg.use_3beam = 'No'
                        elif val == 1: self.Cfg.use_3beam = 'Yes'
                        else: self.Cfg.use_3beam = 'N/a'
                        
                        val = int(self.Cfg.ex[i_ens][7], 2)
                        if val == 0: self.Cfg.map_bins = 'No'
                        elif val == 1: self.Cfg.map_bins = 'Yes'
                        else: self.Cfg.map_bins = 'N/a'
                        
                        self.Cfg.ea_deg[i_ens] = np.fromfile(f, np.int16, count=1)[0] * 0.01
                        self.Cfg.ea_deg[i_ens] = np.fromfile(f, np.uint16, count=1)[0] * 0.01
                        self.Cfg.ez[i_ens] = "{0:08b}".format(np.fromfile(f, np.uint8, count=1)[0])
                        
                        val = int(self.Cfg.ez[i_ens][:2],2)
                        if val == 0: self.Cfg.sos_src[i_ens] = 'Manual EC'
                        elif val == 1: self.Cfg.sos_src[i_ens] = 'Calculated'
                        elif val == 3: self.Cfg.sos_src[i_ens] = 'SVSS Sensor'
                        else: self.Cfg.sos_src[i_ens] = 'N/a'
                        
                        val = int(self.Cfg.ez[i_ens][2], 2)
                        if val == 0: self.Cfg.xdcr_dep_srs[i_ens] = 'Manual ED'
                        if val == 1: self.Cfg.xdcr_dep_srs[i_ens] = 'Sensor'
                        else: self.Cfg.xdcr_dep_srs[i_ens] = 'N/a'
                        
                        val = int(self.Cfg.ez[i_ens][3], 2)
                        if val == 0: self.Cfg.head_src[i_ens] = 'Manual EH'
                        if val == 1: self.Cfg.head_src[i_ens] = 'Int. Sensor'
                        else: self.Cfg.head_src[i_ens] = 'N/a'
                        
                        val = int(self.Cfg.ez[i_ens][4], 2)
                        if val == 0: self.Cfg.pitch_src[i_ens] = 'Manual EP'
                        if val == 1: self.Cfg.pitch_src[i_ens] = 'Int. Sensor'
                        else: self.Cfg.pitch_src[i_ens] = 'N/a'
                        
                        val = int(self.Cfg.ez[i_ens][5], 2)
                        if val == 0: self.Cfg.roll_src[i_ens] = 'Manual ER'
                        if val == 1: self.Cfg.roll_src[i_ens] = 'Int. Sensor'
                        else: self.Cfg.roll_src[i_ens] = 'N/a'
                        
                        val = int(self.Cfg.ez[i_ens][6], 2)
                        if val == 0: self.Cfg.xdcr_dep_srs[i_ens] = 'Manual ES'
                        if val == 1: self.Cfg.xdcr_dep_srs[i_ens] = 'Int. Sensor'
                        else: self.Cfg.xdcr_dep_srs[i_ens] = 'N/a'
                        
                        val = int(self.Cfg.ez[i_ens][7], 2)
                        if val == 0: self.Cfg.temp_src[i_ens] = 'Manual ET'
                        if val == 1: self.Cfg.temp_src[i_ens] = 'Int. Sensor'
                        else: self.Cfg.temp_src[i_ens] = 'N/a'
                        
                        self.Cfg.sensor_avail[i_ens] = "{0:08b}".format(np.fromfile(f, np.uint8, count=1)[0]) 
                        self.Cfg.dist_bin1_cm[i_ens] = np.fromfile(f, np.uint16, count=1)[0]
                        self.Cfg.xmit_pulse_cm[i_ens] = np.fromfile(f, np.uint16, count=1)[0]
                        self.Cfg.ref_lay_str_cell[i_ens] = np.fromfile(f, np.uint8, count=1)[0]
                        self.Cfg.ref_lay_end_cell[i_ens] = np.fromfile(f, np.uint8, count=1)[0]
                        self.Cfg.wa[i_ens] = np.fromfile(f, np.uint8, count=1)[0]
                        self.Cfg.cx[i_ens] = np.fromfile(f, np.uint8, count=1)[0]
                        self.Cfg.lag_cm[i_ens] = np.fromfile(f, np.uint16, count=1)[0]
                        self.Cfg.cpu_ser_no[i_ens] = np.fromfile(f, np.uint8, count=1)[0]
                        self.Cfg.wb[i_ens] = np.fromfile(f, np.uint8, count=1)[0]
                        self.Cfg.cq[i_ens] = np.fromfile(f, np.uint8, count=1)[0]
                        
                        #check if more data types need to be read and position the pointer
                        self.end_reading(f, file_loc, i_data_types, i_ens, bytes_per_ens)
                        
                        
                    elif leader_id == '0x80':
                        #update the data types counter
                        i_data_types += 1
                        #read instrument clock and sensor data
                        self.Sensor.num[i_ens] = np.fromfile(f, np.uint16, count=1)[0]
                        self.Sensor.date_not_y2k[i_ens,:] = np.fromfile(f, np.uint8, count=3)
                        self.Sensor.time[i_ens, :] = np.fromfile(f, np.uint8, count=4)
                        self.Sensor.num_fact[i_ens] = np.fromfile(f, np.uint8, count=1)[0]
                        self.Sensor.num_tot[i_ens] = self.Sensor.num[i_ens] + self.Sensor.num_fact[i_ens]*65535
                        self.Sensor.bit_test[i_ens] = np.fromfile(f, np.uint16, count=1)[0]
                        self.Sensor.sos_mps[i_ens] = np.fromfile(f, np.uint16, count=1)[0]
                        self.Sensor.xdcr_depth_dm[i_ens] = np.fromfile(f, np.uint16, count=1)[0]
                        self.Sensor.heading_deg[i_ens] = np.fromfile(f, np.uint16, count=1)[0] / 100.
                        self.Sensor.pitch_deg[i_ens] = np.fromfile(f, np.int16, count=1)[0] /100.
                        self.Sensor.roll_deg[i_ens] = np.fromfile(f, np.int16, count=1)[0] / 100.
                        self.Sensor.salinity_ppt[i_ens] = np.fromfile(f, np.uint16, count=1)[0]
                        self.Sensor.temperature_deg_c[i_ens] = np.fromfile(f, np.int16, count=1)[0] / 100.
                        self.Sensor.mpt_msc[i_ens,:] = np.fromfile(f, np.uint8, count=3)
                        self.Sensor.heading_std_dev_deg[i_ens] = np.fromfile(f, np.uint8, count=1)[0]
                        self.Sensor.pitch_std_dev_deg[i_ens] = np.fromfile(f, np.uint8, count=1)[0] / 10.
                        self.Sensor.roll_std_dev_deg[i_ens] = np.fromfile(f, np.uint8, count=1) / 10.
                        self.Sensor.xmit_current[i_ens] = np.fromfile(f, np.uint8, count=1)[0]
                        self.Sensor.xmit_voltage[i_ens] = np.fromfile(f, np.uint8, count=1)[0]
                        self.Sensor.ambient_temp[i_ens] = np.fromfile(f, np.uint8, count=1)[0]
                        self.Sensor.pressure_pos[i_ens] = np.fromfile(f, np.uint8, count=1)[0]
                        self.Sensor.pressure_neg[i_ens] = np.fromfile(f, np.uint8, count=1)[0]
                        self.Sensor.attitude_temp[i_ens] = np.fromfile(f, np.uint8, count=1)[0]
                        self.Sensor.attitude[i_ens] = np.fromfile(f, np.uint8, count=1)[0]
                        self.Sensor.contam_sensor[i_ens] = np.fromfile(f, np.uint8, count=1)[0]
                        self.Sensor.error_status_word[i_ens] = ["{0:08b}".format(x) for x in np.fromfile(f, np.uint8, count=4)]
                        f.seek(2,1)
                        self.Sensor.pressure_pascal[i_ens] = np.fromfile(f, np.uint32, count=1)[0]
                        self.Sensor.pressure_var_pascal[i_ens] = np.fromfile(f, np.uint32, count=1)[0]
                        f.seek(1,1)
                        self.Sensor.date_y2k[i_ens,:] = np.fromfile(f, np.uint8, count=4)
                        self.Sensor.time_y2k[i_ens,:] = np.fromfile(f, np.uint8, count=4)
                        self.Sensor.date[i_ens,:] = self.Sensor.date_not_y2k[i_ens,:]
                        self.Sensor.date[i_ens,0] = self.Sensor.date_y2k[i_ens,1] * 100 + self.Sensor.date_y2k[i_ens,2]
                        
                        self.Cfg.lag_near_bottom[i_ens] = np.fromfile(f, np.uint8, count=1)[0]
                        
                        self.end_reading(f, file_loc, i_data_types, i_ens, bytes_per_ens)
                            
                     
                    elif leader_id == '0x100':
                        
                        #update data types counter
                        i_data_types += 1
                        if self.Cfg.wn[i_ens] > self.Wt.vel_mps.shape[1]:
                            self.Wt.vel_mps.resize([self.Wt.vel_mps.shape[0], int(self.Cfg.wn[i_ens]), self.Wt.vel_mps.shape[2]])
                            
                        dummy = np.fromfile(f, np.int16, count=int(self.Cfg.wn[i_ens] *4))
                        dummy = np.reshape(dummy,[int(self.Cfg.wn[i_ens]), n_velocities])
                        self.Wt.vel_mps[:n_velocities,:int(self.Cfg.wn[i_ens]),i_ens] = dummy.T
                        
                        self.end_reading(f, file_loc, i_data_types, i_ens, bytes_per_ens)
                            
                    elif leader_id == '0x200':
                        
                        i_data_types += 1
                        if self.Cfg.wn[i_ens] > self.Wt.corr.shape[1]:
                            self.Wt.corr.resize([self.Wt.corr.shape[0], int(self.Cfg.wn[i_ens]), self.Wt.corr.shape[2]])
                        dummy = np.fromfile(f, np.uint8, count=int(self.Cfg.wn[i_ens]*4))
                        dummy = np.reshape(dummy, [int(self.Cfg.wn[i_ens]), n_velocities])
                        self.Wt.corr[:n_velocities, :int(self.Cfg.wn[i_ens]),i_ens] = dummy.T
                        
                        self.end_reading(f, file_loc, i_data_types, i_ens, bytes_per_ens)
                        
                    elif leader_id == '0x300':
                        
                        i_data_types += 1
                        if self.Cfg.wn[i_ens] > self.Wt.rssi.shape[1]:
                            self.Wt.rssi.resize([self.Wt.rssi.shape[0], int(self.Cfg.wn[i_ens]), self.Wt.rssi.shape[2]])
                        dummy = np.fromfile(f, np.uint8, count=int(self.Cfg.wn[i_ens]*4))
                        dummy = np.reshape(dummy, [int(self.Cfg.wn[i_ens]), n_velocities])
                        self.Wt.rssi[:n_velocities, :int(self.Cfg.wn[i_ens]),i_ens] = dummy.T
                        
                        self.end_reading(f, file_loc, i_data_types, i_ens, bytes_per_ens)
                        
                    elif leader_id == '0x400':
                        
                        i_data_types += 1
                        if self.Cfg.wn[i_ens] > self.Wt.pergd.shape[1]:
                            self.Cfg.wn[i_ens] = self.Wt.pergd.shape[1]
                        dummy = np.fromfile(f, np.uint8, count=int(self.Cfg.wn[i_ens]*4))
                        dummy = np.reshape(dummy, [int(self.Cfg.wn[i_ens]), n_velocities])
                        self.Wt.pergd[:n_velocities, :int(self.Cfg.wn[i_ens]),i_ens] = dummy.T
                        
                        self.end_reading(f, file_loc, i_data_types, i_ens, bytes_per_ens)
                        
                    elif leader_id == '0x600':
                        
                        i_data_types +=1
                        
                        #read bottom track configuration data
                        self.Cfg.bp[i_ens] = np.fromfile(f, np.uint16, count=1)[0]
                        long1 = np.fromfile(f, np.uint16, count=1)[0]
                        self.Cfg.bc[i_ens] = np.fromfile(f, np.uint8, count=1)[0]
                        self.Cfg.ba[i_ens] = np.fromfile(f, np.uint8, count=1)[0]
                        self.Cfg.bg[i_ens] = np.fromfile(f, np.uint8, count=1)[0]
                        self.Cfg.bm[i_ens] = np.fromfile(f, np.uint8, count=1)[0]
                        self.Cfg.be_mmps[i_ens] = np.fromfile(f, np.uint16, count=1)[0]
                        
                        #read winriver 10.06 fromat GPS data
                        self.Gps.lat_deg[i_ens] = (np.fromfile(f, np.int32, count=1)[0]/2)**31 * 180
                        
                        #read the Least Significant Bytes for beam depths
                        dummy = np.fromfile(f, np.uint16, count=4)
                        self.Bt.depth_m[0:4,i_ens] = dummy.T
                        
                        #read bottom-track velocities
                        dummy = np.fromfile(f, np.int16, count=4)
                        self.Bt.vel_mps[0:4, i_ens] = dummy.T
                        
                        #read bottom-track correlations
                        dummy =  np.fromfile(f, np.uint8, count=4)
                        self.Bt.corr[0:4, i_ens] = dummy.T
                        
                        #read bottom-track evaluation amplitude
                        dummy = np.fromfile(f, np.uint8, count=4)
                        self.Bt.eval_amp[0:4, i_ens] = dummy.T
                        
                        #read bottom-track percent good
                        dummy = np.fromfile(f, np.uint8, count=4)
                        self.Bt.pergd[0:4,i_ens] = dummy.T
                        
                        #read WinRiver 10.06 format GPS data
                        dummy = np.fromfile(f, np.uint16, count=1)[0]
                        if dummy != -32768:
                            self.Gps.alt_m[i_ens] = (dummy-32768)/10
                        else:
                            self.Gps.altm[i_ens] = np.nan
                            
                        long2 = np.fromfile(f, np.uint16, count=1)[0]
                        self.Gps.long_deg[i_ens] = ((long1+long2*2**16)/2**31)*180
                        if self.Gps.long_deg[i_ens] > 180:
                            self.Gps.long_deg[i_ens] -= 360
                            
                        self.Bt.ext_depth_cm[i_ens] = np.fromfile(f, np.int16, count=1)[0]
                        dummy=np.fromfile(f,np.int16,count=1)[0]
                        if dummy != -32768:
                            self.Gps.gga_vel_e_mps[i_ens]=dummy*-1 / 1000
                        else:
                            self.Gps.gga_vel_e_mps[i_ens]= np.nan
                            
                        dummy=np.fromfile(f,np.int16,count=1)[0]
                        if dummy != -32768:
                            self.Gps.gga_vel_n_mps[i_ens]=dummy*-1 / 1000
                        else:
                            self.Gps.gga_vel_n_mps[i_ens]= np.nan
                            
                        dummy=np.fromfile(f,np.int16,count=1)[0]
                        if dummy != -32768:
                            self.Gps.vtg_vel_e_mps[i_ens]=dummy*-1 / 1000
                        else:
                            self.Gps.vtg_vel_e_mps[i_ens]= np.nan
                            
                        dummy=np.fromfile(f,np.int16,count=1)[0]
                        if dummy != -32768:
                            self.Gps.vtg_vel_n_mps[i_ens]=dummy*-1 / 1000
                        else:
                            self.Gps.vtg_vel_n_mps[i_ens]= np.nan
                            
                        dummy = np.fromfile(f,np.uint8, count=1)[0]
                        if dummy != 0:
                            self.Gps.gsa_v_dop[i_ens]= dummy
                        dummy = np.fromfile(f,np.uint8, count=1)[0]
                        if dummy != 0:
                            self.Gps.gsa_p_dop[i_ens]= dummy
                        dummy = np.fromfile(f,np.uint8, count=1)[0]
                        if dummy != 0:
                            self.Gps.gga_n_stats[i_ens,0]= dummy
                            
                        f.seek(1,1)
                        self.Gps.gsa_sat[i_ens,4] = np.fromfile(f, np.uint8, count=1)[0]
                        self.Gps.gsa_sat[i_ens,5] = np.fromfile(f, np.uint8, count=1)[0]
                        self.Gps.gga_diff[i_ens] = np.fromfile(f,np.uint8, count=1)[0]
                        
                        dummy = np.fromfile(f,np.uint8, count=1)[0]
                        if dummy != 0:
                            self.Gps.gga_hdop[i_ens]= dummy/10
                        
                        self.Gps.gsa_sat[i_ens,0] = np.fromfile(f,np.uint8, count=1)[0]
                        self.Gps.gsa_sat[i_ens,1] = np.fromfile(f,np.uint8, count=1)[0]
                        self.Gps.gsa_sat[i_ens,2] = np.fromfile(f,np.uint8, count=1)[0]
                        self.Gps.gsa_sat[i_ens,3] = np.fromfile(f,np.uint8, count=1)[0]
                        
                        #read bx configuration setting
                        self.Cfg.bx_dm[i_ens] = np.fromfile(f, np.uint16, count=1)[0]
                        
                        #Read bottom-tracking RSSI
                        self.Bt.rssi[0,i_ens] = np.fromfile(f,np.uint8, count=1)[0]
                        self.Bt.rssi[1, i_ens] = np.fromfile(f,np.uint8, count=1)[0]
                        self.Bt.rssi[2, i_ens] = np.fromfile(f,np.uint8, count=1)[0]
                        self.Bt.rssi[3, i_ens] = np.fromfile(f,np.uint8, count=1)[0]
                        
                        #read wj configuration setting
                        self.Cfg.wj[i_ens] = np.fromfile(f,np.uint8, count=1)[0]
                        
                        #Read most significant byte and compute beam depths
                        dummy = np.fromfile(f, np.uint8, count=1)[0]
                        RR_BtDepth_correction[0:4,i_ens] = dummy.T*2e16/100
                        
                        self.end_reading(f, file_loc, i_data_types, i_ens, bytes_per_ens)
                        
                        
                    elif leader_id == '0x2022':
                        i2022 += 1
                        i_data_types += 1
                        specific_id = np.fromfile(f, np.int16, count=1)[0]
                        msg_size = np.fromfile(f, np.int16, count=1)[0]
                        delta_time = np.fromfile(f, np.double, count=1)[0]
                        
                        if specific_id == 100:
                            j100 += 1
                            
                            if j100>20:
                                self.Gps2.gga_delta_time[:,j100] = np.nan
                                self.Gps2.gga_header[:,j100] = ''
                                self.Gps2.lat_deg[:,j100] = np.nan
                                self.Gps2.utc[:,j100] = np.nan
                                self.Gps2.lat_ref[:][j100] = ''
                                self.Gps2.lon_ref[:][j100] = ''
                                self.Gps2.lon_deg[:,j100] = np.nan
                                self.Gps2.corr_qual[:,j100] = np.nan
                                self.Gps2.num_sats[:,j100] = np.nan
                                self.Gps2.hdop[:,j100] = np.nan
                                self.Gps2.alt[:,j100] = np.nan
                                self.Gps2.alt_unit[:][j100] = ''
                                self.Gps2.geoid[:,j100] = np.nan
                                self.Gps2.geoid_unit[:][j100] = ''
                                self.Gps2.d_gps_age[:,j100] = np.nan
                                self.Gps2.ref_stat_id[:,j100] = np.nan
                                
                            self.Gps2.gga_delta_time[i_ens,j100] = delta_time
                            # this line was incomplete I used later example to add the following line - DSM
                            # self.Gps2.gga_header[i_ens,j100]
                            self.Gps2.gga_header[i_ens][j100]=''.join([chr(x) for x in f.read(10)])
                            temp = re.match('[0-9]+\.[0-9]+', ''.join(chr(x) for x in f.read(10)))
                            try:
                                self.Gps2.utc[i_ens,j100] = float(temp.group(0))
                            except:
                                self.Gps2.utc[i_ens,j100] = np.nan
                                
                            self.Gps2.lat_deg[i_ens,j100] = np.fromfile(f, np.float64, count=1)[0]
                            self.Gps2.lat_ref[i_ens][j100] = chr(f.read(1)[0])
                            self.Gps2.lon_deg[i_ens,j100] = np.fromfile(f, np.float64, count=1)[0]
                            self.Gps2.lon_ref[i_ens][j100] = chr(f.read(1)[0])
                            self.Gps2.corr_qual[i_ens,j100] = np.fromfile(f, np.uint8, count=1)[0]
                            self.Gps2.num_sats[i_ens,j100] = np.fromfile(f, np.uint8, count=1)[0]
                            self.Gps2.hdop[i_ens,j100] = np.fromfile(f, np.float32, count=1)[0]
                            self.Gps2.alt[i_ens,j100] = np.fromfile(f, np.float32, count=1)[0]
                            self.Gps2.alt_unit[i_ens][j100] = chr(f.read(1)[0])
                            self.Gps2.geoid[i_ens,j100] = np.fromfile(f, np.float32, count=1)[0]
                            self.Gps2.geoid_unit[i_ens][j100] = chr(f.read(1)[0])
                            self.Gps2.d_gps_age[i_ens,j100] = np.fromfile(f, np.float32, count=1)[0]
                            self.Gps2.ref_stat_id[i_ens,j100] = np.fromfile(f,np.int16,count=0)[0]
                            
                        elif specific_id == 101:
                            j101 += 1
                            
                            if j101>20:
                                self.Gps2.vtg_delta_time[:,j101] = np.nan
                                self.Gps2.vtg_header[:][j101] = ''
                                self.Gps2.course_true[:,j101] = np.nan
                                self.Gps2.true_indicator[:][j101] = ''
                                self.Gps2.course_mag[:,j101] = np.nan
                                self.Gps2.mag_indicator[:][j101] = ''
                                self.Gps2.speed_knots[:,j101] = np.nan
                                self.Gps2.knots_indicator[:][j101] = ''
                                self.Gps2.speed_k_mph[:,j101] = np.nan
                                self.Gps2.kmph_indicator[:][j101] = ''
                                self.Gps2.mode_indicator[:][j101] = ''
                            
                            self.Gps2.vtg_delta_time[i_ens,j101] = delta_time
                            self.Gps2.vtg_header[i_ens][j101] = ''.join([chr(x) for x in f.read(10)])
                            self.Gps2.course_true[i_ens,j101] = np.fromfile(f, np.float32, count=1)[0]
                            self.Gps2.true_indicator[i_ens][j101] = chr(f.read(1)[0])
                            self.Gps2.course_mag[i_ens,j101] = np.fromfile(f, np.float32, count=1)[0]
                            self.Gps2.mag_indicator[i_ens][j101] = chr(f.read(1)[0])
                            self.Gps2.speed_knots[i_ens,j101] = np.fromfile(f, np.float32, count=1)[0]
                            self.Gps2.kmph_indicator[i_ens][j101] = chr(f.read(1)[0])
                            self.Gps2.mode_indicator[i_ens][j101] = chr(f.read(1)[0])
                        
                        elif specific_id == 102:
                            j102 += 1
                            
                            if j102 > 20:
                                self.Gps2.dbt_delta_time[:,j102] = np.nan
                                self.Gps2.dbt_header[:][j102] = ''
                                self.Gps2.depth_ft[:,j102] = np.nan
                                self.Gps2.ft_indicator[:][j102] = ''
                                self.Gps2.depth_m[:,j102] = np.nan
                                self.Gps2.m_indicator[:][j102] = ''
                                self.Gps2.depth_fath[:,j102] = np.nan
                                self.Gps2.fath_indicator[:][j102] = ''
                                
                            self.Gps2.dbt_delta_time[i_ens,j102] = delta_time
                            self.Gps2.dbt_header[i_ens][j102] = ''.join([chr(x) for x in f.read(10)])
                            self.Gps2.depth_ft[i_ens,j102] = np.fromfile(f, np.float32, count=1)[0]
                            self.Gps2.ft_indicator[i_ens][j102] = chr(f.read(1)[0])
                            self.Gps2.depth_m[i_ens,j102] = np.fromfile(f, np.float32, count=1)[0]
                            self.Gps2.m_indicator[i_ens][j102] = chr(f.read(1)[0])
                            self.Gps2.depth_fath[i_ens, j102] = np.fromfile(f, np.float32, count=1)[0]
                            self.Gps2.fath_indicator[i_ens][j102] = chr(f.read(1)[0])
                            
                        elif specific_id == 103:
                            j103 += 1
                            
                            if j103 > 20:
                                self.Gps2.hdt_delta_time[:, j103] = np.nan
                                self.Gps2.hdt_header[:][j103] = ''
                                self.Gps2.heading_deg[:,j103] = np.nan
                                self.Gps2.h_true_indicator[:][j103] = ''
                            
                            self.Gps2.hdt_delta_time[i_ens,j103] = delta_time
                            self.Gps2.hdt_header[i_ens][j103] = ''.join([chr(x) for x in f.read(10)])
                            self.Gps2.heading_deg[i_ens,j103] = np.fromfile(f, np.double, count=1)[0]
                            self.Gps2.h_true_indicator[i_ens][j103] = chr(f.read(1)[0])
                            
                        elif specific_id == 104:
                            j100 +=1
                            
                            if j100>20:
                                self.Gps2.gga_delta_time[:,j100] = np.nan
                                self.Gps2.gga_header[:][j100] = ''
                                self.Gps2.utc[:,j100] = np.nan
                                self.Gps2.lat_ref[:][j100] = ''
                                self.Gps2.lon_deg[:,j100] = np.nan
                                self.Gps2.lon_ref[:][j100] = ''
                                self.Gps2.corr_qual[:,j100] = np.nan
                                self.Gps2.num_sats[:,j100] = np.nan
                                self.Gps2.hdop[:,j100] = np.nan
                                self.Gps2.alt[:,j100] = np.nan
                                self.Gps2.alt_unit[:][j100] = ''
                                self.Gps2.geoid[:,j100] = np.nan
                                self.Gps2.geoid_unit[:][j100] = ''
                                self.Gps2.d_gps_age[:,j100] = np.nan
                                self.Gps2.ref_stat_id[:,j100] = np.nan
                                
                            self.Gps2.lat_deg[i_ens,j100] = delta_time
                            self.Gps2.gga_header[i_ens] = ''.join([chr(x) for x in f.read(7)])
                            try:
                                temp = ''.join([chr(x) for x in f.read(10)])
                                self.Gps2.utc[i_ens,j100] = float(re.match('[0-9]+\.[0-9]+', temp))
                            except:
                                self.Gps2.utc[i_ens,j100] = np.nan
                            self.Gps2.lat_deg[i_ens,j100] = np.fromfile(f, np.float64, count=1)[0]
                            self.Gps2.lat_ref[i_ens][j100] = chr(f.read(1)[0])
                            self.Gps2.lon_deg[i_ens, j100] = np.fromfile(f, np.float64, count=1)[0]
                            self.Gps2.lon_ref[i_ens][j100] = chr(f.read(1)[0])
                            self.Gps2.corr_qual[i_ens,j100] = np.fromfile(f, np.uint8, count=1)[0]
                            self.Gps2.num_sats[i_ens,j100] = np.fromfile(f, np.uint8, count=1)[0]
                            self.Gps2.hdop[i_ens,j100] = np.fromfile(f, np.float32, count=1)[0]
                            self.Gps2.alt[i_ens,j100] = np.fromfile(f, np.float32, count=1)[0]
                            self.Gps2.alt_unit[i_ens][j100] = chr(f.read(1)[0])
                            self.Gps2.geoid[i_ens,j100] = np.fromfile(f, np.float32, count=1)[0]
                            self.Gps2.geoid_unit[i_ens][j100] = chr(f.read(1)[0])
                            self.Gps2.d_gps_age[i_ens,j100] = np.fromfile(f, np.float32, count=1)[0]
                            self.Gps2.ref_stat_id[i_ens,j100] = np.fromfile(f, np.int16, count=1)[0]
                            
                        elif specific_id == 105:
                            j101 += 1
                            
                            if j101 > 20:
                                self.Gps2.vtg_delta_time[:,j101] = np.nan
                                self.Gps2.vtg_header[:][j101] = ''
                                self.Gps2.course_true[:,j101] = np.nan
                                self.Gps2.true_indicator[:][j101] = ''
                                self.Gps2.course_mag[:,j101] = np.nan
                                self.Gps2.mag_indicator[:][j101] = ''
                                self.Gps2.speed_knots[:,j101] = np.nan
                                self.Gps2.knots_indicator[:][j101] = ''
                                self.Gps2.speed_k_mph[:,j101] = np.nan
                                self.Gps2.kmph_indicator[:][j101] = ''
                                self.Gps2.mode_indicator[:][j101] = ''
                                
                            self.Gps2.vtg_delta_time[i_ens,j101] = delta_time
                            self.Gps2.vtg_header[i_ens][j101] = ''.join([chr(x) for x in f.read(7)])
                            self.Gps2.course_true[i_ens,j101] = np.fromfile(f, np.float32, count=1)[0]
                            self.Gps2.true_indicator[i_ens][j101] = chr(f.read(1)[0])
                            self.Gps2.course_mag[i_ens,j101] = np.fromfile(f, np.float32, count=1)[0]
                            self.Gps2.mag_indicator[i_ens][j101] = chr(f.read(1)[0])
                            self.Gps2.speed_knots[i_ens,j101] = np.fromfile(f, np.float32, count=1)[0]
                            self.Gps2.knots_indicator[i_ens][j101] = chr(f.read(1)[0])
                            self.Gps2.speed_k_mph[i_ens,j101] = np.fromfile(f, np.float32, count=1)[0]
                            self.Gps2.kmph_indicator[i_ens][j101] = chr(f.read(1)[0])
                            self.Gps2.mode_indicator[i_ens][j101] = chr(f.read(1)[0])
                            
                        elif specific_id == 106:
                            j102 += 1
                            
                            if j102 > 20:
                                self.Gps2.dbt_delta_time[:,j102] = np.nan
                                self.Gps2.dbt_header[:][j102] = ''
                                self.Gps2.depth_ft[:,j102] = np.nan
                                self.Gps2.ft_indicator[:][j102] = ''
                                self.Gps2.depth_m[:][j102] = np.nan
                                self.Gps2.m_indicator[:][j102] = ''
                                self.Gps2.depth_fath[:,j102] = np.nan
                                self.Gps2.fath_indicator[:][j102] = ''
                                                               

                            self.Gps2.dbt_delta_time[i_ens,j102] = delta_time
                            self.Gps2.dbt_header[i_ens][j102] = ''.join([chr(x) for x in f.read(7)])
                            self.Gps2.depth_ft[i_ens,j102] = np.fromfile(f, np.float32, count=1)[0]
                            self.Gps2.ft_indicator[i_ens][j102] = chr(f.read(1)[0])
                            self.Gps2.depth_m[i_ens, j102] = np.fromfile(f, np.float32, count=1)[0]
                            self.Gps2.m_indicator[i_ens][j102] =chr(f.read(1)[0])
                            self.Gps2.depth_fath[i_ens,j102] = np.fromfile(f, np.float32, count=1)[0]
                            self.Gps2.fath_indicator[i_ens][j102] = chr(f.read(1)[0])
                                
                        elif specific_id == 107:
                            j103 += 1
                            
                            if j103 > 20:
                                self.Gps2.hdt_delta_time[:,j103] = delta_time
                                self.Gps2.hdt_header[:][j103] = ''
                                self.Gps2.heading_deg[:,j103] = np.nan
                                self.Gps2.h_true_indicator[:][j103] = ''
                                
                            self.Gps2.hdt_delta_time[i_ens,j103] = delta_time
                            self.Gps2.hdt_header[i_ens][j103] = ''.join([chr(x) for x in f.read(7)])
                            self.Gps2.heading_deg[i_ens,j103] = np.fromfile(f, np.double, count=1)[0]
                            self.Gps2.h_true_indicator[i_ens][j103] = chr(f.read(1)[0])
                            
                        elif specific_id == 204:
                            j100 += 1
                            
                            if j100 > 20:
                                self.Gps2.gga_delta_time[:,j100] = np.nan
                                self.Gps2.gga_header[:][j100] = ''
                                self.Gps2.utc[:,j100] = np.nan
                                self.Gps2.lat_deg[:,j100] = np.nan
                                self.Gps2.lat_ref[:][j100] = ''
                                self.Gps2.lon_deg[:,j100] = np.nan
                                self.Gps2.corr_qual[:,j100] = np.nan
                                self.Gps2.num_sats[:,j100] = np.nan
                                self.Gps2.hdop[:,j100] = np.nan
                                self.Gps2.alt[:,j100] = np.nan
                                self.Gp2.alt_unit[:,j100] = ''
                                self.Gps2.geoid[:,j100] = np.nan
                                self.Gps2.geoid_unit[:][j100] = ''
                                self.Gps2.d_gps_age[:,j100] = np.nan
                                self.Gps2.ref_stat_id[:,j100] = np.nan
                                
                            temp = ''.join([chr(x) for x in f.read(msg_size)])
                            self.Gps2.gga_sentence[i_ens][j100] = temp
                            temp_array = np.array(temp.split(','))
                            temp_array[temp_array == '999.9'] = ''
                            
                            try:
                                self.Gps2.gga_delta_time[i_ens,j100] = delta_time
                                self.Gps2.gga_header[i_ens][j100] = temp_array[0]
                                self.Gps2.utc[i_ens,j100] = float(temp_array[1])
                                lat_str = temp_array[2]
                                lat_deg = float(lat_str[0:2])
                                lat_deg = lat_deg+float(lat_str[2:]) / 60
                                self.Gps2.lat_deg[i_ens,j100] = lat_deg
                                self.Gps2.lat_ref[i_ens][j100] = temp_array[3]
                                lon_str = temp_array[4]
                                lon_num = float(lon_str)
                                lon_deg = np.floor(lon_num / 100)
                                self.Gps2.lon_deg[i_ens,j100] = lon_deg
                                self.Gps2.lon_ref[i_ens][j100] = temp_array[5]
                                self.Gps2.corr_qual[i_ens,j100] = float(temp_array[6])
                                self.Gps2.num_sats[i_ens,j100] = float(temp_array[7])
                                self.Gps2.hdop[i_ens,j100] = float(temp_array[8])
                                self.Gps2.alt[i_ens,j100] = float(temp_array[9])
                                self.Gps2.alt_unit[i_ens,j100] = temp_array[10]
                                self.Gps2.geoid[i_ens,j100] = temp_array[11]
                                self.Gps2.geoid_unit[i_ens,j100] = temp_array[12]
                                self.Gps2.d_gps_age[i_ens,j100] = float(temp_array[13])
                                idx_star =  temp_array[14].find('*')
                                self.Gps2.ref_stat_id[i_ens,j100] = float(temp_array[15][:idx_star])
                                
                            except:
                                pass
                            
                        elif specific_id == 205:
                            j101 += 1
                            
                            if j101>20:
                                self.Gps2.vtg_delta_time[:,j101] = np.nan
                                self.Gps2.vtg_header[:][j101] = ''
                                self.Gps2.course_true[:,j101] = np.nan
                                self.Gps2.true_indicator[:][j101] = ''
                                self.Gps2.course_mag[:,j101] = np.nan
                                self.Gps2.mag_indicator[:][j101] = ''
                                self.Gps2.speed_knots[:,j101] = np.nan
                                self.Gps2.knots_indicator[:][j101] = ''
                                self.Gps2.speed_k_mph[:,j101] = np.nan
                                self.Gps2.kmph_indicator[:][j101] = ''
                                self.Gps2.mode_indicator[:][j101] = ''
                                
                            temp = ''.join([chr(x) for x in f.read(msg_size)])
                            self.Gps2.vtg_sentence[i_ens][j100] = temp
                            temp_array = np.array(temp.split(','))
                            temp_array[temp_array == '999.9'] = ''
                            
                            try:
                                self.Gps2.vtg_delta_time[i_ens,j101] = delta_time
                                self.Gps2.vtg_header[i_ens][j101] = temp_array[0]
                                self.Gps2.course_true[i_ens,j101] = float(temp_array[1])
                                self.Gps2.true_indicator[i_ens][j101] = temp_array[2]
                                self.Gps2.course_mag[i_ens,j101] = float(temp_array[3])
                                self.Gps2.mag_indicator[i_ens][j101] = temp_array[4]
                                self.Gps2.speed_knots[i_ens,j101] = float(temp_array[5])
                                self.Gps2.knots_indicator[i_ens][j101] = temp_array[6]
                                self.Gps2.speed_k_mph[i_ens,j101] = float(temp_array[7])
                                self.Gps2.kmph_indicator[i_ens][j101] = temp_array[8]
                                idx_star = temp_array[9].find('*')
                                self.Gps2.mode_indicator[i_ens][j101] = temp_array[9][:idx_star]
                                
                            except:
                                pass
                            
                        elif specific_id == 206:
                            j102 += 1
                            
                            if j102 > 20:
                                self.Gps2.dbt_delta_time[:,j102] = np.nan
                                self.Gps2.dbt_header[:][j102] = ''
                                self.Gps2.depth_ft[:,j102] = np.nan
                                self.Gps2.ft_indicator[:][j102] = ''
                                self.Gps2.depth_m[:,j102] = np.nan
                                self.Gps2.m_indicator[:][j102] = ''
                                self.Gps2.depth_fath[:,j102] = np.nan
                                self.Gps2.fath_indicator[:][j102] = ''
                                
                            temp = ''.join([chr(x) for x in f.read(msg_size)])
                            temp_array = np.array(temp.split(','))
                            temp_array[temp_array == '999.9'] = ''    
                                
                            try:
                                self.Gps2.dbt_delta_time[i_ens,j102] = delta_time
                                self.Gps2.dbt_header[i_ens][j102] = temp_array[0]
                                self.Gps2.depth_ft[i_ens,j102] = float(temp_array[1])
                                self.Gps2.ft_indicator[i_ens][j102] = temp_array[2]
                                self.Gps2.depth_m[i_ens,j102] = float(temp_array[3])
                                self.Gps2.m_indicator[i_ens][j102] = temp_array[4]
                                self.Gps2.depth_fath[i_ens,j102] = float(temp_array[5])
                                idx_star = temp.find('*')
                                self.Gps2.fath_indicator[i_ens][j102] = temp_array[6][:idx_star]
                                
                            except:
                                pass
                            
                        elif specific_id == 207:
                            j103 += 1
                            
                            if j103 > 20:
                                self.Gps2.hdt_delta_time[:,j103] = np.nan
                                self.Gps2.hdt_header[:][j103] = ''
                                self.Gps2.heading_deg[:,j103] = np.nan
                                self.Gps2.h_true_indicator[:][j103] = ''
                                
                            temp = ''.join([chr(x) for x in f.read(msg_size)])
                            temp_array = np.array(temp.split(','))
                            temp_array[temp_array == '999.9'] = '' 
                            
                            try:
                                self.Gps2.hdt_delta_time[i_ens,j103] = delta_time
                                self.Gps2.hdt_header[i_ens][j103] = temp_array[0]
                                self.Gps2.heading_deg[i_ens,j103] = float(temp_array[1])
                                idx_star = temp.find('*')
                                self.Gps2.h_true_indicator[i_ens][j103] = temp_array[2][:idx_star]
                                
                            except:
                                pass
                                
                        self.end_reading(f, file_loc, i_data_types, i_ens, bytes_per_ens)
                           
                    elif leader_id == '0x2100':                
                        
                        #Update data types counter
                        i_data_types += 1 
                        #Reposition file pointer  
                        f.seek(int(self.Hdr.data_offsets[i_ens,i_data_types-1])+file_loc+4,0)  
                        
                        #Determine the number of characters to read    
                        if i_data_types < self.Hdr.n_data_types[i_ens]:
                            num_2_read = self.Hdr.data_offsets[i_ens,i_data_types] - self.Hdr.data_offsets[i_ens,i_data_types - 1] - 4
                        else:
                            num_2_read = bytes_per_ens - self.Hdr.data_offsets[i_ens, i_data_types-1] - 6 
                            
                        #Read DBT sentence
                        self.Nmea.dbt[i_ens] = ''.join([chr(x) for x in f.read(int(num_2_read))])

                        self.end_reading(f, file_loc, i_data_types, i_ens, bytes_per_ens)
                        
                    elif leader_id == '0x2101':
                        #Update data types counter
                        i_data_types += 1 
                        #Reposition file pointer  
                        f.seek(int(self.Hdr.data_offsets[i_ens,i_data_types-1])+file_loc+4,0)   
                        #Determine the number of characters to read  
                         
                        if i_data_types < self.Hdr.n_data_types[i_ens]:
                            num_2_read = self.Hdr.data_offsets[i_ens,i_data_types] - self.Hdr.data_offsets[i_ens,i_data_types - 1] - 4
                        else:
                            num_2_read = bytes_per_ens - self.Hdr.data_offsets[i_ens, i_data_types-1] - 6 
                            
                        #Read GGA sentence
                        self.Nmea.gga[i_ens] = ''.join([chr(x) for x in f.read(int(num_2_read))])  
                              
                        self.end_reading(f, file_loc, i_data_types, i_ens, bytes_per_ens)
                    
                    elif leader_id == '0x2102':
                        #Update data types counter
                        i_data_types += 1 
                        #Reposition file pointer  
                        f.seek(int(self.Hdr.data_offsets[i_ens,i_data_types-1])+file_loc+4,0)   
                        #Determine the number of characters to read  
                         
                        if i_data_types < self.Hdr.n_data_types[i_ens]:
                            num_2_read = self.Hdr.data_offsets[i_ens,i_data_types] - self.Hdr.data_offsets[i_ens,i_data_types - 1] - 4
                        else:
                            num_2_read = bytes_per_ens - self.Hdr.data_offsets[i_ens, i_data_types-1] - 6 
                          
                        #Read VTG sentence
                        self.Nmea.vtg[i_ens] = ''.join([chr(x) for x in f.read(int(num_2_read))])  
                        
                        self.end_reading(f, file_loc, i_data_types, i_ens, bytes_per_ens)
                        
                    elif leader_id == '0x2103':
                        #Update data types counter
                        i_data_types += 1 
                        #Reposition file pointer  
                        f.seek(int(self.Hdr.data_offsets[i_ens,i_data_types-1])+file_loc+4,0)   
                        #Determine the number of characters to read  
                         
                        if i_data_types < self.Hdr.n_data_types[i_ens]:
                            num_2_read = self.Hdr.data_offsets[i_ens,i_data_types] - self.Hdr.data_offsets[i_ens,i_data_types - 1] - 4
                        else:
                            num_2_read = bytes_per_ens - self.Hdr.data_offsets[i_ens, i_data_types-1] - 6 
                            
                        #Read GSA sentence
                        self.Nmea.gsa[i_ens] = ''.join([chr(x) for x in f.read(num_2_read)])  
                        
                        self.end_reading(f, file_loc, i_data_types, i_ens, bytes_per_ens)
                        
                    elif leader_id == '0x10':
                        #Update data types counter
                        i_data_types += 1 
                        
                        self.Surface.no_cells[i_ens] = np.fromfile(f, np.uint8, count=1)[0]
                        self.Surface.cell_size_cm[i_ens] = np.fromfile(f, np.uint16, count=1)[0]
                        self.Surface.dist_bin1_cm[i_ens] = np.fromfile(f, np.uint16, count=1)[0]
                        
                        self.end_reading(f, file_loc, i_data_types, i_ens, bytes_per_ens)
                        
                    #read velocity data
                    elif leader_id == '0x110':
                        #Update data types counter
                        i_data_types += 1 
                        
                        dummy = np.fromfile(f, np.int16, count= int((self.Surface.no_cells[i_ens]*4)))
                        dummy = np.reshape(dummy, [n_velocities,int(self.Surface.no_cells[i_ens])])
                        self.Surface.vel_mps[:n_velocities,:int(self.Surface.no_cells[i_ens]), i_ens] = dummy
                        
                        self.end_reading(f, file_loc, i_data_types, i_ens, bytes_per_ens)
                        
                        
                    #Read Surface Correlation Magnitude   
                    elif leader_id == '0x210':
                        #Update data types counter
                        i_data_types += 1 
                        
                        dummy = np.fromfile(f, np.uint8, count= int((self.Surface.no_cells[i_ens]*4)))
                        dummy = np.reshape(dummy, [n_velocities,int(self.Surface.no_cells[i_ens])])
                        self.Surface.corr[:n_velocities,:int(self.Surface.no_cells[i_ens]), i_ens] = dummy
                        
                        self.end_reading(f, file_loc, i_data_types, i_ens, bytes_per_ens)
                        
                    #Read Surface Echo Intensity
                    elif leader_id == '0x310':
                        #Update data types counter
                        i_data_types += 1
                        
                        dummy = np.fromfile(f, np.uint8, count= int((self.Surface.no_cells[i_ens]*4)))
                        dummy = np.reshape(dummy, [int(self.Surface.no_cells[i_ens]),n_velocities])
                        self.Surface.rssi[:n_velocities,:int(self.Surface.no_cells[i_ens]), i_ens] = dummy.T
                        
                        self.end_reading(f, file_loc, i_data_types, i_ens, bytes_per_ens)
                        
                    elif leader_id == '0x410':
                        #Update data types counter
                        i_data_types += 1
                        
                        dummy = np.fromfile(f, np.uint8, count= int((self.Surface.no_cells[i_ens]*4)))
                        dummy = np.reshape(dummy, [self.Surface.no_cells[i_ens],n_velocities])
                        self.Surface.pergd[:n_velocities,:self.Surface.no_cells[i_ens], i_ens] = dummy.T
                        
                        self.end_reading(f, file_loc, i_data_types, i_ens, bytes_per_ens)
                        
                    elif leader_id == '0x510':
                        #Update data types counter
                        i_data_types += 1
                        
                        self.end_reading(f, file_loc, i_data_types, i_ens, bytes_per_ens)
                        
                    elif leader_id == '0x4401':
                        #Update data types counter
                        i_data_types += 1
                        
                        self.AutoMode.beam_count[i_ens] = np.fromfile(f, np.uint8, count=1)[0]
                        
                        self.AutoMode.Beam1.mode[i_ens] = np.fromfile(f, np.uint8, count=1)[0]
                        self.AutoMode.Beam1.depth_cm[i_ens] = np.fromfile(f, np.uint16, count=1)[0]
                        self.AutoMode.Beam1.ping_count[i_ens] = np.fromfile(f, np.uint8, count=1)[0]
                        self.AutoMode.Beam1.ping_type[i_ens] = np.fromfile(f, np.uint8, count=1)[0]
                        self.AutoMode.Beam1.cell_count[i_ens] = np.fromfile(f, np.uint16, count=1)[0]
                        self.AutoMode.Beam1.cell_size_cm[i_ens] = np.fromfile(f, np.uint16, count=1)[0]
                        self.AutoMode.Beam1.cell_mid_cm[i_ens] = np.fromfile(f, np.uint16, count=1)[0]
                        self.AutoMode.Beam1.code_repeat[i_ens] = np.fromfile(f, np.uint8, count=1)[0]
                        self.AutoMode.Beam1.trans_length_cm[i_ens] = np.fromfile(f, np.uint16, count=1)[0]
                        self.AutoMode.Beam1.lag_length_cm[i_ens] = np.fromfile(f, np.uint16, count=1)[0]
                        self.AutoMode.Beam1.transmit_bw[i_ens] = np.fromfile(f, np.uint8, count=1)[0]
                        self.AutoMode.Beam1.receive_bw[i_ens] = np.fromfile(f, np.uint8, count=1)[0]
                        self.AutoMode.Beam1.ping_interval_ms[i_ens] = np.fromfile(f, np.uint16, count=1)[0]
                       
                        self.AutoMode.Beam2.mode[i_ens] = np.fromfile(f, np.uint8, count=1)[0]
                        self.AutoMode.Beam2.depth_cm[i_ens] = np.fromfile(f, np.uint16, count=1)[0]
                        self.AutoMode.Beam2.ping_count[i_ens] = np.fromfile(f, np.uint8, count=1)[0]
                        self.AutoMode.Beam2.ping_type[i_ens] = np.fromfile(f, np.uint8, count=1)[0]
                        self.AutoMode.Beam2.cell_count[i_ens] = np.fromfile(f, np.uint16, count=1)[0]
                        self.AutoMode.Beam2.cell_size_cm[i_ens] = np.fromfile(f, np.uint16, count=1)[0]
                        self.AutoMode.Beam2.cell_mid_cm[i_ens] = np.fromfile(f, np.uint16, count=1)[0]
                        self.AutoMode.Beam2.code_repeat[i_ens] = np.fromfile(f, np.uint8, count=1)[0]
                        self.AutoMode.Beam2.trans_length_cm[i_ens] = np.fromfile(f, np.uint16, count=1)[0]
                        self.AutoMode.Beam2.lag_length_cm[i_ens] = np.fromfile(f, np.uint16, count=1)[0]
                        self.AutoMode.Beam2.transmit_bw[i_ens] = np.fromfile(f, np.uint8, count=1)[0]
                        self.AutoMode.Beam2.receive_bw[i_ens] = np.fromfile(f, np.uint8, count=1)[0]
                        self.AutoMode.Beam2.ping_interval_ms[i_ens] = np.fromfile(f, np.uint16, count=1)[0]
                      
                        self.AutoMode.Beam3.mode[i_ens] = np.fromfile(f, np.uint8, count=1)[0]
                        self.AutoMode.Beam3.depth_cm[i_ens] = np.fromfile(f, np.uint16, count=1)[0]
                        self.AutoMode.Beam3.ping_count[i_ens] = np.fromfile(f, np.uint8, count=1)[0]
                        self.AutoMode.Beam3.ping_type[i_ens] = np.fromfile(f, np.uint8, count=1)[0]
                        self.AutoMode.Beam3.cell_count[i_ens] = np.fromfile(f, np.uint16, count=1)[0]
                        self.AutoMode.Beam3.cell_size_cm[i_ens] = np.fromfile(f, np.uint16, count=1)[0]
                        self.AutoMode.Beam3.cell_mid_cm[i_ens] = np.fromfile(f, np.uint16, count=1)[0]
                        self.AutoMode.Beam3.code_repeat[i_ens] = np.fromfile(f, np.uint8, count=1)[0]
                        self.AutoMode.Beam3.trans_length_cm[i_ens] = np.fromfile(f, np.uint16, count=1)[0]
                        self.AutoMode.Beam3.lag_length_cm[i_ens] = np.fromfile(f, np.uint16, count=1)[0]
                        self.AutoMode.Beam3.transmit_bw[i_ens] = np.fromfile(f, np.uint8, count=1)[0]
                        self.AutoMode.Beam3.receive_bw[i_ens] = np.fromfile(f, np.uint8, count=1)[0]
                        self.AutoMode.Beam3.ping_interval_ms[i_ens] = np.fromfile(f, np.uint16, count=1)[0]
                      
                        self.AutoMode.Beam4.mode[i_ens] = np.fromfile(f, np.uint8, count=1)[0]
                        self.AutoMode.Beam4.depth_cm[i_ens] = np.fromfile(f, np.uint16, count=1)[0]
                        self.AutoMode.Beam4.ping_count[i_ens] = np.fromfile(f, np.uint8, count=1)[0]
                        self.AutoMode.Beam4.ping_type[i_ens] = np.fromfile(f, np.uint8, count=1)[0]
                        self.AutoMode.Beam4.cell_count[i_ens] = np.fromfile(f, np.uint16, count=1)[0]
                        self.AutoMode.Beam4.cell_size_cm[i_ens] = np.fromfile(f, np.uint16, count=1)[0]
                        self.AutoMode.Beam4.cell_mid_cm[i_ens] = np.fromfile(f, np.uint16, count=1)[0]
                        self.AutoMode.Beam4.code_repeat[i_ens] = np.fromfile(f, np.uint8, count=1)[0]
                        self.AutoMode.Beam4.trans_length_cm[i_ens] = np.fromfile(f, np.uint16, count=1)[0]
                        self.AutoMode.Beam4.lag_length_cm[i_ens] = np.fromfile(f, np.uint16, count=1)[0]
                        self.AutoMode.Beam4.transmit_bw[i_ens] = np.fromfile(f, np.uint8, count=1)[0]
                        self.AutoMode.Beam4.receive_bw[i_ens] = np.fromfile(f, np.uint8, count=1)[0]
                        self.AutoMode.Beam4.ping_interval_ms[i_ens] = np.fromfile(f, np.uint16, count=1)[0]
                        
                        self.AutoMode.Reserved[i_ens] = np.fromfile(f, np.uint8, count=1)[0]
                        
                        self.end_reading(f, file_loc, i_data_types, i_ens, bytes_per_ens)
                    
                    elif leader_id == '0x4100':
                        
                        i_data_types += 1
                        
                        self.Sensor.vert_beam_eval_amp[i_ens] = np.fromfile(f, np.uint8, count=1)[0]
                        self.Sensor.vert_beam_RSSI_amp[i_ens] = np.fromfile(f, np.uint8, count=1)[0]
                        self.Sensor.vert_beam_range_m[i_ens] = np.fromfile(f, np.uint32, count=1)[0] / 1000
                        temp = "{0:08b}".format(np.fromfile(f, np.uint8, count=1)[0])
                        self.Sensor.vert_beam_status[i_ens] = int(temp[6:],2)
                        if (temp[5] == '0'):
                            self.Sensor.vert_beam_gain[i_ens] = 'L'
                        else:
                            self.Sensor.vert_beam_gain[i_ens] = 'H'
                            
                        
                        self.end_reading(f, file_loc, i_data_types, i_ens, bytes_per_ens)
                        
                    elif leader_id == '0x3200':
                        
                        i_data_types += 1
                        
                        self.Inst.t_matrix[0,:] = np.fromfile(f, np.int16, count=4) * .0001
                        self.Inst.t_matrix[1,:] = np.fromfile(f, np.int16, count=4) * .0001
                        self.Inst.t_matrix[2,:] = np.fromfile(f, np.int16, count=4) * .0001
                        self.Inst.t_matrix[3,:] = np.fromfile(f, np.int16, count=4) * .0001
                    
                        self.end_reading(f, file_loc, i_data_types, i_ens, bytes_per_ens)
                        
                    else:
                        
                        self.Hdr.invalid[i_ens] = leader_id
                        i_data_types += 1
                        
                        if (i_data_types + 1) <= self.Hdr.n_data_types[i_ens]:
                            f.seek(int(self.Hdr.data_offsets[i_ens,i_data_types])+file_loc,0)
                        else:
                            if f.tell() < end_file:
                                if i_data_types + 1 > self.Hdr.n_data_types[i_ens] + 1:
                                    current_loc = f.tell()
                                    srch_string = struct.unpack('B'*(end_file-current_loc), f.read(end_file-current_loc))
                                    hex_string = ''.join([hex(x) for x in srch_string])
                                    
                                    next_ens = hex_string.find('0x7f7f')
                                    if next_ens > 0:
                                        next_ens = (next_ens[0]-1)/2
                                        f.seek(current_loc+next_ens,0)
                                        i_data_types = 0
                                    else:
                                        end_file_check = end_file + 1
                                        
                                else:
                                    f.seek(file_loc+bytes_per_ens-2,0)
                                    
                                    
                    #update end_file check to stop iterations
                    if i_ens <= len(self.Hdr.n_data_types):
                        if i_data_types >= self.Hdr.n_data_types[i_ens] and f.tell() <= end_file:
                            
                            try:
                                self.Inst.res_RDI = np.fromfile(f, np.uint16, count=1)[0]
                                checksum = np.fromfile(f, np.uint16, count=1)[0]
                            except:
                                pass
                    else:
                        end_file_check = end_file
                        
                    if end_file_check < end_file:
                        end_file_check = f.tell()   
                    
                # Screen for bad data, and do the unit conversions
                self.Wt.vel_mps[self.Wt.vel_mps == -32768] = np.nan
                self.Wt.vel_mps = self.Wt.vel_mps / 1000
                self.Wt.corr[self.Wt.corr == -32768] = np.nan
                self.Wt.rssi[self.Wt.rssi == -32768] = np.nan
                self.Wt.pergd[self.Wt.pergd == -32768] = np.nan
                
                #Remove bad data, convert units
                self.Bt.depth_m[self.Bt.depth_m == -32768] = np.nan
                self.Bt.depth_m = self.Bt.depth_m / 100
                self.Bt.vel_mps[self.Bt.vel_mps == -32768] = np.nan
                self.Bt.vel_mps = self.Bt.vel_mps / 1000
                self.Bt.corr[self.Bt.corr == -32768] = np.nan
                self.Bt.eval_amp[self.Bt.eval_amp == -32768] = np.nan
                self.Bt.pergd[self.Bt.pergd == -32768] = np.nan
                
                # Correct Bt.depth_m for RiverRay data
                if np.prod(RR_BtDepth_correction == np.nan) == 0:
                    self.RR_BtDepth_correction = RR_BtDepth_correction
                    self.RR_BtDepth_correction[self.RR_BtDepth_correction == (-32768*2e16)/100] = np.nan
                    self.Bt.depth_m += RR_BtDepth_correction
                    
                # Remove bad data from Surface structure (RR), convert where needed
                self.Surface.vel_mps[self.Surface.vel_mps == -32768] = np.nan
                self.Surface.vel_mps = self.Surface.vel_mps / 1000
                self.Surface.corr[self.Surface.corr == -32768] = np.nan
                self.Surface.rssi[self.Surface.rssi == -32768] = np.nan
                self.Surface.pergd[self.Surface.pergd == -32768] = np.nan
                
                if self.WR2:
                    
                    #if vtg data are available compute north and east components
                    if self.Gps2.vtg_header[0,0] == '$':
                        
                        # Find minimum of absolute value of delta time from raw data
                        vtg_delta_time = np.abs(self.Gps2.vtg_delta_time)
                        vtg_min = np.nanmin(vtg_delta_time)
                        
                        # compute the velocity components in m/s
                        for i in range(len(vtg_delta_time)):
                            idx = np.where(vtg_delta_time == vtg_min)[0][0]
                            self.Gps2.vtg_velE_mps[i], self.Gps2.vtg_velN_mps[i] = \
                            self.pol2cart((90 - self.Gps2.course_true[i,idx])*np.pi/180,
                                           self.Gps2.speed_k_mph[i,idx] * 0.2777778)
                            
                    if self.Gps2.gga_header[0,0] == '$':
                        
                        idx = None
                        e_radius = 6378137
                        coeff = e_radius * np.pi / 180
                        ellip = 1 / 298.257223563
                        #Find minimum of absolute value of delta time from raw data
                        gga_delta_time = np.abs(self.Gps2.gga_delta_time)
                        gga_min = np.nanmin(gga_delta_time, axis=1)
                        
                        for i in range(len(self.gga_delta_time)):
                            idx = np.where(gga_delta_time[i:] == gga_min)
                            if idx > 0:
                                L = (self.Gps2.lat_deg[i,idx[i]]+self.Gps2.lat_deg[i-1,idx[i-1]]) / 2
                                sL = sindg(L)
                                RE = coeff * (1+ellip * sL * sL)
                                RN = coeff * (1 - 2*ellip + 3*ellip * sL * sL)
                                dx = RE * (self.Gps2.lon_deg[i, idx[i]] - self.Gps2.lon_deg(i-1,idx(i-1))) * cosdg(L)
                                dy = RN * (self.Gps2.lat_deg[i, idx[i]] - self.Gps2.lat_deg[i-1,idx[i-1]])
                                dt = self.Gps2.utc[i,idx[i]] - self.Gps2.utc[i-1,idx[i-1]]
                                self.Gps2.gga_vel_e_mps[i] = dx / dt
                                self.Gps2.gga_velN_mps[i] = dy / dt
                            else:
                                self.Gps2.gga_velE_mps = np.nan
                                self.Gps2.gga_velN_mps = np.nan
                         
                            
                    
                   
                      
                    
                        
    def number_of_ensembles(self, f, f_info):
        i=0
        leader_id='0000'
            
        #find the first ensemble
        while leader_id != '0x7f7f' and i<f_info:
            f.seek(i,0)
            i=i+1
            leader_id = hex(np.fromfile(f, np.uint16, count=1)[0])
            
            
        #call findeensno to get the first ensemble number
        first_num = self.find_ens_no(f)
        
        i = 0
        leader_id='0000'
        last_num = -1
        
        while last_num < 0:
            while leader_id != '0x7f7f' and i < f_info:
                i=i+1
                f.seek(-i,2)
                
                try:
                    leader_id= hex(np.fromfile(f, np.uint16, count=1)[0])
                except:
                    continue
                    
               
                
                

            last_num = self.find_ens_no(f)
#             print 'last ens no', last_num
            if last_num is None or np.isnan(last_num):
                last_num = -1
            
            leader_id = '0000'
        
        return last_num-first_num+1
    
    def find_ens_no(self, f):
        '''
        This method finds the appropriate ensemble and performs a checksum for QA
        Terminology: bof (Seek from beginning of file x bytes)
                     cof (Seek from current location x bytes)
                     eof (Seek from end of file x bytes)
        '''
        try:
            fileloc = f.tell()-2
            if self.check_sum(f, fileloc):
           
                    f.seek(fileloc+5,0)
                
                    n_data_types = np.fromfile(f, np.uint8, count=1)[0]
                    data_offsets = []
                    for x in range(n_data_types):
                        data_offsets.append(np.fromfile(f, np.uint16, count=1)[0])
                
                    i=0
                    leader_id='0000'
                    
                    while leader_id != '0x80' and i < n_data_types:
                        
                        f.seek(data_offsets[i]+fileloc,0)
                        leader_id = hex(np.fromfile(f, np.uint16, count=1)[0])
                        i = i+1
                    
                #if the correct hexadecimal is found return ensemble number, otherwise -1
                    if leader_id == '0x80':
#                         print 'ftell of found x80', f.tell()
                        return np.fromfile(f, np.uint16, count=1)[0]
    
               
               
            else:
                return -1
        except:
            return np.nan
        
        
    def check_sum(self, f, fileloc, bytes_per_ens=None):
        #get file locatuion and number of bytes in an ensemble
          
        try:
             
            if bytes_per_ens is None:
                bytes_per_ens = np.fromfile(f, np.uint16, count=1)[0] 
            #go to file location from the beginning of file
            f.seek(fileloc,0)
              
            # read in the values for all of the bytes an get a check sum
            Testb = []
            x = f.read(bytes_per_ens)
            for y in x:
                Testb.append(y)
                  
            check_sum = sum(Testb)
            check_h = hex(check_sum)[2:]
              
            #we want a hex that is greater than 4 (including L indicator at the end)
            if len(check_h)>4:
                  
                #seek to locati_hon of check sum and compared to computed
                if check_h[-1] == 'L':
                    check_h = check_h[:-1]
                      
                f.seek(fileloc+bytes_per_ens,0)
                check_sum = np.fromfile(f, np.uint16, count=1)[0]  
                if int('0x'+check_h[1:], 16) == check_sum:
                    return True
                else:
                    return False
            elif len(check_h)>3:
                #seek to locati_hon of check sum and compared to computed
                if check_h[-1] == 'L':
                    check_h = check_h[:-1]
                      
                f.seek(fileloc+bytes_per_ens,0)
                check_sum = np.fromfile(f, np.uint16, count=1)[0]  
                if int('0x'+check_h, 16) == check_sum:
                    return True
                else:
                    return False
            else:
                return False
        except:  
            return False
        
    def bad_check_sum(self, f, file_loc):
        
#         print 'Bad Checksum New Code'
        search_id='    '
        search_loc = file_loc+2
        while search_id != '0x7f7f':
            f.seek(search_loc,0)
            search_loc += 1
              
            try:
                search_id= hex(np.fromfile(f, np.uint16, count=1)[0])
            except:
                continue
        f.seek(search_loc,0)
#         print 'Done with bad check sum'
        
    def end_reading(self, f, file_loc, i_data_types, i_ens, bytes_per_ens):
        if i_data_types + 1 <= self.Hdr.n_data_types[i_ens]:
                            f.seek(int(self.Hdr.data_offsets[i_ens,i_data_types])+file_loc,0)
        else:
            f.seek(file_loc+bytes_per_ens-2,0) 
            
    def pol2cart(self, rho, phi):
        x = rho * np.cos(phi)
        y = rho * np.sin(phi)
        return(x, y)
       
# -------------------------------------------------thus far testing individual pd0 file extraction until comfortable  
if __name__ == '__main__':
#     
    files = [r'C:\Users\gpetrochenkov\Desktop\drive-download-20170522T150040Z-0014\RG_1308000_359\13038000_359_000.PD0',
             r'C:\Users\gpetrochenkov\Desktop\drive-download-20170522T150040Z-0014\RG_1308000_359\13038000_359_001.PD0',
             r'C:\Users\gpetrochenkov\Desktop\drive-download-20170522T150040Z-0014\RG_1308000_359\13038000_359_002.PD0',
             r'C:\Users\gpetrochenkov\Desktop\drive-download-20170522T150040Z-0014\RG_1308000_359\13038000_359_003.PD0',
             r'C:\Users\gpetrochenkov\Desktop\drive-download-20170522T150040Z-0014\RG_1308000_359\13038000_359_004.PD0',
             r'C:\Users\gpetrochenkov\Desktop\drive-download-20170522T150040Z-0014\RG_1308000_359\13038000_359_005.PD0',
             r'C:\Users\gpetrochenkov\Desktop\drive-download-20170522T150040Z-0014\RG_Cal_NoEval\01327750_0_000_14-04-14.PD0',
             r'C:\Users\gpetrochenkov\Desktop\drive-download-20170522T150040Z-0014\RG_Cal_NoEval\01327750_0_000_14-04-14_LBT.PD0',
             r'C:\Users\gpetrochenkov\Desktop\drive-download-20170522T150040Z-0014\RG_Cal_NoEval\01327750_0_001_14-04-14.PD0',
             r'C:\Users\gpetrochenkov\Desktop\drive-download-20170522T150040Z-0014\RG_Cal_NoEval\01327750_0_002_14-04-14.PD0',
             r'C:\Users\gpetrochenkov\Desktop\drive-download-20170522T150040Z-0014\RG_Cal_NoEval\01327750_0_003_14-04-14.PD0',
             r'C:\Users\gpetrochenkov\Desktop\drive-download-20170522T150040Z-0014\RP_02353000_554_New_SysTest\02353000_554_000.PD0',
             r'C:\Users\gpetrochenkov\Desktop\drive-download-20170522T150040Z-0014\RP_02353000_554_New_SysTest\02353000_554_000_LBT.PD0',
             r'C:\Users\gpetrochenkov\Desktop\drive-download-20170522T150040Z-0014\RP_02353000_554_New_SysTest\02353000_554_001.PD0',
             r'C:\Users\gpetrochenkov\Desktop\drive-download-20170522T150040Z-0014\RR_Multi_Cal\05LC004_20140812.AQ1_0_001.PD0',
             r'C:\Users\gpetrochenkov\Desktop\drive-download-20170522T150040Z-0014\SP_10126sp778_3\10126sp778_0_000.PD0',
             r'C:\Users\gpetrochenkov\Desktop\drive-download-20170522T150040Z-0014\SP_10126sp778_3\10126sp778_0_000_LBT.PD0',
             r'C:\Users\gpetrochenkov\Desktop\drive-download-20170522T150040Z-0014\SP_10126sp778_3\10126sp778_0_001.PD0'
            ]
      
    c = [Pd0TRDI(x) for x in files]
    

        
            
                
                
    
            
        
        
            