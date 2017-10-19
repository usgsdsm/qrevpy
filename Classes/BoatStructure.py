'''
Created on Sep 5, 2017

@author: gpetrochenkov
'''
import numpy as np
from Classes.BoatData import BoatData

class BoatStructure(object):
    '''This class organizes the various sources for boat velocity into 
    a single structured class and establishes a selected property that
    contains the select source for velocity and discharge computations'''
    
    def __init__(self):
        
        self.selected = None #Name of BoatData object to be used for discharge computations
        self.bt_vel = None #BoatData object for bottom track velocity
        self.gga_vel = None #BoatData object for gga velocity
        self.vtg_vel = None #BoatData object for vtg velocity
        
        #composite track information is not currently provided by the manufacturers.
        #Future versions may try to determine this setting from SonTek data
        self.composite = 'Off' # Setting for compositir tracks "On", "Off"
        
    def add_boat_object(self, source, vel_in, freq_in, coord_sys_in, nav_ref_in, kargs = None):
        '''Adds a BoatData object to the appropriate property
        
        source: name of manufacturer
        
        vel_in: boat velocity array
        
        freq_in: acoustic frequency, if available
       
        coord_sys_in: coordinate system of data in velIn
       
        nav_ref_in: source of boat
        
        kargs: used for bottom track velocities
        kargs[0]: 3 beam solution setting
        kargs[1]: bottom mode
    
        '''
        
        if nav_ref_in == 'BT':
            self.bt_vel = BoatData()
            self.bt_vel.populate_data(source, vel_in, freq_in, coord_sys_in, nav_ref_in, kargs)
        if nav_ref_in == 'GGA':
            self.gga_vel = BoatData()
            self.gga_vel.populate_data(source, vel_in, freq_in, coord_sys_in, nav_ref_in)
        if nav_ref_in == 'VTG':
            self.vtg_vel = BoatData()
            self.vtg_vel.populate_data(source, vel_in, freq_in, coord_sys_in, nav_ref_in)
            
            
    def set_nav_reference(self, reference):
        '''This function will set the navigation reference property to the specified object reference'''
        
        if reference == 'BT':
            self.selected = 'bt_vel'
        elif reference == 'GGA':
            self.selected = 'gga_vel'
        elif reference == 'VTG':
            self.selected = 'vtg_vel'
            
    def change_nav_reference(self, reference, transect):
        '''This function changes the navigation reference to the specified object reference and recomputes
        the compsite tracks, if necessary'''
        
        
        if reference == 'BT':
            self.selected = 'bt_vel'
        elif reference == 'GGA':
            self.selected = 'gga_vel'
        elif reference == 'VTG':
            self.selected = 'vtg_vel'
        elif reference == 'bt_vel':
            self.selected = 'bt_vel'
        elif reference == 'gga_ve':
            self.selected = 'gga_vel'
        elif reference == 'vtg_vel':
            self.selected = 'vtg_vel'
            
        self.composite_tracks(transect)
        
        
    def change_coord_sys(self, new_coord_sys, sensors, adcp):
        '''This function will change the coordinate system of the boat velocity reference
        
        Input:
        new_coord_sys: specified new coordinate system
        sensors: object of Sensors
        adcp: object of InstrumentData
        '''
        
        #Change coordinate system for all available boat velocity sources
        if self.bt_vel is not None:
            self.bt_vel.change_coord_sys(new_coord_sys, sensors, adcp)
        if self.gga_vel is not None:
            self.gga_vel.change_coord_sys(new_coord_sys, sensors, adcp)
        if self.vtg_vel is not None:
            self.vtg_vel.change_coord_sys(new_coord_sys, sensors, adcp)   
            
            
        
    def composite_tracks(self, transect, kargs=None):
        '''If new composite setting is provided it is used, if not the setting saved in the object is used
        
        Input:
        trasect: object of TransectData
        kargs[0]: new setting for composite tracks
        '''
        
        if kargs is None:
            setting = self.composite
        else:
            #New Setting
            setting = kargs[0]
            self.composite = setting
            
        #Composite depths turned on
        if setting == 'On':
            
            #Prepare Bt data
            if self.bt_vel is not None:
                u_BT = self.bt_vel.__u_processed_mps
                v_BT = self.bt_vel.__v_processed_mps
                #Set to invalid all interpolated velocities
                valid_BT = self.bt_vel.__valid_data[0,:]
                u_BT[valid_BT == False] = np.nan
                v_BT[valid_BT == False] = np.nan
                
            if self.gga_vel is not None:
                #Get gga velocities
                u_GGA = self.gga_vel.__u_processed_mps
                v_GGA = self.gga_vel.__v_processed_mps
                #Set to invalid all interpolated velocities
                valid_GGA = self.gga_vel.__valid_data[0,:]
                u_GGA[valid_GGA == False] = np.nan
                v_GGA[valid_GGA == False] = np.nan
            else:
                u_GGA = np.tile([np.nan], u_BT.shape)
                v_GGA = np.tile([np.nan], v_BT.shape)
                
            if self.vtg_vel is not None:
                #Get vtg velocities
                u_VTG = self.vtg_vel.__u_processed_mps
                v_VTG = self.vtg_vel.__v_processed_mps
                #Set to invalid all interpolated velocities
                valid_VTG = self.vtg_vel.__valid_data[0,:]
                u_VTG[valid_VTG == False] = np.nan
                v_VTG[valid_VTG == False] = np.nan
            else:
                u_VTG = np.tile([np.nan], u_BT.shape)
                v_VTG = np.tile([np.nan], v_BT.shape)
                
            if self.bt_vel is not None:
                
                #Initialize composite source
                comp_source = np.tile(np.nan, u_BT.shape)
                
                #Process u velocity component
                u_comp = u_BT
                comp_source[np.isnan(u_comp) == False] = 1
                
                #If BT data are not valid try VTG and set composite source
                u_comp[np.isnan(u_comp)] = u_VTG[np.isnan(comp_source)]
                comp_source[np.isnan(u_comp) == False & np.isnan(comp_source)] = 3
                
                #If there are still invalid boat velocities, try GGA and set composite source
                u_comp[np.isnan(u_comp)] = u_GGA[np.isnan(u_comp)]
                comp_source[np.isnan(u_comp) == False & np.isnan(comp_source)] = 2
                
                #If there are still invalid boat velocities, use interpolated
                #values if present and set composite source
                u_comp[np.isnan(u_comp)] = self.bt_vel.__u_processed_mps[np.isnan(u_comp)]
                comp_source[np.isnan(u_comp) == False & np.isnan(comp_source)] = 0
                
                #Set composite source to invalid for all remaining invalid boat velocity data
                comp_source[np.isnan(comp_source)] = -1
                
                #Process v velocity component.  Assume that the composite source is the same
                #as the u component
                v_comp = v_BT
                v_comp[np.isnan(v_comp)] = v_VTG[np.isnan(v_comp)]
                v_comp[np.isnan(v_comp)] = v_GGA[np.isnan(v_comp)]
                v_comp[np.isnan(v_comp)] = self.bt_vel.__v_processed_mps[np.isnan(v_comp)]
                
                #Apply the composite settings to the bottom track Boatdata objects
                self.bt_vel = self.apply_composite(u_comp,v_comp,comp_source)
                self.bt_vel = self.interpolate_composite()
                
            if self.gga_vel is not None:
                
                #Initialize the composite source
                comp_source = np.tile([np.nan], u_BT.shape)
                
                #Process the u velocity component
                u_comp = u_GGA
                comp_source[np.isnan(u_comp) == False] = 2
                
                #If GGA data are not valid try VTG and set composite source
                u_comp[np.isnan(u_comp)] = u_VTG[np.isnan(u_comp)]
                comp_source[np.isnan(u_comp) == False & np.isnan(comp_source)] = 3
                
                #If there are still invalid boar velocities, try BT and set composite source
                u_comp[np.isnan(u_comp)] = u_BT[np.isnan(u_comp)]
                comp_source[np.isnan(u_comp) == False & np.isnan(comp_source)] = 1
                
                #If there are still invalid boat velocities, use interpolated values,
                #if present and set composite source
                u_comp[np.isnan(u_comp)] = self.gga_vel.__u_processed_mps[np.isnan(u_comp)]
                comp_source[np.isnan(u_comp) == False & np.isnan(comp_source)] = 0
                
                #Set composite source to invalid for all remaining invalid boat velocity data
                comp_source[np.isnan(comp_source)] = -1
                
                #Process v velocity component.  Assume that the composite source is the 
                #same as the u component
                v_comp = v_GGA
                v_comp[np.isnan(v_comp)] = v_VTG[np.isnan(v_comp)]
                v_comp[np.isnan(v_comp)] = v_BT[np.isnan(v_comp)]
                v_comp[np.isnan(v_comp)] = self.gga_vel.__v_processed_mps[np.isnan(v_comp)]
                
                #Apply the composite settings to the gga BoatData object
                
    
        
        
         
        
        