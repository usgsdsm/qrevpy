"""
Created on Sep 5, 2017

@author: gpetrochenkov

Modified DSM 1/29/2018
    - docstrings to numpy style
    - cleaned up PEP8
    - added code for vtg primary
    - made code consistent with latest version of QRev
    - removed need for karg
"""
import numpy as np
from Classes.BoatData import BoatData


class BoatStructure(object):
    """This class organizes the various sources for boat velocity into
    a single structured class and establishes a selected property that
    contains the select source for velocity and discharge computations.

    Attributes
    ----------
    selected: str
        Name of BoatData object to be used for discharge computations.
    bt_vel: object
        BoatData object for bottom track velocity
    gga_vel: object
        BoatData object for gga velocity
    vtg_vel: object
        BoatData object for vtg velocity
    composite: bool
        Setting to use (True) or not (False) composite tracks.
    """
    
    def __init__(self):
        
        self.selected = None  # Name of BoatData object to be used for discharge computations
        self.bt_vel = None  # BoatData object for bottom track velocity
        self.gga_vel = None  # BoatData object for gga velocity
        self.vtg_vel = None  # BoatData object for vtg velocity
        
        # Composite track information is not currently provided by the manufacturers.
        # Future versions may try to determine this setting from SonTek data
        self.composite = False  # Setting for compositir tracks
        
    def add_boat_object(self, source, vel_in, freq_in=None, coord_sys_in=None, nav_ref_in=None, min_beams=3, bottom_mode='Variable'):
        """Adds a BoatData object to the appropriate property

        Parameters
        ----------
        source: str
            Name of manufacturer.
        vel_in: np.array
            Boat velocity array.
        freq_in: np.array or float
            Acoustic frequency
        coord_sys_in: str
            Coordinate system of velocity data.
        nav_ref_in: str
            Source of boat velocity data
        min_beams: float
            Setting to allow 3 beam solutions or require 4 beam solutions or set to Auto (-1)
        bottom_mode: str
            Bottom mode used
        """
        
        if nav_ref_in == 'BT':
            self.bt_vel = BoatData()
            self.bt_vel.populate_data(source, vel_in, freq_in, coord_sys_in, nav_ref_in, min_beams, bottom_mode)
        if nav_ref_in == 'GGA':
            self.gga_vel = BoatData()
            self.gga_vel.populate_data(source, vel_in, freq_in, coord_sys_in, nav_ref_in)
        if nav_ref_in == 'VTG':
            self.vtg_vel = BoatData()
            self.vtg_vel.populate_data(source, vel_in, freq_in, coord_sys_in, nav_ref_in)
            
    def set_nav_reference(self, reference):
        """This function will set the navigation reference property to the specified object reference.

        Parameters
        ----------
        reference: str
            Navigation reference, BT, GGA, or VTG
        """
        
        if reference == 'BT':
            self.selected = 'bt_vel'
        elif reference == 'GGA':
            self.selected = 'gga_vel'
        elif reference == 'VTG':
            self.selected = 'vtg_vel'
            
    def change_nav_reference(self, reference, transect):
        """This function changes the navigation reference to the specified object reference and recomputes
        the compsite tracks, if necessary.

        Parameters
        ----------
        reference: str
            New navigation reference, BT, GGA, or VTG.
        transect: object
            Object of TransectData.
        """

        if reference == 'BT':
            self.selected = 'bt_vel'
        elif reference == 'GGA':
            self.selected = 'gga_vel'
        elif reference == 'VTG':
            self.selected = 'vtg_vel'
        elif reference == 'bt_vel':
            self.selected = 'bt_vel'
        elif reference == 'gga_vel':
            self.selected = 'gga_vel'
        elif reference == 'vtg_vel':
            self.selected = 'vtg_vel'

        self.composite_tracks(transect)

    def change_coord_sys(self, new_coord_sys, sensors, adcp):
        """This function will change the coordinate system of the boat velocity reference.
        
        Parameters
        ----------
        new_coord_sys: str
            Specified new coordinate system.
        sensors: object
            Object of Sensors.
        adcp: object
            Object of InstrumentData.
        """
        
        # Change coordinate system for all available boat velocity sources
        if self.bt_vel is not None:
            self.bt_vel.change_coord_sys(new_coord_sys, sensors, adcp)
        if self.gga_vel is not None:
            self.gga_vel.change_coord_sys(new_coord_sys, sensors, adcp)
        if self.vtg_vel is not None:
            self.vtg_vel.change_coord_sys(new_coord_sys, sensors, adcp)   

    def composite_tracks(self, transect, setting=None):
        """If new composite setting is provided it is used, if not the setting saved in the object is used
        
        Parameters
        ----------
        transect: object
            Object of TransectData.
        setting: bool
            New setting for composite tracks
        """
        
        if setting is None:
            setting = self.composite
        else:
<<<<<<< HEAD
            # New Setting
=======
            #New Setting
            if type(kargs) is str:
                setting = kargs
            else:
                setting = kargs[0]
>>>>>>> 6ca6c50c231afa610ed3a693864074d7104a5f20
            self.composite = setting
            
        # Composite depths turned on
        if setting:
            
            # Prepare bt data
            if self.bt_vel is not None:
                u_bt = self.bt_vel.u_processed_mps
                v_bt = self.bt_vel.v_processed_mps
                # Set to invalid all interpolated velocities
                valid_bt = self.bt_vel.valid_data[0, :]
                u_bt[valid_bt == False] = np.nan
                v_bt[valid_bt == False] = np.nan

            # Prepare gga data
            if self.gga_vel is not None:
                # Get gga velocities
                u_gga = self.gga_vel.u_processed_mps
                v_gga = self.gga_vel.v_processed_mps
                # Set to invalid all interpolated velocities
                valid_gga = self.gga_vel.valid_data[0, :]
                u_gga[valid_gga == False] = np.nan
                v_gga[valid_gga == False] = np.nan
            elif self.bt_vel is not None:
                u_gga = np.tile([np.nan], u_bt.shape)
                v_gga = np.tile([np.nan], v_bt.shape)

            # Prepare vtg data
            if self.vtg_vel is not None:
                # Get vtg velocities
                u_vtg = self.vtg_vel.u_processed_mps
                v_vtg = self.vtg_vel.v_processed_mps
                # Set to invalid all interpolated velocities
                valid_vtg = self.vtg_vel.valid_data[0, :]
                u_vtg[valid_vtg == False] = np.nan
                v_vtg[valid_vtg == False] = np.nan
            elif self.bt_vel is not None:
                u_vtg = np.tile([np.nan], u_bt.shape)
                v_vtg = np.tile([np.nan], v_bt.shape)

            # Process bt as primary
            if self.bt_vel is not None:
                
                # Initialize composite source
                comp_source = np.tile(np.nan, u_bt.shape)
                
                # Process u velocity component
                u_comp = u_bt
                comp_source[np.isnan(u_comp) == False] = 1
                
                # If BT data are not valid try VTG and set composite source
                u_comp[np.isnan(u_comp)] = u_vtg[np.isnan(comp_source)]
                comp_source[np.isnan(u_comp) == False & np.isnan(comp_source)] = 3
                
                # If there are still invalid boat velocities, try GGA and set composite source
                u_comp[np.isnan(u_comp)] = u_gga[np.isnan(u_comp)]
                comp_source[np.isnan(u_comp) == False & np.isnan(comp_source)] = 2
                
                # If there are still invalid boat velocities, use interpolated
                # values if present and set composite source
                # DSM changed 1/29/2018 u_comp[np.isnan(u_comp)] = self.bt_vel.u_processed_mps[np.isnan(u_comp)]
                comp_source[np.isnan(u_comp) == False & np.isnan(comp_source)] = 0
                
                # Set composite source to invalid for all remaining invalid boat velocity data
                comp_source[np.isnan(comp_source)] = -1
                
                # Process v velocity component.  Assume that the composite source is the same
                # as the u component
                v_comp = v_bt
                v_comp[np.isnan(v_comp)] = v_vtg[np.isnan(v_comp)]
                v_comp[np.isnan(v_comp)] = v_gga[np.isnan(v_comp)]
                v_comp[np.isnan(v_comp)] = self.bt_vel.v_processed_mps[np.isnan(v_comp)]
                
                # Apply the composite settings to the bottom track Boatdata objects
                # DSM changed 1/29/2018 self.bt_vel = self.apply_composite(u_comp,v_comp,comp_source)
                # self.bt_vel = self.interpolate_composite()
                self.bt_vel.apply_composite(u_comp, v_comp, comp_source)
                self.bt_vel.interpolate_composite(transect)

            # Process gga as primary
            if self.gga_vel is not None:
                
                # Initialize the composite source
                comp_source = np.tile([np.nan], u_bt.shape)
                
                # Process the u velocity component
                u_comp = u_gga
                comp_source[np.isnan(u_comp) == False] = 2
                
                # If GGA data are not valid try VTG and set composite source
                u_comp[np.isnan(u_comp)] = u_vtg[np.isnan(u_comp)]
                comp_source[np.isnan(u_comp) == False & np.isnan(comp_source)] = 3
                
                # If there are still invalid boar velocities, try BT and set composite source
                u_comp[np.isnan(u_comp)] = u_bt[np.isnan(u_comp)]
                comp_source[np.isnan(u_comp) == False & np.isnan(comp_source)] = 1
                
                # If there are still invalid boat velocities, use interpolated values,
                # if present and set composite source
                # DSM changed 1/29/2018 u_comp[np.isnan(u_comp)] = self.gga_vel.u_processed_mps[np.isnan(u_comp)]
                comp_source[np.isnan(u_comp) == False & np.isnan(comp_source)] = 0
                
                # Set composite source to invalid for all remaining invalid boat velocity data
                comp_source[np.isnan(comp_source)] = -1
                
                # Process v velocity component.  Assume that the composite source is the
                # same as the u component
                v_comp = v_gga
                v_comp[np.isnan(v_comp)] = v_vtg[np.isnan(v_comp)]
                v_comp[np.isnan(v_comp)] = v_bt[np.isnan(v_comp)]
                v_comp[np.isnan(v_comp)] = self.gga_vel.v_processed_mps[np.isnan(v_comp)]
                
                # Apply the composite settings to the gga BoatData object
                # Apply the composite settings to the bottom track Boatdata objects
                # DSM changed 1/29/2018 self.gga_vel = self.apply_composite(u_comp,v_comp,comp_source)
                # self.gga_vel = self.interpolate_composite()
                self.gga_vel.apply_composite(u_comp, v_comp, comp_source)
                self.gga_vel.interpolate_composite(transect)

            # Process vtg as primary
            if self.vtg_vel is not None:
                # Initialize the composite source
                comp_source = np.tile([np.nan], u_bt.shape)

                # Process the u velocity component
                u_comp = u_vtg
                comp_source[np.isnan(u_comp) == False] = 2

                # If GGA data are not valid try VTG and set composite source
                u_comp[np.isnan(u_comp)] = u_gga[np.isnan(u_comp)]
                comp_source[np.isnan(u_comp) == False & np.isnan(comp_source)] = 2

                # If there are still invalid boar velocities, try BT and set composite source
                u_comp[np.isnan(u_comp)] = u_bt[np.isnan(u_comp)]
                comp_source[np.isnan(u_comp) == False & np.isnan(comp_source)] = 1

                # If there are still invalid boat velocities, use interpolated values,
                # if present and set composite source
                # DSM changed 1/29/2018 u_comp[np.isnan(u_comp)] = self.gga_vel.u_processed_mps[np.isnan(u_comp)]
                comp_source[np.isnan(u_comp) == False & np.isnan(comp_source)] = 0

                # Set composite source to invalid for all remaining invalid boat velocity data
                comp_source[np.isnan(comp_source)] = -1

                # Process v velocity component.  Assume that the composite source is the
                # same as the u component
                v_comp = v_vtg
                # DSM wrong in Matlab version 1/29/2018 v_comp[np.isnan(v_comp)] = v_vtg[np.isnan(v_comp)]
                v_comp[np.isnan(v_comp)] = v_gga[np.isnan(v_comp)]
                v_comp[np.isnan(v_comp)] = v_bt[np.isnan(v_comp)]
                v_comp[np.isnan(v_comp)] = self.vtg_vel.v_processed_mps[np.isnan(v_comp)]

                # Apply the composite settings to the gga BoatData object
                # Apply the composite settings to the bottom track Boatdata objects
                # DSM changed 1/29/2018 self.gga_vel = self.apply_composite(u_comp,v_comp,comp_source)
                # self.gga_vel = self.interpolate_composite()
                self.gga_vel.apply_composite(u_comp, v_comp, comp_source)
                self.gga_vel.interpolate_composite(transect)
            # else:
                # Composite tracks false

                # Use only intepolations for bt
            # DSM Stopped here 1/29/2018


# TODO Finish coding BoatStructure class
