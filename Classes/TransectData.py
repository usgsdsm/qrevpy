import numpy as np
import os

class TransectData(object):
    """Class to hold Transect properties (may be removed on a refactor)
        Analogous to matlab class: clsTransectData
    """

    def __init__(self):
        self.adcp = None            # object of clsInstrument
        self.file_name = None       # filename of transect data file
        self.w_vel = None           # object of clsWaterData
        self.boat_vel = None        # class for various boat velocity references
                                    # (btVel, ggaVel, vtgVel)
        self.gps = None             # object of clsGPSData
        self.sensors = None         # object of clsSensorData
        self.depths = None          # object of clsDepthStructure for depth data including cell depths & ref depths
                                    # (btDepthds, vbDepths, dsDepths)
        self.edges = None           # object of clsEdges
                                    # (left and right object of clsEdgeData)
        self.extrap = None          # object of clsExtrapData
        self.start_edge = None      # starting edge of trasect looking sownstream (Left or Right)
        self.datetime = None
        self.checked = None         #transect was checked for use in mmt file assumed checked for SonTek
        self.in_transect_idx = None # index of ensemble data associated with the moving-boat portion of the transect


    def get_data(self, source, file, args):
        if source == 'TRDI':
            pass

    def TRDI(self, mmt, args=None):
        if args is None or args[0] == 'Q':
            transects = 'transects'
            active_config = 'active_config'


        elif args[0] == 'MB':
            transects = 'mbt_transects'
            activeconfig = 'mbtActiveConfig'
        
        if args is not None and args[1] == 1:
            file_names = [x.Files for x in mmt[transects] if x.Checked == 1]
        else:
            file_names = [x.Files for x in mmt[transects]]

          
        pathname = r"C:/"
        
        #determine if any files are missing
        num_valid_files = 0
        for x in file_names:
            if os.path.exists(x.path):
                num_valid_files += 1
                
        


