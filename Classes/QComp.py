'''
Created on Sep 26, 2017

@author: gpetrochenkov
'''
from Classes.TransectData import TransectData
import numpy as np
from MiscLibs.convenience import cart2pol, pol2cart

class QComp(object):
    '''Computes the discharge for each transect'''
    
    def __init__(self):
        
        self.top = None #Transect total extrapolated top discharge 
        self.middle = None # Transect toal measured middle discharge including interpolations 
        self.bottom = None # ETransect total extrapolated bottom discharge
        self.top_ens = None #Extrapolated top discharge by ensemble
        self.middle_cells = None #Measured middle discharge including interpolation by cell
        self.bottom_ens = None #Extrapolate bottom discharge by ensemble 
        self.left = None #Left edge discharge
        self.left_idx = None #Ensembles used for left edge
        self.right = None #Right edge discharge
        self.right_idx = None #Ensembles used for right edge
        self.total_uncorrected = None #Total discharge for transect uncorrected for moving-bed, if required
        self.total = None #Total discharge with moving-bed correction applied if necessary 
        self.correction_factor = None #Moving-bed correction factor, if required 
        self.int_cells = None #Total discharge computed for invalid depth cells excluding invalid ensembles
        self.int_ens = None #Total discharge computed for invalid ensembles
        
    def populate_data(self, data_in, kargs = None):
        '''Discharge is computed using the data provided to the method.  Water data provided are assumed to
        becorrected for the navigation reference.  If a moving-bed correction is to be applied it is computed
        and applied.  The TRDI method using expanded delta time is applied if the processing method is WR2.
        
        Input:
        data_in: Object TransectData or Measurement
        
        kargs: Allows overrideing the extrapolation methods specified in data in
        kargs[0]: top extrapolation method
        kargs[1]: bottom extrapolation method
        kargs[2]: extrapolation exponent
        '''
        if isinstance(data_in, TransectData):
            trans_data = data_in
            
            #Use bottom track interpolation settings to determine the appropriate algorithms to apply\
            if trans_data.boat_vel.bt_vel.interpolate == 'None':
                processing = 'WR2'
            elif trans_data.boat_vel.bt_vel.interpolate == 'Linear':
                processing = 'QRev'
            else:
                processing = 'RSL'
                
        else:
            
            #If the data in is a Measurement assign variables
            meas = data_in
            trans_data = meas.transects
            processing = meas.processing
            
        correction_flag = False
        
        #Compute cross product
        x_prod = self.cross_product(trans_data)
        
        #Get index of ensembles in moving-boat portion of transect
        in_transect_idx = trans_data.in_transect_idx
        
        if processing == 'WR2':
            #TRDI uses expanded delta time to handle invalid ensembles which can be caused by invalid BT
            #WT, or depth.  QRev by default handles this invalid data through linear interpolation of the
            #invalid data through linear interpolation of the invalid data type.  This if statement and
            #associated code is required to maintain compatibility with WinRiver II discharge computations.
            
            #Determine valid ensembles
            valid_ens = np.any(np.isnan(x_prod) == False) 
            valid_ens = valid_ens[in_transect_idx]
            
            #Compute the ensemble duration using TRDI approach of expanding delta time to compensate
            #for invalid ensembles
            n_ens = len(valid_ens)
            ens_dur = trans_data.date_time.ens_duration_sec[in_transect_idx]
            delta_t = np.tile([np.nan], n_ens)
            cum_dur = 0
            idx = 1
            for j in range(idx,n_ens):
                cum_dur = np.nansum(np.hstack([cum_dur, ens_dur[j]]))
                if valid_ens[j] == True:
                    delta_t[j] = cum_dur
                    cum_dur = 0
                    
        else:
            #For non-WR2 processing use actual ensemble duration
            delta_t = trans_data.date_time.ens_duration_sec[in_transect_idx]
            
        #Compute measured or middle discharge
        q_mid_cells = self.discharge_middle_cells(x_prod, trans_data, delta_t)
        self.middle_cells = q_mid_cells
        self.middle_ens = np.nansum(q_mid_cells, 0)
        self.middle = np.nansum(self.middle_ens, 0)
        
        #Compute the top discharge
        self.top_ens = self.discharge_top(x_prod, trans_data, delta_t, kargs)
        
        self.top = np.nansum(self.top_ens)
        
        #Compute the bottom discharge
        self.bottom_ens = self.discharge_bot(x_prod, trans_data, delta_t, kargs)
        self.bottom = np.nansum(self.bottom_ens)
        
        #Compute interpolated cell and ensemble discharge from computed
        #measured discharge
        self.int_cells, self.int_ens = self.discharge_interpolated(self.top_ens, q_mid_cells, self.bottom_ens, trans_data)
        
        #Compute left edge discharge
        if trans_data.edges.left.type == 'User Q':
            self.right = self.discharge_edge('right', trans_data, kargs)
        else:
            self.right = trans_data.edges.right.__user_Q_cms
            
        #Compute left edge discharge
        if trans_data.edges.left.type == 'User Q':
            self.left = self.discharge_edge('left', trans_data, kargs)
        else:
            self.left = trans_data.edges.left.__user_Q_cms
            
        #Compute moving-bed correction, if applicable.  Two checks are used to account for the
        #way the meas object is created.
        
        #Check to see if the mb_tests property of Measurement exists
        try:
            getattr(data_in, 'mb_tests')
            mb_data = meas.mb_tests
            
            if mb_data is not None:
                
                #Determine if any of the moving-bed tests indicated a moving bed
                mb_valid = mb_data.selected
                if mb_data[mb_valid].__moving_bed == 'Yes':
                    
                    use_2_correct = mb_data[:].__use_2_correct
                    
                    #Determine if a moving-bed test is to be used for correction
                    if np.sum(use_2_correct) > 0:
                        
                        #Make sure bottom track is the navigation reference and composite
                        #tracks are turned off
                        
                        if trans_data.boat_vel.selected == 'bt_vel' and trans_data.boat_vel.composite == 'Off':
                            
                            #Apply appropriate moving-bed test correction method
                            if np.sum(mb_data[use_2_correct].type == 'Stationary') > 0:
                                self.correction_factor = self.stationary_correction(self.top, self.middle, self.bottom, trans_data, mb_data, delta_t)
                            else:
                                self.correction_factor = self.loop_correction(self.top, self.middle, self.bottom, trans_data, mb_data[use_2_correct], delta_t)
                    #
                    else:
                        
                        #Set a flag to generate a warning
                        correction_flag = True
        except:
            pass
        
        self.total_uncorrected = np.nansum(np.hstack([self.top,self.middle,self.bottom,self.left, self.right]))
        
        #Compute final discharge using correction if applicable
        if self.correction_factor is None or self.correction_factor == 1:
            self.total = self.total_uncorrected
        else:
            self.total = self.left + self.right + (self.middle + self.bottom + self.top) * self.correction_factor
            
            
    def cross_product(self, kargs):
        '''Computes the cross product using data from an object TransectData or input
        for water and navigation arrays
        
        Input:
        kargs: variable argument in allows wither an object of TransectData of arrays 
        of water and navigation data to be used
            if TransectData:
                kargs[0] is an object of TransectData that is only input
            if arrays of water and navigation data
                kargs[0] is array of water velocity in the x direction
                kargs[1] is array of water velocity in the y direction
                kargs[2] is vector of navigation velocity in x-direction
                kargs[3] is vector of naviagation velocity in y-direction'''
        
        if len(kargs) < 2:
            #Assign object of TransectData to local Variable
            transect = kargs[0]
            
            #Prepare water track data
            cells_above_sl = np.array(transect.w_vel._WaterData__cells_above_sl).astype(np.double)
            cells_above_sl[cells_above_sl < 0.5] = np.nan
            w_vel_x = transect.w_vel._WaterData__u_processed_mps * cells_above_sl
            w_vel_y = transect.w_vel._WaterData__v_processed_mps * cells_above_sl
            
            #Get navigation data from object properties
            trans_select = getattr(transect.boat_vel, transect.boat_vel.selected)
            if trans_select is not None:
                b_vel_x = trans_select._BoatData__u_processed_mps
                b_vel_y = trans_select._BoatData__v_processed_mps
            else:
                b_vel_x = np.tile([np.nan], transect.boat_vel.bt_vel._BoatData__u_processed_mps.shape)
                b_vel_y = np.tile([np.nan], transect.boat_vel.bt_vel._BoatData__v_processed_mps.shape)
                
            start_edge = transect.start_edge
            
        else:
            
            #Assign data arrays to local variables
            w_vel_x = kargs[0]
            w_vel_y = kargs[1]
            b_vel_x = kargs[2]
            b_vel_y = kargs[3]
            start_edge = kargs[4]
            
        #Compute the cross product
        xprod = np.multiply(w_vel_x, b_vel_y) - np.multiply(w_vel_y, b_vel_x)
        
        if start_edge == 'Right':
            dir = 1
        else:
            dir = -1
        xprod = xprod * dir
        
        return xprod
    
    def discharge_middle_cells(self, xprod, transect, delta_t):
        '''Computes the discharge in the measured or middle portion of the cross section
        
        Input:
        xprod: cross product computed from the cross product method
        transect: object of TransectData
        delta_t: duration of each ensemble computed from QComp
        
        Output:
        q_mid_cells: discharge in each bin or depth cell
        '''
        
        #Assign properties from transect object to local variables
        in_transect_idx = transect.in_transect_idx
        trans_select = getattr(transect.depths, transect.depths.selected)
        cell_size = trans_select.depth_cell_size_m
        
        #Determine is xprod contains edge data and process appropriately
        if len(xprod) > len(in_transect_idx):
            q_mid_cells = np.multiply(xprod[:,in_transect_idx] * cell_size[:,in_transect_idx], delta_t)
        else:
            q_mid_cells = np.multiply(xprod * cell_size[:, in_transect_idx], delta_t)
            
        return q_mid_cells
    
    #Top Extrapolation
    def discharge_top(self, xprod, transect, delta_t, kargs = None):
        '''Coordinates computation of the extrapolated top discharge
        
        Input:
        xprod: cross product computed from the cross product method
        transect: object of TransectData
        delta_t: duration of each ensemble computed from QComp
        
        kargs: allows specifying top and bottom extrapolations
        kargs[0]: top extrapolation method
        kargs[1]: not used in this method
        kargs[2]: exponent
        
        Output:
        q_top: top extrapolated discharge for each ensemble
        '''
        
        if kargs is not None:
            top_method = kargs[0]
            exponent = kargs[2]
        else:
            top_method = transect.extrap.top_method
            exponent = transect.extrap.exponent
            
        #Get index for ensembles in moving-boat portion of transect
        in_transect_idx = transect.in_transect_idx
        
        #Compute top variables
        idx_top, idx_top3, top_rng = self.top_variables(xprod, transect)
        idx_top = idx_top[in_transect_idx]
        idx_top3 = idx_top3[:,in_transect_idx]
        top_rng = top_rng[in_transect_idx]
        
        #Get data from transect object
        trans_select = getattr(transect.depths, transect.depths.selected)
        cell_size = trans_select.depth_cell_size_m[:, in_transect_idx]
        cell_depth = trans_select.depth_cell_depth_m[:, in_transect_idx]
        depth_ens = trans_select.depth_processed_m[in_transect_idx]
        
        #Compute z
        z = np.subtract(depth_ens, cell_depth)
        
        #Use only valid data
        valid_data = np.isnan(xprod[:,in_transect_idx]) == False
        z[valid_data == False] = np.nan
        cell_size[valid_data == False] = np.nan
        cell_depth[valid_data == False] = np.nan
        
        #Compute top discharge
        q_top = self.extrapolate_top(top_method,exponent,idx_top,idx_top3, top_rng,
                                     xprod[:,in_transect_idx], cell_size, cell_depth,
                                     depth_ens, delta_t, z)
        
        return q_top
    
    def extrapolate_top(self, top_method,exponent,idx_top,idx_top_3, top_rng,
                                     component, cell_size, cell_depth,
                                     depth_ens, delta_t, z):
        '''Computes the top extrapolated value of the provided component
        using the specified extrapolation method
        
        Input:
        top_method: top extrapolation method (Power, Constane, 3-Point)
        exponent: exponent for the power extrapolation method
        idx_top: index to the topmost valid depth cell in each ensemble
        idx_top_3: index to the top 3 valid depth cells in each ensemble
        top_rng: range from the water surface to the top of the topmost cell
        component: the variable to be extrapolated
        cell_size: array of cellsizes (n cells x n ensembles)
        cell_depth: depth of each cell (n cells x n ensembles)
        depth_ens: bottom depth for each ensemble
        delta_t: duration of each ensemble compute by QComp
        z: relatice depth from the bottom of each depth cell computed in discharge top method 
        
        Output:
        top_value: total for the specified component integrated over the top range
        '''
        
        #Identify method
        
        #Top power extrapolation
        if top_method == 'Power':
            coef = ((exponent+1) * np.nansum(component * cell_size)) / \
                    np.nansum(((z + 0.5 * cell_size)**(exponent+1)) - \
                    ((z - 0.5 * cell_size)**(exponent+1)))
            top_value = delta_t * (coef / (exponent + 1)) * \
                    (depth_ens**(exponent + 1) - (depth_ens-top_rng)**(exponent + 1))
        
        #Top constant extrapolation
        if top_method == 'Constant':
            n_ensembles = len(delta_t)
            top_value = np.tile([np.nan], (1, n_ensembles))
            for j in range(n_ensembles):
                if idx_top[j] != 0:
                    top_value[j] = delta_t[j] * component(idx_top[j],j) * top_rng[j]
                    
        #Top 3-point extrapolation
        if top_method == '3-Point':
            #Determine number of bins available in each profile
            valid_data = np.isnan(component) == False
            n_bins = np.nansum(valid_data, 0)
            #Determine number of ensembles
            n_ensembles = len(delta_t)
            #Preallocate qtop vector
            top_value = np.tile([np.nan], (1, n_ensembles))
            
            for j in range(n_ensembles):
                
                if n_bins[j] < 6 and n_bins[j] > 0 and idx_top[j] != 0:
                    top_value[j] = delta_t * component[idx_top[j],j] * top_rng[j]
                
                #If 6 or more bins use 3-pg at top
                if n_bins[j] > 5:
                    sumd = np.nansum(cell_depth[idx_top_3[0:3,j]])
                    sumd2 = np.nansum(cell_depth[idx_top_3[0:3, j]]**2)
                    sumq = np.nansum(component[idx_top_3[0:3, j]])
                    sumqd = np.nansum(component[idx_top_3[0:3, j],j])
                    delta = 3 * sumd2 - sumd **2
                    A = (3 * sumqd - sumq * sumd) / delta
                    B = (sumq * sumd2 - sumqd * sumd) / delta
                    #Compute discharge for 3-pt fit
                    Qo = (A * top_rng[j]**2) / (2 + B * top_rng[j])
                    top_value[j] = delta_t[j] * Qo
                    
        return top_value
    
    def top_variables(self, xprod, transect):
        '''Computes the index to the top and top3 valid cells in each ensemble and
        the range from the water surface to the top of the topmost cell
        
        Input:
        xprod: cross product computed from the cross product method
        transect: object of TransectData
        
        Output:
        idx_top: index to the topmost valid depth cell in each ensemble
        idx_top_3: index to the top 3 valid depth cell in each ensemble
        top_rng: range from the water surface to the top of the topmost cell
        '''
        
        #Get data from transect object
        valid_data1 = transect.w_vel.valid_data[:,:,0]
        valid_data2 = np.isnan(xprod) == False
        valid_data = valid_data1 * valid_data2
        trans_select = getattr(transect.depths, transect.depths.selected)
        cell_size = trans_select.depth_cell_size_m
        cell_depth = trans_select.depth_cell_depth_m
        
        #Preallocate variables
        n_ensembles = valid_data.shape[1]
        idx_top = np.zeros(1, valid_data.shape[1])
        idx_top_3 = np.zeros(3, valid_data.shape[1])
        top_rng = np.tile([np.nan], (1, n_ensembles))
        
        #Loop through ensembles
        for n in range(n_ensembles):
            #Identify topmost 1 and 3 valid cells
            idx_temp = np.where(np.isnan(xprod[:,n]))[0]
            if len(idx_temp) > 0:
                idx_top[n] = idx_temp[0]
                if len(idx_temp) > 2:
                    idx_top_3[:,n] = idx_temp    
                #Compute top range
                top_rng[n] = cell_depth[idx_top[n], n] - 0.5 * cell_size[idx_top[n], n]
            else:
                top_rng[n] = 0
                idx_top[n] = 1
        
        return [idx_top, idx_top_3, top_rng]
    
    #Bottom Extrapolation
    def discharge_bot(self, xprod, transect, delta_t, kargs  = None):
        '''Coordinates computation of the extrapolated bottom discharge
        
        Input:
        xprod: cross product computed from the corss product method
        transect: object of TransectData
        delta_t: duration of each ensemble
        kargs:
        kargs[0]: not used in this method
        kargs[1]: bottom extrapolation method
        kargs[2]:exponent
        
        Output:
        q_bot: bottom extrpaolated discharge for each ensemble
        '''
        
        #Determine extrapolation methods and exponent
        if kargs is not None:
            bot_method = kargs[1]
            exponent = kargs[3]
        else:
            bot_method = transect.extrap.__bot_method
            exponent = transect.extrap.exponent
            
        #Get index for ensembles in mocing-boat portion of transect
        in_transect_idx = transect.in_transect_idx
        xprod = xprod[:, in_transect_idx]
        
        #Compute bottom variables
        idx_bot, bot_rng = self.bot_variables(xprod, transect)
        
        #Get data from transect properties
        trans_select = getattr(transect.depths, transect.depths.selected)
        cell_size = trans_select.depth_cell_size_m[:, in_transect_idx]
        cell_depth = trans_select.depth_cell_size_m[:, in_transect_idx]
        depth_ens = trans_select.depth_processed_m[in_transect_idx]
        
        #Compute z
        z = np.subtract(depth_ens, cell_depth)
        valid_data = np.isnan(xprod) == False
        z[valid_data == False] = np.nan
        z[z < 0] = np.nan
        cell_size[valid_data == False] = np.nan
        
        #Computre bottom discharge
        q_bot = self.extrapolate_bottom(bot_method, exponent, idx_bot, bot_rng, xprod, 
                                        cell_size, cell_depth, depth_ens, delta_t, z)
        
    def bot_value(self, bot_method, exponent, idx_bot, bot_rng, component, 
                                        cell_size, cell_depth, depth_ens, delta_t, z):
        '''Computes the bottom extrapolated value of the provided component
        using the specified extrapolation method
        
        Input:
        bot_method: bottom extrapolation method (Power, Constant, 3-Point)
        exponent: exponent for the power
        idx_bot: index to the bottom most valid depth cell in each ensemble
        bot_rng: range from the streambed to the bottom of the bottom most cell
        component: the variable to be extrapolated
        cell_size: array of cellsizes (n cells x n ensembles)
        cell_depth: depth of each cell (n cells x n ensembles)
        depth_ens: bottom depth for each ensemble
        delta_t: duration of each ensemble computed by QComp
        z: relative depth from the bottom of each depth cell comoputed in discharge top
        
        Output:
        bot_value: total for the specified component integrated over the top range
        '''
        
        #Identify the bottom method
        #Bottom power extrapolation
        if bot_method == 'Power':
            coef = ((exponent+1) * np.nansum(component * cell_size)) / \
                np.nansum(((z + 0.5 * cell_size)**(exponent + 1)) - \
                          ((z - 0.5 * cell_size))**(exponent +1))
        
        #Bottom no slip extrapolation
        elif bot_method == 'No Slip':
            #Valid data in the lower 20% of the water column or
            #the last valid depth cell are used to compute the no slip power fit
            cutoff_depth = .8 * depth_ens 
            depth_OK = (cell_depth > np.tile(cutoff_depth, (cell_depth.shape[0], 1)))
            component_OK = np.isnan(component) == False
            use_NS = depth_OK * component_OK
            for j in range(len(delta_t)):
                if idx_bot[j] != 0:
                    use_NS[idx_bot[j],j] = 1
                    
            use_NS[use_NS == 0] = np.nan
            
            #Create cross product and z arrays for the data to be used in
            #no slip computations
            component_NS = use_NS * component
            z_NS = use_NS * z
            coef = ((exponent + 1) * np.nansum(component_NS * cell_size)) / \
                    np.nansum(((z_NS + 0.5 * cell_size)**(exponent + 1)) - \
                              ((z_NS - 0.5 * cell_size)**(exponent + 1)))
                    
            #Compute the bottom discharge of each profile
            bot_value = delta_t * (coef / (exponent + 1)) * (bot_rng**(exponent+1))
            
            return bot_value
        
    def bot_variables(self, x_prod, transect):
        '''Computes the index to the bottom most valid cell in each ensemble and the range from
        the bottom to the bottom of the bottom most cell
        
        Input:
        x_prod: cross product computed from the cross product method
        transect: object of TransectData
        
        Output:
        idx_bot: index to the bottom most valid depth cell in each ensemble
        bot_rng: range from the streambed to the bottom of the bottom most cell
        '''
        
        #Identify valid data
        in_transect_idx = transect.in_transect_idx
        valid_data1 = transect.w_vel.valid_data[:, in_transect_idx]
        valid_data2 = np.isnan(x_prod) == False
        valid_data = valid_data1 * valid_data2
        
        #Assign transect properties to local variables
        trans_selected = getattr(transect.depths, transect.depths.selected)
        cell_size = trans_selected.depth_cell_size_m[:, in_transect_idx]
        cell_depth = trans_selected.depth_cell_depth_m[:, in_transect_idx]
        depth_ens = trans_selected.depth_processed_m[in_transect_idx]
        
        #Preallocate variables
        n_ensembles = valid_data.shape[1]
        idx_bot = np.zeros((1, valid_data.shape[1]))
        bot_rng = np.tile([np.nan], (1, n_ensembles))
        
        for n in range(n_ensembles):
            #Identifying bottom most valid cell
            idx_temp = np.where(valid_data[:,n] == True)[0]
            if len(idx_temp) > 0:
                idx_temp = idx_temp[-1]
                idx_bot[n] = idx_temp
                #Compute bottom range
                bot_rng[n] = depth_ens[n] - cell_depth[idx_bot[n],n] - 0.5 * cell_size[idx_bot]
            else:
                bot_rng[n] = 0
                
        return (idx_bot, bot_rng)
                
    def discharge_edge(self, edge_loc, transect, kargs = None):
        '''Computes edge discharge
        
        Input:
        edge_loc: edge location (left or right)
        transect: object of TransectData
        kargs:
            kargs[0]: top extrapolation method
            kargs[1]: bottom extrapolation method
            kargs[2]: exponent
            
        Output:
        edge computed edge discharge
        '''
        
        #Determine what ensembles to use for edge computation.
        #The method of determining varies by manufacturer
        edge_idx = self.edge_ensembles(edge_loc, transect)
        
        #Average depth for the edge ensembles
        trans_select = getattr(transect.depths, transect.depths.selected)
        depth = trans_select.depth_processed_m[edge_idx]
        depth_avg = np.nanmean(depth)
        
        #Edge distance
        edge_selected = getattr(transect.edges, edge_loc)
        edge_dist = edge_selected.dist_m
        
        #Compute edge velocity and sign
        edge_vel_mag, edge_vel_sign = self.edge_velocity(edge_idx, transect, kargs)

        #Compute edge coefficient
        coef = self.edge_coef(edge_loc, transect)

        #Compute edge discharge
        edge = coef * depth_avg * edge_vel_mag * edge_dist * edge_vel_sign
        if np.isnan(edge):
            edge = 0
            
        return edge
    
    def edge_ensembles(self,edge_loc, transect):
        '''This function computes the starting and ending ensemble numbers for an edge using
        either the method used by TRDI which used the specified number of valid ensembles or SonTek
        which uses the specified number of ensembles prior to screening for valid data
        
        Input:
        edge_loc: edge location (left or right)
        transect: TransectData
        
        Output:
        edge_idx: indices of ensembles used to compute edge discharge
        '''
        
        #Assign number of ensembles in edge to local variable
        edge_select = getattr(transect.edges, edge_loc)
        num_edge_ens = edge_select.num_ens_2_avg
        
        #TRDI method
        if transect.adcp.manufacturer == 'TRDI':
            #Determine the indices of the edge ensembles which contain
            #the specified number of valid ensembles
            valid_ens = self.valid_ens(transect)
            if edge_loc == transect.start_edge:
                edge_idx = np.where(valid_ens == True)[0][0]
            else:
                edge_idx = np.where(valid_ens == True)[0][-1]
        #Sontek Method
        else:
            #Determine the indices of the edge ensembles as collected by RiverSurveyor.  There
            #is no check as to whether the ensembles contain valid data
            if edge_loc == transect.start_edge:
                edge_idx = np.arange(0,num_edge_ens)
            else:
                trans_select = getattr(transect.depths, transect.depths.selected)
                n_ensembles = trans_select.depth_processed_m
                edge_idx = np.arange(n_ensembles - num_edge_ens, n_ensembles)
                
        return edge_idx
    
    def edge_velocity(self, edge_idx, transect, kargs = None):
        '''Coordinates computation of edge velocity
        
        Input:
        edge_idx: indices of ensembles used to compute edge discharge
        transect: object of TransectData
        
        kargs: allows specifying top and bottom extrapolations independent of transect object
            kargs[0]: top extrapolation method
            kargs[1]: bottom extrpolation method
            kargs[2]: exponent
            
        Output:
        edge_vel_mag: magnitude of edge velocity
        edge_vel_sign: sign of edge celocity (discharge)
        '''
        
        #Check to make sure there is edge data
        if edge_idx is not None and len(edge_idx) > 0:
            
            #Compute edge velocity using specified method
            #Used by TRDI
            if transect.edges.vel_method == 'MeasMag':
                edge_vel_mag, edge_vel_sign = self.edge_velocity_TRDI(edge_idx, transect)
                
            #Used by Sontek    
            elif transect.edges.vel_method == 'VectorProf':
                edge_val_mag, edge_vel_sign = self.edge_velocity_SonTek(edge_idx, transect, kargs)
                
            #USGS proposed method
            elif transect.edges.vel_method == 'Profile':
                edge_vel_mag, edge_vel_sign = self.edge_velocity_profile(edge_idx, transect)
        
        #If no data set edge velocity to 0    
        else:
            edge_vel_sign = 1
            edge_vel_mag = 0
        
        return (edge_vel_sign, edge_vel_mag)
    
    def edge_velocity_TRDI(self, edge_idx, transect):
        '''Computes edge velocity magnitude and sign using TRDI's method, which uses only
        the measured data and no extrapolation
        
        Input:
        edge_idx: indices of ensembles used to compute edge discharge
        transect: object of TransectData
        
        Output:
        edge_vel_mage: magnitude of edge velocity
        edge_vel_sign: sign of edge velocity (dicharge)
        '''
        
        #Assign water velocity to local variables
        x_vel = transect.w_vel.__u_processed_mps[:, edge_idx]
        y_vel = transect.w_vel.__v_processed_mps[:, edge_idx]
        
        #Use only valid data
        valid = transect.w_vel.valid_data[:, edge_idx, 1].astype(np.double)
        valid[valid == 0] = np.nan
        x_vel = x_vel * valid
        y_vel = y_vel * valid
        
        #Compute the mean velocity components
        x_vel_avg = np.nanmean(np.nanmean(x_vel, 1))
        y_vel_avg = np.nanmean(np.nanmean(y_vel, 1))
        
        #Compute magnitude and direction
        edge_dir, edge_vel_mag = cart2pol(x_vel_avg, y_vel_avg)
        
        #Compute unit vector to help determine sign
        unit_water_x, unit_water_y = pol2cart(edge_dir, 1)
        if transect.start_edge == 'Right':
            dir = 1
        else:
            dir = -1
            
        #Compute unit boat vector to help determine sign
        ens_delta_time = transect.date_time.ens_duration_sec
        in_transect_idx = transect.in_transect_idx
        trans_selected = getattr(transect.boat_vel, transect.boat_vel.selected)
        if trans_selected is not None:
            b_vel_x = trans_selected.__u_proccesed_mps
            b_vel_y = trans_selected.__v_processed_mps
        else:
            b_vel_x = np.tile([np.nan], transect.boat_vel.__u_processed_mps.shape)
            b_vel_y = np.tile([np.nan], transect.boat_vel.__v_processed_mps.shape)
            
        track_x = np.nancumsum(b_vel_x[in_transect_idx] * ens_delta_time[in_transect_idx])
        track_y = np.nancumsum(b_vel_y[in_transect_idx] * ens_delta_time[in_transect_idx])
        boat_dir, boat_mag = cart2pol(track_x[-1], track_y[-1])
        unit_track_x, unit_track_y = pol2cart(boat_dir, 1)
        unit_x_prod = (unit_water_x * unit_track_y - unit_water_y * unit_track_x) * dir
        edge_vel_sign = np.sign(unit_x_prod)
        
        return (edge_vel_mag, edge_vel_sign)
    
    def edge_coef(self, edge_loc, transect):
        '''Returns the edge coefficient based on the edge settings and transect object
        
        Input:
        edge_loc: edge location (left_or right)
        transect: object of Transect DAta
        
        Output:
        coef: edge coefficient for accounting for velocity distribution and edge shape
        '''
        
        #Process appropriate edge type
        edge_select = getattr(transect.edges, edge_loc)
        if edge_select.type == 'Triangular':
            coef = 0.3535
            
        elif edge_select.type == 'Rectangular':
            #Rectangular edge coefficient dependds on the rec_edge_method.  
            #'Fixed' is compatible with the method used by TRDI.
            #'Variable is compatible with the method used by SonTek
            
            if transect.edges.__rec_edge_method == 'Fixed':
                
                #Fixed Method
                coef = 0.91
                
            else:
                
                #Variable method
                
                #Get edge distance
                dist = edge_select.dist_m
                
                #Get edge ensembles to use
                edge_idx = self.edge_ensembles(edge_loc, transect)
                
                trans_select = getattr(transect.depths, transect.depths.selected)
                #Compute the mean depth for edge
                depth_edge = np.nanmean(trans_select.depth_processed_m[edge_idx])
                
                #Compute coefficient using equation 34 from Principle of River Discharge Measurement, SonTek, 2003
                coef = (1 - ((0.35 / 4) * (depth_edge / dist) * (1 - np.exp(-4 * (dist / depth_edge))))) / \
                    (1 - 0.35 * np.exp(-4 * (dist / depth_edge)))
                    
        elif edge_select.type == 'Custom':
            coef = edge_select.cust_coef
            
        else:
            coef = []
           
        return coef  
    
    def correction_factor(self, top_q, middle_q, bottom_q, trans_data, mb_data, delta_t):
        '''Computes the discharge correction factor from stationary moving-bed tests
        
        Input:
        top_q: top discharge from extrapolation
        middle_q: computed middle discharge
        bottom_q: bottom discharge from extrapolation
        trans_data: object of TransectData
        mb_data: object of MovingBedTests
        delta_t: duration of each ensemble, computed in QComp
        
        Output:
        correction_factor: correction factor to be applied to the discharge to correct for moving-bed tests
        '''
        
        #Find and composite results of all stationary test that are selected for use
        n_mb_tests = len(mb_data)
        n_sta_tests = 0
        
               
        
                
                
            
        
        
        
            