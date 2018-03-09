"""
Created on Aug 1, 2017

@author: gpetrochenkov
"""
import numpy as np

def cosd(angle):
    
    return np.cos(np.pi * angle/180)

def sind(angle):
    
    return np.sin(np.pi * angle/180)

def tand(angle):
    
    return np.tan(np.pi * angle/180)

def arctand(angle):
    
    return np.arctan(angle) * 180/np.pi

def cart2pol(x, y):
    
    rho = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x)
    
    return(phi, rho)

def pol2cart(phi, rho):
    
    x = rho * np.cos(phi)
    y = rho * np.sin(phi)
    
    return(x, y)

def iqr(data):
    q75, q25 = np.nanpercentile(data, [75 ,25])
    iqr = q75 - q25
    
    return iqr

def azdeg2rad(angle):
    direction = np.deg2rad(90-angle)
    idx= np.where(direction < 0)[0]
    if len(idx) > 0:
        direction[idx] = direction[idx]+ 2 * np.pi
        
    return direction

def rad2azdeg(angle):
    if isinstance(angle, float):
        deg = np.rad2deg(angle)
        deg = 90 - deg
        if deg < 0:
            deg += 360
            
        return deg
    else:
        #Multiple values
        deg = np.rad2deg(angle)
        deg = 90 - deg
        sub_zero = np.where(deg < 0)
        deg[sub_zero] = deg[sub_zero] + 360
        
        return deg

def nandiff(values):
    
    final_values = []
    for n in range(len(values) - 1):
        
        if np.isnan(values[n]):
            final_values.append(np.nan)
        else:
            i = n + 1
            while np.isnan(values[i]) and i < len(values) - 1:
                i += 1
                
            
            final_values.append(values[i] - values[n])
        
    return np.array(final_values)


def get_object_values(list_in, item, checked=None):
    if checked is not None:
        working_list = np.array(list_in)[checked]
    else:
        working_list = list_in

    out = []
    for obj in working_list:
        temp = getattr(obj, item)
        out.append(temp)
    return np.array(out)

def sontek_3d_arrange(data_in):
    r1 = np.squeeze(data_in[:, 0, :])
    r2 = np.squeeze(data_in[:, 1, :])
    r3 = np.squeeze(data_in[:, 2, :])
    r4 = np.squeeze(data_in[:, 3, :])
    new_array = np.array([r1, r2, r3, r4])
    return new_array


