'''
Created on Aug 1, 2017

@author: gpetrochenkov
'''
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
    
    return(rho, phi)

def pol2cart(rho, phi):
    
    x = rho * np.cos(phi)
    y = rho * np.sin(phi)
    
    return(x, y)

def iqr(data):
    q75, q25 = np.percentile(data, [75 ,25])
    iqr = q75 - q25
    
    return iqr

def azdeg2rad(angle):
    direction = np.deg2rad(90-angle);
    idx= np.where(direction<0)[0];
    if len(idx) > 0:
        direction[idx] = dir[idx]+ 2 * np.pi;
        
    return direction
