'''
Created on Jun 14, 2017

@author: gpetrochenkov
'''

import binascii
import numpy as np

class Pd0TRDI(object):
    '''Class to extract data from PD0 files
        Analalogous to matlab class: clsPd0TRDI
    '''
    
    def __init__(self):
        '''Constructor initializing properties'''
        
        self.file_name = None
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
        
    def initialize_pd0(self, pathname, files):
        '''Begins the pd0 file read process'''
        
        for x in files:
            self.pd_read(x, None)
            
   
        
    def pd_read(self, fullname, args):
        n_velocities = 4
        max_surface_bins = 5
        
        
                
        initial_pos = f._tell()-2
        
        
# -------------------------------------------------thus far testing individual pd0 file extraction until comfortable  
if __name__ == '__main__':
    files = [r'C:\Users\gpetrochenkov\Documents\Visual Studio 2017\Projects\TestFlask\PythonApplication1\drive-download-20170522T150040Z-001\RG_1308000_359\13038000_359_000.PD0']
    
    def find_ens_no(f):
    
        fileloc = f.tell()-2
        bytes_per_ens = np.fromfile(f, np.uint16, count=1)[0]
        f.seek(fileloc,0)
        Testb = np.fromfile(f,str,bytes_per_ens)
        print 'done'
        
            
    def number_of_ensembles(f):
        i=0
        leader_id='0000'
            
            
        while leader_id != '0x7f7f':
            f.seek(i,0)
            i=i+1
            leader_id = hex(np.fromfile(f,np.uint16, count = 1)[0])
            
        find_ens_no(f)
            
        
       
#     a = Pd0TRDI()
#     a.initialize_pd0('test', files)
    with open(files[0],'r') as f:
        leader_id =  hex(np.fromfile(f,np.uint16, count = 1)[0])
        print leader_id
        if leader_id != '0x7f7f':
            while leader_id != '0x7f7f':
                f.seek(-1,1)
                leader = hex(np.fromfile(f,np.uint16, count = 1)[0])
         
        initialPos = f.tell()-2   
            
        bytespPerEns = np.fromfile(f, dtype=np.uint8, count=1)[0]
        f.seek(1,1)
        nTypes =  np.fromfile(f,np.uint8, count = 1)[0]
        offset = np.fromfile(f,np.uint16, count = 1)[0]
        f.seek(initialPos+offset+8, 0)
        nBeams = np.fromfile(f,np.uint8,count=1)[0]
        nbins = np.fromfile(f, np.uint8,count=1)[0]
        
        number_of_ensembles(f)
        print 'done'
        
            
                
                
    
            
        
        
            