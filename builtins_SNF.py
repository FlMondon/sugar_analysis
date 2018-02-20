# -*- coding: utf-8 -*-
"""
Created on Fri Nov 24 13:24:34 2017

@author: mondon
"""

import sncosmo 
import astropy.units as u
import numpy as np
import pylab as P
import pyfits
from scipy.interpolate import InterpolatedUnivariateSpline as Spline1d
CLIGHT = 2.99792458e18         # [A/s]
HPLANCK = 6.62606896e-27        # [erg s]

def register_SNf_bands():
    
#    wl = np.arange(3000,9500,0.1)
#    transmission_U = np.zeros(len(wl))
#    transmission_B = np.zeros(len(wl))
#    transmission_V = np.zeros(len(wl))
#    transmission_R = np.zeros(len(wl))
#    transmission_I = np.zeros(len(wl))

#    for i in range (len(wl)): 
#        if wl[i] >= 3300 and wl[i] <= 4102:
#            transmission_U[i]= 1.
#        if wl[i] >= 4102 and wl[i] <= 5100:
#            transmission_B[i]= 1. 
#        if wl[i] >= 5200 and wl[i] <= 6289:
#            transmission_V[i]= 1.
#        if wl[i] >= 6289 and wl[i] <= 7607 :
#            transmission_R[i]= 1.
#        if wl[i] >= 7607 and wl[i] <= 9200 :
#            transmission_I[i]= 1.        
#    band_U = sncosmo.Bandpass(wl,transmission_U,wave_unit=u.AA,name='USNf')    
#    sncosmo.registry.register(band_U)
#    band_B = sncosmo.Bandpass(wl,transmission_B,wave_unit=u.AA,name='BSNf')    
#    sncosmo.registry.register(band_B)
#    band_V = sncosmo.Bandpass(wl,transmission_V,wave_unit=u.AA,name='VSNf')    
#    sncosmo.registry.register(band_V)
#    band_R = sncosmo.Bandpass(wl,transmission_R,wave_unit=u.AA,name='RSNf')    
#    sncosmo.registry.register(band_R)
#    band_I = sncosmo.Bandpass(wl,transmission_I,wave_unit=u.AA,name='ISNf')    
#    sncosmo.registry.register(band_I)
            
    band_file_U = 'USNf_3300-4102.dat'
    band_file_B = 'BSNf_4102-5100.dat'
    band_file_V = 'VSNf_5200-6289.dat'
    band_file_R = 'RSNf_6289-7607.dat'
    band_file_I = 'ISNf_7607-9200.dat'
       
        
    filt2 = np.genfromtxt('../sugar_analysis_data/data/Instruments/SNIFS/'+ band_file_U)
    wl_U = filt2[:,0]
    transmission_U = filt2[:,1]

    filt2 = np.genfromtxt('../sugar_analysis_data/data/Instruments/SNIFS/'+ band_file_B)
    wl_B = filt2[:,0]
    transmission_B = filt2[:,1]

    filt2 = np.genfromtxt('../sugar_analysis_data/data/Instruments/SNIFS/'+ band_file_V)
    wl_V = filt2[:,0]
    transmission_V = filt2[:,1]

    filt2 = np.genfromtxt('../sugar_analysis_data/data/Instruments/SNIFS/'+ band_file_R)
    wl_R = filt2[:,0]
    transmission_R = filt2[:,1]

    filt2 = np.genfromtxt('../sugar_analysis_data/data/Instruments/SNIFS/'+ band_file_I)
    wl_I = filt2[:,0]
    transmission_I = filt2[:,1]
    
    
    


#            
#    band_U = sncosmo.Bandpass(wl_U,transmission_U,wave_unit=u.AA,name='USNf',top_hat=True)    
#    sncosmo.registry.register(band_U)
#    band_B = sncosmo.Bandpass(wl_B,transmission_B,wave_unit=u.AA,name='BSNf',top_hat=True)    
#    sncosmo.registry.register(band_B)
#    band_V = sncosmo.Bandpass(wl_V,transmission_V,wave_unit=u.AA,name='VSNf',top_hat=True)    
#    sncosmo.registry.register(band_V)
#    band_R = sncosmo.Bandpass(wl_R,transmission_R,wave_unit=u.AA,name='RSNf',top_hat=True)    
#    sncosmo.registry.register(band_R)
#    band_I = sncosmo.Bandpass(wl_I,transmission_I,wave_unit=u.AA,name='ISNf',top_hat=True)    
#    sncosmo.registry.register(band_I)
#    
    band_U = sncosmo.Bandpass(wl_U,transmission_U,wave_unit=u.AA,name='USNf')    
    sncosmo.registry.register(band_U)
    band_B = sncosmo.Bandpass(wl_B,transmission_B,wave_unit=u.AA,name='BSNf')    
    sncosmo.registry.register(band_B)
    band_V = sncosmo.Bandpass(wl_V,transmission_V,wave_unit=u.AA,name='VSNf')    
    sncosmo.registry.register(band_V)
    band_R = sncosmo.Bandpass(wl_R,transmission_R,wave_unit=u.AA,name='RSNf')    
    sncosmo.registry.register(band_R)
    band_I = sncosmo.Bandpass(wl_I,transmission_I,wave_unit=u.AA,name='ISNf')    
    sncosmo.registry.register(band_I)
    
    
    band_file_U = 'fU.dat'
    band_file_B = 'fB.dat'
    band_file_V = 'fV.dat'
    band_file_R = 'fR.dat'

       
        
    filt2 = np.genfromtxt('../sugar_analysis_data/data/Instruments/Florian/'+ band_file_U)
    wl_U = filt2[:,0]
    transmission_U = filt2[:,1]

    filt2 = np.genfromtxt('../sugar_analysis_data/data/Instruments/Florian/'+ band_file_B)
    wl_B = filt2[:,0]
    transmission_B = filt2[:,1]

    filt2 = np.genfromtxt('../sugar_analysis_data/data/Instruments/Florian/'+ band_file_V)
    wl_V = filt2[:,0]
    transmission_V = filt2[:,1]

    filt2 = np.genfromtxt('../sugar_analysis_data/data/Instruments/Florian/'+ band_file_R)
    wl_R = filt2[:,0]
    transmission_R = filt2[:,1]


            
    band_U = sncosmo.Bandpass(wl_U,transmission_U,wave_unit=u.AA,name='fU')    
    sncosmo.registry.register(band_U)
    band_B = sncosmo.Bandpass(wl_B,transmission_B,wave_unit=u.AA,name='fB')    
    sncosmo.registry.register(band_B)
    band_V = sncosmo.Bandpass(wl_V,transmission_V,wave_unit=u.AA,name='fV')    
    sncosmo.registry.register(band_V)
    band_R = sncosmo.Bandpass(wl_R,transmission_R,wave_unit=u.AA,name='fR')    
    sncosmo.registry.register(band_R)

    band_file_U = 'fU_0.1.dat'
    band_file_B = 'fB_0.1.dat'
    band_file_V = 'fV_0.1.dat'
    band_file_R = 'fR_0.1.dat'
    band_file_I = 'fI_0.1.dat'
       
        
    filt2 = np.genfromtxt('../sugar_analysis_data/data/Instruments/Florian/'+ band_file_U)
    wl_U = filt2[:,0]
    transmission_U = filt2[:,1]

    filt2 = np.genfromtxt('../sugar_analysis_data/data/Instruments/Florian/'+ band_file_B)
    wl_B = filt2[:,0]
    transmission_B = filt2[:,1]

    filt2 = np.genfromtxt('../sugar_analysis_data/data/Instruments/Florian/'+ band_file_V)
    wl_V = filt2[:,0]
    transmission_V = filt2[:,1]

    filt2 = np.genfromtxt('../sugar_analysis_data/data/Instruments/Florian/'+ band_file_R)
    wl_R = filt2[:,0]
    transmission_R = filt2[:,1]

    filt2 = np.genfromtxt('../sugar_analysis_data/data/Instruments/Florian/'+ band_file_I)
    wl_I = filt2[:,0]
    transmission_I = filt2[:,1]
            
    band_U = sncosmo.Bandpass(wl_U,transmission_U,wave_unit=u.AA,name='fU_0.1')    
    sncosmo.registry.register(band_U)
    band_B = sncosmo.Bandpass(wl_B,transmission_B,wave_unit=u.AA,name='fB_0.1')    
    sncosmo.registry.register(band_B)
    band_V = sncosmo.Bandpass(wl_V,transmission_V,wave_unit=u.AA,name='fV_0.1')    
    sncosmo.registry.register(band_V)
    band_R = sncosmo.Bandpass(wl_R,transmission_R,wave_unit=u.AA,name='fR_0.1')    
    sncosmo.registry.register(band_R)    
    band_I = sncosmo.Bandpass(wl_R,transmission_I,wave_unit=u.AA,name='fI_0.1')    
    sncosmo.registry.register(band_I)  
    
    band_file_U = 'fU_1.dat'
    band_file_B = 'fB_1.dat'
    band_file_V = 'fV_1.dat'
    band_file_R = 'fR_1.dat'
    band_file_I = 'fI_1.dat'
       
        
    filt2 = np.genfromtxt('../sugar_analysis_data/data/Instruments/Florian/'+ band_file_U)
    wl_U = filt2[:,0]
    transmission_U = filt2[:,1]

    filt2 = np.genfromtxt('../sugar_analysis_data/data/Instruments/Florian/'+ band_file_B)
    wl_B = filt2[:,0]
    transmission_B = filt2[:,1]

    filt2 = np.genfromtxt('../sugar_analysis_data/data/Instruments/Florian/'+ band_file_V)
    wl_V = filt2[:,0]
    transmission_V = filt2[:,1]

    filt2 = np.genfromtxt('../sugar_analysis_data/data/Instruments/Florian/'+ band_file_R)
    wl_R = filt2[:,0]
    transmission_R = filt2[:,1]

    filt2 = np.genfromtxt('../sugar_analysis_data/data/Instruments/Florian/'+ band_file_I)
    wl_I = filt2[:,0]
    transmission_I = filt2[:,1]
            
    band_U = sncosmo.Bandpass(wl_U,transmission_U,wave_unit=u.AA,name='fU_1')    
    sncosmo.registry.register(band_U)
    band_B = sncosmo.Bandpass(wl_B,transmission_B,wave_unit=u.AA,name='fB_1')    
    sncosmo.registry.register(band_B)
    band_V = sncosmo.Bandpass(wl_V,transmission_V,wave_unit=u.AA,name='fV_1')    
    sncosmo.registry.register(band_V)
    band_R = sncosmo.Bandpass(wl_R,transmission_R,wave_unit=u.AA,name='fR_1')    
    sncosmo.registry.register(band_R)    
    band_I = sncosmo.Bandpass(wl_R,transmission_I,wave_unit=u.AA,name='fI_1')    
    sncosmo.registry.register(band_I)  
    
    band_file_U = 'fU_1.5.dat'
    band_file_B = 'fB_1.5.dat'
    band_file_V = 'fV_1.5.dat'
    band_file_R = 'fR_1.5.dat'
    band_file_I = 'fI_1.5.dat'
        
    filt2 = np.genfromtxt('../sugar_analysis_data/data/Instruments/Florian/'+ band_file_U)
    wl_U = filt2[:,0]
    transmission_U = filt2[:,1]

    filt2 = np.genfromtxt('../sugar_analysis_data/data/Instruments/Florian/'+ band_file_B)
    wl_B = filt2[:,0]
    transmission_B = filt2[:,1]

    filt2 = np.genfromtxt('../sugar_analysis_data/data/Instruments/Florian/'+ band_file_V)
    wl_V = filt2[:,0]
    transmission_V = filt2[:,1]

    filt2 = np.genfromtxt('../sugar_analysis_data/data/Instruments/Florian/'+ band_file_R)
    wl_R = filt2[:,0]
    transmission_R = filt2[:,1]

    filt2 = np.genfromtxt('../sugar_analysis_data/data/Instruments/Florian/'+ band_file_I)
    wl_I = filt2[:,0]
    transmission_I = filt2[:,1]
            

    band_U = sncosmo.Bandpass(wl_U,transmission_U,wave_unit=u.AA,name='fU_1.5')    
    sncosmo.registry.register(band_U)
    band_B = sncosmo.Bandpass(wl_B,transmission_B,wave_unit=u.AA,name='fB_1.5')    
    sncosmo.registry.register(band_B)
    band_V = sncosmo.Bandpass(wl_V,transmission_V,wave_unit=u.AA,name='fV_1.5')    
    sncosmo.registry.register(band_V)
    band_R = sncosmo.Bandpass(wl_R,transmission_R,wave_unit=u.AA,name='fR_1.5')    
    sncosmo.registry.register(band_R)    
    band_I = sncosmo.Bandpass(wl_R,transmission_I,wave_unit=u.AA,name='fI_1.5')    
    sncosmo.registry.register(band_I)
    
    band_file_U = 'fU_5.dat'
    band_file_B = 'fB_5.dat'
    band_file_V = 'fV_5.dat'
    band_file_R = 'fR_5.dat'
    band_file_I = 'fI_5.dat'
       
        
    filt2 = np.genfromtxt('../sugar_analysis_data/data/Instruments/Florian/'+ band_file_U)
    wl_U = filt2[:,0]
    transmission_U = filt2[:,1]

    filt2 = np.genfromtxt('../sugar_analysis_data/data/Instruments/Florian/'+ band_file_B)
    wl_B = filt2[:,0]
    transmission_B = filt2[:,1]

    filt2 = np.genfromtxt('../sugar_analysis_data/data/Instruments/Florian/'+ band_file_V)
    wl_V = filt2[:,0]
    transmission_V = filt2[:,1]

    filt2 = np.genfromtxt('../sugar_analysis_data/data/Instruments/Florian/'+ band_file_R)
    wl_R = filt2[:,0]
    transmission_R = filt2[:,1]

    filt2 = np.genfromtxt('../sugar_analysis_data/data/Instruments/Florian/'+ band_file_I)
    wl_I = filt2[:,0]
    transmission_I = filt2[:,1]
            
    band_U = sncosmo.Bandpass(wl_U,transmission_U,wave_unit=u.AA,name='fU_5')    
    sncosmo.registry.register(band_U)
    band_B = sncosmo.Bandpass(wl_B,transmission_B,wave_unit=u.AA,name='fB_5')    
    sncosmo.registry.register(band_B)
    band_V = sncosmo.Bandpass(wl_V,transmission_V,wave_unit=u.AA,name='fV_5')    
    sncosmo.registry.register(band_V)
    band_R = sncosmo.Bandpass(wl_R,transmission_R,wave_unit=u.AA,name='fR_5')    
    sncosmo.registry.register(band_R)     
    band_I = sncosmo.Bandpass(wl_R,transmission_I,wave_unit=u.AA,name='fI_5')    
    sncosmo.registry.register(band_I)                

    band_file_U = 'fU_5.5.dat'
    band_file_B = 'fB_5.5.dat'
    band_file_V = 'fV_5.5.dat'
    band_file_R = 'fR_5.5.dat'
    band_file_I = 'fI_5.5.dat'
       
        
    filt2 = np.genfromtxt('../sugar_analysis_data/data/Instruments/Florian/'+ band_file_U)
    wl_U = filt2[:,0]
    transmission_U = filt2[:,1]

    filt2 = np.genfromtxt('../sugar_analysis_data/data/Instruments/Florian/'+ band_file_B)
    wl_B = filt2[:,0]
    transmission_B = filt2[:,1]

    filt2 = np.genfromtxt('../sugar_analysis_data/data/Instruments/Florian/'+ band_file_V)
    wl_V = filt2[:,0]
    transmission_V = filt2[:,1]

    filt2 = np.genfromtxt('../sugar_analysis_data/data/Instruments/Florian/'+ band_file_R)
    wl_R = filt2[:,0]
    transmission_R = filt2[:,1]

    filt2 = np.genfromtxt('../sugar_analysis_data/data/Instruments/Florian/'+ band_file_I)
    wl_I = filt2[:,0]
    transmission_I = filt2[:,1]
            
    band_U = sncosmo.Bandpass(wl_U,transmission_U,wave_unit=u.AA,name='fU_5.5')    
    sncosmo.registry.register(band_U)
    band_B = sncosmo.Bandpass(wl_B,transmission_B,wave_unit=u.AA,name='fB_5.5')    
    sncosmo.registry.register(band_B)
    band_V = sncosmo.Bandpass(wl_V,transmission_V,wave_unit=u.AA,name='fV_5.5')    
    sncosmo.registry.register(band_V)
    band_R = sncosmo.Bandpass(wl_R,transmission_R,wave_unit=u.AA,name='fR_5.5')    
    sncosmo.registry.register(band_R)     
    band_I = sncosmo.Bandpass(wl_R,transmission_I,wave_unit=u.AA,name='fI_5.5')    
    sncosmo.registry.register(band_I) 

    band_file_U = 'fU_2.5.dat'
    band_file_B = 'fB_2.5.dat'
    band_file_V = 'fV_2.5.dat'
    band_file_R = 'fR_2.5.dat'
    band_file_I = 'fI_2.5.dat'
       
        
    filt2 = np.genfromtxt('../sugar_analysis_data/data/Instruments/Florian/'+ band_file_U)
    wl_U = filt2[:,0]
    transmission_U = filt2[:,1]

    filt2 = np.genfromtxt('../sugar_analysis_data/data/Instruments/Florian/'+ band_file_B)
    wl_B = filt2[:,0]
    transmission_B = filt2[:,1]

    filt2 = np.genfromtxt('../sugar_analysis_data/data/Instruments/Florian/'+ band_file_V)
    wl_V = filt2[:,0]
    transmission_V = filt2[:,1]

    filt2 = np.genfromtxt('../sugar_analysis_data/data/Instruments/Florian/'+ band_file_R)
    wl_R = filt2[:,0]
    transmission_R = filt2[:,1]

    filt2 = np.genfromtxt('../sugar_analysis_data/data/Instruments/Florian/'+ band_file_I)
    wl_I = filt2[:,0]
    transmission_I = filt2[:,1]
            
    band_U = sncosmo.Bandpass(wl_U,transmission_U,wave_unit=u.AA,name='fU_2.5')    
    sncosmo.registry.register(band_U)
    band_B = sncosmo.Bandpass(wl_B,transmission_B,wave_unit=u.AA,name='fB_2.5')    
    sncosmo.registry.register(band_B)
    band_V = sncosmo.Bandpass(wl_V,transmission_V,wave_unit=u.AA,name='fV_2.5')    
    sncosmo.registry.register(band_V)
    band_R = sncosmo.Bandpass(wl_R,transmission_R,wave_unit=u.AA,name='fR_2.5')    
    sncosmo.registry.register(band_R)     
    band_I = sncosmo.Bandpass(wl_R,transmission_I,wave_unit=u.AA,name='fI_2.5')    
    sncosmo.registry.register(band_I)
    
    band_file_U = 'fU_20.dat'
    band_file_B = 'fB_20.dat'
    band_file_V = 'fV_20.dat'
    band_file_R = 'fR_20.dat'
    band_file_I = 'fI_20.dat'
       
        
    filt2 = np.genfromtxt('../sugar_analysis_data/data/Instruments/Florian/'+ band_file_U)
    wl_U = filt2[:,0]
    transmission_U = filt2[:,1]

    filt2 = np.genfromtxt('../sugar_analysis_data/data/Instruments/Florian/'+ band_file_B)
    wl_B = filt2[:,0]
    transmission_B = filt2[:,1]

    filt2 = np.genfromtxt('../sugar_analysis_data/data/Instruments/Florian/'+ band_file_V)
    wl_V = filt2[:,0]
    transmission_V = filt2[:,1]

    filt2 = np.genfromtxt('../sugar_analysis_data/data/Instruments/Florian/'+ band_file_R)
    wl_R = filt2[:,0]
    transmission_R = filt2[:,1]

    filt2 = np.genfromtxt('../sugar_analysis_data/data/Instruments/Florian/'+ band_file_I)
    wl_I = filt2[:,0]
    transmission_I = filt2[:,1]
            
    band_U = sncosmo.Bandpass(wl_U,transmission_U,wave_unit=u.AA,name='fU_20')    
    sncosmo.registry.register(band_U)
    band_B = sncosmo.Bandpass(wl_B,transmission_B,wave_unit=u.AA,name='fB_20')    
    sncosmo.registry.register(band_B)
    band_V = sncosmo.Bandpass(wl_V,transmission_V,wave_unit=u.AA,name='fV_20')    
    sncosmo.registry.register(band_V)
    band_R = sncosmo.Bandpass(wl_R,transmission_R,wave_unit=u.AA,name='fR_20')    
    sncosmo.registry.register(band_R)     
    band_I = sncosmo.Bandpass(wl_I,transmission_I,wave_unit=u.AA,name='fI_20')    
    sncosmo.registry.register(band_I)        
 
    band_file_U = 'fU_2.dat'
    band_file_B = 'fB_2.dat'
    band_file_V = 'fV_2.dat'
    band_file_R = 'fR_2.dat'
    band_file_I = 'fI_2.dat'
       
        
    filt2 = np.genfromtxt('../sugar_analysis_data/data/Instruments/Florian/'+ band_file_U)
    wl_U = filt2[:,0]
    transmission_U = filt2[:,1]

    filt2 = np.genfromtxt('../sugar_analysis_data/data/Instruments/Florian/'+ band_file_B)
    wl_B = filt2[:,0]
    transmission_B = filt2[:,1]

    filt2 = np.genfromtxt('../sugar_analysis_data/data/Instruments/Florian/'+ band_file_V)
    wl_V = filt2[:,0]
    transmission_V = filt2[:,1]

    filt2 = np.genfromtxt('../sugar_analysis_data/data/Instruments/Florian/'+ band_file_R)
    wl_R = filt2[:,0]
    transmission_R = filt2[:,1]

    filt2 = np.genfromtxt('../sugar_analysis_data/data/Instruments/Florian/'+ band_file_I)
    wl_I = filt2[:,0]
    transmission_I = filt2[:,1]
            
    band_U = sncosmo.Bandpass(wl_U,transmission_U,wave_unit=u.AA,name='fU_2')    
    sncosmo.registry.register(band_U)
    band_B = sncosmo.Bandpass(wl_B,transmission_B,wave_unit=u.AA,name='fB_2')    
    sncosmo.registry.register(band_B)
    band_V = sncosmo.Bandpass(wl_V,transmission_V,wave_unit=u.AA,name='fV_2')    
    sncosmo.registry.register(band_V)
    band_R = sncosmo.Bandpass(wl_R,transmission_R,wave_unit=u.AA,name='fR_2')    
    sncosmo.registry.register(band_R)     
    band_I = sncosmo.Bandpass(wl_I,transmission_I,wave_unit=u.AA,name='fI_2')    
    sncosmo.registry.register(band_I) 
    
    band_file_U = 'fU_100.dat'
    band_file_B = 'fB_100.dat'
    band_file_V = 'fV_100.dat'
    band_file_R = 'fR_100.dat'
    band_file_I = 'fI_100.dat'
       
        
    filt2 = np.genfromtxt('../sugar_analysis_data/data/Instruments/Florian/'+ band_file_U)
    wl_U = filt2[:,0]
    transmission_U = filt2[:,1]

    filt2 = np.genfromtxt('../sugar_analysis_data/data/Instruments/Florian/'+ band_file_B)
    wl_B = filt2[:,0]
    transmission_B = filt2[:,1]

    filt2 = np.genfromtxt('../sugar_analysis_data/data/Instruments/Florian/'+ band_file_V)
    wl_V = filt2[:,0]
    transmission_V = filt2[:,1]

    filt2 = np.genfromtxt('../sugar_analysis_data/data/Instruments/Florian/'+ band_file_R)
    wl_R = filt2[:,0]
    transmission_R = filt2[:,1]

    filt2 = np.genfromtxt('../sugar_analysis_data/data/Instruments/Florian/'+ band_file_I)
    wl_I = filt2[:,0]
    transmission_I = filt2[:,1]
            
    band_U = sncosmo.Bandpass(wl_U,transmission_U,wave_unit=u.AA,name='fU_100')    
    sncosmo.registry.register(band_U)
    band_B = sncosmo.Bandpass(wl_B,transmission_B,wave_unit=u.AA,name='fB_100')    
    sncosmo.registry.register(band_B)
    band_V = sncosmo.Bandpass(wl_V,transmission_V,wave_unit=u.AA,name='fV_100')    
    sncosmo.registry.register(band_V)
    band_R = sncosmo.Bandpass(wl_R,transmission_R,wave_unit=u.AA,name='fR_100')    
    sncosmo.registry.register(band_R)    
    band_I = sncosmo.Bandpass(wl_R,transmission_I,wave_unit=u.AA,name='fI_100')    
    sncosmo.registry.register(band_I)  

    band_file_U = 'fU_10.dat'
    band_file_B = 'fB_10.dat'
    band_file_V = 'fV_10.dat'
    band_file_R = 'fR_10.dat'
    band_file_I = 'fI_10.dat'
       
        
    filt2 = np.genfromtxt('../sugar_analysis_data/data/Instruments/Florian/'+ band_file_U)
    wl_U = filt2[:,0]
    transmission_U = filt2[:,1]

    filt2 = np.genfromtxt('../sugar_analysis_data/data/Instruments/Florian/'+ band_file_B)
    wl_B = filt2[:,0]
    transmission_B = filt2[:,1]

    filt2 = np.genfromtxt('../sugar_analysis_data/data/Instruments/Florian/'+ band_file_V)
    wl_V = filt2[:,0]
    transmission_V = filt2[:,1]

    filt2 = np.genfromtxt('../sugar_analysis_data/data/Instruments/Florian/'+ band_file_R)
    wl_R = filt2[:,0]
    transmission_R = filt2[:,1]

    filt2 = np.genfromtxt('../sugar_analysis_data/data/Instruments/Florian/'+ band_file_I)
    wl_I = filt2[:,0]
    transmission_I = filt2[:,1]
            
    band_U = sncosmo.Bandpass(wl_U,transmission_U,wave_unit=u.AA,name='fU_10')    
    sncosmo.registry.register(band_U)
    band_B = sncosmo.Bandpass(wl_B,transmission_B,wave_unit=u.AA,name='fB_10')    
    sncosmo.registry.register(band_B)
    band_V = sncosmo.Bandpass(wl_V,transmission_V,wave_unit=u.AA,name='fV_10')    
    sncosmo.registry.register(band_V)
    band_R = sncosmo.Bandpass(wl_R,transmission_R,wave_unit=u.AA,name='fR_10')    
    sncosmo.registry.register(band_R)    
    band_I = sncosmo.Bandpass(wl_R,transmission_I,wave_unit=u.AA,name='fI_10')    
    sncosmo.registry.register(band_I)  

    band_file_U = 'fU_3.dat'
    band_file_B = 'fB_3.dat'
    band_file_V = 'fV_3.dat'
    band_file_R = 'fR_3.dat'
    band_file_I = 'fI_3.dat'
       
        
    filt2 = np.genfromtxt('../sugar_analysis_data/data/Instruments/Florian/'+ band_file_U)
    wl_U = filt2[:,0]
    transmission_U = filt2[:,1]

    filt2 = np.genfromtxt('../sugar_analysis_data/data/Instruments/Florian/'+ band_file_B)
    wl_B = filt2[:,0]
    transmission_B = filt2[:,1]

    filt2 = np.genfromtxt('../sugar_analysis_data/data/Instruments/Florian/'+ band_file_V)
    wl_V = filt2[:,0]
    transmission_V = filt2[:,1]

    filt2 = np.genfromtxt('../sugar_analysis_data/data/Instruments/Florian/'+ band_file_R)
    wl_R = filt2[:,0]
    transmission_R = filt2[:,1]

    filt2 = np.genfromtxt('../sugar_analysis_data/data/Instruments/Florian/'+ band_file_I)
    wl_I = filt2[:,0]
    transmission_I = filt2[:,1]
            
    band_U = sncosmo.Bandpass(wl_U,transmission_U,wave_unit=u.AA,name='fU_3')    
    sncosmo.registry.register(band_U)
    band_B = sncosmo.Bandpass(wl_B,transmission_B,wave_unit=u.AA,name='fB_3')    
    sncosmo.registry.register(band_B)
    band_V = sncosmo.Bandpass(wl_V,transmission_V,wave_unit=u.AA,name='fV_3')    
    sncosmo.registry.register(band_V)
    band_R = sncosmo.Bandpass(wl_R,transmission_R,wave_unit=u.AA,name='fR_3')    
    sncosmo.registry.register(band_R)    
    band_I = sncosmo.Bandpass(wl_R,transmission_I,wave_unit=u.AA,name='fI_3')    
    sncosmo.registry.register(band_I)  
    
    band_file_U = 'fU_500.dat'
    band_file_B = 'fB_500.dat'
    band_file_V = 'fV_500.dat'
    band_file_R = 'fR_500.dat'
    band_file_I = 'fI_500.dat'
       
        
    filt2 = np.genfromtxt('../sugar_analysis_data/data/Instruments/Florian/'+ band_file_U)
    wl_U = filt2[:,0]
    transmission_U = filt2[:,1]

    filt2 = np.genfromtxt('../sugar_analysis_data/data/Instruments/Florian/'+ band_file_B)
    wl_B = filt2[:,0]
    transmission_B = filt2[:,1]

    filt2 = np.genfromtxt('../sugar_analysis_data/data/Instruments/Florian/'+ band_file_V)
    wl_V = filt2[:,0]
    transmission_V = filt2[:,1]

    filt2 = np.genfromtxt('../sugar_analysis_data/data/Instruments/Florian/'+ band_file_R)
    wl_R = filt2[:,0]
    transmission_R = filt2[:,1]

    filt2 = np.genfromtxt('../sugar_analysis_data/data/Instruments/Florian/'+ band_file_I)
    wl_I = filt2[:,0]
    transmission_I = filt2[:,1]
            
    band_U = sncosmo.Bandpass(wl_U,transmission_U,wave_unit=u.AA,name='fU_500')    
    sncosmo.registry.register(band_U)
    band_B = sncosmo.Bandpass(wl_B,transmission_B,wave_unit=u.AA,name='fB_500')    
    sncosmo.registry.register(band_B)
    band_V = sncosmo.Bandpass(wl_V,transmission_V,wave_unit=u.AA,name='fV_500')    
    sncosmo.registry.register(band_V)
    band_R = sncosmo.Bandpass(wl_R,transmission_R,wave_unit=u.AA,name='fR_500')    
    sncosmo.registry.register(band_R)    
    band_I = sncosmo.Bandpass(wl_R,transmission_I,wave_unit=u.AA,name='fI_500')    
    sncosmo.registry.register(band_I)

    band_file_U = 'fU_40.dat'
    band_file_B = 'fB_40.dat'
    band_file_V = 'fV_40.dat'
    band_file_R = 'fR_40.dat'
    band_file_I = 'fI_40.dat'
       
        
    filt2 = np.genfromtxt('../sugar_analysis_data/data/Instruments/Florian/'+ band_file_U)
    wl_U = filt2[:,0]
    transmission_U = filt2[:,1]

    filt2 = np.genfromtxt('../sugar_analysis_data/data/Instruments/Florian/'+ band_file_B)
    wl_B = filt2[:,0]
    transmission_B = filt2[:,1]

    filt2 = np.genfromtxt('../sugar_analysis_data/data/Instruments/Florian/'+ band_file_V)
    wl_V = filt2[:,0]
    transmission_V = filt2[:,1]

    filt2 = np.genfromtxt('../sugar_analysis_data/data/Instruments/Florian/'+ band_file_R)
    wl_R = filt2[:,0]
    transmission_R = filt2[:,1]

    filt2 = np.genfromtxt('../sugar_analysis_data/data/Instruments/Florian/'+ band_file_I)
    wl_I = filt2[:,0]
    transmission_I = filt2[:,1]
            
    band_U = sncosmo.Bandpass(wl_U,transmission_U,wave_unit=u.AA,name='fU_40')    
    sncosmo.registry.register(band_U)
    band_B = sncosmo.Bandpass(wl_B,transmission_B,wave_unit=u.AA,name='fB_40')    
    sncosmo.registry.register(band_B)
    band_V = sncosmo.Bandpass(wl_V,transmission_V,wave_unit=u.AA,name='fV_40')    
    sncosmo.registry.register(band_V)
    band_R = sncosmo.Bandpass(wl_R,transmission_R,wave_unit=u.AA,name='fR_40')    
    sncosmo.registry.register(band_R)     
    band_I = sncosmo.Bandpass(wl_R,transmission_I,wave_unit=u.AA,name='fI_40')    
    sncosmo.registry.register(band_I)    
        
def load_spectral_magsys_fits2(relpath, name=None, version=None):
#    fits = pyfits.open(relpath)
#    fit = fits[1]
#    wl = np.zeros(len(fit.data))
#    flux = np.zeros(len(fit.data))
#    for i in range(len(fit.data)):
#         wl[i] = fit.data[i][0]
#         flux[i] = fit.data[i][1]
#    refspectrum = sncosmo.spectrum.Spectrum(wl, flux,
#                           unit=(u.erg / u.s / u.cm**2 / u.AA), wave_unit=u.AA)
    
    data = np.genfromtxt(relpath) 
    dispersion = data[:,0]
    flux_density = data[:,1]
    refspectrum = sncosmo.spectrum.Spectrum(dispersion, flux_density,
                           unit=(u.erg / u.s / u.cm**2 / u.AA), wave_unit=u.AA)
    return sncosmo.magsystems.SpectralMagSystem(refspectrum, name=name)
    

def mag_sys_SNF():
    """
    define magnitude systeme for snf
    """
    
#    sncosmo.registry.register_loader(sncosmo.MagSystem, 'vega_snf_0', load_spectral_magsys_fits2,
#                             args=['../sugar_analysis_data/Vega.fits'],
#                             meta={'description': 'use vega spectrum that come from snf data'}) 
#    bands_snf ={'USNF': ('vega_snf_0', 0.),
#	'BSNF': ('vega_snf_0', 0.),
#	'VSNF': ('vega_snf_0', 0.),
#	'RSNF': ('vega_snf_0', 0.),
#	'ISNF': ('vega_snf_0', 0)}
    
    sncosmo.registry.register_loader(sncosmo.MagSystem, 'vega_snf_0', load_spectral_magsys_fits2,
                             args=['../sugar_analysis_data/data/MagSys/bd_17d4708_stisnic_002.ascii'],
                             meta={'description': 'use bd_17d4708 spectrum that come from snfit'})        

    bands_snf ={'USNf': ('vega_snf_0', 9.787),
	'BSNf': ('vega_snf_0', 9.791),
	'VSNf': ('vega_snf_0', 9.353),
	'RSNf': ('vega_snf_0', 9.011),
	'ISNf': ('vega_snf_0', 8.768),
    'fU': ('vega_snf_0', 9.787),
	'fB': ('vega_snf_0', 9.791),
	'fV': ('vega_snf_0', 9.353),
	'fR': ('vega_snf_0', 9.011),
    'fU_0.1': ('vega_snf_0', 9.787),
	'fB_0.1': ('vega_snf_0', 9.791),
	'fV_0.1': ('vega_snf_0', 9.353),
	'fR_0.1': ('vega_snf_0', 9.011),
    'fI_0.1': ('vega_snf_0', 8.768),
    'fU_1': ('vega_snf_0', 9.787),
	'fB_1': ('vega_snf_0', 9.791),
	'fV_1': ('vega_snf_0', 9.353),
	'fR_1': ('vega_snf_0', 9.011),
    'fI_1': ('vega_snf_0', 8.768),
    'fU_1.5': ('vega_snf_0', 9.787),
    	'fB_1.5': ('vega_snf_0', 9.791),
	'fV_1.5': ('vega_snf_0', 9.353),
	'fR_1.5': ('vega_snf_0', 9.011),
    	'fI_1.5': ('vega_snf_0', 8.768),
	'fU_2': ('vega_snf_0', 9.787),
	'fB_2': ('vega_snf_0', 9.791),
	'fV_2': ('vega_snf_0', 9.353),
	'fR_2': ('vega_snf_0', 9.011),
	'fI_2': ('vega_snf_0', 8.768),
	'fU_3': ('vega_snf_0', 9.787),
	'fB_3': ('vega_snf_0', 9.791),
	'fV_3': ('vega_snf_0', 9.353),
	'fR_3': ('vega_snf_0', 9.011),
	'fI_3': ('vega_snf_0', 8.768),
	'fU_2.5': ('vega_snf_0', 9.787),
	'fB_2.5': ('vega_snf_0', 9.791),
	'fV_2.5': ('vega_snf_0', 9.353),
	'fR_2.5': ('vega_snf_0', 9.011),
	'fI_2.5': ('vega_snf_0', 8.768),  
    'fU_5': ('vega_snf_0', 9.787),
	'fB_5': ('vega_snf_0', 9.791),
	'fV_5': ('vega_snf_0', 9.353),
	'fR_5': ('vega_snf_0', 9.011),
    'fI_5': ('vega_snf_0', 8.768),
	'fU_5.5': ('vega_snf_0', 9.787),
	'fB_5.5': ('vega_snf_0', 9.791),
	'fV_5.5': ('vega_snf_0', 9.353),
	'fR_5.5': ('vega_snf_0', 9.011),
	'fI_5.5': ('vega_snf_0', 8.768),    
    'fU_20': ('vega_snf_0', 9.787),
	'fB_20': ('vega_snf_0', 9.791),
	'fV_20': ('vega_snf_0', 9.353),
	'fR_20': ('vega_snf_0', 9.011),
    'fI_20': ('vega_snf_0', 8.768),
    'fU_10': ('vega_snf_0', 9.787),
	'fB_10': ('vega_snf_0', 9.791),
	'fV_10': ('vega_snf_0', 9.353),
	'fR_10': ('vega_snf_0', 9.011),
    'fI_10': ('vega_snf_0', 8.768),
    'fU_100': ('vega_snf_0', 9.787),
	'fB_100': ('vega_snf_0', 9.791),
	'fV_100': ('vega_snf_0', 9.353),
	'fR_100': ('vega_snf_0', 9.011),
    'fI_100': ('vega_snf_0', 8.768),
    'fU_500': ('vega_snf_0', 9.787),
	'fB_500': ('vega_snf_0', 9.791),
	'fV_500': ('vega_snf_0', 9.353),
	'fR_500': ('vega_snf_0', 9.011),
    'fI_500': ('vega_snf_0', 8.768),
    'fU_40': ('vega_snf_0', 9.787),
	'fB_40': ('vega_snf_0', 9.791),
	'fV_40': ('vega_snf_0', 9.353),
	'fR_40': ('vega_snf_0', 9.011),
    'fI_40': ('vega_snf_0', 8.768),}
    sncosmo.registry.register(sncosmo.CompositeMagSystem(bands=bands_snf),'vega_snf') 
 
def plot_SNf_filters(wl, transmission_U, transmission_V, transmission_B, transmission_R, transmission_I):
    P.plot(wl, transmission_U)
    P.plot(wl, transmission_B)
    P.plot(wl, transmission_V)
    P.plot(wl, transmission_R)
    P.plot(wl, transmission_I)
    P.ylim(0,2)
    P.show()

    