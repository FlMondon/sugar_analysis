# -*- coding: utf-8 -*-
"""
Created on Wed Nov 22 17:50:29 2017
@author: mondon

read data from meta.pkl (Caballo2) to creat AstropyTable for the input of the light curve fitting with sncosmo
"""
import numpy as np
from astropy.table import Table
import sncosmo
from math import isnan
import sys


def mag_to_flux(mag,band,width):
    vega = sncosmo.get_magsystem('vega_snf_'+str(width))
    flux = vega.band_mag_to_flux(mag, band)
    return flux
   

        
def read_meta_SNF(meta,sn_name,filters=['BSNf','VSNf','RSNf'],model='sugar',errorscale=True, width=10):
    """
    """

    vega = sncosmo.get_magsystem('vega_snf_'+str(width))
    time = []
    band = []
    flux = []
    fluxerr = []
    zp = []
    zpsys = []
    errorscale_factor = meta[sn_name]['target.errorscale']
    zhl = meta[sn_name]['host.zhelio']
    zcmb = meta[sn_name]['host.zcmb']
    zerr = meta[sn_name]['host.zhelio.err']
    daymax = meta[sn_name]['salt2.DayMax']
    mwebv = meta[sn_name]['target.mwebv']
    sn_data = meta[sn_name]['spectra']
    for t in sn_data.keys():

        for f in filters:
            if f == 'USNf'or f == 'BSNf':
                try:
                    quality_flag = sn_data[t]['procB.Quality']
                    
                except:
                    quality_flag = 1
                    
            elif f == 'VSNf'or f == 'RSNf' or f == 'ISNf':
                try:
                    quality_flag = sn_data[t]['procR.Quality']
                except:
                    quality_flag = 1
            else:
                print >> sys.stderr, \
                    'filter not recognise' 
                return

            if quality_flag == 1 and not isnan(sn_data[t]['mag.'+f]):
#
                time.append(sn_data[t]['obs.mjd'])
                
                if f == 'USNf':
                    fn = 'fU_'+str(width)
                    band.append(fn)
                if f == 'BSNf':
                    fn = 'fB_'+str(width)
                    band.append(fn)
                if f == 'VSNf':
                    fn = 'fV_'+str(width)
                    band.append(fn)
                if f == 'RSNf':
                    fn = 'fR_'+str(width)
                    band.append(fn)
                if f == 'ISNf':
                    #fn = 'new_fI_'+str(width)
                    fn = 'ISNf'
                    band.append(fn) 

                
                flux.append(mag_to_flux(sn_data[t]['mag.'+f],f,width))
                zp.append(2.5*np.log10(vega.zpbandflux(f)))               
                zpsys.append('vega_snf_'+str(width))
                if errorscale:

                    err_mag = sn_data[t]['mag.'+f+'.err']

                else:  
                    err_mag = sn_data[t]['mag.'+f+'.err']/errorscale_factor

                fluxerr.append(err_mag * mag_to_flux(sn_data[t]['mag.'+f],f,width) / 1.0857362047581294)
    data = Table([time, band, flux, fluxerr, zp, zpsys], names=('time', 'band', 'flux', 'fluxerr', 'zp', 'zpsys'), meta={'name': 'data'})
    return data, zhl, zerr, zcmb, mwebv, daymax
    
    

    
    
    
    
    
    
    
    
    
    