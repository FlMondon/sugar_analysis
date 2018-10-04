#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  1 10:55:09 2018

@author: florian
"""

from sugar_training import extract_idr_data
import numpy as np
import sncosmo
from sugar_analysis import constant as cst
from astropy.table import Table
from scipy.interpolate import InterpolatedUnivariateSpline as Spline1d
from sugar_analysis import builtins
import cPickle as pkl
output_path = '../../'

def extract_spectra():
    
    bsd = extract_idr_data.build_spectral_data(output_path+'sugar_analysis_data/SNF-0203-CABALLO2',
                                               mjd_max=np.inf, day_max=5., redshift_max=0.11)
    bsd.load_spectra()
#    bsd.resampled_spectra(velocity=1500.)
#    bsd.to_ab_mag()
#    bsd.reorder_and_clean()
    bsd.write_pkl(output_path+'sugar_analysis_data/data_input/spectra_snia_all_good_for_florian.pkl')
    return bsd
    
class build_data(object):

    def __init__(self,meta=output_path+'sugar_analysis_data/SNF-0203-CABALLO2/META.pkl'):
        builtins.register_SNf_bands_width(width=10)
        builtins.mag_sys_SNF_width(width=10)
        builtins.register_SUGAR()  
        self.spec = pkl.load(open(output_path+
                                  'sugar_analysis_data/spectra_table/meta_spectra.pkl'))
        self.meta= pkl.load(open(meta))
        self.bsd = extract_spectra()
        
    def integral_to_phot(self, wave, flux, band):
        filt2 = np.genfromtxt(output_path+
                              'sugar_analysis_data/data/Instruments/Florian/'+band+'.dat')
        wlen = filt2[:,0]
        tran = filt2[:,1]
        self.splB = Spline1d(wlen, tran, k=1,ext = 1)    
        #computation of the integral
        dxs = (float(max(wave)-min(wave))/(len(wave)-1))
        inte = np.sum(flux*wave/(cst.CLIGHT*1.0e13*cst.HPLANCK)*self.splB(wave)*dxs)
        return inte
    def build_Astropy_Table(self, sn_name, band_used=['new_fI_10','fB_10','fV_10','fR_10','fU_10']):
        """
        """    
        time = []
        flux = []
        flux_err = []
        band = []
        zp = []
        zpsys = []
        vega = sncosmo.get_magsystem('vega_snf_10')

        for nspec in self.bsd.dico_spectra[sn_name]:
            p = self.bsd.dico_spectra[sn_name][nspec]['days']
            for b in band_used:
                phot_flux = self.integral_to_phot(
                        self.bsd.dico_spectra[sn_name][nspec]['X'],self.bsd.dico_spectra[sn_name][nspec]['Y'],b)  

                flux.append(phot_flux)
                time.append(p)
                band.append(b)
                flux_err.append(None)
                zp.append(2.5*np.log10(vega.zpbandflux(b)))
                zpsys.append('vega_snf_10')
                
        data = Table([time, band, flux, flux_err, zp, zpsys], 
                     names=('time', 'band', 'flux', 'fluxerr', 'zp', 'zpsys'), 
                     meta={'name': 'data'})
        return data


    def build_Astropy_Table_flo(self, sn_name, band_used=['new_fI_10','fB_10','fV_10','fR_10','fU_10']):
        """
        """    
        time = []
        flux = []
        flux_err = []
        band = []
        zp = []
        zpsys = []
        vega = sncosmo.get_magsystem('vega_snf_10')

        for p in self.spec[sn_name]['spectra'].keys():
            for b in band_used:
                phot_flux = self.integral_to_phot(
                        self.spec[sn_name]['spectra'][str(p)]['wave'],
                        self.spec[sn_name]['spectra'][str(p)]['flux'],b)  

                flux.append(phot_flux)
                time.append(p)
                band.append(b)
                flux_err.append(None)
                zp.append(2.5*np.log10(vega.zpbandflux(b)))
                zpsys.append('vega_snf_10')
                
        data = Table([time, band, flux, flux_err, zp, zpsys], 
                     names=('time', 'band', 'flux', 'fluxerr', 'zp', 'zpsys'), 
                     meta={'name': 'data'})
        return data
