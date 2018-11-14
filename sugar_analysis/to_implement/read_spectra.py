#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  9 15:04:14 2018

@author: florian
"""

import pyfits
import cPickle
import numpy as np

class spectrum_table():
    
    def read_spectrum(self, spec_path):
        """
        """
        self._fits = pyfits.open(spec_path)
        self.spec = self._fits[0]
        self._hdr = self.spec.header.copy()       # Spectrum header
        
        self._hdr['CRPIX1'] = self._hdr.get('CRPIX1', 1)  # Make it mandatory
        
        self.npts = self._hdr['NAXIS1']
        self.step = self._hdr['CDELT1']
        self.start= self._hdr['CRVAL1'] - (self._hdr['CRPIX1'] - 1) * self.step
        self.end  = self.start + (self.npts - 1) * self.step
        wave = np.linspace(self.start, self.end, self.npts)  # Wavelength
        flux = self.spec.data.copy()                         # Signal
        var = self._fits['VARIANCE'].data.copy()
        return wave, flux, var

    def table_time_wave_flux(self):
        """
        """
        meta = cPickle.load(open('../sugar_analysis_data/META-CABALLO2.pkl'))
        list_SN = ['SNF20080802-006']
        meta_spectra = {}
        for sn_name in meta.keys():
#        for sn_name in list_SN :
            if meta[sn_name]['idr.subset'] != 'bad' and meta[sn_name]['idr.subset'] != 'auxiliary':
                meta_spectra[sn_name] = {}
                
                zhl = meta[sn_name]['host.zhelio']
                zcmb = meta[sn_name]['host.zcmb']
                daymax = meta[sn_name]['salt2.DayMax']
                mwebv = meta[sn_name]['target.mwebv']
                sn_data = meta[sn_name]['spectra']
                meta_spectra[sn_name]['host.zhelio'] = zhl
                meta_spectra[sn_name]['host.zcmb'] = zcmb
                meta_spectra[sn_name]['salt2.DayMax'] = daymax
                meta_spectra[sn_name]['target.mwebv'] = mwebv
                meta_spectra[sn_name]['spectra'] = {}
                for t in sn_data.keys():
                    try:
                        quality_flag_B = sn_data[t]['procB.Quality']
                    except:
                        quality_flag_B = 1
                    try:
                        quality_flag_R = sn_data[t]['procR.Quality']
                    except:
                        quality_flag_R = 1
                    if quality_flag_B == 1 and quality_flag_R == 1:
                        wave, flux, var = self.read_spectrum('../sugar_analysis_data/spectra_table/'+sn_data[t]['idr.spec_merged'])
                        time = sn_data[t]['obs.mjd']
                        meta_spectra[sn_name]['spectra'][str(time)] = {}
                        meta_spectra[sn_name]['spectra'][str(time)]['wave'] = wave
                        meta_spectra[sn_name]['spectra'][str(time)]['flux'] = flux  
                        meta_spectra[sn_name]['spectra'][str(time)]['var'] = var
                        meta_spectra[sn_name]['spectra'][str(time)]['mag'] = -2.5 *np.log10(flux)  
                        meta_spectra[sn_name]['spectra'][str(time)]['dmag'] = 1.0857362047581294 / flux  * var  # 2.5/log(10) = 1.0857...
        File = open('../sugar_analysis_data/spectra_table/meta_spectra.pkl','w')
        cPickle.dump(meta_spectra, File)
        File.close()
        return meta_spectra
                    
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
