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
    
    def read_spectrum(self, spec_path, sad_path = '../../../'):
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
        self.x = np.linspace(self.start, self.end, self.npts)  # Wavelength
        self.y = self.spec.data.copy()                         # Signal
        self.v = self._fits['VARIANCE'].data.copy()
        self.sad_path = sad_path
        return self.x, self.y ,self.v

    def deredden(self, ebmv, law='OD94', Rv=3.1):
        """Deredden spectrum using E(B-V) and a nextinction law:

        :param float ebmv: E(B-V) value.
        :param int law: Extinction law. Could be CCM89, OD94, FM98 or G08
        :param float Rv: Value for Rv. Default is 3.1"""
        from ToolBox.Astro.Extinction import extinctionFactor

        if hasattr(self, 'zorig'):      # Spectrum has been deredshifted
            raise ValueError, \
                "Dereddening should be done prior to deredshifting."

        # Extinction factor (<1)
        ext = extinctionFactor(self.x, ebmv, Rv=Rv, law=law)
        self.y /= ext
        self.v /= ext**2
            
    def table_time_wave_flux(self):
        """
        """
        meta = cPickle.load(open(self.sad_path+'sugar_analysis_data/SNF-0203-CABALLO2/META.pkl'))
        meta_spectra = {}
        for sn_name in meta.keys():
#        for sn_name in list_SN :
            if meta[sn_name]['idr.subset'] != 'bad' and meta[sn_name]['idr.subset'] != 'auxiliary':
                meta_spectra[sn_name] = {}
                
                zhl = meta[sn_name]['host.zhelio']
                zcmb = meta[sn_name]['host.zcmb']
                daymax = meta[sn_name]['salt2.DayMax']
                self.mwebv = meta[sn_name]['target.mwebv']
                sn_data = meta[sn_name]['spectra']
                meta_spectra[sn_name]['host.zhelio'] = zhl
                meta_spectra[sn_name]['host.zcmb'] = zcmb
                meta_spectra[sn_name]['salt2.DayMax'] = daymax
                meta_spectra[sn_name]['target.mwebv'] = self.mwebv
                meta_spectra[sn_name]['spectra_B'] = {}
                meta_spectra[sn_name]['spectra_R'] = {}
                for t in sn_data.keys():
                    try:
                        quality_flag_B = sn_data[t]['procB.Quality']
                    except:
                        quality_flag_B = 0
                    try:
                        quality_flag_R = sn_data[t]['procR.Quality']
                    except:
                        quality_flag_R = 0
                    if quality_flag_B == 1 and quality_flag_R == 1:
                        for x in ['B','R']:
                            wave, flux, var = self.read_spectrum(self.sad_path+'sugar_analysis_data/SNF-0203-CABALLO2/'+sn_data[t]['idr.spec_'+x])
                            time = sn_data[t]['obs.mjd']
                            meta_spectra[sn_name]['spectra_'+x][str(time)] = {}
                            meta_spectra[sn_name]['spectra_'+x][str(time)]['wave'] = wave
                            meta_spectra[sn_name]['spectra_'+x][str(time)]['flux'] = flux  
                            meta_spectra[sn_name]['spectra_'+x][str(time)]['var'] = var
        File = open('../../../sugar_analysis_data/SNF-0203-CABALLO2/meta_spectra.pkl','w')
        cPickle.dump(meta_spectra, File)
        File.close()
        return meta_spectra
                    
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
