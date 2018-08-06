#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu May 24 13:07:00 2018

@author: florian
"""

import numpy as np
import cPickle as pkl
import sncosmo
import sys
from math import isnan
from matplotlib import pyplot as plt

def mag_to_flux(mag,band,width):
    vega = sncosmo.get_magsystem('vega_snf_'+str(width))
    flux = vega.band_mag_to_flux(mag, band)

    return flux

def distrib():
    width = 10 
    errors_USNf = []
    errors_BSNf = []
    errors_VSNf = []
    errors_RSNf = []
    errors_ISNf = []
    phase = [] 
    filters=['USNf','BSNf','VSNf','RSNf','ISNf']
    meta = pkl.load(open('../sugar_analysis_data/META-CABALLO2.pkl'))
   
    for sn_name in meta.keys():
        if meta[sn_name]['idr.subset'] != 'bad' and meta[sn_name]['idr.subset'] != 'auxiliary':
            daymax = meta[sn_name]['salt2.DayMax']
            
            for t in meta[sn_name]['spectra'].keys():

                phase.append(meta[sn_name]['spectra'][t]['obs.mjd']-daymax)
                for f in filters:
                    if f == 'USNf'or f == 'BSNf':
                        try:
                            quality_flag = meta[sn_name]['spectra'][t]['procB.Quality']
                            
                        except:
                            quality_flag = 1
                            
                    elif f == 'VSNf'or f == 'RSNf' or f == 'ISNf':
                        try:
                            quality_flag = meta[sn_name]['spectra'][t]['procR.Quality']
                        except:
                            quality_flag = 1
                    else:
                        print >> sys.stderr, \
                            'filter not recognise' 
                        return
                        
        #            if isnan(sn_data[t]['mag.'+f]) == False:
                    if quality_flag == 1 and not isnan(meta[sn_name]['spectra'][t]['mag.'+f]):
                
                        if f == 'USNf':
                            fn = 'fU_'+str(width)
                            errors_USNf.append(meta[sn_name]['spectra'][t]['mag.'+f+'.err'] * mag_to_flux(meta[sn_name]['spectra'][t]['mag.'+f],fn,width) / 1.0857362047581294)
                        if f == 'BSNf':
                            fn = 'fB_'+str(width)
                            errors_BSNf.append(meta[sn_name]['spectra'][t]['mag.'+f+'.err'] * mag_to_flux(meta[sn_name]['spectra'][t]['mag.'+f],fn,width) / 1.0857362047581294)
                        if f == 'VSNf':
                            fn = 'fV_'+str(width)
                            errors_VSNf.append(meta[sn_name]['spectra'][t]['mag.'+f+'.err'] * mag_to_flux(meta[sn_name]['spectra'][t]['mag.'+f],fn,width) / 1.0857362047581294)
                        if f == 'RSNf':
                            fn = 'fR_'+str(width)
                            errors_RSNf.append(meta[sn_name]['spectra'][t]['mag.'+f+'.err'] * mag_to_flux(meta[sn_name]['spectra'][t]['mag.'+f],fn,width) / 1.0857362047581294)
                        if f == 'ISNf':
                            fn = 'fI_'+str(width)
                            errors_ISNf.append(meta[sn_name]['spectra'][t]['mag.'+f+'.err'] * mag_to_flux(meta[sn_name]['spectra'][t]['mag.'+f],fn,width) / 1.0857362047581294)
    plt.hist(phase,50)
    plt.show()
    plt.close()
    plt.hist(errors_USNf,50)
    print np.mean(errors_USNf)
    print np.std(errors_USNf)
    plt.show()
    plt.close()
    plt.hist(errors_BSNf,50)
    print np.mean(errors_BSNf)
    print np.std(errors_BSNf)
    plt.show()
    plt.close()
    plt.hist(errors_VSNf,50)
    print np.mean(errors_VSNf)
    print np.std(errors_VSNf)
    plt.show()
    plt.close()
    plt.hist(errors_RSNf,50)
    print np.mean(errors_RSNf)
    print np.std(errors_RSNf)
    plt.show()
    plt.close()
    plt.hist(errors_ISNf,50)
    print np.mean(errors_ISNf)
    print np.std(errors_ISNf)
    plt.show()
    plt.close()