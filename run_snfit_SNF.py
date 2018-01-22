# -*- coding: utf-8 -*-
"""
Created on Mon Dec 18 14:39:25 2017

@author: mondon
"""

import cPickle
import numpy as np
import os 
import sys
import subprocess

def run_snfit_SNF(filters=['BSNf','VSNf','RSNf'], errorscale=True):

    '''
    creat snfit file for each filter (like the file format input in snfit 2.2.2) 
    and run snfit 2.4.0 for each SNF in meta.pkl 
    '''
    meta = cPickle.load(open('../sugar_analysis_data/META-CABALLO2.pkl'))
    outdir = '../sugar_analysis_data/snfit_file_SNF'
    sn_dico = {}
    if errorscale:
        results = open( '../sugar_analysis_data/results/results_snfit.txt','w')
    else:       
        results = open( '../sugar_analysis_data/results/results_snfit_without_errorscale.txt','w')
    
#    sn_list = ['SNF20080323-009','SNF20080510-005']
    for sn_name in meta.keys():
#    for sn_name in sn_list:
        if meta[sn_name]['idr.subset'] != 'bad' and meta[sn_name]['idr.subset'] != 'auxiliary':
            errorscale_factor = meta[sn_name]['target.errorscale']
            redshift = meta[sn_name]['host.zhelio']
            mwebv = meta[sn_name]['target.mwebv']
            lightfile = open(os.path.join(outdir, 'lightfile'), 'w')
            lightfile.write('Name     %s\n' % sn_name)
            lightfile.write('Redshift %f\n' % redshift)
            lightfile.write('MWEBV    %f\n' % mwebv)
            lightfile.close()
            sn_dico[sn_name] = {}
            
            for f in filters:    
                filename = '%s/maglc2fit_%s.dat' % (outdir, f)
                
                
                with open(filename, 'w') as file_name:
                    file_name.write('@BAND %s\n' % f)
                    file_name.write('@INSTRUMENT SNIFS\n')
                    file_name.write('@MAGSYS VEGA\n')
                    file_name.write('# Date :\n')
                    file_name.write('# Mag :\n')
                    file_name.write('# Magerr :\n') 
                    file_name.write('# end :\n')
                    
                    for spectra_name in meta[sn_name]['spectra'].keys():
                        
                        if f == 'USNf'or f == 'BSNf':
                            try:
                                quality_flag = meta[sn_name]['spectra'][spectra_name]['procB.Quality']
                                
                            except:
                                quality_flag = 1
                                
                        elif f == 'VSNf'or f == 'RSNf' or f == 'ISNf':
                            try:
                                quality_flag = meta[sn_name]['spectra'][spectra_name]['procR.Quality']
                            except:
                                quality_flag = 1
                        else:
                            print >> sys.stderr, \
                                'filter not recognise' 
                            return
                          
                        if quality_flag == 1:
                                mjd = meta[sn_name]['spectra'][spectra_name]['obs.mjd']
                                mag = meta[sn_name]['spectra'][spectra_name]['mag.'+f]
                                sn_dico[sn_name]['spectra'] = spectra_name
                                
                                if errorscale:
                                    err_mag = meta[sn_name]['spectra'][spectra_name]['mag.'+f+'.err']
                                
                                
                                else:  
                                    err_mag = float(meta[sn_name]['spectra'][spectra_name]['mag.'+f+'.err']/errorscale_factor)
       
                                if str(mag) != 'nan':
                                    file_name.write('%f %f %f\n' % (mjd, mag, err_mag))
            cmd = 'snfit '
            for f in filters:
                cmd += '%s/maglc2fit_%s.dat ' % (outdir, f)    
    
            cmd += '-o '+outdir+'/results_salt2.dat'
            try:
                snfit_results = subprocess.call(cmd,shell=True,stderr=subprocess.STDOUT)
            except:
                print >> sys.stderr, \
                    'WARNING: snfit error, please check if your SALTPATH is right' 
            results.write('########################################################## \n')
            results.write('sn_name '+sn_name+'\n')
            res_snfit = open(outdir+'/results_salt2.dat') 
            res_snfit_lines = res_snfit.readlines()
            
            for res_snfit_line in res_snfit_lines:
                results.write(res_snfit_line)
                
            res_snfit.close()
    results.close()
    return sn_dico
            
            