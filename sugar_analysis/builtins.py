# -*- coding: utf-8 -*-
"""
Created on Fri Nov 24 13:24:34 2017
@author: mondon
Load SNF Filters and magnitude systeme and SUGAR
"""

import sncosmo 
import astropy.units as u
import numpy as np
from sncosmo.models import _SOURCES


def register_SNf_bands_width(width=10, sad_path = '../../'):
    
 
    band_file_U = 'USNf_3300-4102.dat'
    band_file_B = 'BSNf_4102-5100.dat'
    band_file_V = 'VSNf_5200-6289.dat'
    band_file_R = 'RSNf_6289-7607.dat'
    band_file_I = 'ISNf_7607-9200.dat'
       
        
    filt2 = np.genfromtxt(sad_path+'sugar_analysis_data/data/Instruments/SNIFS/'+ band_file_U)
    wl_U = filt2[:,0]
    transmission_U = filt2[:,1]

    filt2 = np.genfromtxt(sad_path+'sugar_analysis_data/data/Instruments/SNIFS/'+ band_file_B)
    wl_B = filt2[:,0]
    transmission_B = filt2[:,1]

    filt2 = np.genfromtxt(sad_path+'sugar_analysis_data/data/Instruments/SNIFS/'+ band_file_V)
    wl_V = filt2[:,0]
    transmission_V = filt2[:,1]

    filt2 = np.genfromtxt(sad_path+'sugar_analysis_data/data/Instruments/SNIFS/'+ band_file_R)
    wl_R = filt2[:,0]
    transmission_R = filt2[:,1]

    filt2 = np.genfromtxt(sad_path+'sugar_analysis_data/data/Instruments/SNIFS/'+ band_file_I)
    wl_I = filt2[:,0]
    transmission_I = filt2[:,1]
    
    
    band_U = sncosmo.Bandpass(wl_U,transmission_U,wave_unit=u.AA,name='USNf')    
    sncosmo.registry.register(band_U, force=True)
    band_B = sncosmo.Bandpass(wl_B,transmission_B,wave_unit=u.AA,name='BSNf')    
    sncosmo.registry.register(band_B, force=True)
    band_V = sncosmo.Bandpass(wl_V,transmission_V,wave_unit=u.AA,name='VSNf')    
    sncosmo.registry.register(band_V, force=True)
    band_R = sncosmo.Bandpass(wl_R,transmission_R,wave_unit=u.AA,name='RSNf')    
    sncosmo.registry.register(band_R, force=True)
    band_I = sncosmo.Bandpass(wl_I,transmission_I,wave_unit=u.AA,name='ISNf')    
    sncosmo.registry.register(band_I, force=True)
    
    

    
    band_file_U = 'fU_'+str(width)+'.dat'
    band_file_B = 'fB_'+str(width)+'.dat'
    band_file_V = 'fV_'+str(width)+'.dat'
    band_file_R = 'fR_'+str(width)+'.dat'
    band_file_I = 'fI_'+str(width)+'.dat'
    band_file_new_U = 'new_fU_'+str(width)+'.dat'
    band_file_new_I = 'new_fI_'+str(width)+'.dat'
        
    filt2 = np.genfromtxt(sad_path+'sugar_analysis_data/data/Instruments/Florian/'+ band_file_U)
    wl_U = filt2[:,0]
    transmission_U = filt2[:,1]

    filt2 = np.genfromtxt(sad_path+'sugar_analysis_data/data/Instruments/Florian/'+ band_file_B)
    wl_B = filt2[:,0]
    transmission_B = filt2[:,1]

    filt2 = np.genfromtxt(sad_path+'sugar_analysis_data/data/Instruments/Florian/'+ band_file_V)
    wl_V = filt2[:,0]
    transmission_V = filt2[:,1]

    filt2 = np.genfromtxt(sad_path+'sugar_analysis_data/data/Instruments/Florian/'+ band_file_R)
    wl_R = filt2[:,0]
    transmission_R = filt2[:,1]

    filt2 = np.genfromtxt(sad_path+'sugar_analysis_data/data/Instruments/Florian/'+ band_file_I)
    wl_I = filt2[:,0]
    transmission_I = filt2[:,1]

    filt2 = np.genfromtxt(sad_path+'sugar_analysis_data/data/Instruments/Florian/'+ band_file_new_U)
    wl_new_U = filt2[:,0]
    transmission_new_U = filt2[:,1]

    filt2 = np.genfromtxt(sad_path+'sugar_analysis_data/data/Instruments/Florian/'+ band_file_new_I)
    wl_new_I = filt2[:,0]
    transmission_new_I = filt2[:,1]
            
    band_U = sncosmo.Bandpass(wl_U,transmission_U,wave_unit=u.AA,name='fU_'+str(width))    
    sncosmo.registry.register(band_U, force=True)
    band_B = sncosmo.Bandpass(wl_B,transmission_B,wave_unit=u.AA,name='fB_'+str(width))    
    sncosmo.registry.register(band_B, force=True)
    band_V = sncosmo.Bandpass(wl_V,transmission_V,wave_unit=u.AA,name='fV_'+str(width))    
    sncosmo.registry.register(band_V, force=True)
    band_R = sncosmo.Bandpass(wl_R,transmission_R,wave_unit=u.AA,name='fR_'+str(width))    
    sncosmo.registry.register(band_R, force=True)    
    band_I = sncosmo.Bandpass(wl_I,transmission_I,wave_unit=u.AA,name='fI_'+str(width))    
    sncosmo.registry.register(band_I, force=True)  
    band_new_U = sncosmo.Bandpass(wl_new_U,transmission_new_U, wave_unit=u.AA,name='new_fU_'+str(width))    
    sncosmo.registry.register(band_new_U, force=True) 
    band_new_I = sncosmo.Bandpass(wl_new_I,transmission_new_I, wave_unit=u.AA,name='new_fI_'+str(width))    
    sncosmo.registry.register(band_new_I, force=True)
    
def load_spectral_magsys_fits2(relpath, name=None, version=None):
    
    data = np.genfromtxt(relpath)
    dispersion = data[:,0]
    flux_density = data[:,1]
#    import pyfits
#    data = pyfits.getdata(relpath)
#    dispersion = data['wavelength']
#    flux_density = data['flux']

    refspectrum = sncosmo.spectrum.Spectrum(dispersion, flux_density,
                           unit=(u.erg / u.s / u.cm**2 / u.AA), wave_unit=u.AA)
    return sncosmo.magsystems.SpectralMagSystem(refspectrum, name=name)
    

def mag_sys_SNF(sad_path='../../'):
    """
    define magnitude systeme for snf
    """
    sncosmo.registry.register_loader(sncosmo.MagSystem, 'vega_snf_0', load_spectral_magsys_fits2,
                             args=[sad_path+'sugar_analysis_data/data/MagSys/bd_17d4708_stisnic_002.ascii'],
                             meta={'description': 'use bd_17d4708 spectrum that come from snfit'},force=True)        

    bands_snf ={#'USNf': ('vega_snf_0', 9.787),
#	'BSNf': ('vega_snf_0', 9.791),
#	'VSNf': ('vega_snf_0', 9.353),
#	'RSNf': ('vega_snf_0', 9.011),
	'BSNf': ('vega_snf_0', 0),
	'VSNf': ('vega_snf_0', 0),
	'RSNf': ('vega_snf_0', 0)}
#	'ISNf': ('vega_snf_0', 8.768),
#    'fU': ('vega_snf_0', 9.787),
#	'fB': ('vega_snf_0', 9.791),
#	'fV': ('vega_snf_0', 9.353),
#	'fR': ('vega_snf_0', 9.011)}

    sncosmo.registry.register(sncosmo.CompositeMagSystem(bands=bands_snf),'vega_snf', force=True) 


def mag_sys_SNF_width(width=10, sad_path='../../'):
    """
    define magnitude systeme for snf
    """  
    sncosmo.registry.register_loader(sncosmo.MagSystem, 'vega_snf_0_'+str(width), load_spectral_magsys_fits2,
                             args=[sad_path+'sugar_analysis_data/data/MagSys/bd_17d4708_stisnic_002.ascii'],
                             meta={'description': 'use bd_17d4708 spectrum that come from snfit'},
                             force=True)      

    bands_snf ={'fU_'+str(width): ('vega_snf_0_'+str(width), 9.787),
    'new_fU_'+str(width): ('vega_snf_0_'+str(width), 9.807),
    'new_fI_'+str(width): ('vega_snf_0_'+str(width),8.786 ),
	'fB_'+str(width): ('vega_snf_0_'+str(width), 9.791),
	'fV_'+str(width): ('vega_snf_0_'+str(width), 9.353),
	'fR_'+str(width): ('vega_snf_0_'+str(width), 9.011),
    'ISNf': ('vega_snf_0_'+str(width), 9.011),
    'fI_'+str(width): ('vega_snf_0_'+str(width), 8.768)}

    sncosmo.registry.register(sncosmo.CompositeMagSystem(bands=bands_snf),'vega_snf_'+str(width), force=True) 


# =============================================================================
# Bandpasses
def builtins_jla_bandpasses(sad_path='../../'):
    jla_meta = {
        'filterset': 'jla',
        'reference': ('Betoule14', 'http://supernovae.in2p3.fr/salt/doku.php?id=instruments'),
        'description': 'Representation of instruments and magnitude systems used in JLA'}

    bands = [('jla_SDSS::u', 'SDSS/u.dat', jla_meta),
             ('jla_SDSS::g', 'SDSS/g.dat', jla_meta),
             ('jla_SDSS::r', 'SDSS/r.dat', jla_meta),
             ('jla_SDSS::i', 'SDSS/i.dat', jla_meta),
             ('jla_SDSS::z', 'SDSS/z.dat', jla_meta),
            	 ('jla_4SHOOTER2::Us', '4Shooter2/Us_4Shooter2.txt', jla_meta),
            	 ('jla_4SHOOTER2::B', '4Shooter2/B_4Shooter2.txt', jla_meta),
            	 ('jla_4SHOOTER2::V', '4Shooter2/V_4Shooter2.txt', jla_meta),
            	 ('jla_4SHOOTER2::R', '4Shooter2/R_4Shooter2.txt', jla_meta),
            	 ('jla_4SHOOTER2::I', '4Shooter2/I_4Shooter2.txt', jla_meta),
            	 ('jla_KEPLERCAM::Us', 'Keplercam/Us_Keplercam.txt', jla_meta),
            	 ('jla_KEPLERCAM::B', 'Keplercam/B_Keplercam.txt', jla_meta),
            	 ('jla_KEPLERCAM::V', 'Keplercam/V_Keplercam.txt', jla_meta),
            	 ('jla_KEPLERCAM::i', 'Keplercam/i_Keplercam.txt', jla_meta),
            	 ('jla_KEPLERCAM::r', 'Keplercam/r_Keplercam.txt', jla_meta),
            	 ('jla_STANDARD::U', 'SNLS3-Landolt-model/sux-shifted.dat', jla_meta),
            	 ('jla_STANDARD::B', 'SNLS3-Landolt-model/sb-shifted.dat', jla_meta),
            	 ('jla_STANDARD::V', 'SNLS3-Landolt-model/sv-shifted.dat', jla_meta),
            	 ('jla_STANDARD::R', 'SNLS3-Landolt-model/sr-shifted.dat', jla_meta),
            	 ('jla_STANDARD::I', 'SNLS3-Landolt-model/si-shifted.dat', jla_meta),
            	 ('jla_ACSWF::F606W', 'ACSWF/F606W.dat', jla_meta),
            	 ('jla_ACSWF::F625W', 'ACSWF/F625W.dat', jla_meta),
            	 ('jla_ACSWF::F775W', 'ACSWF/F775W.dat', jla_meta),
            	 ('jla_ACSWF::F814W', 'ACSWF/F814W.dat', jla_meta),
            	 ('jla_ACSWF::F850LP', 'ACSWF/F850LP.dat', jla_meta),
            	 ('jla_NICMOS2::F110W', 'NICMOS2/F110W.dat', jla_meta),
            	 ('jla_NICMOS2::F160W', 'NICMOS2/F160W.dat', jla_meta),
            	 ('jla_MEGACAMPSF::u', 'Megacam-PSF/u0.list', jla_meta),
            	 ('jla_MEGACAMPSF::g', 'Megacam-PSF/g0.list', jla_meta),
            	 ('jla_MEGACAMPSF::r', 'Megacam-PSF/r0.list', jla_meta),
            	 ('jla_MEGACAMPSF::i', 'Megacam-PSF/i0.list', jla_meta),
            	 ('jla_MEGACAMPSF::z', 'Megacam-PSF/z0.list', jla_meta),
            	 ('jla_SWOPE2::u', 'Swope2/u_texas_WLcorr_atm.txt', jla_meta),
            	 ('jla_SWOPE2::B', 'Swope2/B_texas_WLcorr_atm.txt', jla_meta),
            	 ('jla_SWOPE2::g', 'Swope2/g_texas_WLcorr_atm.txt', jla_meta),
            	 ('jla_SWOPE2::V', 'Swope2/V_LC3014_texas_WLcorr_atm.txt', jla_meta),
            	 ('jla_SWOPE2::V1', 'Swope2/V_LC3009_texas_WLcorr_atm.txt', jla_meta),
            	 ('jla_SWOPE2::V2', 'Swope2/V_LC9844_texax_WLcorr_atm.txt', jla_meta),
            	 ('jla_SWOPE2::r', 'Swope2/r_texas_WLcorr_atm.txt', jla_meta),
            	 ('jla_SWOPE2::i', 'Swope2/i_texas_WLcorr_atm.txt', jla_meta)]

    for name, fname, meta in bands:
        	gen = (r.encode('utf-8') for r in open(sad_path + 'sugar_analysis_data/data/jla_data/Instruments/' +  fname) if not r[0] in ('@', '#'))
        	data = np.genfromtxt(gen)
        	band = sncosmo.Bandpass(data[:,0],data[:,1],wave_unit=u.AA,name=name)
        	sncosmo.registry.register(band, force=True)

# =============================================================================
# MagSystems
def mag_sys_jla(sad_path = '../../'):
    def load_spectral_magsys_fits2(relpath, name=None, version=None):
        data = np.genfromtxt(relpath)
        dispersion = data[:,0]
        flux_density = data[:,1]
        refspectrum = sncosmo.spectrum.Spectrum(dispersion, flux_density,
                               unit=(u.erg / u.s / u.cm**2 / u.AA), wave_unit=u.AA)
        return sncosmo.magsystems.SpectralMagSystem(refspectrum, name=name)

    website = 'http://supernovae.in2p3.fr/salt/doku.php?id=instruments'
    subclass = '`~sncosmo.SpectralMagSystem`'
    
    VEGAHST_desc = 'Vega0.dat'
    AB_jla_desc = 'B12-AB-off.dat'
    VEGA2_desc = 'BD17-jla1.dat'
    VEGA2_mb_desc = 'BD17-snls3.dat'


    for name, fn, desc in [('jla_VEGAHST', 'alpha_lyr_stis_005.ascii', VEGAHST_desc),
    		       ('jla_AB_B12_0', 'ab-spec.dat', AB_jla_desc),
                           ('jla_VEGA2_0', 'bd_17d4708_stisnic_003.ascii', VEGA2_desc),
    		       ('jla_VEGA2_mb_0', 'bd_17d4708_stisnic_002.ascii', VEGA2_mb_desc)]:
        sncosmo.registry.register_loader(sncosmo.MagSystem, name, load_spectral_magsys_fits2,
                                 args=[sad_path + 'sugar_analysis_data/data/jla_data/MagSys/' + fn],
                                 meta={'subclass': subclass, 'url': website,
                                       'description': desc}, force=True)

    # offsets are in the sense (mag_SDSS - mag_AB) = offset
    # -> for example: a source with AB mag = 0. will have SDSS mag = 0.06791
    bands_ab = {'jla_SDSS::u': ('jla_AB_B12_0',  0.06791),
            'jla_SDSS::g': ('jla_AB_B12_0', -0.02028),
            'jla_SDSS::r': ('jla_AB_B12_0', -0.00493),
            'jla_SDSS::i': ('jla_AB_B12_0', -0.01780),
            'jla_SDSS::z': ('jla_AB_B12_0', -0.01015),
            'jla_MEGACAMPSF::g': ('jla_AB_B12_0', 0),
            'jla_MEGACAMPSF::r': ('jla_AB_B12_0', 0),
            'jla_MEGACAMPSF::i': ('jla_AB_B12_0', 0),
            'jla_MEGACAMPSF::z': ('jla_AB_B12_0', 0)}

    sncosmo.registry.register(sncosmo.CompositeMagSystem(bands=bands_ab),'jla_AB_B12', force=True)
    
    bands_BD17 ={'jla_STANDARD::U': ('jla_VEGA2_0', 9.724),
    	'jla_STANDARD::B': ('jla_VEGA2_0', 9.907),
    	'jla_STANDARD::V': ('jla_VEGA2_0', 9.464),
    	'jla_STANDARD::R': ('jla_VEGA2_0', 9.166),
    	'jla_STANDARD::I': ('jla_VEGA2_0', 8.846),
    	'jla_SDSS::u': ('jla_VEGA2_0', 10.56),
    	'jla_SDSS::g': ('jla_VEGA2_0', 9.64),
    	'jla_SDSS::r': ('jla_VEGA2_0', 9.35),
    	'jla_SDSS::i': ('jla_VEGA2_0', 9.25),
    	'jla_SDSS::z': ('jla_VEGA2_0', 9.23),
    	'jla_KEPLERCAM::Us': ('jla_VEGA2_0', 9.724), # U standard (normalement)
    	'jla_KEPLERCAM::B': ('jla_VEGA2_0', 9.8803),
    	'jla_KEPLERCAM::V': ('jla_VEGA2_0', 9.4722),
    	'jla_KEPLERCAM::r': ('jla_VEGA2_0', 9.3524),
    	'jla_KEPLERCAM::i': ('jla_VEGA2_0', 9.2542),
    	'jla_4SHOOTER2::Us': ('jla_VEGA2_0', 9.724), # U standard (normalement)
    	'jla_4SHOOTER2::B': ('jla_VEGA2_0', 9.8744),
    	'jla_4SHOOTER2::V': ('jla_VEGA2_0', 9.4789),
    	'jla_4SHOOTER2::R': ('jla_VEGA2_0', 9.1554),
    	'jla_4SHOOTER2::I': ('jla_VEGA2_0', 8.8506),
    	'jla_SWOPE2::u': ('jla_VEGA2_0',   10.514),
    	'jla_SWOPE2::g': ('jla_VEGA2_0',   9.64406),
    	'jla_SWOPE2::r': ('jla_VEGA2_0',   9.3516),
    	'jla_SWOPE2::i': ('jla_VEGA2_0',   9.25),
    	'jla_SWOPE2::B': ('jla_VEGA2_0',   9.876433),
    	'jla_SWOPE2::V': ('jla_VEGA2_0',   9.476626),
    	'jla_SWOPE2::V1': ('jla_VEGA2_0',  9.471276),
    	'jla_SWOPE2::V2': ('jla_VEGA2_0',  9.477482)}
    
    sncosmo.registry.register(sncosmo.CompositeMagSystem(bands=bands_BD17),'jla_VEGA2', force=True)
    
    bands_BD17_mb ={'jla_STANDARD::U': ('jla_VEGA2_mb_0', 9.724),
    	'jla_STANDARD::B': ('jla_VEGA2_mb_0', 9.907),
    	'jla_STANDARD::V': ('jla_VEGA2_mb_0', 9.464),
    	'jla_STANDARD::R': ('jla_VEGA2_mb_0', 9.166),
    	'jla_STANDARD::I': ('jla_VEGA2_mb_0', 8.846),
    	'jla_SDSS::u': ('jla_VEGA2_mb_0', 10.56),
    	'jla_SDSS::g': ('jla_VEGA2_mb_0', 9.64),
    	'jla_SDSS::r': ('jla_VEGA2_mb_0', 9.35),
    	'jla_SDSS::i': ('jla_VEGA2_mb_0', 9.25),
    	'jla_SDSS::z': ('jla_VEGA2_mb_0', 9.23),
    	'jla_KEPLERCAM::B': ('jla_VEGA2_mb_0', 9.8803),
    	'jla_KEPLERCAM::V': ('jla_VEGA2_mb_0', 9.4722),
    	'jla_KEPLERCAM::r': ('jla_VEGA2_mb_0', 9.3524),
    	'jla_KEPLERCAM::i': ('jla_VEGA2_mb_0', 9.2542),
    	'jla_4SHOOTER2::B': ('jla_VEGA2_mb_0', 9.8744),
    	'jla_4SHOOTER2::V': ('jla_VEGA2_mb_0', 9.4789),
    	'jla_4SHOOTER2::R': ('jla_VEGA2_mb_0', 9.1554),
    	'jla_4SHOOTER2::I': ('jla_VEGA2_mb_0', 8.8506)}
    
    sncosmo.registry.register(sncosmo.CompositeMagSystem(bands=bands_BD17_mb),'jla_VEGA2_mb', force=True)








   
