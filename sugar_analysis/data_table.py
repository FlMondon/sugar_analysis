#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  1 10:55:09 2018

@author: florian
"""
import os 
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.table import Table
import numpy as np
import sncosmo
try:
   import cPickle as pkl
except ModuleNotFoundError:
   import pickle as pkl
   
from .constant import CLIGHT, HPLANCK
from .cosmo_tools import zhelio2zcmb 
from .builtins import register_SNf_bands_width, mag_sys_SNF_width,  builtins_jla_bandpasses, mag_sys_jla
#from .load_sugar import register_SUGAR


RAD2DEG = 57.295779513082323            # 180/pi
def rec2pol(x, y, deg=False):
    """Conversion of rectangular *(x,y)* to polar *(r,theta)*
    coordinates"""

    r = np.hypot(x,y)
    t = np.arctan2(y,x)
    if deg:                             # Convert to radians
        t *= RAD2DEG

    return r,t    
def radec2gcs(ra, dec, deg=True):
    """Convert *(ra,dec)* equatorial coordinates (J2000, in degrees if
    *deg*) to Galactic Coordinate System coordinates *(lII,bII)* (in
    degrees if *deg*).

    Sources:

    - http://www.dur.ac.uk/physics.astrolab/py_source/conv.py_source
    - Rotation matrix from
      http://www.astro.rug.nl/software/kapteyn/celestialbackground.html

    .. Note:: This routine is only roughly accurate, probably at the
              arcsec level, and therefore not to be used for
              astrometric purposes. For most accurate conversion, use
              dedicated `kapteyn.celestial.sky2sky` routine.

    >>> radec2gal(123.456, 12.3456)
    (210.82842704243518, 23.787110745502183)
    """
    
    if deg:
        ra  =  ra / RAD2DEG
        dec = dec / RAD2DEG

    rmat = np.array([[-0.054875539396, -0.873437104728, -0.48383499177 ],
                    [ 0.494109453628, -0.444829594298,  0.7469822487  ],
                    [-0.867666135683, -0.198076389613,  0.455983794521]])
    cosd = np.cos(dec)
    v1 = np.array([np.cos(ra)*cosd,
                  np.sin(ra)*cosd,
                  np.sin(dec)])
    v2 = np.dot(rmat, v1)
    x,y,z = v2

    c,l = rec2pol(x,y)
    r,b = rec2pol(c,z)
    
    assert np.allclose(r,1), "Precision error"

    if deg:
        l *= RAD2DEG
        b *= RAD2DEG
    
    return l, b
    
class build_data(object):


    def __init__(self,sad_path = '../../', modeldir='../../sugar_model/'):
        self.sad_path = sad_path
        register_SNf_bands_width(width=10, sad_path=self.sad_path)
        mag_sys_SNF_width(width=10 , sad_path=self.sad_path)
#        register_SUGAR(modeldir=modeldir)  

        self.noTH = True
        pkl_file = os.path.join(self.sad_path, 'sugar_analysis_data/SNF-0203-CABALLO2/meta_spectra.pkl')
        self.spec = pkl.load(open(pkl_file, 'rb'), encoding="latin1")
        self.meta= pkl.load(open(self.sad_path+'sugar_analysis_data/SNF-0203-CABALLO2/META.pkl', 'rb'), encoding="latin1")
      
    def integral_to_phot(self, wave, flux, var=None):
        
        self.start = min(wave)
        self.end =  max(wave)
        # Interpolate filter over spectrum (natively regularly sampled)
        f = self.interpolate(wave)
        f *= wave/(CLIGHT*1.0e13*HPLANCK)
        
        #computation of the integral
        dxs = (float(self.end - self.start)/(len(wave)-1))
        phot = np.sum(flux*f*dxs)
        if  isinstance(var, np.ndarray):
            dphot = np.sqrt((var*f**2).sum())*dxs
        else:
            dphot = None
#        for i in range(len(f)):
#            if f[i] > 1:
#                print('%e'%(f[i]))
        return phot, dphot
    
    def interpolate(self, x):
        """
        Interpolate filter over input array x (not necessarily
        regularly sampled).
        Beware that interpolation doesn't go wild *outside* the filter
        definition range, or even *inside*.
        """
        from scipy.interpolate import UnivariateSpline
        filt2 = np.genfromtxt(self.sad_path+'sugar_analysis_data/data/Instruments/Florian/'+self.band+'.dat')
        transp = 0.000
        for i, trans in enumerate(filt2[:,1]):
            if trans > 0.000 and transp == 0.000:
                start = filt2[:,0][i]
            if transp > 0.000 and trans == 0.000:
                end = filt2[:,0][i-1]
            transp = trans
        y = self.pixel_fracs(x, start, end)
        
        if self.noTH :
            inside = y > 0
            spline = UnivariateSpline(filt2[:,0], filt2[:,1], s=0)
            y[inside] = np.maximum(spline(x[inside]), 0)
#            spline = pchip(filt2[:,0], filt2[:,1])
#            y[inside] = spline(x[inside])
        return y

    def pixel_fracs(self, x, xmin, xmax):
        """Return which fraction of each (Voronoi) pixel centered on (not
        necessaryly regularly sampled) x is contained within wavelength
        range [xmin,xmax]."""
    
        x = np.asarray(x)        # Center of pixels
        n = len(x)
    
        # Pixel boundaries
        dx = np.diff(x) / 2        # Half px width
        xb = np.concatenate(([x[0] - dx[0]],    # 1st px left boundary
                            x[1:] - dx,        # Middle pixel left boundaries
                            [x[-1] + dx[-1]]))  # Last px right boundary
    
        # Fractions
        f = np.zeros(n)
        f[(xmin <= x) & (x <= xmax)] = 1.  # Rough approach 1st
    
        # Compute precise fractions at domain boundaries
        imin = xb.searchsorted(xmin)  # xb[imin-1] < xmin < xb[imin]
        imax = xb.searchsorted(xmax)  # xb[imax-1] < xmax < xb[imax]
        if 0 < imin <= n:
            f[imin - 1] = (min(xb[imin], xmax) - xmin) / (xb[imin] - xb[imin - 1])
        if 0 < imax <= n:
            f[imax - 1] = (xmax - max(xb[imax - 1], xmin)) / \
                (xb[imax] - xb[imax - 1])
#        if self.band == 'RSNf':
#            for i in range(len(f)):
#                if f[i]>0.1:
#                    print(x[i], f[i])       
        return f

    def build_Astropy_Table(self, sn_name, band_used=['new_fU_10','fB_10','fV_10','fR_10','new_fI_10']):
        """
        """    
        time = []
        flux = []
        flux_err = []
        band = []
        zp = []
        zpsys = []
        vega = sncosmo.get_magsystem('vega_snf_10')

        
        for b in band_used:
            self.band = b
            if b == 'fB_10' or b == 'new_fU_10' or b == 'BSNf':
                spectra = 'spectra_B'
            elif b == 'fV_10' or  b == 'fR_10' or  b == 'new_fI_10'  or b == 'VSNf'or b == 'RSNf':
                spectra = 'spectra_R'
            for p in self.spec[sn_name][spectra].keys():
                phot_flux, dphot = self.integral_to_phot(
                        self.spec[sn_name][spectra][str(p)]['wave'],
                        self.spec[sn_name][spectra][str(p)]['flux'],
                        var=self.spec[sn_name][spectra][str(p)]['var']) 
                flux.append(phot_flux)
                time.append(float(p))
                band.append(self.band)
                flux_err.append(dphot*self.meta[sn_name]['target.errorscale'])
                zp.append(2.5*np.log10(vega.zpbandflux(self.band)))
                zpsys.append('vega_snf_10')
                
        data = Table([time, band, flux, flux_err, zp, zpsys], 
                     names=('time', 'band', 'flux', 'fluxerr', 'zp', 'zpsys'), 
                     meta={'name': 'data'})
        return data

    
    def dump_pkl(self, path='../../sugar_analysis_data/', filters=['new_fU_10','fB_10','fV_10','fR_10','new_fI_10']):
        pkl_path = path+'phot_SNF_'.join(filters)+'.pkl'
        dic = {}
        for sn_name in self.meta.keys():
            print(sn_name)
            if self.meta[sn_name]['idr.subset'] == 'training' or self.meta[sn_name]['idr.subset'] == 'validation':
                dic[sn_name] = self.build_Astropy_Table(sn_name, band_used=filters)             
        pkl.dump(dic, open(pkl_path, 'wb'))  
        
        
class read_csp(object):
    def __init__(self, sad_path='../../'):
        self.sad_path = sad_path
        self.file_path = self.sad_path+'sugar_analysis_data/DR3/'
    def build_csp_table(self, sn_name,  drop=None):
        """
        """
        csp_filter_dict = {'u' : 'cspu', 'g' : 'cspg', 'r' : 'cspr', 'i' : 'cspi',
                   'B' : 'cspb', 'V0' : 'cspv3014', 'V1' : 'cspv3009',
                   'V' : 'cspv9844', 'Y' : 'cspys', 'H' : 'csphs',
                   'J' : 'cspjs', 'Jrc2' : 'cspjs', 'Ydw' : 'cspyd',
                   'Jdw' : 'cspjd', 'Hdw' : 'csphd'}
        self.infile = open(self.file_path+sn_name,'r')
        self.inline = self.infile.readlines()
        time = []
        flux = []
        flux_err = []
        band = []
        zp = []
        zpsys = []
        vega = sncosmo.get_magsystem('csp')
        zhl = float(self.inline[0].split()[1])
        for line in self.inline :
            line = line.split()
            if line[0] == 'filter':
                f = csp_filter_dict[line[1]]
            elif len(line) == 3:
                if type(drop) == str:
                    if f == drop:
                        continue
                elif type(drop) == list:
                    if f in drop:
                        continue
                time.append(float(line[0]))
                flux.append(vega.band_mag_to_flux(float(line[1]),f))
                flux_err.append(float(line[2])*vega.band_mag_to_flux(float(line[1]),f) / 1.0857362047581294)
                band.append(f)
                zp.append(2.5*np.log10(vega.zpbandflux(f)))  
                zpsys.append('csp')

        data = Table([time, band, flux, flux_err, zp, zpsys], 
                     names=('time', 'band', 'flux', 'fluxerr', 'zp', 'zpsys'), 
                     meta={'name': 'data'})   
        

        return data, zhl

    def get_mwebv(self, ra, dec, zhelio):
        """
        """
        import sfdmap
        from astropy.coordinates import Galactic
        c = SkyCoord(ra+' '+dec, unit=(u.deg, u.deg))
        m = sfdmap.SFDMap(self.sad_path+'sugar_analysis_data/sfddata-master', scaling=1.0)
        galcoords = c.transform_to(Galactic)
        gc = galcoords.galactic.l.deg,  galcoords.galactic.b.deg
        return m.ebv(c.ra.deg, c.dec.deg), zhelio2zcmb(zhelio,gc)

    

         
