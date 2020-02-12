#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 24 16:23:09 2018

@author: florian
"""
import numpy as np
from scipy import integrate

from .constant import CLIGHT, H0
# ------------------ #
#   Cosmology        #
# ------------------ #
#: Pre-defined apices: galactic coordinates [deg] and velocity [km/s]

def int_cosmo(z, omgM, omgK, w):   
    """
    """
    return 1./np.sqrt(omgM*(1+z)**3 + (1-omgM-omgK)*(1+z)**(3*(1+w)) + omgK*(1+z)**2)

def _velproj(coords, apex, deg=True):
    """Compute projection factor of velocity in the direction of
    *apex* [deg|rad].
    """

    lp,bp = coords                      # Input coordinates
    la,ba = apex                        # Apex coordinates
    RAD2DEG = 57.295779513082323            # 180/pi
    if deg:                             # Convert to radians
        la,ba,lp,bp = [ a/RAD2DEG for a in (la,ba,lp,bp) ]

    return np.sin(bp)*np.sin(ba) + np.cos(bp)*np.cos(ba)*np.cos(lp-la)

def v2apex(coords, apex='CMB'):
    """Velocity correction (to be added) from input reference frame
    (*coords* in deg) to *apex* reference frame (*apex* =
    (l[deg],b[deg],v[km/s])).

    By default, `apex='CMB'` is expressed in galactic coordinates, and
    therefore so should be input coordinates. See `APICES` for
    pre-defined (local group, galactocentric or CMB) apices.

    .. Note:: similar to `NED Velocity Correction Calculator
              <http://ned.ipac.caltech.edu/forms/vel_correction.html>`_
              but whith sightly different results due to different
              conversion routines to galactic coordinates and CMB apex
              definition.

    >>> v2apex(radec2gal(123.456, 12.3456))
    245.83192891291446
    """
    APICES = {
        'COBE':(264.14,48.26,371.),  # CMB COBE/Firas (1996ApJ...473..576F) (NED)
        'LG':  ( 93.  ,-4.  ,316.),  # Local Group (1996AJ....111..794K)
        'GSR': ( 87.8 , 1.7 ,232.3), # Galactocentric (1991RC3...M...0000d)
        'WMAP':(263.99,48.26,369.),  # CMB WMAP-5/7yr (2009ApJS..180..225H)
        }
    APICES['CMB'] = APICES['WMAP']   # Use WMAP for CMB
    if apex in APICES:
        apex = APICES[apex]
    else:
        try:
            print("Velocity conversion apex: " \
                  "l=%f deg, b=%f deg, v=%f km/s" % apex)
        except TypeError:
            raise TypeError("Invalid apex '%s': l[deg],b[deg],v[km/s] or %s" %
                            (str(apex),','.join(APICES.keys())))

    la,ba,va = apex      # Apex galactic coords [deg], velocity [km/s]

    return va * _velproj((la,ba), coords, deg=True)


def zhelio2zcmb(zhelio, galcoords):
    """Converts heliocentric redshift *zhelio* to CMB-centric redshift
    *zcmb*, given input galactic coordinates *galcoords* [deg]."""

    return zhelio + v2apex(galcoords, apex='CMB')/CLIGHT

def zcmb2zhelio(zcmb, galcoords):
    """Converts CMB-centric redshift *zcmb* to heliocentric redshift
    *zhelio*, given input galactic coordinates *galcoords* [deg]."""

    return zcmb - v2apex(galcoords, apex='CMB')/CLIGHT
    
def luminosity_distance(omgM, omgK, w, zcmb, zhl):
    """
    """
    """
    """

    if type(zcmb)==np.ndarray:
        integr = np.zeros_like(zcmb)
        for i in range(len(zcmb)):
            integr[i]=integrate.quad(int_cosmo, 0, zcmb[i], args=(omgM, omgK, w))[0]
    else:
        integr = integrate.quad(int_cosmo, 0, zcmb, args=(omgM, omgK, w))[0]
    
    if omgK == 0. :
        return (1+zhl)*(CLIGHT/H0)*integr
    elif omgK < 0.:
        return (1+zhl)*(1/np.sqrt(np.abs(omgK)))*(np.sin(np.sqrt(np.abs(omgK))*integr)) *(CLIGHT/H0)
    elif omgK > 0.:
        return (1+zhl)*(1/np.sqrt(np.abs(omgK)))*(np.sinh(np.sqrt(np.abs(omgK))*integr))*(CLIGHT/H0)
    
def distance_modulus_th(zcmb, zhl, omgM=0.3, omgK=0., w=-1.):
    """
    """
    return 5.*np.log(luminosity_distance(omgM, omgK, w, zcmb, zhl))/np.log(10.)-5.



