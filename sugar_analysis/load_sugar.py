#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 21 17:59:43 2019

@author: florian
"""
import sncosmo
import numpy as np
import os
from astropy.extern import six
from sncosmo.models import _SOURCES
from scipy.interpolate import interp2d

class SUGARSource(sncosmo.Source):
    """
    The SUGAR Type Ia supernova spectral timeseries model.
    
    The spectral flux density of this model is given by
    
    .. math::
       
    F(t, \\lambda) = Xgr 10^{-0.4 (M_0(t, \\lambda) 
                                + q_1 M_1(t, \\lambda)
                                + q_1 M_1(t, \\lambda)
                                + q_2 M_2(t, \\lambda)
                                +q_3 M_3(t, \\lambda) 
                                + A M_4(t, \\lambda))} (10^{-3} c\\lambda^{2})

    where ``Xgr``, ``q_1``, ``q_2``, ``q_3``  and ``A`` are the free parameters
    of the model,``M_0``, ``M_1``, `M_2``, `M_4``, `M_4`` are the zeroth and 
    first components of the model.     

    Parameters
    ----------
    modeldir : str, optional
        Directory path containing model component files. Default is `None`,
        which means that no directory is prepended to filenames when
        determining their path.
        m0file, m1file, m2file, m3file, m4file: str or fileobj, optional
        Filenames of various model components. Defaults are:

        * m0file = 'sugar_template_0.dat' (2-d grid)
        * m1file = 'sugar_template_1.dat' (2-d grid)
        * m2file = 'sugar_template_1.dat' (2-d grid)
        * m3file = 'sugar_template_1.dat' (2-d grid)
        * m4file = 'sugar_template_1.dat' (2-d grid)
        
    Notes
    -----
    The "2-d grid" files have the format ``<phase> <wavelength>
    <value>`` on each line.    
    """
    _param_names = ['Xgr', 'q1', 'q2', 'q3', 'A']
    param_names_latex = ['X_gr', 'q_1', 'q_2', 'q_3', 'A']
    

    def __init__(self, modeldir='../../sugar_model/', m0file='sugar_template_0.dat',
                 m1file='sugar_template_1.dat',
            		m2file='sugar_template_2.dat',
            		m3file='sugar_template_3.dat',
                 m4file='sugar_template_4.dat',
                 mod_errfile='../../sugar_model/model_err_sug.dat',
                 name=None, version=None):
        self.name = name
        self.version = version
        self._model = {}
        self._parameters = np.array([1.0e10-15, 1., 1., 1., 0.])
        self.M_keys = ['M0', 'M1', 'M2', 'M3', 'M4']    
        names_or_objs = {'M0': m0file, 'M1': m1file, 'M2': m2file, 'M3': m3file, 'M4': m4file}

        # Make filenames into full paths.
        if modeldir is not None:
            for k in names_or_objs:
                v = names_or_objs[k]
                if (v is not None and isinstance(v, six.string_types)):
                    names_or_objs[k] = os.path.join(modeldir, v)
                    
        for i, key in enumerate(self.M_keys):
            phase, wave, values = sncosmo.read_griddata_ascii(names_or_objs[key])
            self._model[key] = sncosmo.salt2utils.BicubicInterpolator(phase, wave, values) 
            if key == 'M0':
                # The "native" phases and wavelengths of the model are those
                self._phase = np.array([ -12.,  -9.,  -6.,  -3.,   0.,   3., 
                                        6.,   9.,  12.,  15.,  18., 21.,  24., 
                                        27.,  30.,  33.,  36.,  39.,  42.,
                                        45.,  48.])
                self._wave = wave
        self.mod_errfile = mod_errfile
        phase, wave, values = sncosmo.read_griddata_ascii(mod_errfile)
        self._model['mod_err'] = interp2d(wave, phase, values) 
        
    def _flux(self, phase, wave):

        M_sug = 48.59
        for i, key in enumerate(self.M_keys):
            if i >= 1.: 
                M_sug += self._model[key](phase, wave) * self._parameters[i]
            else: 
                M_sug += self._model[key](phase, wave)
        wave2_factor = (wave ** 2 / 299792458. * 1.e-10)
        return (self._parameters[0] * 10. ** (-0.4 * M_sug) /  wave2_factor)

    def _bandflux_rvar_single(self, band, phase):
        """Model relative variance for a single bandpass."""
        # Raise an exception if bandpass is out of model range.
        if (band.minwave() < self._wave[0] or band.maxwave() > self._wave[-1]):
            raise ValueError('bandpass {0!r:s} [{1:.6g}, .., {2:.6g}] '
                             'outside spectral range [{3:.6g}, .., {4:.6g}]'
                             .format(band.name, band.wave[0], band.wave[-1],
                                     self._wave[0], self._wave[-1]))
        mod_err = self._model['mod_err'](band.wave_eff, phase)[:, 0]
        
        # v is supposed to be variance but can go negative
        # due to interpolation. Correct negative values to some small
        # number. (at present, use prescription of snfit : set
        # negatives to 0.0001)
        mod_err[mod_err < 0.0] = 0.0001
        result = (mod_err * self.bandflux(band, phase)) / 1.0857362047581294
        return result
    
    def bandflux_rcov(self, band, phase):
        """
        model error in comming
        """
        # construct covariance array with relative variance on diagonal
        diagonal = np.zeros(phase.shape, dtype=np.float64)
        for b in set(band):
            mask = band == b
            diagonal[mask] = self._bandflux_rvar_single(b, phase[mask])**2
        result = np.diagflat(diagonal)
        
        return result
    
def register_SUGAR(modeldir='../../sugar_model/',
                   mod_errfile='../../sugar_model/model_err_sug.dat',
                   version='1.0'):
    website = 'http://no'
    PF16ref = ('PF16', 'PF et al. 2016 '
              '<http://arxiv.org/>')
    for topdir, ver, ref in [('SUGAR_model', version, PF16ref)]:
        meta = {'type': 'SN Ia', 'subclass': '`~sncosmo.SUGARSource`',
                'url': website, 'reference': ref}
        print 'ca devrait etre :' , mod_errfile, version

        _SOURCES.register_loader('sugar', lambda relpath, mod_errfile,
        name=None, version=None : SUGARSource(modeldir=relpath,
                                              mod_errfile=mod_errfile,
                                              name=name, version=version)
        , args=([modeldir, mod_errfile]), version=ver, meta=meta, force=True)

        