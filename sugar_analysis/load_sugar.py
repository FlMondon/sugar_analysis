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
    param_names_latex = ['X_r', 'q_1', 'q_2', 'q_3', 'A']
    

    def __init__(self, modeldir='../../sugar_model/', m0file='sugar_template_0.dat',
                 m1file='sugar_template_1.dat',
            		m2file='sugar_template_2.dat',
            		m3file='sugar_template_3.dat',
                 m4file='sugar_template_4.dat',
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
                self._phase = np.array([ -12.,  -9.,  -6.,  -3.,   0.,   3.,   6.,   9.,  12.,  15.,  18.,
                                                21.,  24.,  27.,  30.,  33.,  36.,  39.,  42.,  45.,  48.])
                self._wave = wave


    def _flux(self, phase, wave):
        """
        TO DO
        """
        M_sug = 48.59
        for i, key in enumerate(self.M_keys):
            if i >= 1.: 
                M_sug += self._model[key](phase, wave) * self._parameters[i]
            else: 
                M_sug += self._model[key](phase, wave)
        wave2_factor = (wave ** 2 / 299792458. * 1.e-10)
        return (self._parameters[0] * 10. ** (-0.4 * M_sug) /  wave2_factor)
    

    def bandflux_rcov(self, band, phase):
        """
        model error in comming
        """
        return np.zeros(phase.shape, dtype=np.float64)
    

# Sugar model
def load_sugarmodel(relpath, name=None, version=None):
    return SUGARSource(modeldir=relpath, name=name, version=version)

def register_SUGAR(modeldir='../../sugar_model/'):
    website = 'http://no'
    PF16ref = ('PF16', 'PF et al. 2016 '
              '<http://arxiv.org/>')
    for topdir, ver, ref in [('SUGAR_model', '1.0', PF16ref)]:
        meta = {'type': 'SN Ia', 'subclass': '`~sncosmo.SUGARSource`',
                'url': website, 'reference': ref}
        _SOURCES.register_loader('sugar', load_sugarmodel, args=([modeldir]), version=ver, meta=meta, force=True)