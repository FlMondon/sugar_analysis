#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  5 10:57:36 2019

@author: florian
"""

import sncosmo
import numpy as np
import os
import six
from sncosmo.models import _SOURCES
from scipy.interpolate import interp2d
from .constant import wl_min_sug, wl_max_sug
try:
   import cPickle as pkl
except ModuleNotFoundError:
   import pickle as pkl
   
class SUGAR_extSource(sncosmo.Source):
    """
    The SUGAR Type Ia supernova spectral timeseries model 
    + extension in the UV using SALT2.
    
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
    _param_names = ['Xgr', 'q1',  'q3', 'A']
    param_names_latex = ['X_gr', 'q_1', 'q_3', 'A']
    _SCALE_FACTOR = 1e-12
    def __init__(self, modeldir='../../sugar_model/', 
                 mod_errfile='../../sugar_model/model_err_sug.dat',
                 transmat_path='trans_matrix_init.pkl', name=None, version=None) :
        
        self.name = name
        self.version = version
        self._model = {}
        
        list_file = [f for f in os.listdir(modeldir) if f.startswith('sugar_template_')]
        names_or_objs = {}
        self._parameters = np.zeros(len(list_file)-1)
        self._parameters[0] = 1e-15
        for i in range(len(list_file)):
            names_or_objs['M'+str(i)] = os.path.join(modeldir, 'sugar_template_'+str(i)+'.dat')
        names_or_objs['L0']  = os.path.join(modeldir, 'salt2_template_0.dat')
        names_or_objs['L1']  = os.path.join(modeldir, 'salt2_template_1.dat')
        names_or_objs['clfile'] = os.path.join(modeldir, 'salt2_color_correction.dat')
        wave_list_min = []
        wave_list_max = []
        for i, key in enumerate(names_or_objs.keys()):
            if key != 'clfile':
                phase, wave, values = sncosmo.read_griddata_ascii(names_or_objs[key])
                self._model[key] = sncosmo.salt2utils.BicubicInterpolator(phase, wave, values) 
                if key == 'M0':
                    # The "native" phases and wavelengths of the model are those
                    self._phase = [-15] + np.linspace(-12., 48, 21)
                    wave_m0 = list(wave)
                if key == 'L0':
                    for wl in wave:
                        if wl <= wl_min_sug :
                            wave_list_min.append(wl)
                        elif wl >= wl_max_sug :
                            wave_list_max.append(wl)
        wave_list = wave_list_min + wave_m0 + wave_list_max
        wave_list = np.array(wave_list)     
        argso=np.argsort(wave_list)
        wave_list=wave_list[argso]        
        self._wave = wave_list
        #load transform matrix sugar to salt2
        self.trans_matrix = pkl.load(open(os.path.join(modeldir, transmat_path), 'rb'))
        # Set the colorlaw based on the "color correction" file.
        self._set_colorlaw_from_file(names_or_objs['clfile'])
        self.t_cut = -5

        self.mod_errfile = mod_errfile
        self.M_keys = names_or_objs.keys()
        if mod_errfile is not None:
            phase, wave, values = sncosmo.read_griddata_ascii(mod_errfile)
            self._model['mod_err'] = interp2d(wave, phase, values) 
        self._SCALE_FACTOR_SALT2 = 1e-12

    def _flux(self, phase, wave):
#            wl_min_sug = 3800
#            print(wave[0])
            flux = self._sug_flux(self._parameters, phase, wave)

            
#            if self.t_cut is not None:
#                for i, p in enumerate(phase):
#                    p_trans = 0
#                    if p < self.t_cut + p_trans and p >= -15: 
##                        print(p)
#                        if  p_trans != 0:
#                                
#                            r = (p -self.t_cut) * (1./p_trans)
#                            if r< 0 :
#                                r = 0
#                        else:
#                            r=0
#                        self.salt2_param = self.compute_salt2_params_ext(self._parameters)
#                        m0 = self._model['L0'](p, wave)
#                        m1 = self._model['L1'](p, wave)
#                        flux_salt2 = self.salt2_param[0]* self._SCALE_FACTOR_SALT2 * (m0 +  self.salt2_param[1]* m1)* 10. ** (-0.4 * self.colorlaw(wave) * self.salt2_param[2])
#                        flux[i, :] = r * flux[i, :] + (1-r) * flux_salt2
#                        
                        
#                        print(flux[i, 0], self.salt2_param[0], self.salt2_param[1], self.salt2_param[2])
#            if (wave[0] <= wl_min_sug and wave[len(wave)-1] <= wl_min_sug) or (wave[0] >= wl_max_sug and wave[len(wave)-1] >= wl_max_sug):
#                self.salt2_param = self.compute_salt2_params_ext(self._parameters)
#                m0 = self._model['L0'](phase, wave)
#                m1 = self._model['L1'](phase, wave)
##                flux = self.salt2_param[0]* self._SCALE_FACTOR_SALT2* (m0 + self.salt2_param[1] * m1)*10. ** (-0.4 * self._model['M4'](phase, wave) * self._parameters[3])
#                flux = self.salt2_param[0]* self._SCALE_FACTOR_SALT2 * (m0 +  self.salt2_param[1]* m1)* 10. ** (-0.4 * self.colorlaw(wave) * self.salt2_param[2])
            range_wave_transi = 100
            if wave[0] <= wl_min_sug + range_wave_transi :
                j = 0
                self.salt2_param = self.compute_salt2_params_ext(self._parameters)
                while wave[j] <= wl_min_sug + range_wave_transi  and j <= len(wave)-2:
                    m0 = self._model['L0'](phase, wave)
                    m1 = self._model['L1'](phase, wave)
#                    mc = self._model['M4'](phase, wave)
                    if  range_wave_transi != 0:
                            
                        r = (wave[j] - wl_min_sug) * (1./range_wave_transi )
                        if r< 0 :
                            r = 0
                    else:
                        r=0
#                    flux[:, j] = (self.salt2_param[0]* self._SCALE_FACTOR_SALT2 * (m0[:, j] +  self.salt2_param[1]* m1[:, j])*10. ** (-0.4 * mc[:, j] * self._parameters[3]) )
                    flux_salt2 = self.salt2_param[0]* self._SCALE_FACTOR_SALT2 * (m0[:, j] +  self.salt2_param[1]* m1[:, j])* 10. ** (-0.4 * self.colorlaw(wave[j]) * self.salt2_param[2])
                    flux[:, j] = flux[:, j] * r + (1-r) * flux_salt2
                    j +=1

            return flux

    def compute_salt2_params_ext(self, parameters):
        sug_param = np.ones(len(parameters)+1)
        sug_param[1:] = parameters
        sugar_flux = self._sug_flux(parameters, 0, wl_min_sug)
        
        if sug_param[1] > 0. :
            sug_param[1] = -2.5*np.log10(sug_param[1])
        else:
            sug_param[1] = 1000.
        salt2_param = np.ones(3)
        salt2_param = sug_param.dot(self.trans_matrix)  
        mc0 = self._model['L0'](0, wl_min_sug)
        mc1 = self._model['L1'](0, wl_min_sug)
        mc = self._model['M4'](0, wl_min_sug)
        salt2_param[0] = 10**(-0.4*salt2_param[0])
#        salt2_param[0] = self._sug_flux(parameters, 0, wl_min_sug)/((self._SCALE_FACTOR_SALT2 * (mc0 +  salt2_param[1]* mc1))*10. ** (-0.4 * mc * self._parameters[3]))
        
#        salt2_param[0] = sugar_flux/((self._SCALE_FACTOR_SALT2 * (mc0 +  salt2_param[1]* mc1))* 10. ** (-0.4 * self.colorlaw(wl_min_sug) * salt2_param[2]))
#        print(self._sug_flux(parameters, 0, wl_min_sug), (self._SCALE_FACTOR_SALT2 * (mc0 +  salt2_param[1]* mc1))* 10. ** (-0.4 * self.colorlaw(wl_min_sug) * salt2_param[2]), salt2_param[0])
        return salt2_param
        
        
    def _sug_flux(self, parameters, phase, wave):
            M_sug = 48.59
            for i in range(len(self.M_keys)-4):
                if i == 0: 
                    M_sug += self._model['M'+str(i)](phase, wave)     
                elif i >= 2 :
                    j = i+1
                    M_sug += self._model['M'+str(j)](phase, wave) * parameters[i]
                else:
                    j = i
                    M_sug += self._model['M'+str(j)](phase, wave) * parameters[i]
            wave2_factor = (wave ** 2 / 299792458. * 1.e-10)
            return parameters[0] * 10. ** (-0.4 * M_sug) /  wave2_factor
            
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
        result = mod_err**2
        return result
    
    def bandflux_rcov(self, band, phase):
        """
        model error description in comming
        """
        # construct covariance array with relative variance on diagonal
        diagonal = np.zeros(phase.shape, dtype=np.float64)
        for b in set(band):
            mask = band == b
            diagonal[mask] = self._bandflux_rvar_single(b, phase[mask])
        result = np.diagflat(diagonal)
        
        return result       
        
    def _set_colorlaw_from_file(self, name_or_obj):
        """Read color law file and set the internal colorlaw function."""

        # Read file
        if isinstance(name_or_obj, six.string_types):
            f = open(name_or_obj, 'r')
        else:
            f = name_or_obj
        words = f.read().split()
        f.close()

        # Get colorlaw coeffecients.
        ncoeffs = int(words[0])
        colorlaw_coeffs = [float(word) for word in words[1: 1 + ncoeffs]]

        # If there are more than 1+ncoeffs words in the file, we expect them to
        # be of the form `keyword value`.
        version = 0
        colorlaw_range = [3000., 7000.]
        for i in range(1+ncoeffs, len(words), 2):
            if words[i] == 'Salt2ExtinctionLaw.version':
                version = int(words[i+1])
            elif words[i] == 'Salt2ExtinctionLaw.min_lambda':
                colorlaw_range[0] = float(words[i+1])
            elif words[i] == 'Salt2ExtinctionLaw.max_lambda':
                colorlaw_range[1] = float(words[i+1])
            else:
                raise RuntimeError("Unexpected keyword: {}".format(words[i]))

        # Set extinction function to use.
        if version == 0:
            raise Exception("Salt2ExtinctionLaw.version 0 not supported.")
        elif version == 1:
            self._colorlaw = sncosmo.salt2utils.SALT2ColorLaw(colorlaw_range, colorlaw_coeffs)
        else:
            raise Exception('unrecognized Salt2ExtinctionLaw.version: ' +
                            version)

    def colorlaw(self, wave=None):
        """Return the value of the CL function for the given wavelengths.

        Parameters
        ----------
        wave : float or list_like

        Returns
        -------
        colorlaw : float or `~numpy.ndarray`
            Values of colorlaw function, which can be interpreted as extinction
            in magnitudes.
        """
        if wave is None:
            wave = self._wave
        else:
            wave = np.asarray(wave)
        if wave.ndim == 0: 
            return self._colorlaw(np.ravel(wave))[0]
        else:
            return self._colorlaw(wave)        

    def update_trans_matrix(self, A):
        self.trans_matrix = A
    
    def update_t_cut(self, t_cut):
        self.t_cut = t_cut
        
def register_SUGAR_ext(modeldir='../../sugar_model/',
                   transmat_path='trans_matrix_init.pkl',
                   mod_errfile='../../sugar_model/model_err_sug.dat',
                   version='1.0'):
    website = 'http://no'
    PF16ref = ('PF16', 'PF et al. 2016 '
              '<http://arxiv.org/>')
    for topdir, ver, ref in [('SUGAR_model', version, PF16ref)]:
        meta = {'type': 'SN Ia', 'subclass': '`~sncosmo.SUGARSource`',
                'url': website, 'reference': ref}
        _SOURCES.register_loader('sugar_ext', lambda relpath, mod_errfile, transmat_path, 
        name=None, version=None : SUGAR_extSource(modeldir=relpath,
                                              mod_errfile=mod_errfile,
                                              transmat_path=transmat_path, 
                                              name=name, version=version)
        , args=([modeldir, mod_errfile, transmat_path]), version=ver, meta=meta, force=True)