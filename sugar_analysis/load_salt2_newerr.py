#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 17 17:21:54 2019

@author: florian
"""

import sncosmo
import numpy as np
import os
import six
from sncosmo.models import _SOURCES
from scipy.interpolate import interp2d

class SALT2newerrSource(sncosmo.Source):
    """The SALT2 Type Ia supernova spectral timeseries model.

    The spectral flux density of this model is given by

    .. math::

       F(t, \lambda) = x_0 (M_0(t, \lambda) + x_1 M_1(t, \lambda))
                       \\times 10^{-0.4 CL(\lambda) c}

    where ``x0``, ``x1`` and ``c`` are the free parameters of the model,
    ``M_0``, ``M_1`` are the zeroth and first components of the model, and
    ``CL`` is the colorlaw, which gives the extinction in magnitudes for
    ``c=1``.

    Parameters
    ----------
    modeldir : str, optional
        Directory path containing model component files. Default is `None`,
        which means that no directory is prepended to filenames when
        determining their path.
    m0file, m1file, clfile : str or fileobj, optional
        Filenames of various model components. Defaults are:

        * m0file = 'salt2_template_0.dat' (2-d grid)
        * m1file = 'salt2_template_1.dat' (2-d grid)
        * clfile = 'salt2_color_correction.dat'

    errmod : str or fileobj
        (optional) Filenames of various model components for
        model covariance in synthetic photometry. See
        ``bandflux_rcov`` for details.  Defaults are:
    

    Notes
    -----
    The "2-d grid" files have the format ``<phase> <wavelength>
    <value>`` on each line.

    The phase and wavelength values of the various components don't
    necessarily need to match. (In the most recent salt2 model data,
    they do not all match.) The phase and wavelength values of the
    first model component (in ``m0file``) are taken as the "native"
    sampling of the model, even though these values might require
    interpolation of the other model components.

    """

    _param_names = ['x0', 'x1', 'c']
    param_names_latex = ['x_0', 'x_1', 'c']
    _SCALE_FACTOR = 1e-12

    def __init__(self, modeldir=None,
                 m0file='salt2_template_0.dat',
                 m1file='salt2_template_1.dat',
                 clfile='salt2_color_correction.dat',
                 mod_errfile='../../sugar_analysis/err_mod_training/salt2_newerr/mod_errsalt2.dat',
                 name=None, version=None):
        self.name = name
        self.version = version
        self._model = {}
        self._parameters = np.array([1., 0., 0.])

        names_or_objs = {'M0': m0file, 'M1': m1file, 
                         'mod_err': mod_errfile, 'clfile': clfile}

        # Make filenames into full paths.
        if modeldir is not None:
            for k in names_or_objs:
                v = names_or_objs[k]
                if (v is not None and isinstance(v, six.string_types)):
                    names_or_objs[k] = os.path.join(modeldir, v)

        # model components are interpolated to 2nd order
        for key in ['M0', 'M1']:
            phase, wave, values = sncosmo.read_griddata_ascii(names_or_objs[key])
            values *= self._SCALE_FACTOR
            self._model[key] = sncosmo.salt2utils.BicubicInterpolator(phase, wave, values)

            # The "native" phases and wavelengths of the model are those
            # of the first model component.
            if key == 'M0':
                self._phase = phase
                self._wave = wave
        phase, wave, values = sncosmo.read_griddata_ascii(mod_errfile)
        self._model['mod_err'] = interp2d(wave, phase, values) 

        # Set the colorlaw based on the "color correction" file.
        self._set_colorlaw_from_file(names_or_objs['clfile'])
        
    def _flux(self, phase, wave):
        m0 = self._model['M0'](phase, wave)
        m1 = self._model['M1'](phase, wave)
        return (self._parameters[0] * (m0 + self._parameters[1] * m1) *
                10. ** (-0.4 * self._colorlaw(wave) * self._parameters[2]))

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
        model error in comming
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


def register_salt2_newerr(modeldir='../../sugar_analysis_data/err_mod_training/salt2_newerr/',
                 mod_errfile='mod_errsalt2.dat',
                 version='1.0'):

        _SOURCES.register_loader('salt2_newerr', lambda modeldir, mod_errfile,
        name=None, version=None : SALT2newerrSource(modeldir=modeldir,
                                              mod_errfile=mod_errfile,
                                              name=name, version=version)
        , args=([modeldir, mod_errfile]), version=version, force=True)
      