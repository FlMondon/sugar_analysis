""" test the math toolbox. """

import numpy as np
import sugar_analysis
from timer_test import timer

def generate_hubble_residuals(sigma_int = 0.13, noise = None):
    """
    Generate Hubble residual as a normal law.

    :sigma_int: float. Intrinsic scatter of SNIa
    :noise:     float. Noise in the hubble residuals
                Default : None. 
    """

    np.random.seed(171094)
    n_sn = 10000
    hr = np.random.normal(scale = sigma_int, size = n_sn)
    if noise is not None: 
        hr += np.random.normal(scale = noise, size = n_sn)
        hr_err = np.ones_like(hr) * noise
    else:
        hr_err = None

    return hr, hr_err

@timer
def test_comp_rms():
    """
    test comp_rms from math_toolbox.py.
    """

    sig_int = 0.13
    noise = 0.08

    hr, hr_err = generate_hubble_residuals(sigma_int=sig_int, noise=None)
    rms, rms_err = sugar_analysis.comp_rms(hr, len(hr), err=True, variance=None)
    np.testing.assert_allclose(rms, sig_int, atol=1e-3)

    hr, hr_err = generate_hubble_residuals(sigma_int=sig_int, noise=noise)
    rms, rms_err = sugar_analysis.comp_rms(hr, len(hr), err=True, variance=hr_err**2)
    np.testing.assert_allclose(rms, np.sqrt(sig_int**2 + noise**2), atol=1e-3)


if __name__ == '__main__':

    test_comp_rms()
