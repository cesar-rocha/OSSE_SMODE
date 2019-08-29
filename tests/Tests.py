import numpy as np
import unittest
from SpectralUtils import *

def test_parseval(nt=28,n=256,rtol=1.e-4):
    """
        Test Parseva's relationship.
    """

    # Create a periodic random signal
    A = np.random.randn(nt,n,n)

    # Estimate spectra
    spectra = OneDimSpectra(A, dx=1/n,dy=1/n,
              windowing=False,detrending=False)

    vardensy = spectra['y']['spectrum']
    vardensx = spectra['x']['spectrum']

    # Count Nyquist frequency only once
    vardensy[-1] = vardensy[-1]/2
    vardensx[-1] = vardensx[-1]/2

    # Variance in spectra space
    specvar = np.array([vardensy[1:].sum(), vardensx[1:].sum()])

    # Variance in physical space
    physvar = np.array([A.var(axis=2).mean(),A.var(axis=1).mean()])

    # Test if differences are within rtol
    assert(np.allclose(specvar, physvar, rtol=rtol))

if __name__ == "__main__":
    test_parseval()
