# Utility functions for spectral calculations
# with UCLA's ROMS output.
# Cesar B Rocha
# WHOI, Summer 2019

from scipy.signal import detrend
import numpy as np

def Wavenumbers(nx,dx):
    """
        Calculates the wavenumber given
        the length of an array nx and the
        resolution dx.

        Return
        - wavenumber array
    """
    dk = 1./(nx*dx)

    if nx%2:
        k = dk*np.arange( (nx-1)/2.  + 1 )
    else:
        k = dk*np.arange(nx/2+1)
    return k

def OneDimSpectra(f,dx,dy,windowing=True,detrending=True):
    """
        Calculates one dimensional spectra
        of ROMS model field f.

        Return
        - specdic: a dictionary with one-dimensional
                   spectra, label by axes ('x','y'),
                   with arrays specdic[label]['wavenumber']
                   and specdic[label]['spectrum'].
    """
    # Spectra grid
    nt, ny, nx = f.shape
   
    # Loop over axes (y,x) = (1,2)
    for n,d,axis,label in zip((nx,ny),(dx,dy),(2,1),('x','y')):

        # Get wavenumbers
        k = Wavenumbers(n,d)

        # Spectral window
        if windowing:
            win  = np.hanning(n)
            win *= np.sqrt(n/(win**2).sum())
        else:
            win = np.ones(n)

        if axis == 1:
            win = win[np.newaxis,:,np.newaxis]
        elif axis == 2:
            win = win[np.newaxis,np.newaxis,:]

        # Remove linear trend
        if detrending:
            f = detrend(f,axis=axis, type='linear')

        # Calculate the spectrum
        fh = np.fft.rfft(f*win, axis=axis)
        spec = 2*(np.abs(fh)**2)/(n/d)

        # Average over non-Fourier axes  
        if axis == 1:
            spec = spec.mean(axis=(0,-1))
        elif axis == 2:
            spec = spec.mean(axis=(0,1))

        # Store one-dimensional spectra and waveumber
        # in a dictionary.
        try:
            specdic.update({label: {"wavenumber": k, "spectrum": spec}})
        except:
            specdic = {label: {"wavenumber": k, "spectrum": spec}}

    return specdic
