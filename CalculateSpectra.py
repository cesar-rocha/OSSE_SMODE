from netCDF4 import Dataset
from scipy.signal import detrend
import numpy as np
import matplotlib.pyplot as plt
import gsw

from SpectralUtils import *

plt.close("all")

# Directory with ROMS output provided by 
# Roy Barkan to Eric D'Asaro during LASER.
__datadir__ = "DATA/"

# Surface fields
surface = Dataset(__datadir__+"surface_150.nc")

# Model grid (dx ≈ 450 m, dy ≈ 450 m).
# Vorticity seems to be on the same grid,
# but divergence has one less point in x.
lon, lat = surface['lon_rho'][:], surface['lat_rho'][:]

# Calculate wavenumber spectra.
ntmax = 240         # 10 days
dx, dy = 0.45, 0.45 # km

spectrau = OneDimSpectra(f=surface['u'][:ntmax,0], dx=0.45,dy=0.45)
spectrav = OneDimSpectra(f=surface['v'][:ntmax,0], dx=0.45,dy=0.45)


