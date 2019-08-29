from netCDF4 import Dataset
from scipy.signal import detrend
import numpy as np
import matplotlib.pyplot as plt
import gsw

plt.close("all")

# Directory with ROMS output provided by 
# Roy Barkan to Eric D'Asaro during LASER.
__datadir__ = "DATA/"

# Surface fields and derived quantities.
surface = Dataset(__datadir__+"surface_150.nc")
div     = Dataset(__datadir__+"div150.nc")
vort    = Dataset(__datadir__+"vort150.nc")

# Model grid (dx ≈ 450 m, dy ≈ 450 m).
# Vorticity seems to be on the same grid,
# but divergence has one less point in x.
lon, lat = surface['lon_rho'][:], surface['lat_rho'][:]


