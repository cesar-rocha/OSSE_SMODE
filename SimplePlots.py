from netCDF4 import Dataset
from scipy.signal import detrend
import numpy as np
import matplotlib.pyplot as plt
import cmocean
import gsw

plt.close("all")
__figdir__ = "Figz/"

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
lon, lat   = surface['lon_rho'][:], surface['lat_rho'][:]
lond, latd = (lon[1:]+lon[:-1])/2, (lat[1:]+lat[:-1])/2, 

#  Plotting
snap = 240
figsize = (12,6)

fig = plt.figure(figsize=figsize)

ax1 = fig.add_subplot(131,aspect=1)
speed = np.sqrt( surface['u'][snap,0]**2 + 
                 surface['v'][snap,0]**2 )
imspeed = plt.pcolor(lon,lat,speed,vmin=0.,vmax=1.5, 
                     cmap=cmocean.cm.ice_r)

ax2 = fig.add_subplot(132,aspect=1)
imvort = plt.pcolor(lon,lat,vort['rvort'][snap]/gsw.f(lat),vmin=-5,
                    vmax=5, cmap=cmocean.cm.curl)

ax3 = fig.add_subplot(133, aspect=1)
imdiv = plt.pcolor(lond,latd,div['surface_div'][240]/gsw.f(latd),vmin=-5,
                   vmax=5, cmap=cmocean.cm.balance)

ax1.set_xlabel("Longitude")
ax1.set_ylabel("Latitude")
ax2.set_xlabel("Longitude")
ax2.set_yticks([])
ax3.set_xlabel("Longitude")
ax3.set_yticks([])

fig.subplots_adjust(top=0.8)
cbar_ax = fig.add_axes([.1425, .835, 0.2, 0.02])
fig.colorbar(imspeed, cax=cbar_ax, extend='max', 
             orientation = 'horizontal',label='Speed [m/s]')

cbar_ax = fig.add_axes([.4125, .835, 0.2, 0.02])
fig.colorbar(imvort, cax=cbar_ax, extend='both', 
             orientation = 'horizontal',label=r'$\zeta/f$')

cbar_ax = fig.add_axes([.6915, .835, 0.2, 0.02])
fig.colorbar(imdiv, cax=cbar_ax, extend='both', 
             orientation = 'horizontal',label=r'$\delta/f$')

plt.savefig(__figdir__+"Snapshot"+str(snap)+"_SpeedVortDiv")


fig = plt.figure(figsize=(5,5))

ax1 = fig.add_subplot(121,aspect=1)
imvort = plt.pcolor(lon,lat,vort['rvort'][snap]/gsw.f(lat),vmin=-5,
                    vmax=5, cmap=cmocean.cm.curl)

plt.xlim(-88,-87.75)
plt.ylim(28.75,29.5)

ax2 = fig.add_subplot(122,aspect=1)
imdiv = plt.pcolor(lond,latd,div['surface_div'][240]/gsw.f(latd),vmin=-5,
                   vmax=5, cmap=cmocean.cm.balance)
plt.xlim(-88,-87.75)
plt.ylim(28.75,29.5)

ax1.set_xlabel("Longitude")
ax1.set_ylabel("Latitude")
ax2.set_xlabel("Longitude")
ax2.set_yticks([])
ax3.set_xlabel("Longitude")
ax3.set_yticks([])

plt.savefig(__figdir__+"Snapshot"+str(snap)+"_VortDiv_Zoom")

# Histograms
hist_args = {'density': True, 'bins': 30, 'range': (-7.5,7.5),
             'rwidth': 0.9}

fig = plt.figure(figsize=(8.5,4))

ax1 = fig.add_subplot(121)
plt.hist((div['surface_div'][snap]/gsw.f(latd)).flatten(),**hist_args)

ax2 = fig.add_subplot(122)
plt.hist((vort['rvort'][snap]/gsw.f(lat)).flatten(),**hist_args)

ax1.set_xlabel(r'$\delta/f$')
ax2.set_xlabel(r'$\zeta/f$')
ax1.set_ylabel(r'PDF')
ax2.set_ylabel(r'PDF')

plt.savefig(__figdir__+"PDFs"+str(snap)+"_VortDiv")



