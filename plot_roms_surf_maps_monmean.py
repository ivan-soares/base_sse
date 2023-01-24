###########################################################################################

#       Python program by IDS @ TOC, Florianopolis, Brazil
#                                                                       2021-02-01

###########################################################################################

import os.path, gc, sys
home = os.path.expanduser('~')
sys.path.append(home + '/scripts/python/')

import math
import numpy as np
import scipy.io as sio

from netCDF4 import Dataset
from datetime import datetime, timedelta

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator

from romstools import get_roms_date
from plot_ogcm import plot_contour_and_vectors

###########################################################################################

print(' ')
print(' +++ Starting python program +++ ')
print(' ')

################### *** read input files *** ##############################################

romsfile = str(sys.argv[1])
roms = Dataset(romsfile,'r')

ntime = int(sys.argv[2])
vint = int(sys.argv[3])
vscl = float(sys.argv[4])

W = int(sys.argv[5])
E = int(sys.argv[6])
S = int(sys.argv[7])
N = int(sys.argv[8])

ogcm = 'glorys+era5'

mon = ['jan','feb','mar','apr','may','jun',
       'jul','aug','sep','oct','nov','dec']

MON = mon[ntime]

lon = roms.variables['lon_rho'][:]
lat = roms.variables['lat_rho'][:]

nt, nlat, nlon = roms.variables['zeta'].shape

print(' ')
print(' ... ntimes = ' + str(nt))
print(' ... nlats = ' + str(nlat))
print(' ... nlons = ' + str(nlon))
print(' ')

time = get_roms_date(romsfile,ntime)

zeta = roms.variables['zeta'] [ntime,:,:]
ubar = roms.variables['ubar_eastward'] [ntime,:,:]
vbar = roms.variables['vbar_northward'] [ntime,:,:]
temp = roms.variables['temp'] [ntime,29,:,:]
uvel = roms.variables['u_eastward'][ntime,29,:,:]
vvel = roms.variables['v_northward'][ntime,29,:,:]

############ append 1 column to U and 1 line to V
############ A = numpy.hstack([A, newline])
############ A = numpy.vstack([A, newrow])

#vvel = np.hstack([vvel, vvel[nlat-2,:]])
#uvel = np.vstack([uvel, uvel[:,nlon-2]])

#uvel = np.insert(uvel, nlon-2, values=uvel[:,nlon-2], axis=1)
#vvel = np.insert(vvel, nlat-2, values=vvel[nlat-2,:], axis=0)
#
nlat, nlon = uvel.shape
ndims = str(nlat) + ' x ' + str(nlon)
print(' ... U is shaped ' + ndims)
#
nlat, nlon = vvel.shape
ndims = str(nlat) + ' x ' + str(nlon)
print(' ... V is shaped ' + ndims)

#vel = np.sqrt(uvel**2 + vvel**2)

print(' ')
print(' ... date ' + str(time))
print(' ')

################## *** plot data *** ######################################################

coastfile = "costa_leste.mat"
coast = sio.loadmat(coastfile)
levels = MaxNLocator(nbins=100).tick_values(10., 30.)
cmap = plt.get_cmap('jet')


############## plot B. Campos

wesn = [W, E, S, N]


units = 'deg C'
vartit = ' surface temp on '
varprn = '_sfc_temp_'

levels = np.arange(15,30,0.25)

tit = 'ROMS+' + ogcm + vartit + time
fig = 'roms+' + ogcm + varprn + MON + '.png'

x = lon[::vint, ::vint]
y = lat[::vint, ::vint]
u = uvel[::vint, ::vint]
v = vvel[::vint, ::vint]

plt.figure()
cf = plt.contourf(lon, lat, temp, levels , cmap=cmap, alpha = 0.75)
cbar = plt.colorbar(cf); cbar.ax.set_ylabel(units)
plt.quiver(x, y, u, v,scale=vscl)
plt.title(tit); plt.xlabel('Longitude'); plt.ylabel('Latitude')
plt.fill(coast['lon'],coast['lat'],'g')
plt.axis((wesn[0], wesn[1], wesn[2], wesn[3]))
plt.savefig(fig,dpi=100)
plt.close()

vel = u*u + v*v
vel = np.sqrt(vel)
vartit = ' surface vels on '
varprn = '_sfc_vels_'
units = 'm/s'

tit = 'ROMS+' + ogcm + vartit + time
fig = 'roms+' + ogcm + varprn + MON + '.png'

levels = np.arange(0,1.6,0.1)
cmap = plt.get_cmap('bwr')

plt.figure()
cf = plt.contourf(x, y, vel, levels, cmap=cmap, alpha = 0.75)
cbar = plt.colorbar(cf); cbar.ax.set_ylabel(units)
plt.quiver(x, y, u, v,scale=vscl)
plt.title(tit); plt.xlabel('Longitude'); plt.ylabel('Latitude')
plt.plot(coast['lon'],coast['lat'],'.k')
plt.axis((wesn[0], wesn[1], wesn[2], wesn[3]))
plt.savefig(fig,dpi=100)
plt.close()


units = 'm'
vartit = ' surface elev on '
varprn = '_sfc_elev_'

levels = np.arange(-1.5,1.6,0.1)

tit = 'ROMS+' + ogcm + vartit + time
fig = 'roms+' + ogcm + varprn + MON + '.png'

x = lon[::vint, ::vint]
y = lat[::vint, ::vint]
u = ubar[::vint, ::vint]
v = vbar[::vint, ::vint]

vscl = vscl/2

plt.figure()
cf = plt.contourf(lon, lat, zeta, cmap=cmap, alpha = 0.75)
cbar = plt.colorbar(cf); cbar.ax.set_ylabel(units)
plt.quiver(x, y, u, v,scale=vscl)
plt.title(tit); plt.xlabel('Longitude'); plt.ylabel('Latitude')
plt.plot(coast['lon'],coast['lat'],'.k')
plt.axis((wesn[0], wesn[1], wesn[2], wesn[3]))
plt.savefig(fig,dpi=100)
plt.close()



###########################################################################################

print(' ')
print(' +++ End of python program +++ ')
print(' ')

###########################################################################################


