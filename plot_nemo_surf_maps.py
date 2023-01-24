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

from romstools import get_nemo_date
from plot_ogcm import plot_contour_and_vectors

###########################################################################################

print(' ')
print(' +++ Starting python program +++ ')
print(' ')

################### *** read input files *** ##############################################

nemofile = str(sys.argv[1])
nemo = Dataset(nemofile,'r')

vint = int(sys.argv[2])
vscl = float(sys.argv[3])

W = int(sys.argv[4])
E = int(sys.argv[5])
S = int(sys.argv[6])
N = int(sys.argv[7])

lon = nemo.variables['longitude'][:]
lat = nemo.variables['latitude'][:]

lon, lat = np.meshgrid(lon, lat)

nt, nlat, nlon = nemo.variables['zos'].shape

print(' ')
print(' ... ntimes = ' + str(nt))
print(' ... zeta nlats = ' + str(nlat))
print(' ... zeta nlons = ' + str(nlon))
print(' ')

time = get_nemo_date(nemofile)

zeta = nemo.variables['zos'] [0,:,:]
temp = nemo.variables['thetao'] [0,0,:,:]
uvel = nemo.variables['uo'][0,0,:,:]
vvel = nemo.variables['vo'][0,0,:,:]

#zeta[np.where(zeta > 1.e+30)]=0.0
#u[np.where(u > 1.e+30)]=0.0
#v[np.where(v > 1.e+30)]=0.0

nt, nz, nlat_u, nlon_u = nemo.variables['uo'].shape
nt, nz, nlat_v, nlon_v = nemo.variables['vo'].shape

print(' ')
print(' ... u nlats = ' + str(nlat_u))
print(' ... u nlons = ' + str(nlon_u))
print(' ')


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


############## plot TEMP

wesn = [W, E, S, N]


units = 'deg C'
vartit = ' surface temp on '
varprn = '_sfc_temp_'

levels = np.arange(15.,30.25,0.25)

tit = 'NEMO' + vartit + time
fig = 'nemo' + varprn + time + '_bcampos.png'

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

############ plot SLA

units = 'm'
vartit = ' surface elev on '
varprn = '_sfc_elev_'

cmap = plt.get_cmap('bwr')

zeta[np.where(zeta <= -1.)] = -1.
zeta[np.where(zeta >= +1.)] = +1.

zeta = zeta - np.mean(zeta)

levels = np.arange(-1.,1.01,0.01)

tit = 'ROMS ' + vartit + time
fig = 'nemo+' + varprn + time + '_bcampos.png'

x = lon[::vint, ::vint]
y = lat[::vint, ::vint]
u = uvel[::vint, ::vint]
v = vvel[::vint, ::vint]

#vscl = vscl/2

plt.figure()
cf = plt.contourf(lon, lat, zeta, cmap=cmap, alpha = 0.75)
cbar = plt.colorbar(cf); cbar.ax.set_ylabel(units)
plt.quiver(x, y, u, v,scale=vscl)
plt.title(tit); plt.xlabel('Longitude'); plt.ylabel('Latitude')
plt.fill(coast['lon'],coast['lat'],'g')
plt.axis((wesn[0], wesn[1], wesn[2], wesn[3]))
plt.savefig(fig,dpi=100)
plt.close()



###########################################################################################

print(' ')
print(' +++ End of python program +++ ')
print(' ')

###########################################################################################


