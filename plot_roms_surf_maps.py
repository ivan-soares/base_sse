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

offset_file = str(sys.argv[2])
roms_offset = Dataset(offset_file,'r')

ntime = int(sys.argv[3])
vint = int(sys.argv[4])
vscl = float(sys.argv[5])

lon1 = int(sys.argv[6])
lon2 = int(sys.argv[7])
lat1 = int(sys.argv[8])
lat2 = int(sys.argv[9])

lon = roms.variables['lon_rho'][lat1:lat2,lon1:lon2]
lat = roms.variables['lat_rho'][lat1:lat2,lon1:lon2]

nlat, nlon = lon.shape

print(' ')
print(' ... nlats = ' + str(nlat))
print(' ... nlons = ' + str(nlon))
print(' ')

time = get_roms_date(romsfile,ntime)

zeta = roms.variables['zeta'] [ntime,lat1:lat2,lon1:lon2]
temp = roms.variables['temp'] [ntime,0,lat1:lat2,lon1:lon2]
u = roms.variables['u'][ntime,0,lat1:lat2,lon1:lon2-1]
v = roms.variables['v'][ntime,0,lat1:lat2-1,lon1:lon2]

ang = -1. * roms.variables['angle'] [lat1:lat2,lon1:lon2]
offset = roms_offset.variables['zeta'] [0,lat1:lat2,lon1:lon2]
zeta = zeta - offset

zeta[np.where(zeta > 1.e+30)]=0.0
u[np.where(u > 1.e+30)]=0.0
v[np.where(v > 1.e+30)]=0.0

nlat_u, nlon_u = u.shape
nlat_v, nlon_v = v.shape

print(' ')
print(' ... u nlats = ' + str(nlat_u))
print(' ... u nlons = ' + str(nlon_u))
print(' ')

print(' ')
print(' ... v nlats = ' + str(nlat_v))
print(' ... v nlons = ' + str(nlon_v))
print(' ')

## if nlat_u is not equal to nlat
## then U and V are in a staggered grid
## and must be interpolated to the same
## grid points as zeta and temp

if (nlat_u != nlat) or (nlat_v != nlat):
    tmp_u = 0.5 * (u[:,0:nlon_u-1] + u[:,1:nlon_u])
    tmp_v = 0.5 * (v[0:nlat_v-1,:] + v[1:nlat_v,:])
    col0 = np.zeros((tmp_u.shape[0],1))
    row0 = np.zeros((1,tmp_v.shape[1]))
    tmp_u = np.hstack((col0,tmp_u,col0))
    tmp_v = np.vstack([row0,tmp_v,row0])
    print(str(tmp_u.shape))
    print(str(tmp_v.shape))
else:
    tmp_u = u
    tmp_v = v

## U and V are rotated according to grid rotation angle
## if the grid is not rotated then the angle is zero and the formula below
## has no effect since sin(0) = 0 and cos(0) = 1.

uvel = +tmp_u*np.cos(ang) + tmp_v*np.sin(ang)
vvel = -tmp_u*np.sin(ang) + tmp_v*np.cos(ang)

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

wesn = [-46, -33, -30, -13]


units = 'deg C'
vartit = ' sea surface temperatute on '
varprn = '_sfc_temp_'

levels = np.arange(15.,30.25,0.25)

tit = 'ROMS' + vartit + time
fig = 'roms' + varprn + time + '_bcampos.png'

x = lon[::vint, ::vint]
y = lat[::vint, ::vint]
u = uvel[::vint, ::vint]
v = vvel[::vint, ::vint]

plt.figure()
cf = plt.contourf(lon, lat, temp, levels , cmap=cmap, alpha = 0.75)
cbar = plt.colorbar(cf); cbar.ax.set_ylabel(units)
plt.quiver(x, y, u, v,scale=vscl)
plt.title(tit); plt.xlabel('Longitude'); plt.ylabel('Latitude')
#$wesn = plt.axis()
plt.fill(coast['lon'],coast['lat'],'g')
plt.axis((wesn[0], wesn[1], wesn[2], wesn[3]))
plt.savefig(fig,dpi=100)
plt.close()

############ plot SLA

units = 'm'
vartit = ' surface elevation on '
varprn = '_sfc_elev_'

cmap = plt.get_cmap('bwr')

zeta[np.where(zeta <= -1.)] = -1.
zeta[np.where(zeta >= +1.)] = +1.

print(' ')
print(' ... zeta min val ' + str(np.min(zeta)))
print(' ... zeta max val ' + str(np.max(zeta)))
print(' ... zeta mean val ' + str(np.mean(zeta)))
print(' ')

zeta = zeta - np.mean(zeta)

print(' ')
print(' ... zeta min val ' + str(np.min(zeta)))
print(' ... zeta max val ' + str(np.max(zeta)))
print(' ... zeta mean val ' + str(np.mean(zeta)))
print(' ')

levels = np.arange(-1.,1.01,0.01)

tit = 'ROMS' + vartit + time
fig = 'roms' + varprn + time + '_bcampos.png'

x = lon[::vint, ::vint]
y = lat[::vint, ::vint]
u = uvel[::vint, ::vint]
v = vvel[::vint, ::vint]

plt.figure()
cf = plt.contourf(lon, lat, zeta, cmap=cmap, alpha = 0.75)
cbar = plt.colorbar(cf); cbar.ax.set_ylabel(units)
plt.quiver(x, y, u, v,scale=vscl)
plt.title(tit); plt.xlabel('Longitude'); plt.ylabel('Latitude')
#wesn = plt.axis()
plt.fill(coast['lon'],coast['lat'],'g')
plt.axis((wesn[0], wesn[1], wesn[2], wesn[3]))
plt.savefig(fig,dpi=100)
plt.close()



###########################################################################################

print(' ')
print(' +++ End of python program +++ ')
print(' ')

###########################################################################################


