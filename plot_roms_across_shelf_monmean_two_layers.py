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

from mpl_toolkits.axes_grid1 import make_axes_locatable

from romstools import get_roms_date, scoord, get_sigma_level_depths
from plot_ogcm import plot_contour_and_vectors


###########################################################################################

print(' ')
print(' +++ Starting python program +++ ')
print(' ')

################### *** read input files *** ##############################################

romsfile = str(sys.argv[1])
roms = Dataset(romsfile,'r')

mon = int(sys.argv[2])
year = int(sys.argv[3])

lat0 = int(sys.argv[4])
lon1 = int(sys.argv[5])
lon2 = int(sys.argv[6])

transect_name = str(sys.argv[7])

print(' ')
print(' ... will read roms file ' + romsfile )
print(' ')

MON = '{:02d}'.format(mon)
ntime = mon -1

########### *** read grid vars *** #######################################################

lon = roms.variables['lon_rho'][lat0,lon1:lon2]; nlon = len(lon)
lat = roms.variables['lat_rho'][lat0,lon1:lon2]; nlat = len(lat)
prf = roms.variables['h'][lat0,lon1:lon2]

print(' ')
print(' ... will plot cross shelf velocities @ lat ' + str(lat))
print(' ')

print(' ')
print(' ... will plot cross shelf velocities @ lon ' + str(lon))
print(' ')

################## *** get depths at rho points *** ######################################

class data:
      def __init__(self):
          pass
PARAM = data   # ROMS grid params

PARAM.spheri = 1
PARAM.vtrans = 2 
PARAM.vstret = 4
PARAM.thetas = 4.
PARAM.thetab = 4.
PARAM.tcline = 100
PARAM.hc     = 100
PARAM.ns     = 30

sig = get_sigma_level_depths(prf, PARAM)
prf = sig.z_r[:]

################## *** read roms file *** ################################################

time = get_roms_date(romsfile,ntime)

print(' ')
print(' ... date ' + str(time))
print(' ')

zeta = roms.variables['zeta'] [ntime,lat0,lon1:lon2]
temp = roms.variables['temp'] [ntime,:,lat0,lon1:lon2]
velo = roms.variables['v'][ntime,:,lat0,lon1:lon2]


ns = PARAM.ns
alon = np.zeros((ns,nlon))

for n in range(ns):
    alon[n,:] = lon[:]

print(' ... Lon is shaped ' + str(alon.shape))
print(' ... Prf is shaped ' + str(prf.shape))
print(' ... Vel is shaped ' + str(velo.shape))


############ *** plot a map showing across shelf *** #####################################

coastfile = "costa_leste.mat"
coast = sio.loadmat(coastfile)

lo1 = roms.variables['lon_rho'][:]
la1 = roms.variables['lat_rho'][:]
dep = roms.variables['h'][:] 
LAT = str(round(-lat[0],2))

plt.figure()
plt.contourf(lo1,la1,dep,cmap='jet')
plt.plot(lon,lat,'.r')
plt.plot(coast['lon'],coast['lat'],'g')
plt.axis((-60,-25,-45,-5))
plt.savefig('across_lat_' + LAT + '.png')
plt.close()


################## *** plot data *** ######################################################

levels = MaxNLocator(nbins=100).tick_values(10., 30.)
cmap = plt.get_cmap('bwr')

print(' ... lon1 ' + str(lon[0]))
print(' ... lon2 ' + str(lon[-1]))


############### temp
#
#units = 'deg C'
#vartit = 'Acros shelf Temp @ lat ' + LAT 
#varprn = 'temp_across_lat_' + LAT + '_' + MON + '_' + str(year)
#
#levels = np.arange(0,30,0.25)
#
#tit = vartit + ' ' + time
#fig = varprn + '.png'
#
#plt.figure()
#cf = plt.contourf(alon, prf, temp, levels , cmap=cmap, alpha = 0.75)
#cbar = plt.colorbar(cf); cbar.ax.set_ylabel(units)
#plt.title(tit); plt.xlabel('Longitude'); plt.ylabel('Depth')
#plt.savefig(fig,dpi=100)
#plt.close()


############ v

units = 'm/s'
vartit = 'Across shelf section @ lat ' + LAT
varprn = 'v_across_lat_' + LAT + '_' + str(year) + '_' + MON

#isnan =  roms.variables['v'].__getattribute__('missing_value')
#print(' ... nan is ' + str(isnan))
#np.where(velo > 1.e+30, 0., velo)

velo[velo.mask] = 0.

print(' ')
print(' ... velocities are: ' + str(velo))
print(' ... min velocity is ' + str(min(velo.flatten())))
print(' ... max velocity is ' + str(max(velo.flatten())))

levels = np.arange(min(velo.flatten()),max(velo.flatten()),0.005)
levels = np.arange(-0.5,0.505,0.005)

velo[np.where(velo > 0.5)] = 0.5
velo[np.where(velo < -0.5)] = -0.5

tit_name = vartit + ', ' + time
fig_name = varprn + '.png'

fig, axs = plt.subplots(2)
fig.suptitle(tit_name)

ax1 = axs[0]
ax2 = axs[1]

im = ax1.contourf(alon, prf, velo, levels, cmap=cmap, alpha = 0.95)
#fig.colorbar(im, ax = axs[:], shrink=0.6, label = 'm/s')

ax1.set_ylim((-500., 0.))
ax1.set_ylabel('Depth (m)')
ax1.set_xticklabels([])
ax1.grid()

im = ax2.contourf(alon, prf, velo, levels, cmap=cmap, alpha = 0.95)
fig.colorbar(im, ax = axs[:], shrink=0.6, label = 'm/s')

ax2.set_ylim((-4000., -500.))
ax2.set_xlabel('Longitude')
ax2.grid()

x1, x2, y1, y2 = ax2.axis()
xlabel = x1 + (x2 - x1)/10
ylabel = y1 + (y2 - y1)/10

ax2.text(xlabel, ylabel, transect_name, color = 'k')
fig.savefig(fig_name,dpi=100)

#plt.figure()
#cf = plt.contourf(alon, prf, velo, cmap=cmap, alpha = 0.95)
#cbar = plt.colorbar(cf); cbar.ax.set_ylabel(units)
#plt.title(tit); plt.xlabel('Longitude'); plt.ylabel('Depth')
#plt.savefig(fig,dpi=100)
#plt.close()

###########################################################################################

print(' ')
print(' +++ End of python program +++ ')
print(' ')

###########################################################################################


