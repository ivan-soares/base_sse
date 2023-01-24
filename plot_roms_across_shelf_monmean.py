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
cmap = plt.get_cmap('jet')

print(' ... lon1 ' + str(lon[0]))
print(' ... lon2 ' + str(lon[-1]))


############## temp

units = 'deg C'
vartit = 'Acros shelf Temp @ lat ' + LAT 
varprn = 'temp_across_lat_' + LAT + '_' + MON + '_' + str(year)

levels = np.arange(0,30,0.25)

tit = vartit + ' ' + time
fig = varprn + '.png'

plt.figure()
cf = plt.contourf(alon, prf, temp, levels , cmap=cmap, alpha = 0.75)
cbar = plt.colorbar(cf); cbar.ax.set_ylabel(units)
plt.title(tit); plt.xlabel('Longitude'); plt.ylabel('Depth')
plt.savefig(fig,dpi=100)
plt.close()


############ v

units = 'm/s'
vartit = 'Across shelf V @ lat ' + LAT
varprn = 'v_across_lat_' + LAT + '_' + MON + '_' + str(year)

levels = np.arange(-1.,1.,0.05)

tit = vartit + ' ' + time
fig = varprn + '.png'

plt.figure()
cf = plt.contourf(alon, prf, velo, cmap=cmap, alpha = 0.95)
cbar = plt.colorbar(cf); cbar.ax.set_ylabel(units)
plt.title(tit); plt.xlabel('Longitude'); plt.ylabel('Depth')
plt.savefig(fig,dpi=100)
plt.close()

###########################################################################################

print(' ')
print(' +++ End of python program +++ ')
print(' ')

###########################################################################################


