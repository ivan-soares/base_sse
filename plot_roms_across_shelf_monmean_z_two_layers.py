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
from math import cos, sin

#import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt

from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
from scipy.interpolate import interp1d, interp2d, griddata
from romstools import get_roms_date2, scoord, get_sigma_level_depths
from plot_ogcm import plot_contour_and_vectors



###########################################################################################

'''------------------------------ begin aux_function -----------------------------------'''

def transp(vel,lon,prf,lat0):
    dlon = np.diff(lon, axis = 1)
    dprf = np.diff(prf, axis = 0)
    #print(dlon)
    #print(dprf)
    twopir = 2.0 * np.pi * 6.371e+06
    dx = dlon * twopir * np.cos(lat0 * np.pi / 180.) / 360.
    dz = dprf
    V = vel * dx * dz
    if lat0 == -11.:
       V = sum(V[np.where(V > 0.)]) * 1.e-06
    else:
       V = sum(V[np.where(V < 0.)]) * 1.e-06
    return V


'''-------------------------------- end aux_function -----------------------------------'''

###########################################################################################

print(' ')
print(' +++ Starting python program +++ ')
print(' ')

################### *** argument input vars *** ###########################################

romsfile = str(sys.argv[1])
roms = Dataset(romsfile,'r')

mon = int(sys.argv[2])
year = int(sys.argv[3])

lat0 = float(sys.argv[4])
lon1 = float(sys.argv[5])
lon2 = float(sys.argv[6])

maxlon = float(sys.argv[7])
maxprf = float(sys.argv[8])
transect_name = str(sys.argv[9])

print(' ')
print(' ... will read roms file ' + romsfile )
print(' ')

MON = '{:02d}'.format(mon)
ntime = mon -1

################## *** read roms velocity *** ############################################

lon = roms.variables['lon'][:]
lat = roms.variables['lat'][:]

time = get_roms_date2(romsfile,ntime)

roms_prf = roms.variables['depth'][:]
roms_u3d = roms.variables['u'][ntime,:,:,:]
roms_v3d = roms.variables['v'][ntime,:,:,:]
roms_nan = roms.variables['v'].__getattribute__('missing_value')


roms_u3d[np.where(roms_u3d < -9000.)] = np.nan
roms_v3d[np.where(roms_v3d < -9000.)] = np.nan

ang = 36. *np.pi/180.

####### formula to rotate u & v
#       uvel[k,:,:] = +u*cos(ang) + v*sin(ang)
#       vvel[k,:,:] = -u*sin(ang) + v*cos(ang)

nz, nlat, nlon = roms_u3d.shape
roms_vel = np.zeros((nz,nlat,nlon))

for k in range(nz):
    roms_vel[k,:,:] = -roms_u3d[k,:,:]*sin(ang) + roms_v3d[k,:,:]*cos(ang)

roms_vel[np.where(roms_vel < -9000.)] = np.nan

print(' ')
print(' ... roms NaN ' + str(roms_nan) )
print('                               ')
print(' ...      n. lats ' + str(nlat) )
print(' ...      n. lons ' + str(nlon) )
print(' ')

################## *** set across shelf vars *** #########################################

across_lon = np.arange(lon1,lon2,0.05)
nlon_across = len(across_lon)

roms_across_lon = np.zeros((nz, nlon_across))
roms_across_prf = np.zeros((nz, nlon_across))
roms_across_vel = np.zeros((nz, nlon_across))

for i in range(nlon_across):
    roms_across_prf[:,i] = roms_prf[:]

for n in range(nz):
    roms_across_lon[n,:] = across_lon[:]


for n in range(nz):
    roms_across_vel[n, :] = griddata(np.column_stack((lon.ravel(), lat.ravel())),
                            roms_vel[n, :, :].ravel(),(across_lon.ravel(), lat0))    

print(' ')
print(' ... across shelf transect lat ' + str(lat0))
print('                                          ')
print(' ... across shelf transect lon ' + str(across_lon[np.where(across_lon <= maxlon)]))
print(' ')
print(' ... across shelf velocities ' + str(roms_across_vel))
print(' ')
print(' ... roms depths ' + str(roms_prf[np.where(roms_prf <= maxprf)]))
print(' ')
print(' ')

############## *** compute transport *** #################################################

I = across_lon[np.where(across_lon <= maxlon)]
J = roms_prf[np.where(roms_prf <= maxprf)]

NI = len(I) - 1
NJ = len(J) - 1

bc_vel = roms_across_vel[0:NJ, 0:NI]
bc_lon = roms_across_lon[0:NJ, 0:NI+1]
bc_prf = roms_across_prf[0:NJ+1, 0:NI]

BC_transp = transp(bc_vel, bc_lon, bc_prf, lat0)
BC = '{:.2f}'.format(BC_transp)

print(' ')
print(' ... BC transport on date ' + str(time) + ' is ' + str(BC_transp))
print(' ')

################# *** plot data *** ######################################################

cmap = plt.get_cmap('bwr')

units = 'm/s'
vartit = 'Across shelf section Lat ' + str(lat0) 
varprn = 'v-northward_across_lat_' + str(-lat0) + '_' + str(year) + '_' + MON

isnan = np.isnan(roms_across_vel)

roms_across_vel[isnan] = 0.
minvel = min(roms_across_vel.flatten())
maxvel = max(roms_across_vel.flatten())
print(' ... min vel is ' + str(minvel))
print(' ... max vel is ' + str(maxvel))
roms_across_vel[isnan] = np.nan

levels = np.arange(minvel,maxvel,0.005)

roms_vel[np.where(roms_vel == -0.222)] = np.nan

tit_name = vartit + ', ' + time
fig_name = varprn + '.png'

fig, axs= plt.subplots(2)
fig.suptitle(tit_name)

ax1 = axs[0]
ax2 = axs[1]

im = ax1.contourf(roms_across_lon, -roms_across_prf, roms_across_vel, levels, cmap=cmap)

ax1.set_ylim((-500., 0.))
ax1.set_ylabel('Depth (m)')
ax1.set_xticklabels([])
ax1.grid()

im = ax2.contourf(roms_across_lon, -roms_across_prf, roms_across_vel, levels, cmap=cmap)
fig.colorbar(im, ax = axs[:], shrink=0.6, label = 'm/s')

ax2.set_ylim((-4000., -500.))
ax2.set_xlabel('Longitude')
ax2.grid()

x1, x2, y1, y2 = ax2.axis()
xlabel = x1 + (x2 - x1)/10
ylabel = y1 + (y2 - y1)/10

if lat0 == -11.:
    ax2.text(xlabel, ylabel, 'CNB transp ' + str(BC) + ' Sv', color = 'k')
else:
    ax2.text(xlabel, ylabel, 'CB transp ' + str(BC) + ' Sv', color = 'k')

fig.savefig(fig_name,dpi=100)
#fig.close()


###########################################################################################

print(' ')
print(' +++ End of python program +++ ')
print(' ')

###########################################################################################


