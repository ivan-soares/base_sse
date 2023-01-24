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

gridfile = str(sys.argv[2])
grid = Dataset(gridfile,'r')

mon = int(sys.argv[3])
year = int(sys.argv[4])

lat1 = -23.18
lat2 = -20.72
lon1 = -41.66
lon2 = -29.33

print(' ')
print(' ... will read roms file ' + romsfile )
print(' ')

MON = '{:02d}'.format(mon)
ntime = mon -1

minlon = -41.5
maxlon = -39.0
maxprf = 550.0

################## *** read roms velocity *** ############################################

lon = roms.variables['lon'][320:530,:]
lat = roms.variables['lat'][320:530,:]

time = get_roms_date2(romsfile,ntime)

roms_prf = roms.variables['depth'][:]
roms_u3d = roms.variables['u'][ntime,:,320:530,:]
roms_v3d = roms.variables['v'][ntime,:,320:530,:]
roms_nan = roms.variables['v'].__getattribute__('missing_value')

roms_u3d[np.where(roms_u3d < -9000.)] = np.nan
roms_v3d[np.where(roms_v3d < -9000.)] = np.nan

ang = 60. *np.pi/180.

####### formula to rotate u & v
#       uvel[k,:,:] = +u*cos(ang) + v*sin(ang)
#       vvel[k,:,:] = -u*sin(ang) + v*cos(ang)


nz, nlat, nlon = roms_v3d.shape
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

dy = (lat2 - lat1)/100.
dx = (lon2 - lon1)/100.

across_lat = np.arange(lat1,lat2,dy)
across_lon = np.arange(lon1,lon2,dx)
nlat_across = len(across_lat)
nlon_across = len(across_lon)

print(' ... nlat across ' + str(nlat_across))
print(' ... nlon across ' + str(nlon_across))

############ plot a map

coastfile = "costa_leste.mat"
coast = sio.loadmat(coastfile)

lo1 = grid.variables['lon_rho'][:]
la1 = grid.variables['lat_rho'][:]
prf = grid.variables['h'][:]

plt.figure()
plt.contourf(lo1,la1,prf,cmap='jet')
plt.fill(coast['lon'],coast['lat'],'g')
plt.plot(across_lon,across_lat,'*c')
plt.title('AX97 XBT line')
plt.xlabel('Longitude'); plt.ylabel('Latitude')
plt.axis([-55., -25., -42., 0.]); plt.grid()
plt.savefig('ax97.png',dpi=100)
plt.close()


roms_across_lon = np.zeros((nz, nlon_across))
roms_across_prf = np.zeros((nz, nlon_across))
roms_across_vel = np.zeros((nz, nlon_across))

for i in range(nlon_across):
    roms_across_prf[:,i] = roms_prf[:]

for n in range(nz):
    roms_across_lon[n,:] = across_lon[:]

for n in range(nz):
    roms_across_vel[n, :] = griddata(np.column_stack((lon.ravel(), lat.ravel())),
                            roms_vel[n, :, :].ravel(),(across_lon.ravel(), across_lat.ravel()))    

print(' ')
print(' ... across shelf transect lat ' + str(across_lat))
print('                                          ')
print(' ... across shelf transect lon ' + str(across_lon[np.where((across_lon >= minlon) & (across_lon <= maxlon))]))
print(' ')
print(' ... across shelf velocities ' + str(roms_across_vel))
print(' ')
print(' ... roms depths ' + str(roms_prf[np.where(roms_prf <= maxprf)]))
print(' ')
print(' ')

############## *** compute transport *** #################################################

I = across_lon[np.where((across_lon >= minlon) & (across_lon <= maxlon))]
J = roms_prf[np.where(roms_prf <= maxprf)]

NI = len(I) - 1
NJ = len(J) - 1

bc_vel = roms_across_vel[0:NJ, 0:NI]
bc_lon = roms_across_lon[0:NJ, 0:NI+1]
bc_prf = roms_across_prf[0:NJ+1, 0:NI]

BC_transp = transp(bc_vel, bc_lon, bc_prf,-22.)
BC = '{:.2f}'.format(BC_transp)

print(' ')
print(' ... BC transport on date ' + str(time) + ' is ' + str(BC_transp))
print(' ')

################# *** plot data *** ######################################################

cmap = plt.get_cmap('bwr')

units = 'm/s'
vartit = 'AX97 ' 
varprn = 'v-northward_across_ax97_' + str(year) + '_' + MON

isnan = np.isnan(roms_across_vel)

roms_across_vel[isnan] = 0.
minvel = min(roms_across_vel.flatten())
maxvel = max(roms_across_vel.flatten())
print(' ... min vel is ' + str(minvel))
print(' ... max vel is ' + str(maxvel))
roms_across_vel[isnan] = np.nan

levels = np.arange(minvel,maxvel,0.005)
levels = np.arange(-.6,.7,.1)

roms_across_vel[np.where(roms_across_vel >= +0.6)] = +0.6
roms_across_vel[np.where(roms_across_vel <= -0.6)] = -0.6


tit_name = vartit + time
fig_name = varprn + '.png'

fig, ax= plt.subplots(1)
fig.suptitle(tit_name)

im = ax.contourf(roms_across_lon, -roms_across_prf, roms_across_vel, levels, cmap=cmap)
fig.colorbar(im, shrink=0.6, label = 'm/s')

ax.set_ylim((-700., 0.))
ax.set_ylabel('Depth (m)')
ax.set_xlabel('Longitude')
ax.grid()

x1, x2, y1, y2 = ax.axis()
xlabel = x1 + (x2 - x1)/10
ylabel = y1 + (y2 - y1)/10

ax.text(xlabel, ylabel, 'CB transp ' + str(BC) + ' Sv', color = 'k')

fig.savefig(fig_name,dpi=100)
#fig.close()


###########################################################################################

print(' ')
print(' +++ End of python program +++ ')
print(' ')

###########################################################################################


