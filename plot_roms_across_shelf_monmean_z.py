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

from romstools import get_roms_date2, scoord, get_sigma_level_depths
from plot_ogcm import plot_contour_and_vectors

###########################################################################################

'''------------------------------ begin aux_function -----------------------------------'''

def transp(romsfile, lat0, lon1, lon2, nt):

    roms = Dataset(romsfile,'r')

    vel = roms.variables['v'][nt,:,lat0,lon1:lon2]
    lat = roms.variables['lat'][lat0,lon1]

    dlon = 0.03999
    twopir = 2.0 * np.pi * 6.371e+06
    dx = dlon * twopir * np.cos(lat * np.pi / 180.) / 360.
    dz = 50.

    print(' ... lat is ' + str(lat))
    print(' ...  dx is ' + str(dx))
    print(' ')

    V = vel * dx * dz

    print(' ... V TRANSP is shape ' + str(V.shape))
    print(' ')

    if lat > -13.:
       V[np.where(V < 0.)] = 0.
       T = sum(sum(V)) * 1.e-06
    else:
       V[np.where(V > 0.)] = 0. 
       T = sum(sum(V)) * 1.e-06

    print(' ... T = ' + str(T))
    print(' ')   

    return V, T


'''-------------------------------- end aux_function -----------------------------------'''

###########################################################################################

print(' ')
print(' +++ Starting python program +++ ')
print(' ')


################### *** read input files *** ##############################################

romsfile = str(sys.argv[1])
gridfile = str(sys.argv[2])

roms = Dataset(romsfile,'r')
grid = Dataset(gridfile,'r')

mon = int(sys.argv[3])
year = int(sys.argv[4])

lat0 = int(sys.argv[5])
lon1 = int(sys.argv[6])
lon2 = int(sys.argv[7])

transect_name = str(sys.argv[8])

print(' ')
print(' ... will read roms file ' + romsfile )
print(' ... will read grid file ' + gridfile )
print(' ')

MON = '{:02d}'.format(mon)
ntime = mon -1

########### *** read grid vars *** #######################################################

lon = roms.variables['lon'][lat0,lon1:lon2]; nlon = len(lon)
lat = roms.variables['lat'][lat0,lon1:lon2]; nlat = len(lat)
prf = roms.variables['depth'][:]; ns = len(prf)

print(' ')
print(' ... will compute BC transport @ lat ' + str(lat[0]))
print(' ')

################## *** get depths at rho points *** ######################################


################## *** read roms file *** ################################################

time = get_roms_date2(romsfile,ntime)

print(' ')
print(' ... date ' + str(time))
print(' ')

zeta = roms.variables['zeta'] [ntime,lat0,lon1:lon2]
temp = roms.variables['temp'] [ntime,:,lat0,lon1:lon2]
velo = roms.variables['v'][ntime,:,lat0,lon1:lon2]

alon = np.zeros((ns,nlon))
aprf = np.zeros((ns,nlon))

for n in range(ns):
    alon[n,:] = lon[:]

for n in range(nlon):
    aprf[:,n] = prf[:]

print(' ... Lon is shaped ' + str(alon.shape))
print(' ... Prf is shaped ' + str(aprf.shape))
print(' ... Vel is shaped ' + str(velo.shape))


############# *** compute transport *** ##################################################

BC_Vzao, BC_transp = transp(romsfile, lat0, lon1, lon2, ntime)

BC = '{0:.2f}'.format(BC_transp)
     
print(' ')
print(' ... BC transport on date ' + str(time) + ' is ' + str(BC_transp))
print(' ')

############ *** plot a map showing across shelf *** #####################################

coastfile = "costa_leste.mat"
coast = sio.loadmat(coastfile)

lo1 = roms.variables['lon'][:]
la1 = roms.variables['lat'][:]
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

print(' ')
print(' ... lon1 ' + str(lon[0]))
print(' ... lon2 ' + str(lon[-1]))
print(' ')

label_xpos = lon[0] + (lon[-1] - lon[0] )/10.
label_ypos = -(np.max(prf) - np.max(prf)/10.)
label_ypos2 = -(np.max(prf) - np.max(prf)/25.)

print(' ')
print(' ... label xpos ' + str(label_xpos))
print(' ... label ypos ' + str(label_ypos))
print(' ... label ypos2 ' + str(label_ypos2))
print(' ')

############## temp

units = 'deg C'
vartit = 'Acros shelf Temp @ lat ' + LAT 
varprn = 'temp_across_lat_' + LAT + '_' + MON + '_' + str(year)

levels = np.arange(0,30,0.25)

tit = vartit + ' ' + time
fig = varprn + '.png'

plt.figure()
cf = plt.contourf(alon, -aprf, temp, levels , cmap=cmap, alpha = 0.75)
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
cf = plt.contourf(alon, -aprf, velo, cmap=cmap, alpha = 0.95)
cbar = plt.colorbar(cf); cbar.ax.set_ylabel(units)
plt.title(tit); plt.xlabel('Longitude'); plt.ylabel('Depth')

if lat[0] < -13.:
    plt.text(label_xpos, label_ypos,  transect_name, color = 'k')
    plt.text(label_xpos, label_ypos2, 'CB transp ' + str(BC) + ' Sv', color = 'k')
else:
    plt.text(label_xpos, label_ypos,  transect_name, color = 'k')
    plt.text(label_xpos, label_ypos2, 'CNB transp ' + str(BC) + ' Sv', color = 'k')

#plt.tight_layout()


plt.savefig(fig,dpi=100)
plt.close()

############ vzao

units = 'm³/s'
vartit = 'Across shelf Vzao @ lat ' + LAT
varprn = 'vzao_across_lat_' + LAT + '_' + MON + '_' + str(year)

levels = np.arange(-1.,1.,0.05)

tit = vartit + ' ' + time
fig = varprn + '.png'

plt.figure()
cf = plt.contourf(alon, -aprf, BC_Vzao, cmap=cmap, alpha = 0.95)
cbar = plt.colorbar(cf); cbar.ax.set_ylabel(units)
plt.title(tit); plt.xlabel('Longitude'); plt.ylabel('Depth')

if lat0 < -13.:
    plt.text(label_xpos, label_ypos,  transect_name, color = 'k')
    plt.text(label_xpos, label_ypos2, 'CB transp ' + str(BC) + ' Sv', color = 'k')
else:
    plt.text(label_xpos, label_ypos,  transect_name, color = 'k')
    plt.text(label_xpos, label_ypos2, 'CNB transp ' + str(BC) + ' Sv', color = 'k')

#plt.tight_layout()


plt.savefig(fig,dpi=100)
plt.close()


###########################################################################################

print(' ')
print(' +++ End of python program +++ ')
print(' ')

###########################################################################################


