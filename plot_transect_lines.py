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

gridfile = str(sys.argv[1])

grid = Dataset(gridfile,'r')

lat0 = [-22.0, -24.0, -26.0, -28.0, -30.0, -11.0]
lon1 = [-42.5, -46.5, -49.5, -50.0, -50.5, -37.0]
lon2 = [-33.0, -34.0, -35.0, -36.0, -37.0, -33.0]

N = len(lat0)

############ *** plot a map showing across shelf *** #####################################

coastfile = "costa_leste.mat"
coast = sio.loadmat(coastfile)

lo1 = grid.variables['lon_rho'][:]
la1 = grid.variables['lat_rho'][:]
prf = grid.variables['h'][:] 

plt.figure()
plt.contourf(lo1,la1,prf,cmap='jet')
plt.plot(coast['lon'],coast['lat'],'g')
plt.title('Seções para avaliação do transporte da Corrente do Brasil')
plt.xlabel('Longitude'); plt.ylabel('Latitude')
plt.grid()

for n in range(N):
    print(' ')
    print(' ... transect ' + str(-lat0[n]) + ' S')
    print(' ')
    lon = np.arange(lon1[n],lon2[n],0.05)
    lat = lat0[n] * np.ones((len(lon)))
    plt.plot(lon,lat,'.r')

plt.axis((-59,-25,-45,-5))
plt.savefig('transects.png',dpi=100)
plt.close()


###########################################################################################

print(' ')
print(' +++ End of python program +++ ')
print(' ')

###########################################################################################


