#-------Plot Chandra photometry data-------#
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from astropy.table import Table
import scipy, pylab
import math
import astropy
import pyregion
from matplotlib.lines import Line2D
from astropy.io import fits
import matplotlib.cm as cmx
import matplotlib.colors as colors

wittmann = Table.read('wittmann-2017.csv')
netcounts11714 = Table.read('netcounts/11714_netcounts.csv')
netcounts11715 = Table.read('netcounts/11715_netcounts.csv')
netcounts11713 = Table.read('netcounts/11713_netcounts.csv')

netcounts = np.append(netcounts11713['net_counts'], netcounts11714['net_counts']) 
netcounts = np.append(netcounts, netcounts11715['net_counts'])                                      # all netcount info

yerr = np.append(netcounts11713['net_err'], netcounts11714['net_err'])
yerr = np.append(yerr, netcounts11715['net_err'])

netcounts_300 = netcounts[0::2]     # separate by energy level (300kev-2000kev, 2000kev-7000kev)
netcounts_2000 = netcounts[1::2]
yerr_300 = yerr[0::2]
yerr_2000 = yerr[1::2]

idarray = np.append(netcounts11713['ID'], (netcounts11714['ID']))
idarray = np.append(idarray, (netcounts11715['ID']))

idarray
idnums_300 = []
for i in range(len(idarray)):
    idnum = str(idarray[i])
    num = idnum.split('.')[0] 
    idnums_300.append(num)             # get all of the investigated object IDs

idnums_300 = idnums_300[0::2]              # take every other id number (doubles of each)
idnums_300 = [int(n) for n in idnums_300]  # convert to integers


    
ra = np.zeros(len(idnums_300))         # initialize the ra and dec arrays for desired objects
dec = np.zeros(len(idnums_300))

wittra = wittmann['ra']
wittdec = wittmann['dec']
wittID = wittmann['ID']


for i in range(len(wittID)):          # get the associated ra and dec for each investigated object from the wittmann data
    for j in range(len(idnums_300)):
        if wittID[i] == idnums_300[j]:
            ra[j] = wittra[i]
            dec[j] = wittdec[i]
            
#-----compute distances from NGC1275-----
def degreedist(ra, dec):
    raNGC = math.radians(49.95041666)
    decNGC = math.radians(41.51138889)
 
    ra = math.radians(ra)
    dec = math.radians(dec)
    cosdist = math.sin(decNGC)*math.sin(dec) + math.cos(decNGC)*math.cos(dec)*math.cos(raNGC - ra)
    return math.degrees(math.acos(cosdist))

distarray = np.zeros(len(ra))
for i in range(len(ra)):
    distarray[i] = degreedist(ra[i], dec[i])



#------regression fitting-----

def line_model(pars, x):
	return pars[0]*np.array(x) + pars[1]

def weighted_squared_deviation(pars, x, y, y_err):
    chi = (y - line_model(pars, x)) / y_err
    return np.sum(chi**2)

_pars = [0, 0]
x = distarray
y_300 = netcounts_300
y_2000 = netcounts_2000


#print(".3-2 kev chi2 =", weighted_squared_deviation(_pars, x, y_300, yerr_300))
#print("2-7 kev chi2 =", weighted_squared_deviation(_pars, x, y_2000, yerr_2000))
chi2_300 = str(weighted_squared_deviation(_pars, x, y_300, yerr_300))
chi2_2000 = str(weighted_squared_deviation(_pars, x, y_2000, yerr_2000))
#------make the plot-----
ax1 = pylab.subplot(111)

ax1.scatter(distarray, netcounts_300, s=80, facecolors='none', edgecolors='r', label='$0.3-2.0 keV$')
ax1.scatter(distarray, netcounts_2000, s=60, c = 'b', marker = "s", label = '$2.0-7.0 keV$')

ax1.errorbar(distarray, netcounts_300, yerr=yerr_300, linestyle='none', ecolor='r')
ax1.errorbar(distarray, netcounts_2000, yerr=yerr_2000, linestyle='none', ecolor='b')

ax1.text(2, 6, r'$X^2_300=$chi2_300', fontsize=15)
#pylab.ylim([40.8, 42.0])
#pylab.xlim([])
plt.title('Chandra Photometry Data of Wittmann UDG objects')
plt.xlabel('distance from NGC1275 \n [deg]')
plt.ylabel('net counts \n [counts]')
plt.legend(loc='best')
plt.tight_layout()

plt.show()
    
