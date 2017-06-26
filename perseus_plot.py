# perseus cluster plot to incorporate CHANDRA panel data

import numpy as np
import math
import matplotlib.pyplot as plt
import pylab
from csv import reader
import astropy
import pyregion

def create_panel(fname):
    r1 = pyregion.open(fname) # returns shapelist object
    patch_list , artist_list = r1.get_mpl_patches_texts()
    return patch_list,artist_list

def plotperseus(ds9fname):
    
    from astropy.table import Table
    import scipy, pylab

    songhigh = Table.read('song-high-sb.csv')
    songlow = Table.read('song-low-sb.csv')
    wittmann = Table.read('wittmann-2017.csv')

    ax1 = pylab.subplot(111)

    ax1.scatter(songlow['ra'], songlow['dec'], s=80, facecolors='none', edgecolors='r', label='Song_low')
    ax1.scatter(songhigh['ra'], songhigh['dec'], s=80, facecolors='none', edgecolors='r', label='Song_high')
    ax1.scatter(wittmann['ra'], wittmann['dec'], s=10, c = 'b', marker = "o", label = 'Wittmann')

                    
    plt.legend(loc='upper right');
                    

    pylab.ylim([40.8, 42.0])
    pylab.xlim([49.0, 51.0])
    plt.xlabel('$ra [deg]$')
    plt.ylabel('$dec [deg]$')
    plt.tight_layout()
    
    import matplotlib.patches as patches
    
    plist,artlist = create_panel(ds9fname)
    for p in plist:
        ax1.add_patch(p)
    for t in artlist:
        ax1.add_artist(t)

    plt.show()
plotperseus("regfiles/12037.reg")