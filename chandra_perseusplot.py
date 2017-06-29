import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from astropy.table import Table
import scipy, pylab
import math
import pylab
import astropy
import pyregion
from matplotlib.lines import Line2D
from astropy.io import fits
import matplotlib.cm as cmx
import matplotlib.colors as colors

def image_file(imagenum):  # takes each file's number and returns a string filename
    fname = 'regfiles/%s.reg' %imagenum

    return fname


def perseusplot(ax1):

    songhigh = Table.read('song-high-sb.csv')              # Takes in Whittmann's and Song's data
    songlow = Table.read('song-low-sb.csv')
    wittmann = Table.read('wittmann-2017.csv')
    songra = np.append(songhigh['ra'], songlow['ra'])      # add the Song files together
    songdec = np.append(songhigh['dec'], songlow['dec'])   
    #ax1 = pylab.subplot(111)
    ax1.scatter(songra, songdec, s=80, facecolors='none', edgecolors='r', label='Song')
    ax1.scatter(wittmann['ra'], wittmann['dec'], s=10, c = 'b', marker = "o", label = 'Wittmann')

    plt.legend(loc='upper right');

    pylab.ylim([40.8, 42.0])
    pylab.xlim([49.0, 51.0])
    plt.xlabel('$ra [deg]$')
    plt.ylabel('$dec [deg]$')
    plt.tight_layout()

    #plt.show()

    return ax1

def panelplot(axis, ds9number, panelarray, colornum):

    ax1 = axis
    imagenum = ds9number                                   # which file are we looking at?
    
    paneltags = []                                         # list of panel tags
    for j in range(len(panelarray)):
        paneltags.append("tag={%s}" %panelarray[j])        # creates a reference tag for each panel number in the input

        
    r = pyregion.open(image_file(imagenum))                # Use pyregion.open to get the information from the region file

    xyarrays = []                                          # stores the desireable panels' xy-coordinates from r[index]
    commentarray = []
    

    for i in range(len(r)):                                # choose the panels we want to plot
        for j in range(len(paneltags)):
            if r[i].comment == paneltags[j]:
                xyarrays.append(r[i].coord_list)
                commentarray.append(r[i].comment)          #  Save each comment from selected panel (saves panel name)


    # upperlimit = 100
    # values = range(upperlimit)
    # jet = cm = plt.get_cmap('jet')
    # cNorm = colors.Normalize(vmin=0, vmax=values[-1])
    # scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)
    # print(scalarMap.get_clim())

    # colorVal = scalarMap.to_rgba(colornum)



    colorarray = np.zeros(len(commentarray))               #  Associate a color with each panel name
    for c in range(len(commentarray)):
        if "I3" in commentarray[c]:
            colorarray[c] = colornum
        elif "I2" in commentarray[c]:
            colorarray[c] = colornum
        elif "I1" in commentarray[c]:
            colorarray[c] = colornum
        elif "I0" in commentarray[c]:
            colorarray[c] = colornum
        elif "S3" in commentarray[c]:
            colorarray[c] = colornum

       
    patches = []                                            # initialize list of patches
    #edgecolors = [matplotlib.cm.jet(x) for x in colorarray]
    for i in range(len(xyarrays)):                          # hard-code the polynomials using four calculated panel points
        xy = xyarrays[i]
        panel = Polygon([[xy[0],xy[1]], [xy[2],xy[3]], [xy[4],xy[5]], [xy[6],xy[7]]], fill=False)  
        patches.append(panel)


    #p = PatchCollection(patches, cmap=matplotlib.cm.jet, alpha=0.2, facecolor=edgecolors, edgecolors=edgecolors) # collect all the patches together
    #p = PatchCollection(patches)
    #colors = colorarray                                     
    #p.set_array(np.array(colors))                           # set all the colors

    
    #ax1.add_collection(p)                                   # add all our patches to the axes
    #plt.title(str(imagenum))                                # add what number we're looking at
    #plt.show()
    return patches

def fits_file(imagenum):  # takes each file's number and returns a string fits filename
    if int(imagenum) < 1000:
        return 'FITSfiles/acisf00%s_repro_evt2.fits' %imagenum
    if int(imagenum) < 10000:
        return 'FITSfiles/acisf0%s_repro_evt2.fits' %imagenum
    else:
        return 'FITSfiles/acisf%s_repro_evt2.fits' %imagenum

chandrapanels = ['I3', 'S3', 'I2', 'I1', 'I0']
nums = [502, 4949, 17274]

nums = [int(n.split('.')[0]) for n in os.listdir('regfiles/') if n[-1]=='g']  # Plot all .reg fules in regfiles/
#print(nums)


fig, ax1 = plt.subplots()
perseusplot(ax1)
timearray = []

for i in range(len(nums)):                               # plot all panel files specified
    fname = fits_file(nums[i])
    hdulist = fits.open(fname)
    hduheader = hdulist[1].header

    exposuretime = hduheader['EXPOSURE']                  # get the exposure time from the FITS header
    colornum = exposuretime                        # adjust exposure time to ksecs
                                                          # also adjusts for a [0,100] color map
    timearray.append(colornum)

colorindex = np.array(timearray)
index_norm = colorindex / colorindex.max()

for i in range(len(nums)):
    #print("exposure time of " + str(nums[i]) + " = " + str(exposuretime))
    rightcolor = index_norm[i]
    patches = panelplot(ax1, nums[i], chandrapanels, rightcolor)
    edgecolors = [matplotlib.cm.jet(rightcolor)]

    p = PatchCollection(patches, cmap=matplotlib.cm.jet, alpha=0.2, facecolor=edgecolors, edgecolors=edgecolors) # collect all the patches together
    #p.set_array()
    ax1.add_collection(p) 


    #plt.title(str(nums[loop]))




plt.show()


    



