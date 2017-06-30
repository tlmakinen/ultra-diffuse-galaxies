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

# takes each file's number and returns a string filename #
def image_file(imagenum):                                  
    fname = 'regfiles/%s.reg' %imagenum

    return fname

# plot all the Wittmann and Song Perseus Cluster objects #
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
    x1275 = [49.95041666]
    y1275 = [41.51138889]
    ax1.plot(x1275,y1275, '+r')                             # plot NGC1275 (center of Perseus Cluster)

    return ax1

# Creates a list of patches for each Chandra observation #
def panelplot(axis, ds9number, panelarray, colornum):       
    ax1 = axis
    imagenum = ds9number                                    # which file are we looking at? 
    paneltags = []                                          # list of panel tags
    for j in range(len(panelarray)):
        paneltags.append("tag={%s}" %panelarray[j])         # creates a reference tag for each panel number in the input    
    
    r = pyregion.open(image_file(imagenum))                 # Use pyregion.open to get the information from the region file
    xyarrays = []                                           # stores the desireable panels' xy-coordinates from r[index]
    commentarray = []
    
    for i in range(len(r)):                                 # choose the panels we want to plot
        for j in range(len(paneltags)):
            if r[i].comment == paneltags[j]:
                xyarrays.append(r[i].coord_list)
                commentarray.append(r[i].comment)           #  Save each comment from selected panel (saves panel name)

    # colorarray = np.zeros(len(commentarray))               #  Associate a color with each panel name
    # for c in range(len(commentarray)):
    #     if "I3" in commentarray[c]:
    #         colorarray[c] = colornum
    #     elif "I2" in commentarray[c]:
    #         colorarray[c] = colornum
    #     elif "I1" in commentarray[c]:
    #         colorarray[c] = colornum
    #     elif "I0" in commentarray[c]:
    #         colorarray[c] = colornum
    #     elif "S3" in commentarray[c]:
    #         colorarray[c] = colornum                                                               
    patches = []                                            # initialize list of patches
    for i in range(len(xyarrays)):                          # hard-code the polynomials using four calculated panel points
        xy = xyarrays[i]
        panel = Polygon([[xy[0],xy[1]], [xy[2],xy[3]], [xy[4],xy[5]], [xy[6],xy[7]]], fill=False)  
        patches.append(panel)
    #if len(xyarrays) > 0:
    #    print(imagenum)
    return patches

 # takes each file's number and returns a string FITS filename #
def fits_file(imagenum):                                    
    if int(imagenum) < 1000:
        return 'FITSfiles/acisf00%s_repro_evt2.fits' %imagenum
    if int(imagenum) < 10000:
        return 'FITSfiles/acisf0%s_repro_evt2.fits' %imagenum
    else:
        return 'FITSfiles/acisf%s_repro_evt2.fits' %imagenum


# Create plot of all panels by exposure time #

chandrapanels = ['I1', 'I0', 'I2', 'I3', 'S3']
nums = [11713, 11715, 11714]
#nums = [int(n.split('.')[0]) for n in os.listdir('regfiles/') if n[-1]=='g']  # Plot all .reg fules in regfiles/
print(nums)


fig, ax1 = plt.subplots()
perseusplot(ax1)

figtitle = 'Chandra panels over Perseus data by exposure time [ksecs]'
t = fig.text(0.45, 0.97, figtitle, horizontalalignment='center')
                                                            
timearray = []                                              # list to store exposure times from each fits file
newnums = []
for i in range(len(nums)):                                
    fname = fits_file(nums[i])                              # open the designated fits file header
    hdulist = fits.open(fname) 
    hduheader = hdulist[1].header                           # get the second header for the data required

    exposuretime = hduheader['EXPOSURE']                    # get the exposure time from the FITS header
    colornum = exposuretime/1000                            # adjust exposure time to ksecs
    upperlimit = 140
    lowerlimit = 50                                          # set desired upper limit to exposure time to plot only panels in threshold
    if colornum > lowerlimit and colornum < upperlimit:                               
        newnums.append(nums[i])                             # use only obs IDs that fall under upper limit                            
        timearray.append(colornum)                          # append the color index array

colorindex = np.array(timearray)
index_norm = colorindex / colorindex.max()                  # Normalize color range to work with the colormap

for i in range(len(index_norm)):                            # plot all panel files specified
    print("exposure time of " + str(newnums[i]) + " = " + str(timearray[i]))
    rightcolor = index_norm[i]                                      # use normalized color index for our plots
    patches = panelplot(ax1, newnums[i], chandrapanels, rightcolor)    # call panelplot for each observation number
    edgecolors = [matplotlib.cm.jet(rightcolor)]                       
    p = PatchCollection(patches, cmap=matplotlib.cm.jet, alpha=0.2, facecolor=edgecolors, edgecolors=edgecolors) # collect all the patches together
    ax1.add_collection(p)                                    # add collection of panels to the plot

p.set_array(np.array(colorindex))                            # Set the colorbar, using the last "p" patch collection limits for the scale (same for all collections)
#p.set_clim([np.ma.min(colorindex), np.ma.max(colorindex)])   # Set the colorbar's min and max values (in ksecs)
p.set_clim([lowerlimit, upperlimit])
plt.colorbar(p, shrink=0.5)

plt.show()


    



