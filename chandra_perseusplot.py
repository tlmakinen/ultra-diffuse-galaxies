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

def image_file(imagenum):  # takes each file's number and returns a string filename
    fname = 'regfiles/%s.reg' %imagenum
    return fname

def perseus_panelplot(ds9number, panelarray):

    songhigh = Table.read('song-high-sb.csv')              # Takes in Whittmann's and Song's data
    songlow = Table.read('song-low-sb.csv')
    wittmann = Table.read('wittmann-2017.csv')
    songra = np.append(songhigh['ra'], songlow['ra'])      # add the Song files together
    songdec = np.append(songhigh['dec'], songlow['dec'])   
    ax1 = pylab.subplot(111)
    ax1.scatter(songra, songdec, s=80, facecolors='none', edgecolors='r', label='Song')
    ax1.scatter(wittmann['ra'], wittmann['dec'], s=10, c = 'b', marker = "o", label = 'Wittmann')

    plt.legend(loc='upper right');

    pylab.ylim([40.8, 42.0])
    pylab.xlim([49.0, 51.0])
    plt.xlabel('$ra [deg]$')
    plt.ylabel('$dec [deg]$')
    plt.tight_layout()


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

    colorarray = np.zeros(len(commentarray))               #  Associate a color with each panel name
    for c in range(len(commentarray)):
        if "I3" in commentarray[c]:
            colorarray[c] = 25
        elif "I2" in commentarray[c]:
            colorarray[c] = 30
        elif "I1" in commentarray[c]:
            colorarray[c] = 45
        elif "I0" in commentarray[c]:
            colorarray[c] = 20
        elif "S3" in commentarray[c]:
            colorarray[c] = 10
       
    patches = []                                            # initialize list of patches

    for i in range(len(xyarrays)):                          # hard-code the polynomials using four calculated panel points
        xy = xyarrays[i]
        panel = Polygon([[xy[0],xy[1]], [xy[2],xy[3]], [xy[4],xy[5]], [xy[6],xy[7]]], True)  
        patches.append(panel)


    p = PatchCollection(patches, cmap=matplotlib.cm.jet, alpha=0.4)  # collect all the patches together

    colors = colorarray                                     
    p.set_array(np.array(colors))                           # set all the colors

    ax1.add_collection(p)                                   # add all our patches to the axes
    plt.show()

chandrapanels = ['I3', 'S3', 'I2', 'I0', 'I1']
nums = [17274, 502, 4949]

for f in range(len(nums)):
    perseus_panelplot(nums[f], chandrapanels)
