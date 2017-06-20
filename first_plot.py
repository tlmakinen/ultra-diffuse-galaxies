import numpy as np
import math
import matplotlib.pyplot as plt
import csv


def read_data(objnum):
    fname = '%s.csv' %objnum
    
    data = np.loadtxt(fname)

with open(filename) as f:
    reader = csv.reader(f)
    columns = next(reader)
    colmap = dict(zip(columns, range(len(columns))))
    
testarr = np.array(np.loadtxt(song-high-sb, delimeter=",", skiprows=1))
testarr

