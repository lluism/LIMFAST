#!/usr/bin/env python
# SIMPLEST USAGE: python sliceplot.py -i file1 file2 file3...
#
#  Default slice is at a fixed y = DIM/2
#
# More complex options:
USAGE = "USAGE: python sliceplot.py [--savearray] [--zindex=<z-index of slice>] [--delzindex=<offset for subtracting image>] [--filter=<smoothing sigma>] [--filterx=<x-axis smoothing sigma>] [--filtery=<y-axis smoothing sigma>] [--filterz=<z-axis smoothing sigma>] [--min=<min of plot>] [--max=<max of plot>] -i <filename1> <filename2>..."

###### LIST OF OPTIONAL ARGUMENTS:
###  --savearray is a flag indicating you want to also save the 2D slice as a .npy file.
###  --zindex= lets you specify which array index cut through the z axis (usualy the LOS axis).
###  --delzindex= if this is specified, then we will plot the difference between slices at array[:,:,zindex] - array[:,:,zindex+delzindex]
###  --filter= allows you to smooth the array with a Gaussian filter with the specified standard deviation (in units of array cells).  DEFAULT is no smoothing.
###  --filterx= smooth only along the horizontal axis.  DEFAULT is no smoothing.
###  --filtery= smooth only along the vertical axis.  DEFAULT is no smoothing.
###  --filterz= smooth only along the line of sight axis.  DEFAULT is no smoothing.
###  --min= specify a minimum value for the plotted range
###  --max= specify a maximum value for the plotted range

import matplotlib
matplotlib.use('Agg')
from matplotlib.mlab import griddata
import matplotlib.pyplot as plt
import numpy as np
from pylab import *
import scipy
from scipy import ndimage
from matplotlib.ticker import *
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.colors import Normalize
from os.path import basename
import os
import sys, argparse


#To normalize the midpoint of the colorbar
class MidpointNormalize(Normalize):
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y))

def load_binary_data(filename, dtype=np.float32):
     """
     We assume that the data was written
     with write_binary_data() (little endian).
     """
     f = open(filename, "rb")
     data = f.read()
     f.close()
     _data = np.fromstring(data, dtype)
     if sys.byteorder == 'big':
       _data = _data.byteswap()
     return _data

# Parse the command line options
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', help='Input filenames', nargs='+', required=True)
parser.add_argument("-f", "--filter", type=float, default=-1, help="smooth the array with a Gaussian filter with the specified standard deviation (in units of array cells).  DEFAULT is no smoothing.")
parser.add_argument("-x", "--filterx", type=float, default=-1, help="smooth only along the horizontal axis.  DEFAULT is no smoothing.")
parser.add_argument("-y", "--filtery", type=float, default=-1, help="smooth only along the vertical axis.  DEFAULT is no smoothing.")
parser.add_argument("-z", "--filterz", type=float, default=-1, help="smooth only along the LOS (z) axis.  DEFAULT is no smoothing.")
parser.add_argument("--zindex", type=int, default=-1, help="specify which array index cut through the z axis (usualy the LOS axis).  DEFAULT is the midpoint, i.e. DIM/2.")
parser.add_argument("--delzindex", type=int, default=-1, help="if this is specified, then we will plot the difference between slices at array[:,:,zindex] - array[:,:,zindex+delzindex]")
parser.add_argument("--min", type=float, default=1e5, help="specify a minimum value for the plotting range.")
parser.add_argument("--max", type=float, default=-1e5, help="specify a maximum value for the plotting range.")
parser.add_argument("--savearray", help="flag indicating you want to also save the 2D slice as a .npy file.", action="store_true")

args = parser.parse_args()

z_pts = []
I_pts = []

# go through list of files and process each one
for path in args.input:

    print 'Processing input file:'
    print '  '+path
    filename="" + path.split("/")[-1]

    # lightcone?
    if basename(filename)[-11:]=='lighttravel':
        DIM = int("" + path.split("_")[-3])
        label=str("" + path.split("_")[-2])

    else:
        DIM = int("" + path.split("_")[-2])
        label=str("" + path.split("_")[-1])
        length = int(label[0:3])

    if args.zindex >= 0:
        z_index = args.zindex

    # read in the data cube located in 21cmFast/Boxes/delta_T*
    data1 = load_binary_data(path)
    data1.shape = (DIM, DIM, DIM)
    data1 = data1.reshape((DIM, DIM, DIM), order='F')

    # smooth the field?
    if args.filter >= 0:
        iso_sigma = args.filter
        print "Smoothing the entire cube with a Gassian filter of width="+str(iso_sigma)
        data1 = scipy.ndimage.filters.gaussian_filter(data1, sigma=iso_sigma)
    else:
        if args.filterx >= 0:
            x_sigma = args.filterx
            print "Smoothing along the x (horizontal) axis with a Gassian filter of width="+str(x_sigma)
            data1 = scipy.ndimage.filters.gaussian_filter1d(data1, sigma=x_sigma, axis=1)
        if args.filtery >= 0:
            y_sigma = args.filtery
            print "Smoothing along the y (vertical) axis with a Gassian filter of width="+str(y_sigma)
            data1 = scipy.ndimage.filters.gaussian_filter1d(data1, sigma=y_sigma, axis=0)
        if args.filterz >= 0:
            z_sigma = args.filterz
            print "Smoothing along the z (line of sight) axis with a Gassian filter of width="+str(z_sigma)
            data1 = scipy.ndimage.filters.gaussian_filter1d(data1, sigma=z_sigma, axis=2)


    fig = plt.figure(dpi=72)
    sub_fig = fig.add_subplot(111)

    # extract a slice from the 3D cube
    y_index = DIM/2
    if args.zindex < 0:
        print "Taking an xz slice at y index="+str(y_index)
        slice = data1[:,y_index,:]
        endstr = '_yindex'+str(y_index)
    else:
        print "Taking an xy slice at z index="+str(z_index)
        slice = data1[:,:,z_index]
        endstr = '_zindex'+str(z_index)
        if args.delzindex >= 0: #difference image is wanted
            del_z_index = args.delzindex
            other_z_index = int(z_index+del_z_index)
            print "Subtracting the slice at index="+str(other_z_index)
            slice = slice - data1[:,:,other_z_index]
            endstr = '_zindex'+str(z_index)+'-'+str(z_index+del_z_index)

    # check if it is an Ha box -GS
    if basename(filename)[8:15]=='HI6563A':
        if args.min > 1e4:
            minrange = -2 #to set later
        if args.max < -1e4:
            maxrange = 0 #to set later

        # Renormalize data
        #slice *= 3.67e10

        i_z = basename(filename).find('A_z') + 3
        z_read = float(basename(filename)[i_z:i_z+5])
        print "z_read:", z_read

        slice = np.log10(slice)
        print "Min: "+str(np.min(slice))
        print "Max: "+str(np.max(slice))

        print "data1.shape is:", data1.shape
        np.savez('./Ha_Box_Coarse_%.1f.npz'%z_read, **{'box': data1})
        mean_intensity = np.sum(data1)/np.prod(data1.shape)/1.0e-23
        print "Mean intensity:", mean_intensity, np.log10(mean_intensity)
        z_pts.append(z_read)
        I_pts.append(mean_intensity)
    # check if it is an Ha box -GS
    elif basename(filename)[8:17]=='OIII5007A':
        if args.min > 1e4:
            minrange = -2 #to set later
        if args.max < -1e4:
            maxrange = 0 #to set later

        # Renormalize data
        #slice *= 3.67e10

        i_z = basename(filename).find('A_z') + 3
        z_read = float(basename(filename)[i_z:i_z+5])
        print "z_read:", z_read

        slice = np.log10(slice)
        print "Min: "+str(np.min(slice))
        print "Max: "+str(np.max(slice))

        mean_intensity = np.sum(data1)/np.prod(data1.shape)/1.0e-23
        print "Mean intensity:", mean_intensity, np.log10(mean_intensity)
        z_pts.append(z_read)
        I_pts.append(mean_intensity)
    # check if it is an Ha box -GS
    elif basename(filename)[8:16]=='OII3727A':
        if args.min > 1e4:
            minrange = -2 #to set later
        if args.max < -1e4:
            maxrange = 0 #to set later

        # Renormalize data
        #slice *= 3.67e10

        i_z = basename(filename).find('A_z') + 3
        z_read = float(basename(filename)[i_z:i_z+5])
        print "z_read:", z_read

        slice = np.log10(slice)
        print "Min: "+str(np.min(slice))
        print "Max: "+str(np.max(slice))

        mean_intensity = np.sum(data1)/np.prod(data1.shape)/1.0e-23
        print "Mean intensity:", mean_intensity, np.log10(mean_intensity)
        z_pts.append(z_read)
        I_pts.append(mean_intensity)
    # check if it is an Ha box -GS
    elif basename(filename)[0:7]=='delta_T':
        if args.min > 1e4:
            minrange = -2 #to set later
        if args.max < -1e4:
            maxrange = 0 #to set later

        # Renormalize data
        #slice *= 3.67e10

        i_z = basename(filename).find('T_z') + 3
        z_read = float(basename(filename)[i_z:i_z+5])
        print "z_read:", z_read

        slice = np.log10(slice)
        print "Min: "+str(np.min(slice))
        print "Max: "+str(np.max(slice))

        np.savez('./deltaT_Box_Coarse_%.1f.npz'%z_read, **{'box': data1})
        mean_intensity = np.sum(data1)/np.prod(data1.shape)/1.0e-23
        print "Mean intensity:", mean_intensity, np.log10(mean_intensity)
        z_pts.append(z_read)
        I_pts.append(mean_intensity)
    # check if it is an SFRD box -GS
    elif basename(filename)[0:3]=='SFR':
        if args.min > 1e4:
            minrange = -2 #to set later
        if args.max < -1e4:
            maxrange = 0 #to set later

        # Renormalize data
        #slice *= 3.67e10

        i_z = basename(filename).find('D_z') + 3
        z_read = float(basename(filename)[i_z:i_z+5])
        print "z_read:", z_read

        slice = np.log10(slice)
        print "Min: "+str(np.min(slice))
        print "Max: "+str(np.max(slice))

        print "data1.shape is:", data1.shape
        np.savez('./SFRD_Box_Coarse_%.1f.npz'%z_read, **{'box': data1})
        mean_intensity = np.sum(data1)/np.prod(data1.shape)
        print "Mean SFRD:", mean_intensity, np.log10(mean_intensity)
        z_pts.append(z_read)
        I_pts.append(mean_intensity)
    # check if it is an SFRD box -GS
    elif basename(filename)[0:3]=='SMD':
        if args.min > 1e4:
            minrange = -2 #to set later
        if args.max < -1e4:
            maxrange = 0 #to set later

        # Renormalize data
        #slice *= 3.67e10

        i_z = basename(filename).find('D_z') + 3
        z_read = float(basename(filename)[i_z:i_z+5])
        print "z_read:", z_read

        slice = np.log10(slice)
        print "Min: "+str(np.min(slice))
        print "Max: "+str(np.max(slice))

        mean_intensity = np.sum(data1)/np.prod(data1.shape)
        print "Mean SMD:", mean_intensity, np.log10(mean_intensity)
        z_pts.append(z_read)
        I_pts.append(mean_intensity)
    # check if it is an SFRD box -GS
    elif basename(filename)[0:8]=='starmass':
        if args.min > 1e4:
            minrange = -2 #to set later
        if args.max < -1e4:
            maxrange = 0 #to set later

        # Renormalize data
        #slice *= 3.67e10

        i_z = basename(filename).find('s_z') + 3
        z_read = float(basename(filename)[i_z:i_z+5])
        print "z_read:", z_read

        slice = np.log10(slice)
        print "Min: "+str(np.min(slice))
        print "Max: "+str(np.max(slice))

        print "data1.shape is:", data1.shape
        np.savez('./Starmass_Box_Coarse.npz', **{'box': data1})
        mean_intensity = np.sum(data1)/np.prod(data1.shape)
        print "Mean Starmass:", mean_intensity, np.log10(mean_intensity)
        z_pts.append(z_read)
        I_pts.append(mean_intensity)
    # check if it is an SFRD box -GS
    elif basename(filename)[0:6]=='Fcoll_':
        if args.min > 1e4:
            minrange = -2 #to set later
        if args.max < -1e4:
            maxrange = 0 #to set later

        # Renormalize data
        #slice *= 3.67e10

        i_z = basename(filename).find('l_z') + 3
        z_read = float(basename(filename)[i_z:i_z+5])
        print "z_read:", z_read

        slice = np.log10(slice)
        print "Min: "+str(np.min(slice))
        print "Max: "+str(np.max(slice))

        print "data1.shape is:", data1.shape
        np.savez('./Fcoll_Box_Coarse_%.1f.npz'%z_read, **{'box': data1})
        mean_intensity = np.sum(data1)/np.prod(data1.shape)
        print "Mean Fcoll:", mean_intensity, np.log10(mean_intensity)
        z_pts.append(z_read)
        I_pts.append(mean_intensity)
    # check if it is an SFRD box -GS
    elif basename(filename)[0:8]=='Fcollfs_':
        if args.min > 1e4:
            minrange = -2 #to set later
        if args.max < -1e4:
            maxrange = 0 #to set later

        # Renormalize data
        #slice *= 3.67e10

        i_z = basename(filename).find('s_z') + 3
        z_read = float(basename(filename)[i_z:i_z+5])
        print "z_read:", z_read

        slice = np.log10(slice)
        print "Min: "+str(np.min(slice))
        print "Max: "+str(np.max(slice))

        print "data1.shape is:", data1.shape
        np.savez('./Fcollfs_Box_Coarse.npz', **{'box': data1})
        mean_intensity = np.sum(data1)/np.prod(data1.shape)
        print "Mean Fcollfs:", mean_intensity, np.log10(mean_intensity)
        z_pts.append(z_read)
        I_pts.append(mean_intensity)
    # check if it is an SFRD box -GS
    elif basename(filename)[0:3]=='met':
        if args.min > 1e4:
            minrange = -2 #to set later
        if args.max < -1e4:
            maxrange = 0 #to set later

        # Renormalize data
        #slice *= 3.67e10

        #try:
        #    i_z = basename(filename).find('y_z') + 3
        #except:
        #i_z = basename(filename).find('cm_z') + 4
        i_z = basename(filename).find('y_z') + 3
        #print 'i_z:', i_z
        z_read = float(basename(filename)[i_z:i_z+5])
        print "z_read:", z_read

        slice = np.log10(slice/0.014)
        print "Min: "+str(np.min(slice))
        print "Max: "+str(np.max(slice))

        mean_Z = np.sum(data1)/np.prod(data1.shape)
        print "Mean Z:", mean_Z/0.014, np.log10(mean_Z/0.014)
        z_pts.append(z_read)
        I_pts.append(mean_Z/0.014)
    else:
        raise NotImplementedError("Oops!")

print "We end up with z_pts: %s"%str(z_pts).format('%.1f')
print "We end up with I_pts: %s"%str(I_pts).format('%.3e')
print "DONE!"
