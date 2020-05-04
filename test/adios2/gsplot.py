#!/usr/bin/env python3
from __future__ import absolute_import, division, print_function, unicode_literals
import adios2
import argparse
from mpi4py import MPI
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.cm as cm
#import decomp
import time
import os


def SetupArgs():
    parser = argparse.ArgumentParser()
    parser.add_argument("--instream", "-i", help="Name of the input stream", required=True)
    parser.add_argument("--outfile", "-o", help="Name of the output file", default="screen")
    #parser.add_argument("--varname", "-v", help="Name of variable read", default="tune")
    parser.add_argument("--displaysec", "-dsec", help="Float representing gap between plot window refresh", default=0.1)
    parser.add_argument("--nx", "-nx", help="Integer representing process decomposition in the x direction",default=1)
    parser.add_argument("--ny", "-ny", help="Integer representing process decomposition in the y direction",default=1)
    parser.add_argument("--nz", "-nz", help="Integer representing process decomposition in the z direction",default=1)
    parser.add_argument("--plane", "-plane", help="The 2D plane to be displayed/stored xy/yz/xz/all", default='yz')
    args = parser.parse_args()

    args.displaysec = float(args.displaysec)
    args.nx = int(args.nx)
    args.ny = int(args.ny)
    args.nz = int(args.nz)

    if args.plane not in ('xz', 'yz', 'xy', 'all'):
        raise "Input argument --plane must be one of xz/yz/xy/all"

    return args



def read_data(args, fr, start_coord, size_dims):
    var1 = args.varname
    data= fr.read(var1, start_coord, size_dims )
    data = np.squeeze(data)
    return data

def plotMe(t1, t2, step):
    i = step % N;
    color = cmap(float(i)/(N-1))
    line, = axes.plot(t1, t2, 'o')
    line.set_color(color)
    if (step == 0):
        axes.set_ylim(min(t2)*0.99, max(t2)*1.01)
        axes.set_xlim(min(t1)*0.99, max(t1)*1.01)
    plt.title(str(i))
    plt.draw()
    plt.show();
    plt.pause(2.1)



if __name__ == "__main__":
    # fontsize on plot
    fontsize = 24

    args = SetupArgs()

##@
    plt.ion()
    fig = plt.figure(figsize=(8,4))
    axes = fig.add_subplot(111)
    data_plot=plt.plot(0,0)
    #cmap = cm.get_cmap('viridis')
    cmap = cm.get_cmap('tab20')
    N = 20 # every 20 steps color will reoccur
##@

    # Read the data from this object
    fr = adios2.open(args.instream, "r", MPI.COMM_WORLD, "adios2_config.xml", "TuneFoot")
#    vars_info = fr.availablevariables()

    # Read through the steps, one at a time
    plot_step = 0
    for fr_step in fr:
#        if fr_step.current_step()
        cur_step= fr_step.current_step()
        vars_info = fr.available_variables()
#        print (vars_info)        
        tune_name1 = "tune1"; 
        tune_name2 = "tune2"
        numPtls = vars_info[tune_name1]["Shape"].split(',')
        print("check numPts=", numPtls)
        
        tune1_data= fr.read(tune_name1)
        tune2_data= fr.read(tune_name2)
        print (tune1_data[0:3], tune2_data[0:3])
        plotMe(tune1_data, tune2_data, plot_step);
        #if args.plane in ('xy', 'all'):
        #    data = read_data (args, fr_step, [0,0,int(shape3[2]/2)], [shape3[0],shape3[1],1])
        #    Plot2D ('xy', data, args, fullshape, sim_step, fontsize)

        plot_step = plot_step + 1;

    fr.close()

