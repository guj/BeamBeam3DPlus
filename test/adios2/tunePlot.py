#!/usr/bin/env python3
from __future__ import absolute_import, division, print_function, unicode_literals
import adios2
import argparse
from mpi4py import MPI
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.cm as cm
import collections
#import decomp
import time
import os
from math  import gcd

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

def plotHistogram(t1, t2):
    d = zip(t1,t2)
    counter =  collections.Counter(d)

    #numColors =  max(counter.values())+1
    #colormap = plt.cm.rainbow
    #colorst = [colormap(i) for i in np.linspace(0, 0.9, numColors)]

    cc=[]
    for x,y in zip(t1,t2):
        c = counter [(x,y)]
        cc.append(c)
    #axesHist.scatter(t1,t2, c=cc, s=ss, cmap=plt.cm.rainbow)
    pic=axesHist.scatter(t1,t2, s=2, c=cc, cmap=plt.cm.plasma_r)

    
    #cb = plt.colorbar(pic);

    cb = fig.colorbar(pic, ax=axesHist);
    return cb

def plotDiff(t1, t2, step):
    d = zip(t1,t2)
    counter =  collections.Counter(d)

    cc=[0]*len(t1)

    if (step > 0):
        delta1 = (t1-pTune1)**2;
        delta2 = (t2-pTune2)**2;
        cc = delta1 + delta2;

    pic=axesDiff.scatter(t1,t2, s=2, c=cc, cmap=plt.cm.plasma_r)


    pTune1[:] = t1
    pTune2[:] = t2

    return fig.colorbar(pic, ax=axesDiff);


def plotMeScatter(t1,t2,step):
    if (step == 0):
        axesHist.set_xlim(xmin_max[0], xmin_max[1]);
        axesHist.set_ylim(ymin_max[0], ymin_max[1]);

        axesDiff.set_ylim(ymin_max[0], ymin_max[1]);
        axesDiff.set_xlim(xmin_max[0], xmin_max[1]);
        plotLines(a1,b1,a2,b2)

    cb1 = plotHistogram(t1,t2)
    cb2 = plotDiff(t1,t2, step)
    plt.savefig("figs/"+str(step)+".png")
    plt.show();
    plt.pause(5)
    cb1.remove()
    cb2.remove()

def toMovie():
    os.system("ffmpeg -r 1 -i figs/%01d.png -vcodec mpeg4 -y figs/movie.mp4")


def index(m, n):
    return m+n

def findTextTip(m, n, c, xlim, ylim):
    if (m == 0):
        return [xlim[0], 0, 0]
    plotSlope = (ylim[1] - ylim[0])/(xlim[1] - xlim[0])
    if (m/n > plotSlope):
        xTop = (c - n*ylim[1])/m
        xTop += (xlim[1] - xTop)*0.1
        xStep = (xTop - xlim[0])*0.1
        return [xTop-xStep, (c-m*xTop)/n, 30];
    else:
        yTop = (c - m*xlim[0])/n
        yStep = (yTop - ylim[0])*0.1
        yTop -= yStep
        return [(c - n*yTop)/m, yTop-yStep, 0]


# input points are (a1/b1, a2/b2) 
# so mb1b2*x + nb1b2*y = ma1b2 + na2b
    #xmin_max = [0.98*a1/b1, 1.03*a1/b1]
    #ymin_max = [a2/b2, 1.067*a2/b2]
    #axesHist.set_xlim(xmin_max[0], xmin_max[1]);
    #axesHist.set_ylim(ymin_max[0], ymin_max[1]);

def plotLines(a1, b1, a2, b2):
    result = [];
    myInts = range(0,15);
    for m in myInts:
        for n in range(0, 15-m):
            if (m + n == 0):
                continue
            constant = m*a1*b2 + n*a2*b1
            tmp1 = gcd(m*b1*b2, constant)
            tmp2 = gcd(n*b1*b2, constant)
            tt = gcd(tmp1, tmp2)
            idx = index(m*b1*b2/tt, n*b1*b2/tt);
            if (idx <= 15):
                result.append((m*b1*b2/tt, n*b1*b2/tt))

    final = list(set(result))
    #for (m,n) in final:
    #    if (n == 0):
    #        print(m,n);
    #    else:
    #        print (m, n, np.arctan(m/n)*180/3.1415926)
    # start drawing:
    x = np.linspace(-1,1,30)

    x0 = a1/b1
    y0 = a2/b2
    for (m,n) in final:
        if (n != 0):
            y = (m*x0 + n*y0 - m*x)/n
            slope = -m*1.0/n*1.0
            equation = str(m)+"x + "+str(n)+"y="+str(m*x0+n*y0)
            angle = np.arctan(slope)*180/3.1415926
            #print ("m=",m, "n=",n, "=>", equation, angle)
            axesHist.plot(x,y, '-r')
            xmin_max[1] = a1/b1
            tip=findTextTip(m,n, m*x0+n*y0, xmin_max, ymin_max)
            text=str(int(m))+"x+"+str(int(n))+"y"
            axesHist.text(tip[0], tip[1], text, fontsize=7, rotation=tip[2])
        else:
            axesHist.axvline(x=x0, c='r')

                                                                                                           
a1 = 3
b1 = 10
a2 = 3
b2 = 10

xmin_max = [0.98*a1/b1, 1.03*a1/b1]
ymin_max = [a2/b2, 1.067*a2/b2]

#previous step
pTune1 = []
pTune2 = []

if __name__ == "__main__":
    # fontsize on plot
    fontsize = 24

    args = SetupArgs()

##@
    plt.ion()
    fig = plt.figure(figsize=(10,6))
    #axesHist = fig.add_subplot(111)
    axesHist = fig.add_subplot(121)
    axesDiff = fig.add_subplot(122)
    data_plot=plt.plot(0,0)
    #cmap = cm.get_cmap('viridis')
    cmap = cm.get_cmap('tab20')
    N = 20 # every 20 steps color will reoccur
##@

    # Read the data from this object
    fr = adios2.open(args.instream, "r", MPI.COMM_WORLD, "adios2_config.xml", "TuneFoot")
#    vars_info = fr.availablevariables()

    tStart = -1;
    tCount = -1;
    # Read through the steps, one at a time
    plot_step = 0
    for fr_step in fr:
#        if fr_step.current_step()
        cur_step= fr_step.current_step()
        vars_info = fr.available_variables()
        varNames = list(vars_info.keys())

        # only draws the first two for 2D image display
        # can extend to 3D later if needed
        tune_name1 = varNames[0]
        tune_name2 = varNames[1];

        numPtls = vars_info[tune_name1]["Shape"].split(',')
        if (plot_step == 0):
            print("Number of points=", numPtls)
        
        tune1_data= fr.read(tune_name1)
        tune2_data= fr.read(tune_name2)

        if (tStart < 0):
            tStart =  fr_step.read_attribute("TurnStart")
            tCount =  fr_step.read_attribute("TurnSize")
                
        currStart =  tStart + plot_step*tCount;
        currEnd   = currStart + tCount;
        titleStr  =  "Turn: " + str(currStart) + " to " + str(currEnd)
        plt.xlabel(tune_name1);
        plt.ylabel(tune_name2);
        plt.title(titleStr)
        #plt.suptitle(numPtls[0]+" particles", ha='left')
        plotMeScatter(tune1_data, tune2_data, plot_step);
        
        plot_step = plot_step + 1;

    fr.close()
    toMovie()

