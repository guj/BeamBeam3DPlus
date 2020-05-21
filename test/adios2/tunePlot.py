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
import sys
from math  import gcd
from fractions import Fraction

def SetupArgs():
    parser = argparse.ArgumentParser()
    parser.add_argument("--instream", "-i", help="Name of the input stream", required=True)

    parser.add_argument('--attr', '-a',  help='attr names (default x y)', nargs='+', type=str, default=['x', 'y'])
    parser.add_argument('--point', '-p', help='tune point (default 0 0) (default does not draw)', nargs='+', type=str, default=['0','0'])

    parser.add_argument("--refreshSecond", "-s", help="refresh after (default: 1) seconds. 0 means save to file directly", default=1)
    
    args = parser.parse_args()

    args.refreshSecond = float(args.refreshSecond)
    
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
        cc = np.sqrt(delta1 + delta2);

    pic=axesDiff.scatter(t1,t2, s=2, c=cc, cmap=plt.cm.plasma_r)


    pTune1[:] = t1
    pTune2[:] = t2

    return fig.colorbar(pic, ax=axesDiff);

#
# adjust min/max according to actual data
#
def adjust(xx, yy):
    xmin  = xx[0]
    xmax  = xx[1]

    ymin  = yy[0]
    ymax  = yy[1]

    xmin_max[0] = min(0.99*xmin, xmin_max[0]);
    xmin_max[1] = max(1.01*xmax, xmin_max[1]);
    
    ymin_max[0] = min(0.99*ymin, ymin_max[0])
    ymin_max[1] = max(1.01*ymax, ymin_max[1])

def plotMeScatter(t1,t2,step):
    if (step == 0):
        m1=[min(t1), max(t1)]
        m2=[min(t2), max(t2)]

        adjust(m1, m2);
        axesHist.set_xlim(xmin_max[0], xmin_max[1]);
        axesHist.set_ylim(ymin_max[0], ymin_max[1]);

        axesDiff.set_ylim(ymin_max[0], ymin_max[1]);
        axesDiff.set_xlim(xmin_max[0], xmin_max[1]);
        if ((a1 > 0) | (a2 > 0)) :
            plotLines(a1,b1,a2,b2)

    cb1 = plotHistogram(t1,t2)
    cb2 = plotDiff(t1,t2, step)
    plt.savefig("figs/"+str(step)+".png")
    if (args.refreshSecond > 0):
        plt.show();
        plt.pause(args.refreshSecond)
    cb1.remove()
    cb2.remove()

def toMovie():
    os.system("ffmpeg -r 1 -i figs/%01d.png -vcodec mpeg4 -y figs/movie.mp4 > /dev/null 2>&1")


def index(m, n):
    return m+n

## draw on the Upper half of the plane
def findTextTipUH(m, n, c, xlim, ylim):
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


def findTextTipLH(m, n, c, xlim, ylim):
    if (m == 0):
        return [xlim[0], 0, 0]

    plotSlope = (ylim[1] - ylim[0])/(xlim[1] - xlim[0])
    if (m/n > plotSlope):
        xBottom = (c - n*ylim[0])/m

        xStep = (xBottom - xlim[0])*0.15
        xBottom -= xStep;
        return [xBottom, (c-m*xBottom)/n, 40];
    else:
        yBottom = (c - m*xlim[1])/n

        yStep = (yBottom - ylim[0])*0.15
        yBottom += yStep

        return [(c - n*yBottom)/m, yBottom, 0]


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
    #x = np.linspace(xmin_max[0], a1/b1,30)

    x0 = a1/b1
    y0 = a2/b2
    for (m,n) in final:
        if (n != 0):
            y = (m*x0 + n*y0 - m*x)/n
            slope = -m*1.0/n*1.0
            equation = str(m)+"x + "+str(n)+"y="+str(m*x0+n*y0)
            angle = np.arctan(slope)*180/3.1415926
            axesHist.plot(x,y, '-r', linewidth=0.8)

            if (xmin_max[0] < a1/b1):
                xmin_max[1] = a1/b1
                tip=findTextTipUH(m,n, m*x0+n*y0, xmin_max, ymin_max)
            else:
                ymin_max[1] = a2/b2;
                tip=findTextTipLH(m,n, m*x0+n*y0, xmin_max, ymin_max)

            text=str(int(m))+tune_name1+"+"+str(int(n))+tune_name2
            axesHist.text(tip[0], tip[1], text, fontsize=7, rotation=tip[2])
        else:
            axesHist.axvline(x=x0, c='r', linewidth=0.8)


  
                                                                                                         
#a1 = 3
#b1 = 10
#a2 = 3
#b2 = 10

#previous step
pTune1 = []
pTune2 = []

if __name__ == "__main__":
    # fontsize on plot
    fontsize = 24

    args = SetupArgs()

    if (len(args.attr) < 2):
        print("Please input 2 attr names ")
        sys.exit()
    if (len(args.point) < 2):
        print("Please input 2 D point ")
        sys.exit()
    if (args.attr[0] is args.attr[1]):
        print("Both attributes are the same. No work to be done. Bye!")
        sys.exit()

    #tune points
    a1 = Fraction(args.point[0]).limit_denominator().numerator
    b1 = Fraction(args.point[0]).limit_denominator().denominator
    a2 = Fraction(args.point[1]).limit_denominator().numerator
    b2 = Fraction(args.point[1]).limit_denominator().denominator

    print ("- Input file: ", args.instream)
    print ("- Attributes: ", args.attr, "Tune Point: ", args.point)
    print ("- Outputs pngs will be in subdir figs/\n")

    if ((a1 == 0) & (a2 == 0)):
        xmin_max=[10,-10]
        ymin_max=[10,-10]
    else:
        xmin_max = [a1/b1, a1/b1]
        ymin_max = [a2/b2, a2/b2]

    plt.ion()
    fig = plt.figure(figsize=(16,6))

    axesHist = fig.add_subplot(121)
    axesDiff = fig.add_subplot(122)
    data_plot=plt.plot(0,0)

    #cmap = cm.get_cmap('viridis')
    cmap = cm.get_cmap('tab20')


    # Read the data from this object
    fr = adios2.open(args.instream, "r", MPI.COMM_WORLD, "adios2_config.xml", "TuneFoot")
#    vars_info = fr.availablevariables()

    tStart = -1;
    tCount = -1;
    # Read through the steps, one at a time
    plot_step = 0

    tune_name1 = args.attr[0]
    tune_name2 = args.attr[1]
    for fr_step in fr:
#        if fr_step.current_step()
        cur_step= fr_step.current_step()
        vars_info = fr.available_variables()
        varNames = list(vars_info.keys())

        # only draws the first two for 2D image display
        # can extend to 3D later if needed

        if ((tune_name1 not in vars_info) or (tune_name2 not in vars_info)):
            print("Please make sure attrs are in the input file. Current known attrs:", list(vars_info.keys()))
            sys.exit()

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

        axesHist.set_xlabel(tune_name1);
        axesHist.set_ylabel(tune_name2);
        axesHist.set_title(titleStr)

        axesDiff.set_xlabel(tune_name1);
        axesDiff.set_ylabel(tune_name2);
        axesDiff.set_title(titleStr)

        #plt.suptitle(numPtls[0]+" particles", ha='left')
        plotMeScatter(tune1_data, tune2_data, plot_step);
        
        plot_step = plot_step + 1;

    fr.close()
    toMovie()

