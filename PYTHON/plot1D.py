#!/usr/bin/env python
import argparse
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import matplotlib

def getUserArgs():
    """Read command line arguments and return as an argparse type""" 
    parser = argparse.ArgumentParser(description='Plot 2D contour plot from file specified by user.')
    parser.add_argument('file', metavar='file', nargs=1, help='the name of the file containing the data')
    parser.add_argument('--columns', '-cols', metavar='x y', type=int, nargs=2, help='which columns in file are x and y data')
    parser.add_argument('--xbounds', '-xb', metavar='x0 x1', type=float, nargs=2, help='lower and upper x range')
    parser.add_argument('--ybounds', '-yb', metavar='y0 y1', type=float, nargs=2, help='lower and upper y range')
    parser.add_argument('--xlabel', '-xl', metavar='x label', type=str, nargs=1, help='x label')
    parser.add_argument('--ylabel', '-yl', metavar='y label', type=str, nargs=1, help='y label')
    parser.add_argument('--output', '-o', metavar='output file', type=str, nargs=1, help='output file name')
    return parser
 
args = getUserArgs()
fn_in = args.parse_args().file[0]
cols = args.parse_args().columns if args.parse_args().columns else [0,1]

print("Open file {}.".format( fn_in ))
data=np.loadtxt(fn_in)
print("Done reading file.")

for col in cols:
    if col > len(data[0,:]):
        raise ValueError('Column {} is larger than the number of columns in {}.'.format(col,fn_in))

print ("Plotting {}.".format(fn_in))
# normalize data to max abs z value
maxv=data[:,cols[1]].max()
minv=data[:,cols[1]].min()
maxa=max( abs(maxv), abs(minv) )
data[:,cols[1]] /= maxa
fontsize=8
lw=1. 

Xdata=data[:,cols[0]]
Ydata=data[:,cols[1]]

fig, ax = plt.subplots(figsize=(3.3,2.6))
ax.plot( Xdata, Ydata )
if ( args.parse_args().xbounds ):
    xlim = args.parse_args().xbounds
    ax.set_xlim(xlim)
if ( args.parse_args().ybounds ):
    ylim = args.parse_args().ybounds
    ax.set_ylim(ylim)
if ( args.parse_args().xlabel ):
    ax.set_xlabel(args.parse_args().xlabel[0], fontsize=fontsize)
if ( args.parse_args().ylabel ):
    ax.set_ylabel(args.parse_args().ylabel[0], fontsize=fontsize)

ax.minorticks_on()
ax.tick_params(axis='both',which='both', labelsize=fontsize, direction='in', top='on',bottom='on', left='on', right='on')

plt.tight_layout(pad=0.1)
if ( args.parse_args().output ):
    plt.savefig(args.parse_args().output[0])
else:
    plt.show()
