#!/usr/bin/env python

# Copyright 2016-2018 Juliane Mai - juliane.mai(at)uwaterloo.ca
#
# License
# This file is part of Juliane Mai's personal code library.
#
# Juliane Mai's personal code library is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Juliane Mai's personal code library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.

# You should have received a copy of the GNU Lesser General Public License
# along with Juliane Mai's personal code library.  If not, see <http://www.gnu.org/licenses/>.
#

from __future__ import print_function

# -------------------------------------------------------------------------
# General settings
#
dobw      = False # True: black & white
docomp    = True  # True: Print classification on top of modules
dosig     = False # True: add signature to plot
dolegend  = False # True: add legend to each subplot
doabc     = True  # True: add subpanel numbering
dotitle   = True  # True: add catchment titles to subpanels

# -------------------------------------------------------------------------
# Command line arguments
#

import argparse

pngbase   = ''
pdffile   = 'test.pdf'
usetex    = False
inputpath = '../../data'
inputfile = ['../../data/monthly_posterior_distributions_Smith_Gronewold/nbs_lake_erie_summary-stats.csv',
             '../../data/monthly_posterior_distributions_Smith_Gronewold/nbs_lake_st-clair_summary-stats.csv',
             '../../data/monthly_posterior_distributions_Smith_Gronewold/inflow_lake_erie_summary-stats.csv',
             '../../data/monthly_posterior_distributions_Smith_Gronewold/inflow_lake_st-clair_summary-stats.csv']

parser  = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                  description='''Plot basin shape.''')
parser.add_argument('-g', '--pngbase', action='store',
                    default=pngbase, dest='pngbase', metavar='pngbase',
                    help='Name basis for png output files (default: open screen window).')
parser.add_argument('-p', '--pdffile', action='store',
                    default=pdffile, dest='pdffile', metavar='pdffile',
                    help='Name of pdf output file (default: open screen window).')
parser.add_argument('-t', '--usetex', action='store_true', default=usetex, dest="usetex",
                    help="Use LaTeX to render text in pdf.")
parser.add_argument('-i', '--inputfile', action='store', default=inputfile, dest='inputfile',
                    help="Name of CSV file containing summary statistics. File with smaples will automatically added based on filename base.")
parser.add_argument('inputpath', nargs='?', default=inputpath, metavar='path of input files',
                    help='Path with all the input files (Default: ../data/)')

args      = parser.parse_args()
pngbase   = args.pngbase
pdffile   = args.pdffile
usetex    = args.usetex
inputpath = args.inputpath
inputfile = args.inputfile

if (pdffile != '') & (pngbase != ''):
    print('\nError: PDF and PNG are mutually exclusive. Only either -p or -g possible.\n')
    parser.print_usage()
    import sys
    sys.exit()

del parser, args

# -----------------------
# add subolder scripts/lib to search path
# -----------------------
import sys
import os 
dir_path = os.path.dirname(os.path.realpath(__file__))
sys.path.append(dir_path+'/../lib')

import color                            # in lib/
from position      import position      # in lib/
from abc2plot      import abc2plot      # in lib/
from brewer        import get_brewer    # in lib/
from autostring    import astr          # in lib/
from str2tex       import str2tex       # in lib/
from fread         import fread         # in lib/

# import fiona          # some shapefile coordinate stuff
import numpy as np
import xarray as xr
import pandas as pd
import copy                       # deep copy objects, arrays etc
import time
import os
import glob
import shapefile
from matplotlib.patches import Polygon, Ellipse
import matplotlib.dates as mdates
import datetime
import calendar
t1 = time.time()

# [ Year,Month,Median,2.5 Percentile,97.5 Percentile ]
data_summary = [ [] for ifile in inputfile]
for ii,ifile in enumerate(inputfile):
    print("Reading file: ", ifile)
    data_summary[ii] = fread(ifile,skip=1,separator=',')
    

times = np.array([ datetime.date(*[int(data_summary[0][ii,0]),int(data_summary[0][ii,1]),1]) for ii in range(np.shape(data_summary[0])[0]) ])

data_summary = [ data_summary[ii][:,2:]  for ii in range(len(inputfile)) ]   # only Median,2.5 Percentile,97.5 Percentile

# [ ncols=ntime, nrows=nens ]
inputensemblefile = [ '_'.join(ifile.split('_')[0:-1])+'_samples.csv' for ifile in inputfile]
data_ens = [ [] for efile in inputensemblefile]
for ii,efile in enumerate(inputensemblefile):

    print("Reading file: ", efile)
    data_ens[ii] = fread(efile,separator=',')

# -------------------------------------------------------------------------
# Customize plots
#

if (pdffile == ''):
    if (pngbase == ''):
        outtype = 'x'
    else:
        outtype = 'png'
else:
    outtype = 'pdf'
    
# Main plot
nrow        = 6           # # of rows of subplots per figure
ncol        = 2           # # of columns of subplots per figure
hspace      = 0.12        # x-space between subplots
vspace      = 0.03        # y-space between subplots
right       = 0.9         # right space on page
textsize    = 10          # standard text size
dxabc       = 0.99        # % of (max-min) shift to the right from left y-axis for a,b,c,... labels
dyabc       = 0.02        # % of (max-min) shift up from lower x-axis for a,b,c,... labels
dxsig       = 1.23        # % of (max-min) shift to the right from left y-axis for signature
dysig       = -0.05       # % of (max-min) shift up from lower x-axis for signature
dxtit       = 0           # % of (max-min) shift to the right from left y-axis for title
dytit       = 1.3         # % of (max-min) shift up from lower x-axis for title

lwidth      = 1.5         # linewidth
elwidth     = 1.0         # errorbar line width
alwidth     = 1.0         # axis line width
glwidth     = 0.5         # grid line width
msize       = 3.0         # marker size
mwidth      = 1.0         # marker edge width
mcol1       = '0.0'       # primary marker colour
mcol2       = '0.4'       # secondary
mcol3       = '0.0'       # third
mcols       = ['0.0', '0.4', '0.4', '0.7', '0.7', '1.0']
lcol1       = color.colours('blue')   # primary line colour
lcol2       = '0.0'
lcol3       = '0.0'
lcols       = ['None', 'None', 'None', 'None', 'None', '0.0']
hatches     = [None, None, None, None, None, '//']

# Legend
llxbbox     = 0.85        # x-anchor legend bounding box
llybbox     = -0.2        # y-anchor legend bounding box
llrspace    = 0.          # spacing between rows in legend
llcspace    = 1.0         # spacing between columns in legend
llhtextpad  = 0.4         # the pad between the legend handle and text
llhlength   = 1.5         # the length of the legend handles
frameon     = False       # if True, draw a frame around the legend. If None, use rc

# PNG
dpi         = 300         # 150 for testing
transparent = False
bbox_inches = 'tight'
pad_inches  = 0.035

# Clock options
ymax = 0.6
ntextsize   = 'medium'       # normal textsize
# modules
bmod        = 0.5            # fraction of ymax from center to start module colours
alphamod    = 0.7            # alpha channel for modules
fwm         = 0.05           # module width to remove at sides
ylabel1     = 1.15           # position of module names
ylabel2     = 1.35           # position of class names
mtextsize   = 'large'        # 1.3*textsize # textsize of module labels
# bars
bpar        = 0.4            # fraction of ymax from center to start with parameter bars
fwb         = [0.7,0.4,0.3]  # width of bars
plwidth     = 0.5
# parameters in centre
bplabel     = 0.1            # fractional distance of ymax of param numbers in centre from 0-line
ptextsize   = 'medium'       # 'small' # 0.8*textsize # textsize of param numbers in centre
# yaxis
space4yaxis = 2              # space for y-axis (integer)
ytextsize   = 'medium'       # 'small' # 0.8*textsize # textsize of y-axis
sig         = 'J Mai' # sign the plot

import matplotlib as mpl
if (outtype == 'pdf'):
    mpl.use('PDF') # set directly after import matplotlib
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    # Customize: http://matplotlib.sourceforge.net/users/customizing.html
    mpl.rc('ps', papersize='a4', usedistiller='xpdf') # ps2pdf
    mpl.rc('figure', figsize=(8.27,11.69)) # a4 portrait
    if usetex:
        mpl.rc('text', usetex=True)
    else:
        mpl.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
        #mpl.rc('font',**{'family':'serif','serif':['times']})
    mpl.rc('text.latex', unicode=True)
elif (outtype == 'png'):
    mpl.use('TkAgg') # set directly after import matplotlib
    import matplotlib.pyplot as plt
    mpl.rc('figure', figsize=(8.27,11.69)) # a4 portrait
    if usetex:
        mpl.rc('text', usetex=True)
    else:
        mpl.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
        #mpl.rc('font',**{'family':'serif','serif':['times']})
    mpl.rc('text.latex', unicode=True)
    mpl.rc('savefig', dpi=dpi, format='png')
else:
    import matplotlib.pyplot as plt
    mpl.rc('figure', figsize=(4./5.*8.27,4./5.*11.69)) # a4 portrait
mpl.rc('font', size=textsize)
mpl.rc('lines', linewidth=lwidth, color='black')
mpl.rc('axes', linewidth=alwidth, labelcolor='black')
mpl.rc('path', simplify=False) # do not remove

from matplotlib.patches import Rectangle, Circle, Polygon
from mpl_toolkits.basemap import Basemap


# colors
if dobw:
    c = np.linspace(0.2, 0.85, nmod)
    c = np.ones(nmod)*0.7
    c = [ str(i) for i in c ]
    ocean_color = '0.1'
else:
    # c = [(165./255.,  0./255., 38./255.), # interception
    #      (215./255., 48./255., 39./255.), # snow
    #      (244./255.,109./255., 67./255.), # soil moisture
    #      (244./255.,109./255., 67./255.), # soil moisture
    #      (253./255.,174./255., 97./255.), # direct runoff
    #      (254./255.,224./255.,144./255.), # Evapotranspiration
    #      (171./255.,217./255.,233./255.), # interflow
    #      (116./255.,173./255.,209./255.), # percolation
    #      ( 69./255.,117./255.,180./255.), # routing
    #      ( 49./255., 54./255.,149./255.)] # geology
    c = get_brewer('rdylbu11', rgb=True)
    tmp = c.pop(5)   # rm yellow
    #c.insert(2,c[2]) # same colour for both soil moistures
    ocean_color = (151/256., 183/256., 224/256.)

    cc = color.get_brewer('dark_rainbow_256', rgb=True)
    cc = cc[::-1] # reverse colors
    cmap = mpl.colors.ListedColormap(cc)


    
# -------------------------------------------------------------------------
# Plot
#

if (outtype == 'pdf'):
    print('Plot PDF ', pdffile)
    pdf_pages = PdfPages(pdffile)
elif (outtype == 'png'):
    print('Plot PNG ', pngbase)
else:
    print('Plot X')
# figsize = mpl.rcParams['figure.figsize']

ifig = 0

# -------------------------------------------------------------------------
# Fig 1 - summary statistics
#
ifig += 1
iplot = 0
print('Plot - Fig ', ifig)
fig = plt.figure(ifig)

for ii,ifile in enumerate(inputfile):
    
    # -----------------------------------------
    # Band of Median and 5% and 95% percentiles of current file
    # ----------------------------------------
    iplot += 1
    sub = fig.add_axes(position(nrow,ncol,iplot,hspace=hspace,vspace=vspace)) #, axisbg='none')

    variable = ifile.split('/')[-1].split('.')[0].split('_')[0]  # inflow or nbs
    lake     = ifile.split('/')[-1].split('.')[0].split('_')[2]  # st-clair or erie
    datatype = ifile.split('/')[-1].split('.')[0].split('_')[3]  # summary or sample

    fig.suptitle(datatype.replace('-',' ').title(), fontsize=textsize+3, x=0.5, y=0.93)

    if variable == 'nbs':
        variable = 'NBS'
    else:
        variable = variable.title()

    if usetex:
        xlab   = r'$\mathrm{Date time}'
        ylab   = r'$\mathrm{'+variable+'} \; [\mathrm{m}^{3} \mathrm{s}^{-1}]$'
    else:
        xlab   = r'Date time'
        ylab   = r''+variable+' [$m^3 s^{-1}$]'
    if iplot == 1:
        sub.set_title('')
    if iplot == 5:
        plt.setp(sub, xlabel=xlab) # axis labels
    if iplot < 6:
        plt.setp(sub, ylabel=ylab)

    # plot only after 2010
    idx = np.where(times>datetime.date(*[2010,1,1]))[0]
        
    sub.plot(         times[idx], data_summary[ii][idx,0], linewidth=0.75, color=cc[0], label='Median')
    sub.fill_between( times[idx], data_summary[ii][idx,1],data_summary[ii][idx,2], linewidth=0.75, color='gray', alpha=0.5, label='Percentiles [$p_{5}$,$p_{95}$]')

    if iplot == 1:
        ll = sub.legend(frameon=frameon, ncol=1,
                    fontsize=textsize-3,
                    labelspacing=llrspace, handletextpad=llhtextpad, handlelength=llhlength,
                    loc='upper right', bbox_to_anchor=(1.0,1.0), scatterpoints=1, numpoints=1)

    # text for lake
    textbox_x = 0.05
    textbox_y = 0.9
    sub.text(textbox_x, textbox_y, lake.replace('-','. ').title(), transform=sub.transAxes,
                         rotation=0, fontsize=textsize,
                         horizontalalignment='left', verticalalignment='center')

    # subplot numbering
    abc2plot(sub, dxabc, dyabc, iplot, lower=True, parenthesis='close',
                 bold=True, large=True,
                 mathrm=True, usetex=usetex,
                 horizontalalignment='right', verticalalignment='bottom')

if (outtype == 'pdf'):
    pdf_pages.savefig(fig)
    plt.close(fig)
elif (outtype == 'png'):
    pngfile = pngbase+"{0:04d}".format(ifig)+".png"
    fig.savefig(pngfile, transparent=transparent, bbox_inches=bbox_inches, pad_inches=pad_inches)
    plt.close(fig)


# -------------------------------------------------------------------------
# Fig 2 - Ensemble members
#
ifig += 1
iplot = 0
print('Plot - Fig ', ifig)
fig = plt.figure(ifig)

for ii,efile in enumerate(inputensemblefile):
    
    # -----------------------------------------
    # Band of Median and 5% and 95% percentiles of current file
    # ----------------------------------------
    iplot += 1
    sub = fig.add_axes(position(nrow,ncol,iplot,hspace=hspace,vspace=vspace)) #, axisbg='none')

    variable = efile.split('/')[-1].split('.')[0].split('_')[0]  # inflow or nbs
    lake     = efile.split('/')[-1].split('.')[0].split('_')[2]  # st-clair or erie
    datatype = efile.split('/')[-1].split('.')[0].split('_')[3]  # summary or sample

    fig.suptitle(datatype.replace('-',' ').title(), fontsize=textsize+3, x=0.5, y=0.93)

    if variable == 'nbs':
        variable = 'NBS'
    else:
        variable = variable.title()
        
    if usetex:
        xlab   = r'$\mathrm{Date time}'
        ylab   = r'$\mathrm{'+variable+'} \; [\mathrm{m}^{3} \mathrm{s}^{-1}]$'
    else:
        xlab   = r'Date time'
        ylab   = r''+variable+' [$m^3 s^{-1}$]'
    if iplot == 1:
        sub.set_title('')
    if iplot == 5:
        plt.setp(sub, xlabel=xlab) # axis labels
    if iplot < 6:
        plt.setp(sub, ylabel=ylab)

    # plot only after 2010
    idx = np.where(times>datetime.date(*[2010,1,1]))[0]

    # only to get the label
    sub.plot(         times[idx], np.transpose(data_ens[ii][0,idx]), linewidth=0.75, color='gray', alpha=0.5, label='Samples')
    sub.plot(         times[idx], np.transpose(data_ens[ii][:,idx]), linewidth=0.75, color='gray', alpha=0.5)
    #sub.fill_between( times[idx], data_ens[ii][idx,1],data_ens[ii][idx,2], linewidth=0.75, color='gray', alpha=0.5, label='Percentiles [$p_{5}$,$p_{95}$]')

    if iplot == 1:
        ll = sub.legend(frameon=frameon, ncol=1,
                    fontsize=textsize-3,
                    labelspacing=llrspace, handletextpad=llhtextpad, handlelength=llhlength,
                    loc='upper right', bbox_to_anchor=(1.0,1.0), scatterpoints=1, numpoints=1)

    # text for lake
    textbox_x = 0.05
    textbox_y = 0.9
    sub.text(textbox_x, textbox_y, lake.replace('-','. ').title(), transform=sub.transAxes,
                         rotation=0, fontsize=textsize,
                         horizontalalignment='left', verticalalignment='center')

    # subplot numbering
    abc2plot(sub, dxabc, dyabc, iplot, lower=True, parenthesis='close',
                 bold=True, large=True,
                 mathrm=True, usetex=usetex,
                 horizontalalignment='right', verticalalignment='bottom')

if (outtype == 'pdf'):
    pdf_pages.savefig(fig)
    plt.close(fig)
elif (outtype == 'png'):
    pngfile = pngbase+"{0:04d}".format(ifig)+".png"
    fig.savefig(pngfile, transparent=transparent, bbox_inches=bbox_inches, pad_inches=pad_inches)
    plt.close(fig)

# -------------------------------------------------------------------------
# Finished
#

if (outtype == 'pdf'):
    pdf_pages.close()
elif (outtype == 'png'):
    pass
else:
    plt.show()

t2    = time.time()
strin = '[m]: '+astr((t2-t1)/60.,1) if (t2-t1)>60. else '[s]: '+astr(t2-t1,0)
print('Time ', strin)




    
