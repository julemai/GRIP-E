#!/usr/bin/env python
from __future__ import print_function

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

# run:
#

#    run barchart.py -p barchart.pdf

# -----------------------
# add subolder scripts/lib to search path
# -----------------------
import sys
import os 
dir_path = os.path.dirname(os.path.realpath(__file__))
sys.path.append(dir_path+'/../lib')

import color                      # in lib/
from position   import position   # in lib/
from str2tex    import str2tex    # in lib/
from readnetcdf import readnetcdf # in lib/
from autostring import astr       # in lib/

import argparse
import collections
import textwrap                # nicer formatting of help text in parser
import numpy as np             # to perform numerics
import shutil                  # file operations
import copy                    # deep copy objects, arrays etc
import datetime                # converting dates
import errormeasures           # error measures like NSE
import xarray as xr            # xarray
import pandas as pd            # pandas
import netCDF4 as nc           # NetCDF library
import matplotlib.dates as mdates

pngbase     = ''
pdffile     = 'test.pdf'
usetex      = False

parser      = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
              description='''Just plots median NSEs for all models (calibration) as barchart. All numbers hard-coded!!!''')
parser.add_argument('-g', '--pngbase', action='store',
                    default=pngbase, dest='pngbase', metavar='pngbase',
                    help='Name basis for png output files (default: open screen window).')
parser.add_argument('-p', '--pdffile', action='store',
                    default=pdffile, dest='pdffile', metavar='pdffile',
                    help='Name of pdf output file (default: open screen window).')
parser.add_argument('-t', '--usetex', action='store_true', default=usetex, dest="usetex",
                    help="Use LaTeX to render text in pdf.")

args        = parser.parse_args()
pngbase     = args.pngbase
pdffile     = args.pdffile
usetex      = args.usetex

del parser, args


dict_results = collections.OrderedDict()
#                              calibration      validation
#                              obj 1  obj 2     obj 1  obj 2
dict_results['ML-LSTM']      = [0.79, 0.79,     0.41,  0.41 ]
dict_results['ML-XGBoost']   = [0.36, 0.26,     0.17,  0.17 ]
dict_results['LBRM']         = [0.66, 0.72,     0.70,  0.70 ]
dict_results['GR4J-lp']      = [0.63, 0.67,     0.50,  0.50 ]
dict_results['GR4J-sd']      = [0.64, 0.67,     0.44,  0.44 ]
dict_results['HYMOD2-DS']    = [0.74, 0.73,     0.59,  0.59 ]    
dict_results['SWAT-EPA']     = [0.19, np.nan,   np.nan,np.nan ]
dict_results['SWAT-Guelph']  = [0.55, 0.59,     0.26,  0.23 ]
dict_results['mHM-Waterloo'] = [0.76, 0.78,     0.68,  0.68 ]
dict_results['mHM-UFZ']      = [0.66, 0.67,     0.64,  0.60 ]
dict_results['HYPE']         = [0.52, 0.48,     0.41,  0.41 ]
dict_results['VIC']          = [0.41, 0.43,     0.53,  0.51 ]
dict_results['VIC-GRU']      = [0.42, 0.43,     0.51,  0.51 ]
dict_results['GEM-Hydro']    = [0.51, 0.44,     0.54,  0.54 ]
dict_results['MESH-SVS']     = [0.48, 0.46,     0.58,  0.58 ]
dict_results['MESH-CLASS']   = [0.34, 0.40,     0.51,  0.51 ]
dict_results['WATFLOOD']     = [0.33, 0.32,     0.05,  0.05 ]     

models = dict_results.keys()




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
ncol        = 1           # # of columns of subplots per figure
nrow        = 4           # # of rows of subplots per figure
hspace      = 0.05        # x-space between subplots
vspace      = 0.06        # y-space between subplots
right       = 0.9         # right space on page
textsize    = 10           # standard text size
textsize_clock = 0.6*textsize        # standard text size
dxabc       = 0.95        # % of (max-min) shift to the right from left y-axis for a,b,c,... labels
dyabc       = 0.9         # % of (max-min) shift up from lower x-axis for a,b,c,... labels
dxabc_clock = 2.3         # % of (max-min) shift to the right from left y-axis for a,b,c,... labels
dyabc_clock = 0           # % of (max-min) shift up from lower x-axis for a,b,c,... labels
dxsig       = 1.23        # % of (max-min) shift to the right from left y-axis for a,b,c,... labels
dysig       = -0.05       # % of (max-min) shift up from lower x-axis for a,b,c,... labels

lwidth      = 0.5         # linewidth
elwidth     = 1.0         # errorbar line width
alwidth     = 0.5         # axis line width
glwidth     = 0.5         # grid line width
msize       = 1.0         # marker size
mwidth      = 1.0         # marker edge width
mcol1       = '0.0'       # primary marker colour
mcol2       = '0.4'                     # secondary
mcol3       = '0.0' # third
mcols       = ['0.0', '0.4', '0.4', '0.7', '0.7', '1.0']
lcol1       = color.colours('blue')   # primary line colour
lcol2       = '0.4'
lcol3       = '0.0'
lcols       = color.colours(['darkgray','blue','green','yellow','orange','red','darkblue','black','darkgreen','gray'])
hatches     = [None, None, None, None, None, '//']
dobw      = False # True: black & white

# Legend
llxbbox       = 0.0         # x-anchor legend bounding box
llybbox       = 1.0        # y-anchor legend bounding box
llxbbox_clock = -0.2        # x-anchor legend bounding box clock_plot
llybbox_clock = 1.08        # y-anchor legend bounding box clock_plot
llrspace      = 0.          # spacing between rows in legend
llcspace      = 1.0         # spacing between columns in legend
llhtextpad    = 0.4         # the pad between the legend handle and text
llhlength     = 1.5         # the length of the legend handles
frameon       = False       # if True, draw a frame around the legend. If None, use rc

# PNG
dpi           = 600
transparent   = False
bbox_inches   = 'tight'
pad_inches    = 0.035

# Clock options
ymax          = 0.8
ntextsize     = textsize_clock # 'medium'   # normal textsize
# modules
bmod          = 0.5         # fraction of ymax from center to start module colours
alphamod      = 0.7         # alpha channel for modules
fwm           = 0.05         # module width to remove at sides
ylabel1       = 1.15        # position of module names
ylabel2       = 1.35        # position of class names
mtextsize     = 1.3*textsize_clock #'large' # 1.3*textsize # textsize of module labels
# bars
bpar          = 0.4         # fraction of ymax from center to start with parameter bars
fwb           = [0.7,0.4,0.3]  # width of bars
plwidth       = 0.5
# parameters in centre
bplabel       = 0.1        # fractional distance of ymax of param numbers in centre from 0-line
ptextsize     = textsize_clock/1.3 #'small' # 0.8*textsize # textsize of param numbers in centre
# yaxis
ytextsize     = textsize_clock/1.3 #'small' # 0.8*textsize # textsize of y-axis

import matplotlib as mpl
if (outtype == 'pdf'):
    mpl.use('PDF') # set directly after import matplotlib
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    # Customize: http://matplotlib.sourceforge.net/users/customizing.html
    mpl.rc('ps', papersize='a4', usedistiller='xpdf') # ps2pdf
    mpl.rc('figure', figsize=(10.97,11.69)) # a4 portrait
    if usetex:
        mpl.rc('text', usetex=True)
    else:
        #mpl.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
        mpl.rc('font',**{'family':'serif','serif':['times']})
    mpl.rc('text.latex', unicode=True)
elif (outtype == 'png'):
    mpl.use('Agg') # set directly after import matplotlib
    import matplotlib.pyplot as plt
    mpl.rc('figure', figsize=(10.97,11.69)) # a4 portrait
    if usetex:
        mpl.rc('text', usetex=True)
    else:
        #mpl.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
        mpl.rc('font',**{'family':'serif','serif':['times']})
    mpl.rc('text.latex', unicode=True)
    mpl.rc('savefig', dpi=dpi, format='png')
else:
    import matplotlib.pyplot as plt
    mpl.rc('figure', figsize=(4./5.*10.97,4./5.*11.69)) # a4 portrait
mpl.rc('font', size=textsize)
mpl.rc('lines', linewidth=lwidth, color='black')
mpl.rc('axes', linewidth=alwidth, labelcolor='black')
mpl.rc('path', simplify=False) # do not remove
mpl.rcParams['hatch.linewidth'] = lwidth

from matplotlib.patches import Rectangle, Circle, Polygon, Patch
from matplotlib.lines import Line2D
from mpl_toolkits.basemap import Basemap

# colors
if dobw:
    c = np.linspace(0.2, 0.85, nmod)
    c = np.ones(nmod)*0.7
    c = [ str(i) for i in c ]
    ocean_color = '0.1'
else:
    c = color.get_brewer('dark_rainbow_256', rgb=True)
    c = color.get_brewer('BlueWhiteOrangeRed', rgb=True)
    c = c[2:] # drop first two
    c = c[::-1] # reverse colors

    c = color.get_brewer('rdylbu11', rgb=True)

    
    ocean_color = (151/256., 183/256., 224/256.)

cmap = mpl.colors.ListedColormap(c)

# use a colormap from Brewer module
cc        = color.get_brewer('gsdtol', rgb=True)[10:]
cmap_gray = mpl.colors.ListedColormap(cc)

cols_para = color.get_brewer('Paired8', rgb=True)   # need to be at least 9 colors

# -------------------------------------------------------------------------
# Plot
# -------------------------------------------------------------------------

if (outtype == 'pdf'):
    print('Plot PDF ', '.'.join(pdffile.split('.')[0:-1])+'_1.pdf')
    pdf_pages = PdfPages('.'.join(pdffile.split('.')[0:-1])+'_1.pdf')
elif (outtype == 'png'):
    print('Plot PNG ', pngbase)
else:
    print('Plot X')
# figsize = mpl.rcParams['figure.figsize']




# ----------------------------------
# Figure NSE
# ----------------------------------
ifig = 0

ifig += 1
iplot = 0
print('Plot - Fig ', ifig, ' ::  barchart NSE')
fig = plt.figure(ifig)

# ---------------------------------------
# calibration
# ---------------------------------------
iplot += 1
sub    = fig.add_axes(position(nrow,ncol,iplot,hspace=hspace,vspace=vspace), frame_on=False )

results = np.array([ dict_results[imodel][0:2] for imodel in models ])

x = np.arange(len(models))
width = 0.35  # the width of the bars

imodels = np.arange(0,2)
rects1 = sub.bar(x[imodels] - width/2, results[imodels,0], width, label='Objective 1 (Machine Learning)', color=[cols_para[0],cols_para[0]])
rects2 = sub.bar(x[imodels] + width/2, results[imodels,1], width, label='Objective 2 (Machine Learning)', color=[cols_para[1],cols_para[1]])
imodels = np.arange(2,9)
rects3 = sub.bar(x[imodels] - width/2, results[imodels,0], width, label='Objective 1 (single-basin calibr.)', color=[cols_para[2],cols_para[2],cols_para[2],cols_para[2],cols_para[2],cols_para[2],cols_para[2]])
rects4 = sub.bar(x[imodels] + width/2, results[imodels,1], width, label='Objective 2 (single-basin calibr.)', color=[cols_para[3],cols_para[3],cols_para[3],cols_para[3],cols_para[3],cols_para[3],cols_para[3]])
imodels = np.arange(9,17)
rects5 = sub.bar(x[imodels] - width/2, results[imodels,0], width, label='Objective 1 (global calibration)', color=[cols_para[6],cols_para[6],cols_para[6],cols_para[6],cols_para[6],cols_para[6],cols_para[6],cols_para[6]])
rects6 = sub.bar(x[imodels] + width/2, results[imodels,1], width, label='Objective 2 (global calibration)', color=[cols_para[7],cols_para[7],cols_para[7],cols_para[7],cols_para[7],cols_para[7],cols_para[7],cols_para[7]])

# Add some text for labels, title and custom x-axis tick labels, etc.
sub.set_ylabel('Median NSE [-]')
sub.set_title('GRIP-E results (calibration)')
sub.set_xticks(x)
sub.set_xticklabels([ dd.replace('-','-\n') for dd in dict_results.keys() ],fontsize='x-small')

sub.set_yticks([])
sub.set_yticklabels([""],fontsize='x-small')

# set labels on top of bars
def autolabel(rects):
    """Attach a text label above each bar in *rects*, displaying its height."""
    for rect in rects:
        height = rect.get_height()
        sub.annotate('{}'.format(height),
                    xy=(rect.get_x() + rect.get_width() / 2, height),
                    xytext=(0, 3),  # 3 points vertical offset
                    textcoords="offset points",
                    fontsize='x-small',
                    rotation=90,
                    ha='center', va='bottom')


autolabel(rects1)
autolabel(rects2)
autolabel(rects3)
autolabel(rects4)
autolabel(rects5)
autolabel(rects6)

# # legend
# ll = sub.legend(frameon=frameon, ncol=3, labelspacing=llrspace, handletextpad=llhtextpad, handlelength=llhlength,
#                      loc='upper center', bbox_to_anchor=(0.5,-0.15), scatterpoints=1, numpoints=1)
# plt.setp(ll.get_texts(), fontsize='x-small')

plt.setp(sub,xlim=[-0.5,len(models)-0.5])
#plt.setp(sub,ylim=[0.0,1.0])
plt.setp(sub,ylim=[0.0,np.nanmax(results)+0.1])

# ---------------------------------------
# validation
# ---------------------------------------
iplot += 1
sub    = fig.add_axes(position(nrow,ncol,iplot,hspace=hspace,vspace=vspace), frame_on=False )

results = np.array([ dict_results[imodel][2:4] for imodel in models ])

x = np.arange(len(models))
width = 0.35  # the width of the bars

patterns = ["////","","","","","",""]

imodels = np.arange(0,2)
rects1  = sub.bar(x[imodels] - width/2, results[imodels,0], width, label='Objective 1 (Machine Learning)', color=[cols_para[0],cols_para[0]])
rects2  = sub.bar(x[imodels] + width/2, results[imodels,1], width, label='Objective 2 (Machine Learning)', color=[cols_para[1],cols_para[1]])
imodels = np.arange(2,3)
rects3  = sub.bar(x[imodels] - width/2, results[imodels,0], width, label='', color=[cols_para[2]], hatch=patterns[0])
rects4  = sub.bar(x[imodels] + width/2, results[imodels,1], width, label='', color=[cols_para[3]], hatch=patterns[0])
imodels = np.arange(3,9)
rects5  = sub.bar(x[imodels] - width/2, results[imodels,0], width, label='Objective 1 (single-basin calibr.)', color=[cols_para[2],cols_para[2],cols_para[2],cols_para[2],cols_para[2],cols_para[2],cols_para[2]])
rects6  = sub.bar(x[imodels] + width/2, results[imodels,1], width, label='Objective 2 (single-basin calibr.)', color=[cols_para[3],cols_para[3],cols_para[3],cols_para[3],cols_para[3],cols_para[3],cols_para[3]])
imodels = [9,10,11,12,13,14,16] # np.arange(9,17)
rects7  = sub.bar(x[imodels] - width/2, results[imodels,0], width, label='Objective 1 (global calibration)', color=[cols_para[6],cols_para[6],cols_para[6],cols_para[6],cols_para[6],cols_para[6],cols_para[6]])
rects8  = sub.bar(x[imodels] + width/2, results[imodels,1], width, label='Objective 2 (global calibration)', color=[cols_para[7],cols_para[7],cols_para[7],cols_para[7],cols_para[7],cols_para[7],cols_para[7]])
imodels = [15] # np.arange(9,17)
rects9  = sub.bar(x[imodels] - width/2, results[imodels,0], width, label='', color=[cols_para[6]], hatch=patterns[0])
rects10 = sub.bar(x[imodels] + width/2, results[imodels,1], width, label='', color=[cols_para[7]], hatch=patterns[0])


rectsdummy = sub.bar(len(models)+3, 0.0, width, label='some stations used for calibration', color=['white'], hatch=patterns[0], linewidth=lwidth, edgecolor='black' )

# Add some text for labels, title and custom x-axis tick labels, etc.
sub.set_ylabel('Median NSE [-]')
sub.set_title('GRIP-E results (validation)')
sub.set_xticks(x)
sub.set_xticklabels([ dd.replace('-','-\n') for dd in dict_results.keys() ],fontsize='x-small')

sub.set_yticks([])
sub.set_yticklabels([""],fontsize='x-small')

# set labels on top of bars
def autolabel(rects):
    """Attach a text label above each bar in *rects*, displaying its height."""
    for rect in rects:
        height = rect.get_height()
        sub.annotate('{}'.format(height),
                    xy=(rect.get_x() + rect.get_width() / 2, height),
                    xytext=(0, 3),  # 3 points vertical offset
                    textcoords="offset points",
                    fontsize='x-small',
                    rotation=90,
                    ha='center', va='bottom')


autolabel(rects1)
autolabel(rects2)
autolabel(rects3)
autolabel(rects4)
autolabel(rects5)
autolabel(rects6)
autolabel(rects7)
autolabel(rects8)
autolabel(rects9)
autolabel(rects10)

# legend
legend_elements = [ Patch(facecolor=cols_para[0], edgecolor=cols_para[0], label='Objective 1 (Machine Learning)'),
                    Patch(facecolor=cols_para[1], edgecolor=cols_para[1], label='Objective 2 (Machine Learning)'),
                    Patch(facecolor='white', edgecolor='white', linewidth=lwidth, hatch=patterns[0], label=''),
                    Patch(facecolor=cols_para[2], edgecolor=cols_para[2], label='Objective 1 (single-basin calibr.)'),
                    Patch(facecolor=cols_para[3], edgecolor=cols_para[3], label='Objective 2 (single-basin calibr.)'),
                    Patch(facecolor='white', edgecolor='black', linewidth=lwidth, hatch=patterns[0], label='Some stations used for calibration'),
                    Patch(facecolor=cols_para[6], edgecolor=cols_para[6], label='Objective 1 (global calibration)'),
                    Patch(facecolor=cols_para[7], edgecolor=cols_para[7], label='Objective 2 (global calibration)'),
                        ]
ll = sub.legend(handles=legend_elements,
                    frameon=frameon, ncol=3, labelspacing=llrspace, handletextpad=llhtextpad, handlelength=llhlength,
                    loc='upper center', bbox_to_anchor=(0.5,-0.15), scatterpoints=1, numpoints=1)
plt.setp(ll.get_texts(), fontsize='x-small')

# # legend
# ll = sub.legend(frameon=frameon, ncol=3, labelspacing=llrspace, handletextpad=llhtextpad, handlelength=llhlength,
#                      loc='upper center', bbox_to_anchor=(0.5,-0.15), scatterpoints=1, numpoints=1)
# plt.setp(ll.get_texts(), fontsize='x-small')

plt.setp(sub,xlim=[-0.5,len(models)-0.5])
plt.setp(sub,ylim=[0.0,np.nanmax(results)+0.1])


if (outtype == 'pdf'):
    pdf_pages.savefig(fig)
    plt.close(fig)
elif (outtype == 'png'):
    pngfile = pngbase+"{0:04d}".format(ifig)+".png"
    fig.savefig(pngfile, transparent=transparent, bbox_inches=bbox_inches, pad_inches=pad_inches)
    plt.close(fig)


# -------------------------------------------------------------------------
# Finished
# -------------------------------------------------------------------------

if (outtype == 'pdf'):
    pdf_pages.close()
elif (outtype == 'png'):
    pass
else:
    plt.show()


if (outtype == 'pdf'):
    print('Plot PDF ', '.'.join(pdffile.split('.')[0:-1])+'_2.pdf')
    pdf_pages = PdfPages('.'.join(pdffile.split('.')[0:-1])+'_2.pdf')
elif (outtype == 'png'):
    print('Plot PNG ', pngbase)
else:
    print('Plot X')
# figsize = mpl.rcParams['figure.figsize']

# ----------------------------------
# Figure NSE
# ----------------------------------
ifig = 0

ifig += 1
iplot = 0
print('Plot - Fig ', ifig, ' ::  barchart NSE')
fig = plt.figure(ifig)

iplot += 1
sub    = fig.add_axes(position(nrow,ncol,iplot,hspace=hspace,vspace=vspace), frame_on=False )

results = np.array([ dict_results[imodel] for imodel in models ])

x = np.arange(len(models))
width = 0.35  # the width of the bars

imodels = np.arange(0,2)
rects1 = sub.bar(x[imodels], results[imodels,0], width, label='Machine Learning', color=[cols_para[1],cols_para[1]])
imodels = np.arange(2,9)
rects3 = sub.bar(x[imodels], results[imodels,0], width, label='Models with single-basin calibr.', color=[cols_para[3],cols_para[3],cols_para[3],cols_para[3],cols_para[3],cols_para[3],cols_para[3]])
imodels = np.arange(9,17)
rects5 = sub.bar(x[imodels], results[imodels,0], width, label='Models with global calibration', color=[cols_para[7],cols_para[7],cols_para[7],cols_para[7],cols_para[7],cols_para[7],cols_para[7],cols_para[7]])

# Add some text for labels, title and custom x-axis tick labels, etc.
sub.set_ylabel('Median NSE [-]')
sub.set_title('GRIP-E results (calibration; objective 1)')
sub.set_xticks(x)
sub.set_xticklabels([ dd.replace('-','-\n') for dd in dict_results.keys() ],fontsize='x-small')

sub.set_yticks([])
sub.set_yticklabels([""],fontsize='x-small')

# set labels on top of bars
def autolabel(rects):
    """Attach a text label above each bar in *rects*, displaying its height."""
    for rect in rects:
        height = rect.get_height()
        sub.annotate('{}'.format(height),
                    xy=(rect.get_x() + rect.get_width() / 2, height),
                    xytext=(0, 3),  # 3 points vertical offset
                    textcoords="offset points",
                    fontsize='x-small',
                    rotation=90,
                    ha='center', va='bottom')


autolabel(rects1)
autolabel(rects3)
autolabel(rects5)

# legend
ll = sub.legend(frameon=frameon, ncol=3, labelspacing=llrspace, handletextpad=llhtextpad, handlelength=llhlength,
                     loc='upper center', bbox_to_anchor=(0.5,-0.15), scatterpoints=1, numpoints=1)
plt.setp(ll.get_texts(), fontsize='x-small')

plt.setp(sub,xlim=[-0.5,len(models)-0.5])
plt.setp(sub,ylim=[0.0,np.nanmax(results)+0.1])


if (outtype == 'pdf'):
    pdf_pages.savefig(fig)
    plt.close(fig)
elif (outtype == 'png'):
    pngfile = pngbase+"{0:04d}".format(ifig)+".png"
    fig.savefig(pngfile, transparent=transparent, bbox_inches=bbox_inches, pad_inches=pad_inches)
    plt.close(fig)


# -------------------------------------------------------------------------
# Finished
# -------------------------------------------------------------------------

if (outtype == 'pdf'):
    pdf_pages.close()
elif (outtype == 'png'):
    pass
else:
    plt.show()
    

if (outtype == 'pdf'):
    print('Plot PDF ', '.'.join(pdffile.split('.')[0:-1])+'_3.pdf')
    pdf_pages = PdfPages('.'.join(pdffile.split('.')[0:-1])+'_3.pdf')
elif (outtype == 'png'):
    print('Plot PNG ', pngbase)
else:
    print('Plot X')
# figsize = mpl.rcParams['figure.figsize']

# ----------------------------------
# Figure NSE
# ----------------------------------
ifig = 0

ifig += 1
iplot = 0
print('Plot - Fig ', ifig, ' ::  barchart NSE')
fig = plt.figure(ifig)

iplot += 1
sub    = fig.add_axes(position(nrow,ncol,iplot,hspace=hspace,vspace=vspace), frame_on=False )

results = np.array([ dict_results[imodel] for imodel in models ])

x = np.arange(len(models))
width = 0.35  # the width of the bars

imodels = np.arange(0,2)
rects2 = sub.bar(x[imodels], results[imodels,1], width, label='Machine Learning', color=[cols_para[1],cols_para[1]])
imodels = np.arange(2,9)
rects4 = sub.bar(x[imodels], results[imodels,1], width, label='Models with single-basin calibr.', color=[cols_para[3],cols_para[3],cols_para[3],cols_para[3],cols_para[3],cols_para[3],cols_para[3]])
imodels = np.arange(9,17)
rects6 = sub.bar(x[imodels], results[imodels,1], width, label='Models with global calibration', color=[cols_para[7],cols_para[7],cols_para[7],cols_para[7],cols_para[7],cols_para[7],cols_para[7],cols_para[7]])

# Add some text for labels, title and custom x-axis tick labels, etc.
sub.set_ylabel('Median NSE [-]')
sub.set_title('GRIP-E results (calibration; objective 2)')
sub.set_xticks(x)
sub.set_xticklabels([ dd.replace('-','-\n') for dd in dict_results.keys() ],fontsize='x-small')

sub.set_yticks([])
sub.set_yticklabels([""],fontsize='x-small')

# set labels on top of bars
def autolabel(rects):
    """Attach a text label above each bar in *rects*, displaying its height."""
    for rect in rects:
        height = rect.get_height()
        sub.annotate('{}'.format(height),
                    xy=(rect.get_x() + rect.get_width() / 2, height),
                    xytext=(0, 3),  # 3 points vertical offset
                    textcoords="offset points",
                    fontsize='x-small',
                    rotation=90,
                    ha='center', va='bottom')


autolabel(rects2)
autolabel(rects4)
autolabel(rects6)

# legend
ll = sub.legend(frameon=frameon, ncol=3, labelspacing=llrspace, handletextpad=llhtextpad, handlelength=llhlength,
                     loc='upper center', bbox_to_anchor=(0.5,-0.15), scatterpoints=1, numpoints=1)
plt.setp(ll.get_texts(), fontsize='x-small')

plt.setp(sub,xlim=[-0.5,len(models)-0.5])
plt.setp(sub,ylim=[0.0,np.nanmax(results)+0.1])


if (outtype == 'pdf'):
    pdf_pages.savefig(fig)
    plt.close(fig)
elif (outtype == 'png'):
    pngfile = pngbase+"{0:04d}".format(ifig)+".png"
    fig.savefig(pngfile, transparent=transparent, bbox_inches=bbox_inches, pad_inches=pad_inches)
    plt.close(fig)



    

# -------------------------------------------------------------------------
# Finished
# -------------------------------------------------------------------------

if (outtype == 'pdf'):
    pdf_pages.close()
elif (outtype == 'png'):
    pass
else:
    plt.show()


