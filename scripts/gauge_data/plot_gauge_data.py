#!/usr/bin/env python
from __future__ import print_function
"""

Plots gauge streamflow data

Run with::
      run plot_gauge_data.py -i ../../data/objective_1/netcdf/all_gauges.nc
      run plot_gauge_data.py -i ../../data/objective_2/netcdf/all_gauges.nc
      
"""

# -------------------------------------------------------------------------
# Command line arguments
#

import argparse

pngbase   = ''
inputfile = ''
variable  = 'Q'
pdffile   = 'test.pdf'
usetex    = False

parser  = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                  description='''Plots for streamflow gauging data.''')
parser.add_argument('-g', '--pngbase', action='store',
                    default=pngbase, dest='pngbase', metavar='pngbase',
                    help='Name basis for png output files (default: open screen window).')
parser.add_argument('-p', '--pdffile', action='store',
                    default=pdffile, dest='pdffile', metavar='pdffile',
                    help='Name of pdf output file (default: open screen window).')
parser.add_argument('-t', '--usetex', action='store_true', default=usetex, dest="usetex",
                    help="Use LaTeX to render text in pdf.")
parser.add_argument('-v', '--variable', action='store', default=variable, dest='variable',
                    help="Name of variable which will be plotted.")
parser.add_argument('-i', '--inputfile', action='store', default=inputfile, dest='inputfile',
                    help="Name of NC file containing data.")

args            = parser.parse_args()
pngbase         = args.pngbase
pdffile         = args.pdffile
usetex          = args.usetex
variable        = args.variable
inputfile       = args.inputfile

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

# import packages after help so that help with command line -h is fast
import numpy as np
import datetime
import matplotlib.dates as mdates

import color                      # in lib/
from fread      import fread      # in lib/
from position   import position   # in lib/
from str2tex    import str2tex    # in lib/
from readnetcdf import readnetcdf # in lib/


# Latlon
fname = inputfile
lon          = readnetcdf(fname,var='lon')             # 1D field
lat          = readnetcdf(fname,var='lat')             # 1D field
var          = readnetcdf(fname,var=variable)          # 2D field
station_id   = readnetcdf(fname,var='station_id')      # 1D field
station_info = readnetcdf(fname,var='station_info')    # 1D field
times        = readnetcdf(fname,var='time')            # 1D field

# some variable properties
varnames      = readnetcdf(fname,variables=True)
var_idx       = varnames.index(variable)
var_unit      = readnetcdf(fname,units=True)[var_idx]
var_longname  = readnetcdf(fname,longnames=True)[var_idx]
time_idx      = varnames.index('time')
time_unit     = readnetcdf(fname,units=True)[time_idx]
time_longname = readnetcdf(fname,longnames=True)[time_idx]

ref_delta = time_unit.split(' ')[0]
ref_date  = time_unit.split(' ')[2]
ref_hour  = time_unit.split(' ')[3]

reference_day = datetime.datetime(    int(ref_date.split('-')[0]),
                                      int(ref_date.split('-')[1]),
                                      int(ref_date.split('-')[2]),
                                      int(ref_hour.split(':')[0]),
                                      int(ref_hour.split(':')[1]))

if ( ref_delta == 'minutes' ):

    times = np.array([ reference_day + datetime.timedelta(minutes=tt) for tt in times ])

elif ( ref_delta == 'seconds' ):

    times = np.array([ reference_day + datetime.timedelta(seconds=tt) for tt in times ])

elif ( ref_delta == 'hours' ):    

    times = np.array([ reference_day + datetime.timedelta(hours=tt) for tt in times ])

elif ( ref_delta == 'days' ):    

    times = np.array([ reference_day + datetime.timedelta(days=tt) for tt in times ])

else:

    print('')
    print('Time units present:    ',ref_delta)
    print('Time units implemented:','[seconds, minutes,hours,days]')
    raise ValueError( 'The time unit is not implemented yet!' )

nstations = np.shape(station_id)[0]

# plot data only within period
idx = np.where((times>=datetime.datetime(2010,1,1,0,0)) & (times<=datetime.datetime(2015,1,1,0,0)))[0]

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
ncol        = 4           # # of columns of subplots per figure
nrow        = np.int(np.ceil(nstations/4.))           # # of rows of subplots per figure
hspace      = 0.05        # x-space between subplots
vspace      = 0.3/nrow    # y-space between subplots
right       = 0.9         # right space on page
textsize    = 7           # standard text size
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
    mpl.rc('figure', figsize=(8.27,11.69)) # a4 portrait
    if usetex:
        mpl.rc('text', usetex=True)
    else:
        #mpl.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
        mpl.rc('font',**{'family':'serif','serif':['times']})
    mpl.rc('text.latex', unicode=True)
elif (outtype == 'png'):
    mpl.use('Agg') # set directly after import matplotlib
    import matplotlib.pyplot as plt
    mpl.rc('figure', figsize=(8.27,11.69)) # a4 portrait
    if usetex:
        mpl.rc('text', usetex=True)
    else:
        #mpl.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
        mpl.rc('font',**{'family':'serif','serif':['times']})
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
    c = color.get_brewer('dark_rainbow_256', rgb=True)
    c = c[::-1] # reverse colors
    ocean_color = (151/256., 183/256., 224/256.)

cmap = mpl.colors.ListedColormap(c)

# use a colormap from Brewer module
cc        = color.get_brewer('gsdtol', rgb=True)[10:]
cmap_gray = mpl.colors.ListedColormap(cc)

# -------------------------------------------------------------------------
# Plot
# -------------------------------------------------------------------------

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
# Fig 1 :: map with color bar (whole domain)
# -------------------------------------------------------------------------
ifig += 1
iplot = 0
print('Plot - Fig ', ifig, ' ::  ')
fig = plt.figure(ifig)

for istation in range(nstations):
    
    # -------------------------------------------------------------------------
    # (1) streamflow data
    # -------------------------------------------------------------------------
    iplot += 1
    sub    = fig.add_axes(position(nrow,ncol,iplot,hspace=hspace,vspace=vspace) )
    # line1 = sub.plot(times[idx],var[istation,idx])
    line1 = sub.semilogy(times[idx],var[istation,idx])
    
    plt.setp(line1, linestyle='-',
                    linewidth=lwidth, color=lcols[0])

    if usetex:
        xlab   = r'$\mathrm{time}'
        ylab   = r'$\mathrm{'+var_longname+'} \; ['+var_unit.replace('**','^').replace('-1','{-1}')+']$'
    else:
        xlab   = r'time'
        ylab   = r''+var_longname+' [$'+var_unit.replace('**','^').replace('-1','{-1}')+'$]'
        
    if iplot == 1:
        sub.set_title('')
    if (iplot-1)//ncol+1 == nrow:        # last row
        plt.setp(sub, xlabel=xlab)
    else:
        sub.tick_params(axis='x',labelbottom='off')

    if nstations/ncol > 8:
        if iplot%ncol == 1  and ((iplot-1)//ncol+1)%3 == 0:                  # first column, every third row
            plt.setp(sub, ylabel=ylab)
    else:
        if iplot%ncol == 1:                  # first column
            plt.setp(sub, ylabel=ylab)

    # text for gauge ID
    textbox_x = 0.5
    textbox_y = 1.00
    sub.text(textbox_x, textbox_y, station_id[istation], transform=sub.transAxes,
                         rotation=0, fontsize=textsize,
                         horizontalalignment='center', verticalalignment='bottom')
    textbox_x = 0.5
    textbox_y = 1.14

    # text for gauge name
    # insert line break if station_info (without spaces) is longer than ~20 charachters
    if nstations/ncol <= 8:
        tmp = station_info[istation].split(':')[2].strip().replace('"','').split(' ')
        if len(' '.join(tmp)) > 25:
            try:
                cut = len(' '.join(tmp))//2
                nn = np.where( np.cumsum(np.array([ len(ii) for ii in tmp ])) > cut )[0][0]
                tmp2 = str(' '.join(tmp[0:nn])+' \n '+' '.join(tmp[nn:]))
            except:
                tmp2 = str(' '.join(tmp))
        else:
            tmp2 = ' '.join(tmp)
        sub.text(textbox_x, textbox_y, tmp2, transform=sub.transAxes,
                             rotation=0, fontsize=textsize,
                             horizontalalignment='center', verticalalignment='bottom')

    # limits
    plt.setp(sub,xlim=[mdates.date2num(datetime.datetime(2010,1,1,0,0)),mdates.date2num(datetime.datetime(2015,1,1,0,0))])
    plt.setp(sub,ylim=[0.1,np.max(var[:,idx])*1.05])

    # rotate x-ticks
    plt.xticks(rotation=45)

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


