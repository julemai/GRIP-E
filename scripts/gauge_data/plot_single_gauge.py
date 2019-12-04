#!/usr/bin/env python
from __future__ import print_function
"""

Plots gauge streamflow data

Run with::
      run plot_single_gauge.py -i ../../data/objective_1/great-lakes/calibration/netcdf/all_gauges.nc -s 02HC030 -p 02HC030.pdf -v Q -a '2001-01-01:2017-01-01'
      
"""

# -------------------------------------------------------------------------
# Command line arguments
#

import argparse

pngbase          = ''
inputfile        = ''
variable         = 'Q'
pdffile          = 'test.pdf'
usetex           = False
station          = ''
time_period      = ''
model_qsim_files = 'None'

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
parser.add_argument('-s', '--station', action='store', default=station, dest='station',
                    help="Name of station to plot.")
parser.add_argument('-a', '--time_period', action='store',
                    default=time_period, dest='time_period', metavar='time_period', nargs=1,
                    help='Time period to plot and calculate performance of (format: YYYY-MM-DD:YYYY-MM-DD; default: all time points).')
parser.add_argument('-m', '--model_qsim_files', '--model_qsim_files', action='store',
                    default=model_qsim_files, dest='model_qsim_files', metavar='model_qsim_files', nargs=1,
                    help='Name of NetCDF file containing model outputs.')

args             = parser.parse_args()
pngbase          = args.pngbase
pdffile          = args.pdffile
usetex           = args.usetex
variable         = args.variable
inputfile        = args.inputfile
station          = args.station
time_period      = args.time_period[0]
model_qsim_files = (args.model_qsim_files[0]).split(' ')

if (pdffile != '') & (pngbase != ''):
    print('\nError: PDF and PNG are mutually exclusive. Only either -p or -g possible.\n')
    parser.print_usage()
    import sys
    sys.exit()

if (station == ''):
    raise ValueError("Station must be set (option -s)!")

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
import xarray as xr
import pandas as pd            # pandas

import errormeasures              # in lib/ :: error measures like NSE
import color                      # in lib/
from fread      import fread      # in lib/
from position   import position   # in lib/
from str2tex    import str2tex    # in lib/
from readnetcdf import readnetcdf # in lib/


# -------------------------------------
# Read modelled values
# -------------------------------------
        
if model_qsim_files != 'None':

    dicts_nse     = {}
    dicts_rmse    = {}
    dicts_pbias   = {}
    dicts_lognse  = {}
    dicts_sqrtnse = {}
    dicts_dates   = {}
    dicts_qobs    = {}
    dicts_qsim    = {}
    dicts_info    = {}  # station info, i.e. station long name

    for imodel_qsim_file,model_qsim_file in enumerate(model_qsim_files):
        
        # simulated streamflow
        data_sim = xr.open_dataset(model_qsim_file)

        # observed streamflow
        data_obs = xr.open_dataset(os.path.dirname(model_qsim_file)+'/../../netcdf/all_gauges.nc')

        # sort all according to observations
        idx_station = [ list(data_sim.station_id).index(ii) for ii in data_obs.station_id ]
        idx_time    = [ np.where(data_obs.time==sim_tt)[0][0] for sim_tt in data_sim.time ]

        # cut out time period requested
        if time_period != '':
            start = time_period.split(':')[0]
            end   = time_period.split(':')[1]
            idx_period = np.where((data_sim.time>=np.datetime64(start)).data & (data_sim.time<=np.datetime64(end)).data)[0]

        # some stuff for labels in plots
        var_longname = data_obs['Q'].long_name
        var_unit     = data_obs['Q'].units

        # calculate objectives
        dict_dates   = {}
        dict_qsim    = {}
        dict_info    = {}
        nstations    = len(data_obs.station_id)
        idx_station_sim = np.where(data_sim.station_id==station)[0]
        
        for istation in [idx_station_sim]:
            stat_id = str(data_sim.station_id[istation][0].data)

            Qsim    = data_sim.Q[idx_station_sim,:]
            dates   = data_sim.time.data
            
            if time_period != '':
                Qsim  = Qsim[0,idx_period]
                dates = dates[idx_period]
            
            dict_qsim[stat_id]    = Qsim.data
            # convert from numpy.datetime64 to datetime.datetime
            dict_dates[stat_id]   = np.array([ pd.Timestamp(itime).to_pydatetime() for itime in dates ])
            # station long name
            dict_info[stat_id]    = data_obs['station_info'][istation]

        # save everything in dictionaries
        dicts_dates[model_qsim_file.split('/')[-1].split('_')[0]]   = dict_dates
        dicts_qsim[model_qsim_file.split('/')[-1].split('_')[0]]    = dict_qsim
        dicts_info[model_qsim_file.split('/')[-1].split('_')[0]]    = dict_info

        data_obs.close()
        data_sim.close()

    # print("dicts_nse  = ",dicts_nse)

# -------------------------------------
# Read observations
# -------------------------------------

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
# idx   = np.where((data_sim.time>=np.datetime64(start)).data & (data_sim.time<=np.datetime64(end)).data)[0]

# cut out time period requested
if time_period != '':
    start = time_period.split(':')[0]
    start = datetime.datetime(np.int(start.split('-')[0]),np.int(start.split('-')[1]),np.int(start.split('-')[2]),0,0)
    end   = time_period.split(':')[1]
    end   = datetime.datetime(np.int(end.split('-')[0]),np.int(end.split('-')[1]),np.int(end.split('-')[2]),0,0)
    idx_period = np.where((times>=start) & (times<=end))[0]
else:
    start = times[0]
    end   = times[-1]
    idx_period = np.arange(len(times))

# find station
idx_station = np.where(station_id==station)[0]
if len(idx_station) == 0:
    print('Stations found: ',station_id)
    raise ValueError('Station not found in input file!')
else:
    idx_station = idx_station[0]


stop


times = times[idx_period]
Qobs = var[idx_station,idx_period]
    

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
nrow        = 4          # # of rows of subplots per figure
hspace      = 0.05        # x-space between subplots
vspace      = 0.35/nrow    # y-space between subplots
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


    
# -------------------------------------------------------------------------
# (1) streamflow data
# -------------------------------------------------------------------------

if ( model_qsim_files == 'None'):
    # --------
    # only observations
    # --------
    iplot += 1
    sub    = fig.add_axes(position(nrow,ncol,iplot,hspace=hspace,vspace=vspace) )
    line1 = sub.plot(times,Qobs,
                        color = 'gray',
                        linestyle='-',
                        linewidth=lwidth,
                        label=str2tex('$Q_{obs}$'))

    if usetex:
        xlab   = r'$\mathrm{time}'
        ylab   = r'$\mathrm{'+var_longname+'} \; ['+var_unit.replace('**','^').replace('-1','{-1}')+']$'
    else:
        xlab   = r'time'
        ylab   = r''+var_longname+' [$'+var_unit.replace('**','^').replace('-1','{-1}')+'$]'
        
    plt.setp(sub, xlabel=xlab)
    plt.setp(sub, ylabel=ylab)

    # text for gauge ID
    textbox_x = 0.5
    textbox_y = 1.00
    sub.text(textbox_x, textbox_y, station_id[idx_station], transform=sub.transAxes,
                         rotation=0, fontsize=textsize,
                         horizontalalignment='center', verticalalignment='bottom')
    textbox_x = 0.5
    textbox_y = 1.14

    # text for gauge name
    # insert line break if station_info (without spaces) is longer than ~20 charachters
    if nstations/ncol <= 8:
        tmp = station_info[idx_station].split(':')[2].strip().replace('"','').split(' ')
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
    plt.setp(sub,xlim=[mdates.date2num(start),mdates.date2num(end)])
    plt.setp(sub,ylim=[0.1,np.max(var[idx_station,idx_period])*1.05])

    # legend
    if (iplot > 0):
        illxbbox     =  0.0        # x-anchor legend bounding box
        illybbox     =  1.0        # y-anchor legend bounding box
        locat       = 'upper left' #'upper left'   
        sub.legend(frameon=frameon,
                       labelspacing=llrspace, handletextpad=llhtextpad, handlelength=llhlength,
                       loc=locat, bbox_to_anchor=(illxbbox,illybbox), scatterpoints=1, numpoints=1,fontsize=textsize)

    # rotate x-ticks
    plt.xticks(rotation=45)
        
else:
    # --------
    # observations and simulations
    # --------
    for model in dicts_qsim.keys():

        iplot += 1
        sub    = fig.add_axes(position(nrow,ncol,iplot,hspace=hspace,vspace=vspace) )
    
        line1 = sub.plot(times,Qobs,
                        color = 'gray',
                        linestyle='-',
                        linewidth=lwidth,
                        label=str2tex('$Q_{obs}$'))
                        #linestyle='None',
                        #linewidth=lwidth, 
                        #marker='o', markeredgecolor=lcols[0], markerfacecolor='None',
                        #markersize=msize*2.5, markeredgewidth=lwidth)
        Qsim = dicts_qsim[model][station]
        line1 = sub.plot(dicts_dates[model][station],Qsim,
                        linestyle='-',
                        linewidth=lwidth,
                        label=str2tex('$Q_{sim}^{'+model+'}$'))

        # period 2001-2010:
        sstart = datetime.datetime(2001, 1, 1,0,0)
        eend   = datetime.datetime(2010,12,31,0,0)
        tt_obs = np.where((times>=sstart) & (times<=eend))[0]
        tt_sim = np.where((dicts_dates[model][station]>=sstart) & (dicts_dates[model][station]<=eend))[0]
        # check that timesteps are in both time series
        ttt_obs = np.array([ tt for tt in tt_obs if times[tt] in dicts_dates[model][station] ])
        ttt_sim = np.array([ tt for tt in tt_sim if dicts_dates[model][station][tt] in times[ttt_obs] ])
        # calculate NSE
        performance = float(errormeasures.nse(Qobs[ttt_obs],Qsim[ttt_sim]))
        performance = "NSE(Q) = {0:6.3f} (2001-2010)".format(performance)
        # text for performance
        textbox_x = 0.98
        textbox_y = 0.98
        sub.text(textbox_x, textbox_y, performance, transform=sub.transAxes,
                 rotation=0, fontsize=textsize-2,
                 horizontalalignment='right', verticalalignment='top')
        
        # period 2011-2016:
        sstart = datetime.datetime(2011, 1, 1,0,0)
        eend   = datetime.datetime(2016,12,31,0,0)
        tt_obs = np.where((times>=sstart) & (times<=eend))[0]
        tt_sim = np.where((dicts_dates[model][station]>=sstart) & (dicts_dates[model][station]<=eend))[0]
        # check that timesteps are in both time series
        ttt_obs = np.array([ tt for tt in tt_obs if times[tt] in dicts_dates[model][station] ])
        ttt_sim = np.array([ tt for tt in tt_sim if dicts_dates[model][station][tt] in times[ttt_obs] ])
        # calculate NSE
        performance = float(errormeasures.nse(Qobs[ttt_obs],Qsim[ttt_sim]))
        performance = "NSE(Q) = {0:6.3f} (2011-2016)".format(performance)

        # text for performance
        textbox_x = 0.98
        textbox_y = 0.9
        sub.text(textbox_x, textbox_y, performance, transform=sub.transAxes,
                 rotation=0, fontsize=textsize-2,
                 horizontalalignment='right', verticalalignment='top')

        if usetex:
            xlab   = r'$\mathrm{time}'
            ylab   = r'$\mathrm{'+var_longname+'} \; ['+var_unit.replace('**','^').replace('-1','{-1}')+']$'
        else:
            xlab   = r'time'
            ylab   = r''+var_longname+' [$'+var_unit.replace('**','^').replace('-1','{-1}')+'$]'
            
        plt.setp(sub, xlabel=xlab)
        plt.setp(sub, ylabel=ylab)

        # text for gauge ID
        textbox_x = 0.5
        textbox_y = 1.00
        sub.text(textbox_x, textbox_y, station_id[idx_station], transform=sub.transAxes,
                             rotation=0, fontsize=textsize,
                             horizontalalignment='center', verticalalignment='bottom')
        textbox_x = 0.5
        textbox_y = 1.14

        # text for gauge name
        # insert line break if station_info (without spaces) is longer than ~20 charachters
        if nstations/ncol <= 8:
            tmp = station_info[idx_station].split(':')[2].strip().replace('"','').split(' ')
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
        plt.setp(sub,xlim=[mdates.date2num(start),mdates.date2num(end)])
        plt.setp(sub,ylim=[0.1,np.max(var[idx_station,idx_period])*1.05])

        # legend
        if (iplot > 0):
            illxbbox     =  0.0        # x-anchor legend bounding box
            illybbox     =  1.0        # y-anchor legend bounding box
            locat       = 'upper left' #'upper left'   
            sub.legend(frameon=frameon,
                           labelspacing=llrspace, handletextpad=llhtextpad, handlelength=llhlength,
                           loc=locat, bbox_to_anchor=(illxbbox,illybbox), scatterpoints=1, numpoints=1,fontsize=textsize)

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


