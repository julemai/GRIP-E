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
#    python plot_nc_model_output.py -i '../../data/objective_1/model/LBRM/lbrm_phase_0_objective_1.nc ../../data/objective_1/model/GEM-Hydro/gem-hydro_phase_0_objective_1.nc ../../data/objective_1/model/VIC-GRU/vic-gru_phase_0_objective_1.nc' -p test.pdf
#
#    python plot_nc_model_output.py -i '../../data/objective_1/model/VIC-GRU/vic-gru_phase_0_objective_1.nc' -a '2011-01-01:2014-12-31' -p test.pdf

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

input_files = ['vic-gru_phase_0_objective_1.nc', 'gem-hydro_phase_0_objective_1.nc']
pngbase     = ''
pdffile     = 'test.pdf'
usetex      = False
time_period = ''

parser      = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
              description='''Convert LBRM raw streamflow model outputs into NetDF format (consistent across all models in GRIP-E).''')
parser.add_argument('-i', '--input_files', action='store',
                    default=input_files, dest='input_files', metavar='input_files', nargs=1,
                    help='Name of NetCDF file containing standard model outputs.')
parser.add_argument('-a', '--time_period', action='store',
                    default=time_period, dest='time_period', metavar='time_period', nargs=1,
                    help='Time period to plot and calculate performance of (format: YYYY-MM-DD:YYYY-MM-DD; default: all time points).')
parser.add_argument('-g', '--pngbase', action='store',
                    default=pngbase, dest='pngbase', metavar='pngbase',
                    help='Name basis for png output files (default: open screen window).')
parser.add_argument('-p', '--pdffile', action='store',
                    default=pdffile, dest='pdffile', metavar='pdffile',
                    help='Name of pdf output file (default: open screen window).')
parser.add_argument('-t', '--usetex', action='store_true', default=usetex, dest="usetex",
                    help="Use LaTeX to render text in pdf.")

args        = parser.parse_args()
input_files = (args.input_files[0]).split(' ')
pngbase     = args.pngbase
pdffile     = args.pdffile
usetex      = args.usetex
time_period = args.time_period[0]

del parser, args

dicts_nse     = {}
dicts_pbias   = {}
dicts_lognse  = {}
dicts_sqrtnse = {}
dicts_dates   = {}
dicts_qobs    = {}
dicts_qsim    = {}
dicts_info    = {}  # station info, i.e. station long name

for iinput_file,input_file in enumerate(input_files):
    
    # simulated streamflow
    data_sim = xr.open_dataset(input_file)

    # observed streamflow
    data_obs = xr.open_dataset(os.path.dirname(input_file)+'/../../netcdf/all_gauges.nc')

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
    dict_nse     = {}
    dict_pbias   = {}
    dict_lognse  = {}
    dict_sqrtnse = {}
    dict_dates   = {}
    dict_qobs    = {}
    dict_qsim    = {}
    dict_info    = {}
    nstations    = len(data_obs.station_id)
    for istation in range(nstations):
        stat_id = str(data_obs.station_id[istation].data)
        Qobs    = data_obs.Q[istation,idx_time]
        Qsim    = data_sim.Q[idx_station[istation],:]
        dates   = data_sim.time.data
        
        if time_period != '':
            Qobs  = Qobs[idx_period]
            Qsim  = Qsim[idx_period]
            dates = dates[idx_period]
        nse     = float(errormeasures.nse(Qobs,Qsim).data)
        pbias   = float(errormeasures.pbias(Qobs,Qsim).data)
        ttidx = np.where((Qsim > 0.0) & (Qobs > 0.0))[0]
        lognse  = float(errormeasures.nse(np.log(Qobs[ttidx]),np.log(Qsim[ttidx])).data)
        sqrtnse = float(errormeasures.nse(np.sqrt(Qobs),np.sqrt(Qsim)).data)
        
        dict_nse[stat_id]     = nse
        dict_pbias[stat_id]   = pbias
        dict_lognse[stat_id]  = lognse
        dict_sqrtnse[stat_id] = sqrtnse
        dict_qobs[stat_id]    = Qobs.data
        dict_qsim[stat_id]    = Qsim.data
        # convert from numpy.datetime64 to datetime.datetime
        dict_dates[stat_id]   = np.array([ pd.Timestamp(itime).to_pydatetime() for itime in dates ])
        # station long name
        dict_info[stat_id]    = data_obs['station_info'][istation]

    # save everything in dictionaries
    dicts_nse[input_file.split('/')[-1].split('_')[0]]     = dict_nse
    dicts_pbias[input_file.split('/')[-1].split('_')[0]]   = dict_pbias
    dicts_lognse[input_file.split('/')[-1].split('_')[0]]  = dict_lognse
    dicts_sqrtnse[input_file.split('/')[-1].split('_')[0]] = dict_sqrtnse
    dicts_dates[input_file.split('/')[-1].split('_')[0]]   = dict_dates
    dicts_qobs[input_file.split('/')[-1].split('_')[0]]    = dict_qobs
    dicts_qsim[input_file.split('/')[-1].split('_')[0]]    = dict_qsim
    dicts_info[input_file.split('/')[-1].split('_')[0]]    = dict_info

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
nrow        = 8           # # of rows of subplots per figure
hspace      = 0.05        # x-space between subplots
vspace      = 0.04        # y-space between subplots
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

for imodel in [ ii.split('/')[-2] for ii in input_files ]:
    
    # -------------------------------------------------------------------------
    # Fig 1 :: plot selected stations for each model individually
    # -------------------------------------------------------------------------
    ifig += 1
    iplot = 0
    print('Plot - Fig ', ifig, ' ::  ',imodel)
    fig = plt.figure(ifig)

    imodel_lc = imodel.lower()

    sel_stations = ['02GG002','02GG003','04196800']
    max_obs = np.max(np.array([ np.nanmax(dicts_qobs[imodel_lc][key]) for key in sel_stations ]))
    max_sim = np.max(np.array([ np.nanmax(dicts_qsim[imodel_lc][key]) for key in sel_stations ]))

    for iistation,istation in enumerate(sel_stations):
        
        # -------------------------------------------------------------------------
        # (1) streamflow data
        # -------------------------------------------------------------------------
        iplot += 1
        sub    = fig.add_axes(position(nrow,ncol,iplot,hspace=hspace,vspace=vspace) )
        line1 = sub.plot(dicts_dates[imodel_lc][istation],dicts_qsim[imodel_lc][istation])    
        plt.setp(line1, linestyle='-',
                        linewidth=lwidth*2, color=lcols[0], label=str2tex('$Q_{sim}^{'+imodel+'}$'))

        # line2 = sub.plot(dicts_dates[imodel_lc][istation],dicts_qobs[imodel_lc][istation])    
        # plt.setp(line2, linestyle='-',
        #                 linewidth=lwidth*2, color=lcols[1], label=str2tex('$Q_{obs}$'))

        mark1 = sub.plot(dicts_dates[imodel_lc][istation][::5],dicts_qobs[imodel_lc][istation][::5])
        plt.setp(mark1, linestyle='None',
                        marker='o', markeredgecolor=lcols[0], markerfacecolor='None',
                        markersize=msize*2.5, markeredgewidth=lwidth, label=str2tex('$Q_{obs}$'))

        if usetex:
            xlab   = r'' #r'$\mathrm{time}'
            ylab   = r'$\mathrm{'+var_longname+'} \; ['+var_unit.replace('**','^').replace('-1','{-1}')+']$'
        else:
            xlab   = r'' #r'time'
            ylab   = r''+var_longname+' [$'+var_unit.replace('**','^').replace('-1','{-1}')+'$]'
            
        if iplot == 1:
            sub.set_title('')
        if (iplot-1)//ncol+1 == len(sel_stations):        # last row
            plt.setp(sub, xlabel=xlab)
        else:
            sub.tick_params(axis='x',labelbottom='off')
        if iplot%ncol == 1:                  # first column
            plt.setp(sub, ylabel=ylab)

        # legend
        if (iplot > 0):
            illxbbox     =  1.0        # x-anchor legend bounding box
            illybbox     =  1.0        # y-anchor legend bounding box
            locat       = 'upper right' #'upper left'   
            sub.legend(frameon=frameon,
                           labelspacing=llrspace, handletextpad=llhtextpad, handlelength=llhlength,
                           loc=locat, bbox_to_anchor=(illxbbox,illybbox), scatterpoints=1, numpoints=1,fontsize=textsize)

        # text for performance metrics
        textbox_x = 0.5
        textbox_y = 1.00
        nse     = "NSE(Q) = {0:6.3f}".format(dicts_nse[imodel_lc][istation]).strip()
        pbias   = "PBIAS(Q) = {0:6.3f}".format(dicts_pbias[imodel_lc][istation]).strip()
        lognse  = "NSE(logQ) = {0:6.3f}".format(dicts_lognse[imodel_lc][istation]).strip()
        sqrtnse = "NSE(sqrtQ) = {0:6.3f}".format(dicts_sqrtnse[imodel_lc][istation]).strip()
        performance = nse+'  '+lognse+'  '+sqrtnse+'  '+pbias
        
        sub.text(textbox_x, textbox_y, performance, transform=sub.transAxes,
                             rotation=0, fontsize=textsize,
                             horizontalalignment='center', verticalalignment='bottom')

        textbox_x = 0.5
        textbox_y = 1.18
        # text for gauge name
        # insert line break if station_info (without spaces) is longer than ~100 charachters
        tmp = str(dicts_info[imodel_lc][istation].data).strip().replace('"','').split(' ')
        if len(' '.join(tmp)) > 100:
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
        plt.setp(sub,ylim=[0.1,max(max_obs,max_sim)*1.05])

        # rotate x-ticks
        plt.xticks(rotation=0)

    if (outtype == 'pdf'):
        pdf_pages.savefig(fig)
        plt.close(fig)
    elif (outtype == 'png'):
        pngfile = pngbase+"{0:04d}".format(ifig)+".png"
        fig.savefig(pngfile, transparent=transparent, bbox_inches=bbox_inches, pad_inches=pad_inches)
        plt.close(fig)


    # -------------------------------------------------------------------------
    # Fig 2 :: statistics as boxplot
    # -------------------------------------------------------------------------
    ifig += 1
    iplot = 0
    print('Plot - Fig ', ifig, ' ::  ',imodel)
    fig = plt.figure(ifig)

    nncol = 2
    vvspace = 0.02

    basins = dicts_nse[imodel_lc].keys()

    # -------------
    # NSE
    # -------------
    data   = [ dicts_nse[imodel_lc][kk] for kk in basins ]
    dmin   = np.min([ [ dicts_nse[mm][kk] for kk in dicts_nse[mm].keys() ] for mm in dicts_nse.keys() ])
    dmax   = np.max([ [ dicts_nse[mm][kk] for kk in dicts_nse[mm].keys() ] for mm in dicts_nse.keys() ])

    iplot += 1
    pos    = position(nrow,nncol,iplot,hspace=hspace,vspace=vvspace)
    pos[2] = 0.6 # width
    sub    = fig.add_axes(pos)
    line1  = sub.plot(basins,data,linestyle='None',
                        marker='o', markeredgecolor=lcols[0], markerfacecolor='None')
    sub.tick_params(axis='x',labelbottom='off')
    plt.setp(sub, ylabel='NSE(Q)')
    plt.setp(sub,ylim=[dmin,max(dmax,1.0)])

    # title
    sub.text(0.5, 1.02, imodel, transform=sub.transAxes,
                             rotation=0, fontsize=textsize+2,
                             horizontalalignment='center', verticalalignment='bottom')

    iplot += 1
    pos    = position(nrow,nncol,iplot,hspace=hspace,vspace=vvspace)
    pos[0] = 0.75  # left
    pos[2] = 0.2  # width
    sub    = fig.add_axes(pos)
    line1  = sub.boxplot(data)
    edge_color = 'black'
    for element in ['boxes', 'whiskers', 'fliers', 'means', 'medians', 'caps']:
        plt.setp(line1[element], color=edge_color)
    median = astr(np.median(data), prec=2)
    sub.text(1.1, np.median(data), median, #transform=sub.transAxes,
                             rotation=0, fontsize=textsize-2,
                             horizontalalignment='left', verticalalignment='center')
    sub.tick_params(axis='x',labelbottom='off')
    sub.tick_params(axis='y',labelleft='off')
    plt.setp(sub,ylim=[dmin,max(dmax,1.0)])

    # -------------
    # LOG NSE
    # -------------
    data   = [ dicts_lognse[imodel_lc][kk] for kk in basins ]
    dmin   = np.min([ [ dicts_lognse[mm][kk] for kk in dicts_lognse[mm].keys() ] for mm in dicts_lognse.keys() ])
    dmax   = np.max([ [ dicts_lognse[mm][kk] for kk in dicts_lognse[mm].keys() ] for mm in dicts_lognse.keys() ])

    iplot += 1
    pos    = position(nrow,nncol,iplot,hspace=hspace,vspace=vvspace)
    pos[2] = 0.6 # width
    sub    = fig.add_axes(pos)
    line1  = sub.plot(basins,data,linestyle='None',
                        marker='o', markeredgecolor=lcols[0], markerfacecolor='None')
    sub.tick_params(axis='x',labelbottom='off')
    plt.setp(sub, ylabel='NSE(log[Q])')
    plt.setp(sub,ylim=[dmin,max(dmax,1.0)])

    iplot += 1
    pos    = position(nrow,nncol,iplot,hspace=hspace,vspace=vvspace)
    pos[0] = 0.75  # left
    pos[2] = 0.2  # width
    sub    = fig.add_axes(pos)
    line1  = sub.boxplot(data)
    edge_color = 'black'
    for element in ['boxes', 'whiskers', 'fliers', 'means', 'medians', 'caps']:
        plt.setp(line1[element], color=edge_color)
    median = astr(np.median(data), prec=2)
    sub.text(1.1, np.median(data), median, #transform=sub.transAxes,
                             rotation=0, fontsize=textsize-2,
                             horizontalalignment='left', verticalalignment='center')
    sub.tick_params(axis='x',labelbottom='off')
    sub.tick_params(axis='y',labelleft='off')
    plt.setp(sub,ylim=[dmin,max(dmax,1.0)])

    # -------------
    # SQRT NSE
    # -------------
    data   = [ dicts_sqrtnse[imodel_lc][kk] for kk in basins ]
    dmin   = np.min([ [ dicts_lognse[mm][kk] for kk in dicts_sqrtnse[mm].keys() ] for mm in dicts_sqrtnse.keys() ])
    dmax   = np.max([ [ dicts_lognse[mm][kk] for kk in dicts_sqrtnse[mm].keys() ] for mm in dicts_sqrtnse.keys() ])

    iplot += 1
    pos    = position(nrow,nncol,iplot,hspace=hspace,vspace=vvspace)
    pos[2] = 0.6 # width
    sub    = fig.add_axes(pos)
    line1  = sub.plot(basins,data,linestyle='None',
                        marker='o', markeredgecolor=lcols[0], markerfacecolor='None')
    sub.tick_params(axis='x',labelbottom='off')
    plt.setp(sub, ylabel='NSE(sqrt[Q])')
    plt.setp(sub,ylim=[dmin,max(dmax,1.0)])

    iplot += 1
    # [left, bottom, width, height)
    pos    = position(nrow,nncol,iplot,hspace=hspace,vspace=vvspace)
    pos[0] = 0.75  # left
    pos[2] = 0.2  # width
    sub    = fig.add_axes(pos)
    line1  = sub.boxplot(data)
    edge_color = 'black'
    for element in ['boxes', 'whiskers', 'fliers', 'means', 'medians', 'caps']:
        plt.setp(line1[element], color=edge_color)
    median = astr(np.median(data), prec=2)
    sub.text(1.1, np.median(data), median, #transform=sub.transAxes,
                             rotation=0, fontsize=textsize-2,
                             horizontalalignment='left', verticalalignment='center')
    sub.tick_params(axis='x',labelbottom='off')
    sub.tick_params(axis='y',labelleft='off')
    plt.setp(sub,ylim=[dmin,max(dmax,1.0)])

    # -------------
    # PBIAS
    # -------------
    data   = [ dicts_pbias[imodel_lc][kk] for kk in basins ]
    dmin   = np.min([ [ dicts_pbias[mm][kk] for kk in dicts_pbias[mm].keys() ] for mm in dicts_pbias.keys() ])
    dmax   = np.max([ [ dicts_pbias[mm][kk] for kk in dicts_pbias[mm].keys() ] for mm in dicts_pbias.keys() ])

    iplot += 1
    pos    = position(nrow,nncol,iplot,hspace=hspace,vspace=vvspace)
    pos[2] = 0.6 # width
    sub    = fig.add_axes(pos)
    line1  = sub.plot(basins,data,linestyle='None',
                        marker='o', markeredgecolor=lcols[0], markerfacecolor='None')
    sub.tick_params(axis='x',labelbottom='on')
    plt.setp(sub, ylabel='PBIAS')
    plt.setp(sub,ylim=[dmin,dmax*1.02])
    # rotate x-ticks
    plt.xticks(rotation=90)

    iplot += 1
    pos    = position(nrow,nncol,iplot,hspace=hspace,vspace=vvspace)
    pos[0] = 0.75  # left
    pos[2] = 0.2  # width
    sub    = fig.add_axes(pos)
    line1  = sub.boxplot(data)
    edge_color = 'black'
    for element in ['boxes', 'whiskers', 'fliers', 'means', 'medians', 'caps']:
        plt.setp(line1[element], color=edge_color)
    median = astr(np.median(data), prec=2)
    sub.text(1.1, np.median(data), median, #transform=sub.transAxes,
                             rotation=0, fontsize=textsize-2,
                             horizontalalignment='left', verticalalignment='center')
    sub.tick_params(axis='x',labelbottom='off')
    sub.tick_params(axis='y',labelleft='off')
    plt.setp(sub,ylim=[dmin,dmax*1.02])


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


# # plot a selected station: "02GG003" but all model outputs
# sub.plot(dicts_dates[imodel_lc]['02GG003'],dicts_qobs[imodel_lc]['02GG003']) 
# sub.plot(dicts_dates[imodel_lc]['02GG003'],dicts_qsim[imodel_lc]['02GG003']) 

# # plot performance at all stations
# sub.plot(dicts_nse['vic-gru'].keys(),dicts_nse['vic-gru'].values(),linewidth=0.0,marker='o')

# # plot 3 stations per model individually
# # 


    



