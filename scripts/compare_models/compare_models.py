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

#    run compare_models.py -i '../../data/objective_1/lake-erie//model/GEM-Hydro/gem-hydro_phase_1_objective_1.nc ../../data/objective_1/lake-erie//model/GR4J-Raven-lp/gr4j-raven-lp_phase_1_objective_1.nc ../../data/objective_1/lake-erie//model/GR4J-Raven-sd/gr4j-raven-sd_phase_1_objective_1.nc ../../data/objective_1/lake-erie//model/HYPE/hype_phase_1_objective_1.nc ../../data/objective_1/lake-erie//model/LBRM/lbrm_phase_1_objective_1.nc ../../data/objective_1/lake-erie//model/MESH-CLASS/mesh-class_phase_1_objective_1.nc ../../data/objective_1/lake-erie//model/MESH-SVS/mesh-svs_phase_1_objective_1.nc ../../data/objective_1/lake-erie//model/ML-ConvLSTM-DEM/ml-convlstm-dem_phase_1_objective_1.nc ../../data/objective_1/lake-erie//model/ML-ConvLSTM-LC-DEM/ml-convlstm-lc-dem_phase_1_objective_1.nc ../../data/objective_1/lake-erie//model/ML-ConvLSTM-LC/ml-convlstm-lc_phase_1_objective_1.nc ../../data/objective_1/lake-erie//model/ML-ConvLSTM/ml-convlstm_phase_1_objective_1.nc ../../data/objective_1/lake-erie//model/ML-LinReg/ml-linreg_phase_1_objective_1.nc ../../data/objective_1/lake-erie//model/ML-XGBoost/ml-xgboost_phase_1_objective_1.nc ../../data/objective_1/lake-erie//model/SWAT/swat_phase_1_objective_1.nc ../../data/objective_1/lake-erie//model/VIC-GRU/vic-gru_phase_1_objective_1.nc ../../data/objective_1/lake-erie//model/VIC/vic_phase_1_objective_1.nc ../../data/objective_1/lake-erie//model/WATFLOOD/watflood_phase_1_objective_1.nc' -a 2013-01-01:2014-12-31 -p compare_models_phase_1_objective_1_lake-erie_2013-01-01_2014-12-31.pdf -y

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
nosorty     = False

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
parser.add_argument('-y', '--nosorty', action='store_true', default=nosorty, dest="nosorty",
                    help="If given, y-axis will not be sorted.")

args        = parser.parse_args()
input_files = (args.input_files[0]).split(' ')
pngbase     = args.pngbase
pdffile     = args.pdffile
usetex      = args.usetex
time_period = args.time_period[0]
nosorty     = args.nosorty

del parser, args

dicts_nse     = {}
dicts_rmse    = {}
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
    dict_rmse    = {}
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
        rmse    = float(errormeasures.rmse(Qobs,Qsim).data)
        pbias   = float(errormeasures.pbias(Qobs,Qsim).data)
        ttidx = np.where((Qsim > 0.0) & (Qobs > 0.0))[0]
        lognse  = float(errormeasures.nse(np.log(Qobs[ttidx]),np.log(Qsim[ttidx])).data)
        sqrtnse = float(errormeasures.nse(np.sqrt(Qobs),np.sqrt(Qsim)).data)
        
        dict_nse[stat_id]     = nse
        dict_rmse[stat_id]    = rmse
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
    dicts_rmse[input_file.split('/')[-1].split('_')[0]]    = dict_rmse
    dicts_pbias[input_file.split('/')[-1].split('_')[0]]   = dict_pbias
    dicts_lognse[input_file.split('/')[-1].split('_')[0]]  = dict_lognse
    dicts_sqrtnse[input_file.split('/')[-1].split('_')[0]] = dict_sqrtnse
    dicts_dates[input_file.split('/')[-1].split('_')[0]]   = dict_dates
    dicts_qobs[input_file.split('/')[-1].split('_')[0]]    = dict_qobs
    dicts_qsim[input_file.split('/')[-1].split('_')[0]]    = dict_qsim
    dicts_info[input_file.split('/')[-1].split('_')[0]]    = dict_info

# print("dicts_nse  = ",dicts_nse)




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
nrow        = 6           # # of rows of subplots per figure
hspace      = 0.05        # x-space between subplots
vspace      = 0.04        # y-space between subplots
right       = 0.9         # right space on page
textsize    = 14           # standard text size
textsize_clock = 0.6*textsize        # standard text size
dxabc       = 0.95        # % of (max-min) shift to the right from left y-axis for a,b,c,... labels
dyabc       = 0.9         # % of (max-min) shift up from lower x-axis for a,b,c,... labels
dxabc_clock = 2.3         # % of (max-min) shift to the right from left y-axis for a,b,c,... labels
dyabc_clock = 0           # % of (max-min) shift up from lower x-axis for a,b,c,... labels
dxsig       = 1.23        # % of (max-min) shift to the right from left y-axis for a,b,c,... labels
dysig       = -0.05       # % of (max-min) shift up from lower x-axis for a,b,c,... labels

lwidth      = 1.0         # linewidth
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
    c = color.get_brewer('BlueWhiteOrangeRed', rgb=True)
    c = c[2:] # drop first two
    c = c[::-1] # reverse colors

    c = color.get_brewer('rdylbu11', rgb=True)

    
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

import scipy
import scipy.cluster.hierarchy as sch
# from scipy.spatial.distance import squareform
import scipy.spatial.distance as dist



model_order = ['ml-linreg', 'ml-convlstm', 'ml-convlstm-dem', 'ml-convlstm-lc', 'ml-convlstm-lc-dem', 'ml-lstm', 'ml-ea-lstm', 'ml-xgboost', 'lbrm', 'gr4j-raven-lp', 'gr4j-raven-sd', 'hmets-raven-lp', 'swat', 'mhm-waterloo', 'mhm-ufz', 'hype', 'vic', 'vic-gru', 'gem-hydro', 'mesh-svs', 'mesh-class', 'watflood', 'wrf-hydro']

# lines will appear after all models of each group
models_group = [['ml-linreg', 'ml-convlstm', 'ml-convlstm-dem', 'ml-convlstm-lc', 'ml-convlstm-lc-dem', 'ml-lstm', 'ml-ea-lstm', 'ml-xgboost'],
                ['lbrm', 'gr4j-raven-lp', 'gr4j-raven-sd', 'hmets-raven-lp', 'swat', 'mhm-waterloo', 'mhm-ufz']]

models = np.sort(dicts_nse.keys())
nmodels = np.shape(models)[0]
# sort models in model_order and append models not existing in model_order at the end
# (these should be only new models where I forgot to add them to model_order yet)
models = np.array([ imodel for imodel in model_order if imodel in models]+[ imodel for imodel in models if not(imodel in model_order) ])
gauges = np.sort(dicts_nse[models[0]].keys())

# only models in groups that actually are available
models_group = [ [ list(models).index(imodel) for imodel in igroup if imodel in models ] for igroup in models_group ]
# find last model in group --> to plot line after
lines_after_model = [ np.max(igroup) for igroup in models_group if len(igroup) > 0 ]

ngauges = np.shape(gauges)[0]



# ----------------------------------
# Figure NSE
# ----------------------------------
ifig = 0

ifig += 1
iplot = 0
print('Plot - Fig ', ifig, ' ::  heatmap NSE')
fig = plt.figure(ifig)

nse_results = np.array([ [ dicts_nse[imodel][igauge] for igauge in gauges ] for imodel in models ])


# median NSE
median_NSE = np.array([ np.median(nse_results[imodel,:]) for imodel in np.arange(nmodels) ])

# truncation of really bad results
nse_results_truncated = copy.deepcopy(nse_results)
min_nse = 0.0
max_nse = 0.2 * (np.int(np.max(nse_results)*5.0)+1)   # closest to [..., 0.6, 0.8, 1.0] to spread colorbar a bit
max_nse = 1.0
# print('max_nse = ',max_nse)
nse_results_truncated[np.where(nse_results_truncated < min_nse)] = min_nse  # NSE will be truncated to min_nse


[ax0_x, ax0_y, ax0_w, ax0_h] = [0.2,0.1,0.01,0.01]

# Cross-Validation: Clustering of columns
d2   = dist.pdist(nse_results_truncated.T)
D2   = dist.squareform(d2)
ax2  = fig.add_axes([ax0_x, ax0_y, ax0_w, ax0_h], frame_on=False) # [x,y,w,h]
Y2   = sch.linkage(D2, method='single', metric='euclidean')
Z2   = sch.dendrogram(Y2,color_threshold=np.inf)
ind2 = sch.fcluster(Y2,0.7*max(Y2[:,2]),'distance') ### This is the default behavior of dendrogram
ax2.set_xticks([]) ### Hides ticks
ax2.set_yticks([]) ### Hides ticks

[ax0_x, ax0_y, ax0_w, ax0_h] = [0.2,0.1,0.01,0.01]

# Cross-Validation: Clustering of rows
d1   = dist.pdist(nse_results_truncated)
D1   = dist.squareform(d1)
ax1  = fig.add_axes([ax0_x, ax0_y, ax0_w, ax0_h], frame_on=False) # [x,y,w,h]
Y1   = sch.linkage(D1, method='average', metric='euclidean') 
Z1   = sch.dendrogram(Y1,color_threshold=np.inf)
ind1 = sch.fcluster(Y1,0.7*max(Y1[:,2]),'distance') ### This is the default behavior of dendrogram
ax1.set_xticks([]) ### Hides ticks
ax1.set_yticks([]) ### Hides ticks


# Cross-Validation: Sort the matrix new
nse_results_sort = copy.deepcopy(nse_results_truncated)
# Sort columns
idx2             = Z2['leaves']     ### apply the clustering for the array-dendrograms to the actual matrix data
nse_results_sort = nse_results_sort[:,idx2]
ind2             = ind2[idx2]       ### JMJM reorder the flat cluster to match the order of the leaves the dendrogram

# Sort rows
if nosorty:
    idx1 = np.arange(nmodels)[::-1]
else:
    idx1 = Z1['leaves']     ### apply the clustering for the gene-dendrograms to the actual matrix data
nse_results_sort = nse_results_sort[idx1,:]   # xt is transformed x
median_NSE       = median_NSE[idx1]
ind1             = ind1[idx1]       ### JMJM  reorder the flat cluster to match the order of the leaves the dendrogram

# Plot distance matrix.
[ax1_x, ax1_y, ax1_w, ax1_h] = [0.2,0.1,0.7,0.3]
sub  = fig.add_axes([ax1_x, ax1_y, ax1_w, ax1_h], axisbg='none',frame_on=False)
im = sub.matshow(nse_results_sort, aspect='auto', origin='lower', cmap=cmap)
im.set_clim(min_nse,max_nse)

# ticklabels
sub.set_xticks(np.arange(ngauges))
if ngauges < 32:
    sub.set_xticklabels(gauges[idx2],rotation=90)
elif ngauges < 70:
    sub.set_xticklabels(gauges[idx2],rotation=90,fontsize='x-small')
else:
    sub.set_xticklabels(gauges[idx2],rotation=90,fontsize='xx-small')
sub.set_yticks(np.arange(nmodels))
sub.set_yticklabels(models[idx1])


# label text for median NSE
sub.text(  ngauges+0.5, nmodels-0.5, 'median NSE',
                        ha = 'left', va = 'bottom', rotation=90,
                        fontsize=textsize )
    
# text for median NSE values
for imodel in np.arange(nmodels):
    # print("y: ",imodel,"   NSE: ",astr(median_NSE[imodel],prec=2))
    sub.text(  ngauges, imodel, astr(median_NSE[imodel],prec=2),
                        ha = 'left', va = 'center',
                        fontsize=textsize )

# draw line after model groups
for iline in lines_after_model:
    sub.plot([-0.5,ngauges-0.5],[nmodels-iline-1.5,nmodels-iline-1.5],linewidth=lwidth,color='black')


# Legend
[ax4_x, ax4_y, ax4_w, ax4_h] = [0.2,0.07,0.7,0.02]
axcb    = fig.add_axes([ax4_x, ax4_y, ax4_w, ax4_h],axisbg='none',frame_on=False)
nticks  = 5.0
ticks   = [ min_nse + itick/(nticks-1)*(max_nse-min_nse) for itick in np.arange(nticks) ]
# print('ticks = ',ticks)
cb      = mpl.colorbar.ColorbarBase(axcb, cmap=cmap, alpha=1.0, orientation='horizontal',ticks=ticks,extend='min',
                                        norm=mpl.colors.Normalize(vmin=min_nse, vmax=max_nse))
cb.set_label('NSE')
labels = [ astr(itick,prec=2) for itick in ticks ]
cb.set_ticklabels(labels)

if (outtype == 'pdf'):
    pdf_pages.savefig(fig)
    plt.close(fig)
elif (outtype == 'png'):
    pngfile = pngbase+"{0:04d}".format(ifig)+".png"
    fig.savefig(pngfile, transparent=transparent, bbox_inches=bbox_inches, pad_inches=pad_inches)
    plt.close(fig)

# ----------------------------------
# Figure PBIAS
# ----------------------------------
ifig = 0

ifig += 1
iplot = 0
print('Plot - Fig ', ifig, ' ::  heatmap PBIAS')
fig = plt.figure(ifig)

pbias_results = np.array([ [ dicts_pbias[imodel][igauge] for igauge in gauges ] for imodel in models ])


# median PBIAS
median_PBIAS = np.array([ np.median(pbias_results[imodel,:]) for imodel in np.arange(nmodels) ])
median_abs_PBIAS = np.array([ np.median(np.abs(pbias_results[imodel,:])) for imodel in np.arange(nmodels) ])

# truncation of really bad results
pbias_results_truncated = copy.deepcopy(pbias_results)
min_pbias = - np.percentile(np.abs(pbias_results),90)
max_pbias =   np.percentile(np.abs(pbias_results),90)
min_pbias = - 30.0
max_pbias =   30.0
# print('max_pbias = ',max_pbias)
pbias_results_truncated[np.where(pbias_results_truncated < min_pbias)] = min_pbias  # PBIAS will be truncated to min_pbias
pbias_results_truncated[np.where(pbias_results_truncated > max_pbias)] = max_pbias  # PBIAS will be truncated to max_pbias

# take sorting from NSE ....

# [ax0_x, ax0_y, ax0_w, ax0_h] = [0.2,0.1,0.01,0.01]

# # Cross-Validation: Clustering of columns
# d2   = dist.pdist(pbias_results_truncated.T)
# D2   = dist.squareform(d2)
# ax2  = fig.add_axes([ax0_x, ax0_y, ax0_w, ax0_h], frame_on=False) # [x,y,w,h]
# Y2   = sch.linkage(D2, method='single', metric='euclidean')
# Z2   = sch.dendrogram(Y2,color_threshold=np.inf)
# ind2 = sch.fcluster(Y2,0.7*max(Y2[:,2]),'distance') ### This is the default behavior of dendrogram
# ax2.set_xticks([]) ### Hides ticks
# ax2.set_yticks([]) ### Hides ticks

# [ax0_x, ax0_y, ax0_w, ax0_h] = [0.2,0.1,0.01,0.01]

# # Cross-Validation: Clustering of rows
# d1   = dist.pdist(pbias_results_truncated)
# D1   = dist.squareform(d1)
# ax1  = fig.add_axes([ax0_x, ax0_y, ax0_w, ax0_h], frame_on=False) # [x,y,w,h]
# Y1   = sch.linkage(D1, method='average', metric='euclidean') 
# Z1   = sch.dendrogram(Y1,color_threshold=np.inf)
# ind1 = sch.fcluster(Y1,0.7*max(Y1[:,2]),'distance') ### This is the default behavior of dendrogram
# ax1.set_xticks([]) ### Hides ticks
# ax1.set_yticks([]) ### Hides ticks


# Cross-Validation: Sort the matrix new
pbias_results_sort = copy.deepcopy(pbias_results_truncated)
# Sort columns
idx2             = Z2['leaves']     ### apply the clustering for the array-dendrograms to the actual matrix data
pbias_results_sort = pbias_results_sort[:,idx2]
ind2             = ind2[idx2]       ### JMJM reorder the flat cluster to match the order of the leaves the dendrogram

# Sort rows
if nosorty:
    idx1 = np.arange(nmodels)[::-1]
else:
    idx1 = Z1['leaves']     ### apply the clustering for the gene-dendrograms to the actual matrix data
pbias_results_sort = pbias_results_sort[idx1,:]   # xt is transformed x
median_PBIAS       = median_PBIAS[idx1]
median_abs_PBIAS   = median_abs_PBIAS[idx1]
ind1               = ind1[idx1]       ### JMJM  reorder the flat cluster to match the order of the leaves the dendrogram

# Plot distance matrix.
[ax1_x, ax1_y, ax1_w, ax1_h] = [0.2,0.1,0.7,0.3]
sub  = fig.add_axes([ax1_x, ax1_y, ax1_w, ax1_h], axisbg='none',frame_on=False)
im = sub.matshow(pbias_results_sort, aspect='auto', origin='lower', cmap=cmap)
im.set_clim(min_pbias,max_pbias)

# ticklabels
sub.set_xticks(np.arange(ngauges))
if ngauges < 32:
    sub.set_xticklabels(gauges[idx2],rotation=90)
elif ngauges < 70:
    sub.set_xticklabels(gauges[idx2],rotation=90,fontsize='x-small')
else:
    sub.set_xticklabels(gauges[idx2],rotation=90,fontsize='xx-small')
sub.set_yticks(np.arange(nmodels))
sub.set_yticklabels(models[idx1])


# label text for median PBIAS
sub.text(  ngauges*(1+0.025/0.8), nmodels-0.5, 'median |PBIAS|',
                        ha = 'center', va = 'bottom', rotation=90,
                        fontsize=textsize )
    
# text for median of absolute PBIAS values
for imodel in np.arange(nmodels):
    sub.text( ngauges*(1+0.05/0.8), imodel, astr(median_abs_PBIAS[imodel],prec=1),
                        ha = 'right', va = 'center',
                        fontsize=textsize )

# draw line after model groups
for iline in lines_after_model:
    sub.plot([-0.5,ngauges-0.5],[nmodels-iline-1.5,nmodels-iline-1.5],linewidth=lwidth,color='black')


# Legend
[ax4_x, ax4_y, ax4_w, ax4_h] = [0.2,0.07,0.7,0.02]
axcb    = fig.add_axes([ax4_x, ax4_y, ax4_w, ax4_h],axisbg='none',frame_on=False)
nticks  = 5.0
ticks   = [ min_pbias + itick/(nticks-1)*(max_pbias-min_pbias) for itick in np.arange(nticks) ]
# print('ticks = ',ticks)
cb      = mpl.colorbar.ColorbarBase(axcb, cmap=cmap, alpha=1.0, orientation='horizontal',ticks=ticks,extend='both',
                                        norm=mpl.colors.Normalize(vmin=min_pbias, vmax=max_pbias))
cb.set_label('PBIAS [%]')
labels = [ astr(itick,prec=1) for itick in ticks ]
cb.set_ticklabels(labels)

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


