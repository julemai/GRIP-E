#!/usr/bin/env python
from __future__ import print_function
"""

Plots gauge streamflow data

Run with::
      run plot_gauge_data_grouped.py -i ../../data/objective_1/great-lakes/gauge_info.csv ../../data/objective_2/great-lakes/gauge_info.csv -p test.pdf -t -v Q
      
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
years     = '2010:2014'

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
parser.add_argument('-i', '--inputfile', action='store', default=inputfile, dest='inputfile', nargs='+',
                    help="Gauge info file(s) conatining 'Watershed_NO', 'Calibration/Validation', and 'ID'.")
parser.add_argument('-y', '--years', action='store', default=years, dest='years', 
                    help="Years to plot in format YYYY:YYYY.")

args            = parser.parse_args()
pngbase         = args.pngbase
pdffile         = args.pdffile
usetex          = args.usetex
variable        = args.variable
inputfile       = args.inputfile
years           = args.years

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
import glob

import color                      # in lib/
from fread      import fread      # in lib/
from position   import position   # in lib/
from str2tex    import str2tex    # in lib/
from readnetcdf import readnetcdf # in lib/

# read CSV files
gaugeinfo = []
for gaugeinfo_file in inputfile:

    gaugeinfo_all = []
    
    with open(gaugeinfo_file) as f:
        content = f.readlines()

    idx_1 = content[0].strip().split(',').index('ID')
    idx_2 = content[0].strip().split(',').index('Name')
    idx_3 = content[0].strip().split(',').index('Calibration/Validation')
    idx_4 = content[0].strip().split(',').index('Watershed_NO')
    idx_5 = content[0].strip().split(',').index('Objective_info')
    idx = np.array([idx_1,idx_2,idx_3,idx_4,idx_5])

    for cc in content[1:]:
        ccc = ' '.join(cc.strip().split())

        # split at "," but not when in quotes
        import re
        PATTERN = re.compile(r'''((?:[^,"']|"[^"]*"|'[^']*')+)''')
        # print( PATTERN.split(ccc)[1::2] )
        ccc = PATTERN.split(ccc)[1::2]
        
        gaugeinfo_all.append(list(np.array(ccc)[idx]))

        
    # slim down gauge information to only gauges present (in right order)
    #station_id_all = [ ii[0] for ii in gaugeinfo_all ]
    #idx_stations   = [ station_id_all.index(ss) for ss in station_id ]

    if gaugeinfo != []:
        ids_new = [ ii[0] for ii in gaugeinfo_all ]  
        ids_old = [ ii[0] for ii in gaugeinfo ]  
        to_add = [ not(ii in ids_old) for ii in ids_new ]

        for ii,iadd in enumerate(to_add):

            if iadd:
                gaugeinfo += [ gaugeinfo_all[ii] ]

    else:
        idx_stations = range(len(gaugeinfo_all))
        gaugeinfo += [ gaugeinfo_all[idx] for idx in idx_stations ]   # gauge information of only the requested gauges: [ID,Name,Lat,Lon,Country,Drainage_area]
    
print("Number of (unique) gauges given in all inputfiles = ",len(gaugeinfo))
    
gaugeinfo = np.array(gaugeinfo)
groups    = np.unique(gaugeinfo[:,3])

lon_all           = []
lat_all           = []
var_all           = []
station_id_all    = []
station_info_all  = []
times_all         = []
varnames_all      = []
var_unit_all      = []
var_longname_all  = []
time_unit_all     = []
time_longname_all = []
reference_day_all = []
nstations_all     = []
idx_all           = []
for gaugeinfo_file in inputfile:

    datafiles = glob.glob('/'.join(gaugeinfo_file.split('/')[0:-1])+'/*/*/*.nc')

    for fname in datafiles:

        # Latlon
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
        year_start = np.int(years.split(':')[0])
        year_end   = np.int(years.split(':')[1])
        idx = np.where((times>=datetime.datetime(year_start,1,1,0,0)) & (times<=datetime.datetime(year_end+1,1,1,0,0)))[0]

        lon_all.append(           lon          )
        lat_all.append(           lat          )
        var_all.append(           var          )
        station_id_all.append(    station_id   )
        station_info_all.append(  station_info )
        times_all.append(         times        )
        varnames_all.append(      varnames     )
        var_unit_all.append(      var_unit     )
        var_longname_all.append(  var_longname )
        time_unit_all.append(     time_unit    )
        time_longname_all.append( time_longname)
        reference_day_all.append( reference_day)
        nstations_all.append(     nstations    )
        idx_all.append(           idx          )

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
vspace      = 0.03         # y-space between subplots
right       = 0.9         # right space on page
textsize    = 12           # standard text size
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
llybbox       = -0.30        # y-anchor legend bounding box
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
# Fig 1 :: for each group plot discharge (multiple gauge stations
# -------------------------------------------------------------------------
for group in groups:   # groups[-1:]:

    ifig += 1
    iplot = 0
    print('Plot - Fig ', ifig, ' ::  Group ',group)
    fig = plt.figure(ifig)

    iplot += 1
    sub    = fig.add_axes(position(nrow,ncol,iplot,hspace=hspace,vspace=vspace) )

    # which gauge IDs in this group
    gaugeinfo_group = gaugeinfo[np.where(gaugeinfo[:,3]==group)]
    gaugeID_group   = gaugeinfo_group[:,0]

    for iigauge,igauge in enumerate(gaugeID_group):

        # find where this gauge is available
        # '02JC008', '02KF011', '02JB013', '02JE027', '02KC018', '02KD002', '02KF005'
        idx_gauge = [ list(np.where(station_id_all[ifile]==igauge)[0]) for ifile in range(len(station_id_all)) ]      # e.g. [[22], [], [], []]

        # find data for this gauge
        found = False
        ifile = 0
        while not(found):
            if len(idx_gauge[ifile]) > 0:
                iidx = idx_gauge[ifile][0]   # take first appearance

                istation = iidx
                idx   = idx_all[ifile]          # timestep indexes for specific period
                times = times_all[ifile]        # time steps in this file
                var   = var_all[ifile]          # actual discharge values of this file
                var_longname = var_longname_all[ifile]
                var_unit     = var_unit_all[ifile]
                station_id   = station_id_all[ifile]
                
                found = True
            else:
                ifile += 1
                if (ifile == len(station_id_all)):
                    print("Station: ",igauge)
                    raise ValueError("Station does not seem to be available anywhere.")

        # -------------------------------------------------------------------------
        # (1) streamflow data
        # -------------------------------------------------------------------------

        # line1 = sub.plot(times[idx],var[istation,idx])
        line1 = sub.semilogy(times[idx],var[istation,idx])
    
        plt.setp(line1, linestyle='-',
                        linewidth=lwidth, color=cmap.colors[int(iigauge*256/len(gaugeID_group))],
                        label=str2tex(igauge+" : "+gaugeinfo_group[iigauge][1]+"  ("+gaugeinfo_group[iigauge][2]+", "+gaugeinfo_group[iigauge][4]+")",usetex=usetex))

        if usetex:
            xlab   = r'$\mathrm{time}'
            ylab   = r'$\mathrm{'+var_longname+'} \; ['+var_unit.replace('**','^').replace('-1','{-1}')+']$'
        else:
            xlab   = r'time'
            ylab   = r''+var_longname+' [$'+var_unit.replace('**','^').replace('-1','{-1}')+'$]'
            
        if iplot == 1:
            sub.set_title('')
        if (iplot-1)//ncol+1 == 1:        # first row
            plt.setp(sub, xlabel=xlab)
        else:
            sub.tick_params(axis='x',labelbottom='off')

        if nstations/ncol > 8:
            if iplot%ncol == 0  and ((iplot-1)//ncol+1)%1 == 0:                  # first column, every row
                plt.setp(sub, ylabel=ylab)
        else:
            if iplot%ncol == 0:                  # first column
                plt.setp(sub, ylabel=ylab)

        # text for group ID
        if iigauge == 1:
            textbox_x = 0.5
            textbox_y = 1.00
            sub.text(textbox_x, textbox_y, str2tex(group,usetex=usetex), transform=sub.transAxes,
                         rotation=0, fontsize=textsize,
                         horizontalalignment='center', verticalalignment='bottom')

        # # text for gauge ID
        # textbox_x = 0.5
        # textbox_y = 1.00
        # sub.text(textbox_x, textbox_y, station_id[istation], transform=sub.transAxes,
        #                      rotation=0, fontsize=textsize,
        #                      horizontalalignment='center', verticalalignment='bottom')

        # text for gauge name
        # insert line break if station_info (without spaces) is longer than ~20 charachters
        textbox_x = 0.5
        textbox_y = 1.14
        if nstations/ncol <= 8:
            tmp = gaugeinfo_group[igauge].split(':')[2].strip().replace('"','').split(' ')
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
        plt.setp(sub,xlim=[mdates.date2num(datetime.datetime(year_start,1,1,0,0)),mdates.date2num(datetime.datetime(year_end+1,1,1,0,0))])
        plt.setp(sub,ylim=[0.1,np.max(var[:,idx])*1.05])

        # rotate x-ticks
        plt.xticks(rotation=45)

        # legend
        ll = plt.legend(frameon=frameon, ncol=1,bbox_to_anchor=(llxbbox,llybbox), loc='upper left',
                    scatterpoints=1, numpoints=1,
                    labelspacing=llrspace, columnspacing=llcspace, handletextpad=llhtextpad, handlelength=llhlength)
        plt.setp(ll.get_texts(), fontsize='small')
        

    # done with this page
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


