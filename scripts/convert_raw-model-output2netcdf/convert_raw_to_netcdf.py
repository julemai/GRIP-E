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
#    ------------
#    LBRM
#    ------------
#    python convert_raw_to_netcdf.py -m LBRM -i ../../data/objective_1/model/LBRM/lbrm_phase_0_objective_1.csv -o ../../data/objective_1/model/LBRM/lbrm_phase_0_objective_1.nc -a ../../data/objective_1/gauge_info.csv
#
#    ------------
#    VIC
#    ------------
#    python convert_raw_to_netcdf.py -m VIC -i ../../data/objective_1/model/VIC/vic_phase_0_objective_1.csv -o ../../data/objective_1/model/VIC/vic_phase_0_objective_1.nc -a ../../data/objective_1/gauge_info.csv -b ../../data/objective_1/model/VIC/subid2gauge.csv

#    ------------
#    VIC-GRU
#    ------------
#    python convert_raw_to_netcdf.py -m VIC-GRU -i ../../data/objective_1/model/VIC-GRU/vic-gru_phase_0_objective_1.csv -o ../../data/objective_1/model/VIC-GRU/vic-gru_phase_0_objective_1.nc -a ../../data/objective_1/gauge_info.csv -b ../../data/objective_1/model/VIC-GRU/subid2gauge.csv

#    ------------
#    GEM-Hydro
#    ------------
#    python convert_raw_to_netcdf.py -m GEM-Hydro -i ../../data/objective_1/model/GEM-Hydro/gem-hydro_phase_0_objective_1.csv -o ../../data/objective_1/model/GEM-Hydro/gem-hydro_phase_0_objective_1.nc -a ../../data/objective_1/gauge_info.csv

#    ------------
#    HYPE
#    ------------
#    python convert_raw_to_netcdf.py -m HYPE -i ../../data/objective_1/model/HYPE/hype_phase_0_objective_1_ -o ../../data/objective_1/model/HYPE/hype_phase_0_objective_1.nc -a ../../data/objective_1/gauge_info.csv

# -----------------------
# add subolder scripts/lib to search path
# -----------------------
import sys
import os 
dir_path = os.path.dirname(os.path.realpath(__file__))
sys.path.append(dir_path+'/../lib')

import argparse
import textwrap                # nicer formatting of help text in parser
import numpy as np             # to perform numerics
import shutil                  # file operations
import copy                    # deep copy objects, arrays etc
import datetime                # converting dates
import netCDF4 as nc           # NetCDF writing
import glob                    # for file listing

#import netcdf4     as     nc4         # in lib/
from fread         import fread        # in lib/
from fsread        import fsread       # in lib/
from writenetcdf   import writenetcdf  # in lib/

model                      = ['LBRM']
input_file                 = ['lbrm_phase_0_objective_1.csv']
output_file                = ['lbrm_phase_0_objective_1.nc']
gaugeinfo_file             = ['gauge_info.csv']
mapping_subbasinID_gaugeID = ['']
parser      = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
              description='''Convert LBRM raw streamflow model outputs into NetDF format (consistent across all models in GRIP-E).''')
parser.add_argument('-i', '--input_file', action='store',
                    default=input_file, dest='input_file', metavar='input_file', nargs=1,
                    help='Name of input file (raw model output file; for HYPE basename (filename to be expected to be basename<gaugeId>.txt).')
parser.add_argument('-o', '--output_file', action='store',
                    default=output_file, dest='output_file', metavar='output_file', nargs=1,
                    help='Name of output file (NetCDF file).')
parser.add_argument('-m', '--model', action='store',
                    default=model, dest='model', metavar='model', nargs=1,
                    help='Model (e.g., LBRM, VIC, VIC-GRU, HYPE, GEM-Hydro).')
parser.add_argument('-a', '--gaugeinfo_file', action='store',
                    default=gaugeinfo_file, dest='gaugeinfo_file', metavar='gaugeinfo_file', nargs=1,
                    help='File containing additional information about gauges (e.g., darinage area, lat/lon).')
parser.add_argument('-b', '--mapping_subbasinID_gaugeID', action='store',
                    default=mapping_subbasinID_gaugeID, dest='mapping_subbasinID_gaugeID', metavar='mapping_subbasinID_gaugeID',nargs=1,
                    help='File containing mapping of subbasin ID (col 1) to gauge ID (col 2). All other columns are ignored. One header line. Only required for VIC and VIC-GRU.')

args                       = parser.parse_args()
model                      = args.model[0]
input_file                 = args.input_file[0]
output_file                = args.output_file[0]
gaugeinfo_file             = args.gaugeinfo_file[0]
mapping_subbasinID_gaugeID = args.mapping_subbasinID_gaugeID[0]

del parser, args

if (model != 'LBRM') and (model != 'VIC') and (model != 'VIC-GRU') and (model != 'GEM-Hydro') and (model != 'HYPE'):
    raise ValueError('This model is not supported yet!')

if ((model == 'VIC-GRU') and (mapping_subbasinID_gaugeID == '')) or ((model == 'VIC') and (mapping_subbasinID_gaugeID == '')):
    raise ValueError('For VIC model CSV file containing the mapping of subbasin ID (col 1) to gauge ID (col 2) needs to be provided. All other columns in that file will be ignored. Exactly one header line needs to be provided.')

# read model output file
if (model == 'HYPE'):
    # ---------------
    # read model outputs
    # - every gauge is in a separate file
    # ---------------
    input_files    = glob.glob(input_file+"*")
    model_stations = [ ii.split(input_file)[1].split('.')[0] for ii in input_files ]
    model_data     = [ [] for ii in input_files ]
    model_dates    = None
    for ii,iinput_file in enumerate(input_files):
        
        # find column containing discharge "cout"
        head = fread(iinput_file,skip=2,cskip=1,header=True)
        idx = head[0].index('cout')
        
        model_data[ii]  = fread(iinput_file,skip=2,cskip=1,header=False)[:,idx]

        # make sure all model dates are same in all files
        if model_dates is None:
            # save first file's dates
            model_dates = fsread(iinput_file,skip=2,snc=1)
            model_dates = [ datetime.datetime( int(str(iii[0])[0:4]),int(str(iii[0])[5:7]),int(str(iii[0])[8:10]),0,0 ) for iii in model_dates ]
        else:
            # check if dates are same as already saved
            tmp_dates = fsread(iinput_file,skip=2,snc=1)
            tmp_dates = [ datetime.datetime( int(str(iii[0])[0:4]),int(str(iii[0])[5:7]),int(str(iii[0])[8:10]),0,0 ) for iii in tmp_dates ]
            if not np.all(tmp_dates == model_dates):
                print('Time steps first file: ',input_files[0])
                print('     ',model_dates)
                print('Time steps current file: ',input_files[ii])
                print('     ',tmp_dates)
                raise ValueError('Time step in files must be all the same!')
            

    model_data  = np.transpose(np.array(model_data))
    model_dates = np.transpose(np.array(model_dates))
    
if (model == 'LBRM'):
    # ---------------
    # read model outputs
    # ---------------
    model_stations = fread(input_file,skip=1,cskip=1,header=True)
    model_data     = fread(input_file,skip=1,cskip=1,header=False)
    model_data     = np.array(model_data,dtype=np.float32)
    model_dates    = fsread(input_file,skip=1,snc=1)
    model_dates    = [ datetime.datetime( int(str(ii[0])[0:4]),int(str(ii[0])[5:7]),int(str(ii[0])[8:10]),0,0 ) for ii in model_dates ]

if (model == 'GEM-Hydro'):
    # ---------------
    # read model outputs
    # - model outputs already pre-processes using "strf_graphs_scores.py"
    # ---------------
    model_stations = fread(input_file,skip=1,cskip=2,header=True)
    model_data     = fread(input_file,skip=1,cskip=2,header=False)
    model_data     = np.array(model_data,dtype=np.float32)
    model_dates    = fsread(input_file,skip=1,snc=2)
    model_dates    = [ datetime.datetime( int(str(ii[0])[0:4]),int(str(ii[0])[5:7]),int(str(ii[0])[8:10]),int(str(ii[1])[0:2]),int(str(ii[1])[3:5]) ) for ii in model_dates ]

if (model == 'VIC-GRU' or model == 'VIC'):
    # ---------------
    # read model outputs
    # - model outputs contain subbasin ID and not gauge ID --> need to be remapped
    # ---------------
    model_stations = fread(input_file,skip=1,cskip=4,header=True)   # this is subbasin IDs    
    model_data     = fread(input_file,skip=1,cskip=4,header=False)
    model_data     = np.array(model_data,dtype=np.float32)
    model_dates    = fsread(input_file,skip=1,cskip=1,snc=2)
    model_dates    = [ datetime.datetime( int(str(ii[0])[0:4]),int(str(ii[0])[5:7]),int(str(ii[0])[8:10]),int(str(ii[1])[0:2]),int(str(ii[1])[3:5]) ) for ii in model_dates ]

    # ---------------
    # read mapping info subbasin ID --> gauge station ID
    # ---------------
    mapping = fsread(mapping_subbasinID_gaugeID,skip=0,snc=2)
    mapping = np.array(mapping)
    
    # ---------------
    # map subbasin ID's to gauging stations IDs
    # ---------------
    for ii,isubbasin in enumerate(model_stations):  # they look like "sub676 [m3/s]" --> "676"

        subbasin_ID = isubbasin.split(' ')[0].split('sub')[1]
        idx = np.where(mapping[:,0]==subbasin_ID)[0][0]
        gauge_id = mapping[idx,1]

        model_stations[ii] = gauge_id

# get gauge station file
gaugeinfo_header = fsread(gaugeinfo_file,comment='#',separator=',',skip=1,header=True)
ncols = len(gaugeinfo_header)
gaugeinfo_data   = fsread(gaugeinfo_file,comment='#',separator=',',skip=1,snc=ncols,header=False)

# read gauge information from CSV file
# [NO,ID,Name,Lat,Lon,Country,Drainage_area,...]
gaugeinfo_all  = []
station_id     = model_stations                                   # station IDs
with open(gaugeinfo_file) as f:
    content = f.readlines()

for cc in content[1:]:
    ccc = ' '.join(cc.strip().split())

    # split at "," but not when in quotes
    import re
    PATTERN = re.compile(r'''((?:[^,"']|"[^"]*"|'[^']*')+)''')
    # print( PATTERN.split(ccc)[1::2] )
    ccc = PATTERN.split(ccc)[1::2]
    
    gaugeinfo_all.append(ccc[1:7])

# slim down gauge information to only gauges present (in right order)
station_id_all = [ ii[0] for ii in gaugeinfo_all ]

# ------------------------
# only dump stations to NetCDF that are in gauge_info.csv
# ------------------------
# check if data read from model are only stations required (see "station_id_all"); if not, delete data from:
#      model_stations
#      model_data
nn = len(model_stations)
ii = 0
while ii < nn:
    if not(model_stations[ii] in station_id_all):
        # print('not in: ',model_stations[ii])
        model_stations.pop(ii)
        model_data = np.delete(model_data,ii,axis=1)
        nn -= 1
    else:
        # print('in:     ',model_stations[ii])
        ii += 1

# ------------------------
# do not dump stations multiple times
# ------------------------
unique_stations = list(np.unique(model_stations))
for uu in unique_stations:

    idx = np.sort(np.where(np.array(model_stations) == uu)[0])
    for ii in idx[:0:-1]:
        # delete that column from
        #      model_stations
        #      model_data
        model_stations.pop(ii)
        model_data = np.delete(model_data,ii,axis=1)

idx_stations   = [ station_id_all.index(ss) for ss in station_id ]
gaugeinfo      = [ gaugeinfo_all[idx] for idx in idx_stations ]   # gauge information of only the requested gauges: [ID,Name,Lat,Lon,Country,Drainage_area]
station_info   = [[] for istation in range(len(gaugeinfo))]
for istation in range(len(gaugeinfo)):
    station_info[istation] = gaugeinfo[istation][0]+' : '+gaugeinfo[istation][1]+' ('+gaugeinfo[istation][4]+')'

# times
start_day            = model_dates[0]
ref_date             = 'days since '+str(start_day)
times_since_in_days = [ ((tt - start_day).total_seconds())/(60.*60.*24.) for tt in model_dates ]

# nodata
nodata = -9999.0

# ----------------------------------------------
# write NetCDF
# ----------------------------------------------

# open netcdf file and add some general information
print("Write '"+output_file+"' ...")
fh = nc.Dataset( output_file, 'w', 'NETCDF4' )
# File attributes
FiAtt   = ([['description', 'Gauging station file created from '+input_file],
            ['history'    , 'Created by Juliane Mai'],
            ['Conventions', 'CF-1.6'],
            ['featureType', 'timeSeries']])
handle  = writenetcdf(fh, fileattributes=FiAtt)

# dummy arrays for dimensions and dimension lengths
dim_name     = np.array([ 'dimmmmmmmm_'+str(ii) for ii in range(100) ])  # ['dim_0', 'dim_1', 'dim_2', ...., 'dim_99']
dims         = np.zeros(np.shape(dim_name),dtype=int)
dd = 0

# ------------------
# station IDs
# ------------------
nstations = len(station_id)

varName      = 'station_id'
varAtt       = ([['long_name', 'station ID'],
                 ['units',     '1'],
                 ['cf_role',   'timeseries_id']])
dims[dd]     = nstations
dim_name[dd] = 'nstations'

# create variable
dh = fh.createDimension(dim_name[dd], dims[dd])
vh = fh.createVariable(varName, str, tuple(['nstations']), zlib=True)
for ii in range(len(varAtt)):
    vh.setncattr(varAtt[ii][0], varAtt[ii][1])
    
# set variable values
vh[:] = np.array(station_id)

dd += 1

# ------------------
# station info
# ------------------
nstations = len(station_id)

varName      = 'station_info'
varAtt       = ([['long_name', 'station long information'],
                 ['units',     '1']])
dims[dd]     = nstations
dim_name[dd] = 'nstations'

# create variable
vh = fh.createVariable(varName, str, tuple(['nstations']), zlib=True)
for ii in range(len(varAtt)):
    vh.setncattr(varAtt[ii][0], varAtt[ii][1])
    
# set variable values
vh[:] = np.array(station_info)

dd += 1

# only if additional info about gauges is given
if (not(gaugeinfo_file is None)):
    
    # ------------------
    # latitudes
    # ------------------
    varName    = 'lat'
    attributes = {"long_name":     "latitude",
                  "standard_name": "latitude",
                  "units":         "degrees_north"}
        
    # set variable values
    # gaugeinfo = [ID,Name,Lat,Lon,Country,Drainage_area]
    lat   = np.array([ np.float(gg[2]) for gg in gaugeinfo ], dtype=np.float32)

    arr_shape = list([len(lat)])
    idx = [ np.where(dims==arr_shape[ddd])[0][0] for ddd in np.arange(len(arr_shape)) ]

    # write values to dimension
    vh   = writenetcdf(fh, name=varName, dims=dim_name[idx],
                       comp=True,
                       attributes=attributes, 
                       vartype=lat.dtype)
    writenetcdf( fh, vh, var = lat)

    # ------------------
    # longitudes
    # ------------------
    varName    = "lon"
    attributes = {"long_name":     "longitude",
                  "standard_name": "longitude",
                  "units":         "degrees_east"}
        
    # set variable values
    # gaugeinfo = [ID,Name,Lat,Lon,Country,Drainage_area]
    lon   = np.array([ np.float(gg[3]) for gg in gaugeinfo ], dtype=np.float32)

    arr_shape = list([len(lon)])
    idx = [ np.where(dims==arr_shape[ddd])[0][0] for ddd in np.arange(len(arr_shape)) ]

    # write values to dimension
    vh   = writenetcdf(fh, name=varName, dims=dim_name[idx],
                       comp=True,
                       attributes=attributes, 
                       vartype=lon.dtype)
    writenetcdf( fh, vh, var = lon)
    
    # ------------------
    # drainage area
    # ------------------
    varName    = "drainage_area"
    attributes = {"long_name":     "drainage area",
                  "units":         "km**2"}
        
    # set variable values
    # gaugeinfo = [ID,Name,Lat,Lon,Country,Drainage_area]
    area   = np.array([ np.float(gg[5]) for gg in gaugeinfo ], dtype=np.float32)

    arr_shape = list([len(area)])
    idx = [ np.where(dims==arr_shape[ddd])[0][0] for ddd in np.arange(len(arr_shape)) ]

    # write values to dimension
    vh   = writenetcdf(fh, name=varName, dims=dim_name[idx],
                       comp=True,
                       attributes=attributes, 
                       vartype=area.dtype)
    writenetcdf( fh, vh, var = area)

# ------------------
# time
# ------------------
# set attributes
varName     = 'time'
varAtt  = ([['units',         ref_date],
            ['calendar',      'gregorian'],
            ['standard_name', 'time']])
ntime        = np.shape(times_since_in_days)[0]
dims[dd]     = ntime
dim_name[dd] = 'time'

# create dimension; dims=None makes it unlimited
th = writenetcdf(fh, name=dim_name[dd], dims=None, attributes=varAtt, isdim=True)
# write values to dimension
writenetcdf(fh, th, time=list(range(ntime)), var=times_since_in_days)

dd += 1

# ------------------
# the discharge data itself
# ------------------

# set attributes
varName     = 'Q'
attributes = {'long_name': "discharge",
              'units':     'm**3 s**-1',
              '_FillValue': np.float32(nodata)}

model_data = np.array(model_data,np.float32)
model_data = np.transpose(model_data) # to match order of dimensions of observations

arr_shape = np.shape(model_data)
idx = [ np.where(dims==arr_shape[ddd])[0][0] for ddd in np.arange(len(arr_shape)) ]

# write values to dimension
vh   = writenetcdf(fh, name=varName, dims=dim_name[idx],
                       comp=True,
                       attributes=attributes, 
                       vartype=model_data.dtype)
writenetcdf( fh, vh, var = model_data)

# set global attributes
fh.setncattr('product', 'raster')
fh.setncattr('Conventions', 'CF-1.6')
fh.setncattr('License', 'These data are provided by the GWF funded IMPC project A5 - Hydrologic model inter-comparison and multi-model analysis for improved prediction. The data are under a GWF licence.')
fh.setncattr('Remarks', "This data were converted from "+input_file+" using convert_raw_to_netcdf.py")

# close netcdf
fh.close() 

