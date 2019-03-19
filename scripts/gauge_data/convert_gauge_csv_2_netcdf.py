from __future__ import print_function

# This function is to convert individual gauging station csv files into a single NetCDF

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
# run convert_gauge_csv_2_netcdf.py --filetype 'WSC' --input_files "../../data/objective_1/csv/02GA010.csv" --output_file "test.nc"

# -----------------------
# add subolder scripts/lib to search path
# -----------------------
import sys
import os 
dir_path = os.path.dirname(os.path.realpath(__file__))
sys.path.append(dir_path+'/../lib')

import argparse
import numpy        as np      # to perform numerics
import datetime
import netCDF4 as nc

#import netcdf4     as nc                                         # in lib/
from   writenetcdf import writenetcdf                             # in lib/
from   fsread      import fsread                                  # in lib/
from   date2dec    import date2dec                                # in lib/
from   dec2date    import dec2date                                # in lib/

verbose          = False
input_files      = None
output_file      = 'output.nc'
filetype         = None
gaugeinfo_file   = None
parser  = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                  description='''Reads gauging data from CSV file and converts them into NetCDF. If multiple input files are given, all are merged into one single NetCDF.''')
parser.add_argument('-i', '--input_files', action='store',
                    default=input_files, dest='input_files', metavar='input_files', nargs=1,
                    help='Name of input CSV file. Multiple can be given; separated with blank in quotes. E.g. "file1 file2 file3"')
parser.add_argument('-o', '--output_file', action='store',
                    default=output_file, dest='output_file', metavar='output_file', nargs=1,
                    help='Name of output NetCDF file. Default: output.nc')
parser.add_argument('-g', '--gaugeinfo_file', action='store',
                    default=gaugeinfo_file, dest='gaugeinfo_file', metavar='gaugeinfo_file', nargs=1,
                    help='Name of file containing gauge infos of at least [NO,ID,Name,Lat,Lon,Country,Drainage_area]. If not given most station information will not be written. Default: None')
parser.add_argument('-f', '--filetype', action='store',
                    default=filetype, dest='filetype', metavar='filetype', nargs=1,
                    help='File type. Either USGS (default) or WSC.')

args           = parser.parse_args()
input_files    = args.input_files[0].split(" ")
output_file    = args.output_file[0]
filetype       = args.filetype[0].split(" ")

try:
    gaugeinfo_file = args.gaugeinfo_file[0]
except:
    gaugeinfo_file = None


if (input_files is None):

    raise ValueError("convert_gauge_csv_2_netcdf.py: input file (option -i) needs to be given")
    stop

if filetype is None:
    filetype = [ 'USGS' for ifile in input_files ]

if len(input_files) != len(filetype):
    print("")
    print("# of input files:  ",len(input_files))
    print("# of file types:   ",len(filetype))
    raise ValueError("Number of input files must match number of file types")

# read all data
content_slim    = [ [] for ifile in input_files ]      # only content that contains the station ID
station_id      = []                                   # station IDs
times_all_files = [ [] for ifile in input_files ]      # time points of all files
for ifile,input_file in enumerate(input_files):

    # station_id is expected to be the filename
    station_id.append( input_file.split('.')[-2].split('/')[-1] )

    # read data from CSV file
    with open(input_file) as f:
        content = f.readlines()

    # merge multiple blanks and make separator " "
    for cc in content:
        if station_id[ifile] in cc:
            ccc = ' '.join(cc.strip().split())
            if filetype[ifile] == 'USGS':
                content_slim[ifile].append(ccc)
            elif filetype[ifile] == 'WSC':
                if ccc.split(',')[1] == '1':
                    content_slim[ifile].append(' '.join(ccc.split(',')))

            else:
                raise ValueError("This filetype is not known!")

    # find time period covering all data
    if filetype[ifile] == 'USGS':
        times_all_files[ifile] = [ ii.split(' ')[2] for ii in content_slim[ifile][2:] ]
    elif filetype[ifile] == 'WSC':
        times_all_files[ifile] = [ ii.split(' ')[2] for ii in content_slim[ifile] ]
    else:
        raise ValueError("This filetype[ifile] is not known!")
    
    times_all_files[ifile] = [ datetime.datetime( int(str(ii)[0:4]),int(str(ii)[5:7]),int(str(ii)[8:10]),0,0 ) for ii in times_all_files[ifile] ]

times_all_files_unique = np.unique(np.array([item for sublist in times_all_files for item in sublist]))
deltat                 = np.min([ (times_all_files_unique[ii+1]-times_all_files_unique[ii]).total_seconds() for ii in range(len(times_all_files_unique)-1) ])
start_day              = times_all_files_unique[0]
end_day                = times_all_files_unique[-1]

ntime = int(np.ceil( (end_day-start_day).total_seconds()/deltat + 1))
times = [ start_day + ii*datetime.timedelta(seconds=deltat) for ii in range(int(np.ceil((end_day-start_day).total_seconds()/deltat+1))) ]

# NODATA value
nodata = -9999.

# all possible quality codes
dict_quality_codes={}
dict_quality_codes["USGS:A"] = "Approved for publication -- Processing and review completed."
dict_quality_codes["USGS:P"] = "Provisional data subject to revision."
dict_quality_codes["USGS:e"] = "Value has been estimated."
dict_quality_codes["USGS:M"] = "Value is missing."
dict_quality_codes["WSC:A"]  = "partial day"
dict_quality_codes["WSC:B"]  = "ice conditions"
dict_quality_codes["WSC:D"]  = "dry"
dict_quality_codes["WSC:E"]  = "estimated"
dict_quality_codes["WSC:S"]  = "sample(s) collected this day"
dict_quality_codes["WSC:M"]  = "Value is missing."
nquality_codes = len(dict_quality_codes)

# how to handle missing values in data array
# [agency, station_id, date, data point, quality]
def fill_missing(string,filetype):

    if filetype == 'USGS':
        ss = cc.split(" ")
        if len(ss) == 3:
            ss+=[str(nodata),"M"]
            
        try:
            np.float(ss[3])
        except:
            ss[3] = str(nodata)
            ss[4] = "M"

    elif filetype == 'WSC':
        ss = cc.split(" ")
        if len(ss) == 3:
            ss+=[str(nodata),"M"]
            
        try:
            np.float(ss[3])
        except:
            ss[3] = str(nodata)
            ss[4] = "M"
            
    else:
        raise ValueError("This filetype is not known!")
        
    return ss

# read gauge information from CSV file
# [NO,ID,Name,Lat,Lon,Country,Drainage_area,...]
gaugeinfo_all = []
if (not(gaugeinfo_file is None)):
    
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
    idx_stations   = [ station_id_all.index(ss) for ss in station_id ]
    gaugeinfo = [ gaugeinfo_all[idx] for idx in idx_stations ]   # gauge information of only the requested gauges: [ID,Name,Lat,Lon,Country,Drainage_area]
    

# gather all data we need for writing to NetCDF
station_info = [[] for ifile in input_files]
data         = [[] for ifile in input_files]
quality_code = [[] for ifile in input_files]
time_idx     = [[] for ifile in input_files]
for ifile,input_file in enumerate(input_files):

    if (gaugeinfo_file is None):
        if filetype[ifile] == 'USGS':
            station_info[ifile] = ('#'.join(content_slim[ifile][0].split('#')[1:]).strip())
            print('Data for station: '+station_info[ifile])
        elif filetype[ifile] == 'WSC':
            station_info[ifile] = station_id[ifile]   
            print('Data for station: '+station_info[ifile])
        else:
            raise ValueError("This filetype[ifile] is not known!")
    else:
        # gaugeinfo = [ID,Name,Lat,Lon,Country,Drainage_area]
        station_info[ifile] = filetype[ifile]+' : '+gaugeinfo[ifile][0]+' : '+gaugeinfo[ifile][1]+' ('+gaugeinfo[ifile][4]+')'

    if filetype[ifile] == 'USGS':
        # [agency, station_id, date, data point, quality]
        # skip first two lines
        data[ifile] = np.array([ fill_missing(cc,filetype[ifile]) for cc in content_slim[ifile][2:] ])
    elif filetype[ifile] == 'WSC':
        # [ station_id, PARAM (1 or 2), Date, Value, SYM ]
        #     PARAM = 1: Daily Discharge (m3/s)
        #     PARAM = 2: Daily Water Level (m)
        data[ifile] = np.array([ fill_missing(cc,filetype[ifile]) for cc in content_slim[ifile] ])
    else:
        raise ValueError("This filetype is not known!")

    # check if all quality codes are known and setup matrix of quality codes per data point
    quality_code[ifile] = np.array( [ [ 0 for ii in range(nquality_codes) ] for jj in range(np.shape(data[ifile])[0]) ], dtype=np.int32 )
    for idata,dd in enumerate(data[ifile][:,4]):
        if dd != '':
            if (not( all([ filetype[ifile]+':'+ddd in dict_quality_codes.keys() for ddd in dd.split(':') ]) )):
                print('Quality code "'+str(dd.split(':'))+'" found in CSV file not known. Should be added to dictionary "quality_codes".')
                
            idx = [ np.where(filetype[ifile]+':'+ddd==np.array(dict_quality_codes.keys()))[0][0] for ddd in dd.split(':') ]
            
            quality_code[ifile][idata,idx] = 1

    # times of current file
    itimes  = [ datetime.datetime( int(str(tt)[0:4]),
                                  int(str(tt)[5:7]),
                                  int(str(tt)[8:10]),
                                  0,0 ) for tt in data[ifile][:,2] ]
    # indexes of this time series in overall time step array (times)
    time_idx[ifile] = [ int((ii-start_day).total_seconds()/deltat) for ii in itimes ]

ref_date             = 'minutes since '+str(start_day)
times_since_in_hours = [ ((tt - start_day).total_seconds())/(60.) for tt in times ]

# expand all time dependent variables to overall number of time steps
qc_all_ts   = [[] for ifile in input_files]
qobs_all_ts = [[] for ifile in input_files]
for ifile,input_file in enumerate(input_files):                

    # ------------------
    # data-value qualification codes
    # ------------------
    qc                                = quality_code[ifile]
    qc_all_ts[ifile]                  = np.array( [ [ np.int(nodata) for ii in range(nquality_codes) ] for jj in range(ntime) ], dtype=np.int32 )
    qc_all_ts[ifile][time_idx[ifile]] = qc

    # ------------------
    # the discharge data itself
    # ------------------
    qobs                                = np.array( [np.float(idata) for idata in data[ifile][:,3]] )
    qobs_all_ts[ifile]                  = np.array([nodata for ii in range(ntime)], dtype=np.float32)
    qobs_all_ts[ifile][time_idx[ifile]] = qobs

    if filetype[ifile] == 'USGS':
        # convert data from ft^3/s to [m^3/s]
        idx                     = np.where(qobs_all_ts[ifile] != nodata)
        qobs_all_ts[ifile][idx] = qobs_all_ts[ifile][idx] * 0.3048**3
    
qc_all_ts   = np.array(qc_all_ts)
qobs_all_ts = np.array(qobs_all_ts)
    

# ----------------------------------------------
# write NetCDF
# ----------------------------------------------

# open netcdf file and add some general information
print("Write '"+output_file+"' ...")
fh = nc.Dataset( output_file, 'w', 'NETCDF4' )
# File attributes
FiAtt   = ([['description', 'Gauging station file created from '+', '.join(input_files)+' (filetype: '+', '.join(filetype)+')'],
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
ntime        = np.shape(times_since_in_hours)[0]
dims[dd]     = ntime
dim_name[dd] = 'time'

# create dimension; dims=None makes it unlimited
th = writenetcdf(fh, name=dim_name[dd], dims=None, attributes=varAtt, isdim=True)
# write values to dimension
writenetcdf(fh, th, time=list(range(ntime)), var=times_since_in_hours)

dd += 1


# ------------------
# dimension for quality codes
# ------------------
dims[dd]     = nquality_codes
dim_name[dd] = 'nquality_codes'

# create dimension; dims=None makes it unlimited
dh = fh.createDimension(dim_name[dd], dims[dd])
dd += 1


# ------------------
# qualification code short info, e.g. 'A', 'e', etc.
# ------------------
# set attributes
varName     = 'quality_code_info_short'
attributes = {'long_name': "short ID for qualification code, e.g. 'A', 'e', etc.",
              'units':     '1'}

arr_shape = list([len(dict_quality_codes)])
idx = [ np.where(dims==arr_shape[ddd])[0][0] for ddd in np.arange(len(arr_shape)) ]

# write values to dimension
vh   = writenetcdf(fh, name=varName, dims=dim_name[idx],
                       comp=True,
                       attributes=attributes, 
                       vartype=np.array(dict_quality_codes.keys()).dtype)
writenetcdf( fh, vh, var = np.array(dict_quality_codes.keys()))

# ------------------
# qualification code long info, e.g. 'Value has been estimated', 'missing', etc.
# ------------------
# set attributes
varName     = 'quality_code_info_long'
attributes = {'long_name': "decription of short ID for qualification code, e.g. 'e' = Value has been estimated.",
              'units':     '1'}

arr_shape = list([len(dict_quality_codes)])
idx = [ np.where(dims==arr_shape[ddd])[0][0] for ddd in np.arange(len(arr_shape)) ]

# write values to dimension
vh   = writenetcdf(fh, name=varName, dims=dim_name[idx],
                       comp=True,
                       attributes=attributes, 
                       vartype=np.array(dict_quality_codes.values()).dtype)
writenetcdf( fh, vh, var = np.array(dict_quality_codes.values()))

# ------------------
# data-value qualification codes
# ------------------

# set attributes
varName     = 'quality_code'
attributes = {'long_name': 'boolean indicating quality flags (see quality_code_info for quality code description)',
              'units':     '1',
              '_FillValue': np.int32(nodata)}

arr_shape = np.shape(qc_all_ts)
idx = [ np.where(dims==arr_shape[ddd])[0][0] for ddd in np.arange(len(arr_shape)) ]

# write values to dimension
vh   = writenetcdf(fh, name=varName, dims=dim_name[idx],
                       comp=True,
                       attributes=attributes, 
                       vartype=qc_all_ts.dtype)
writenetcdf( fh, vh, var = qc_all_ts)

# ------------------
# the discharge data itself
# ------------------

# set attributes
varName     = 'Q'
attributes = {'long_name': "discharge",
              'units':     'm**3 s**-1',
              '_FillValue': np.float32(nodata)}

arr_shape = np.shape(qobs_all_ts)
idx = [ np.where(dims==arr_shape[ddd])[0][0] for ddd in np.arange(len(arr_shape)) ]

# write values to dimension
vh   = writenetcdf(fh, name=varName, dims=dim_name[idx],
                       comp=True,
                       attributes=attributes, 
                       vartype=qobs_all_ts.dtype)
writenetcdf( fh, vh, var = qobs_all_ts)



# close netcdf
fh.close() 

