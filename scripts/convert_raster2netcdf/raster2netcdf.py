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
#    python raster2netcdf.py -i test-arcgis-raster.txt -o test-arcgis-raster.nc -v test,int32,1,test variable -a addinfo.csv

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

import netcdf4     as     nc4      # in lib/
from fread         import fread    # in lib/
from fsread        import fsread   # in lib/

input_file  = 'test-arcgis-raster.dat'
output_file = 'test-arcgis-raster.nc'
varname     = 'test,int32,1,test variable'
addinfo     = None
parser      = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
              description='''Convert files from ArcGIS raster format into NetDF file usable in CaSPAr.''')
parser.add_argument('-i', '--input_file', action='store',
                    default=input_file, dest='input_file', metavar='input_file', nargs=1,
                    help='Name of input file (raster file).')
parser.add_argument('-o', '--output_file', action='store',
                    default=output_file, dest='output_file', metavar='output_file', nargs=1,
                    help='Name of output file (NetCDF file).')
parser.add_argument('-v', '--varname', action='store',
                    default=varname, dest='varname', metavar='varname', nargs=1,
                    help='Variable name, datatype, unit, longname.')
parser.add_argument('-a', '--addinfo', action='store',
                    default=addinfo, dest='addinfo', metavar='addinfo',
                    help='File containing additional information.')


args          = parser.parse_args()
input_file    = args.input_file[0]
output_file   = args.output_file[0]
varname       = np.array(args.varname[0].split(','))
addinfo       = args.addinfo

del parser, args

# read raster file
header = fread(input_file,skip=6,nc=2,header=True)
#
if ( header[0][0].lower() == 'ncols'):
    ncols = int(header[0][1])
else:
    raise ValueError('Argument in header of raster data not understood: '+header[0][0])
#
if ( header[1][0].lower() == 'nrows'):
    nrows = int(header[1][1])
else:
    raise ValueError('Argument in header of raster data not understood: '+header[1][0])
#
if ( header[4][0].lower() == 'cellsize'):
    cellsize = np.float(header[4][1])
else:
    raise ValueError('Argument in header of raster data not understood: '+header[4][0])
#
if ( header[2][0].lower() == 'xllcorner'):
    xllcorner = np.float(header[2][1])
    xllcenter = xllcorner + cellsize/2.
elif ( header[2][0].lower() == 'xllcenter'):
    xllcenter = np.float(header[2][1])
else:
    raise ValueError('Argument in header of raster data not understood: '+header[2][0])
#
if ( header[3][0].lower() == 'yllcorner'):
    yllcorner = np.float(header[3][1])
    yllcenter = yllcorner + cellsize/2.    
elif ( header[3][0].lower() == 'yllcenter'):
    yllcenter = np.float(header[3][1])
else:
    raise ValueError('Argument in header of raster data not understood: '+header[3][0])
#
if ( header[5][0].lower() == 'nodata_value'):
    nodata_value = np.float(header[5][1])
else:
    raise ValueError('Argument in header of raster data not understood: '+header[5][0])

print('   Reading raster data ...')
# longitudes are x-values
lon = [ xllcenter + ii*cellsize for ii in range(ncols) ]
lon = np.tile(lon, (nrows,1))
nlon = ncols

# latitudes are y-values
lat = [ [yllcenter + ii*cellsize] for ii in range(nrows)[::-1] ] 
lat = np.tile(lat, (ncols))
nlat = nrows

# the variable itself applying mask where nodata value found
var = fread(input_file,skip=6)
var = var[np.newaxis,:,:]    # adding the time axis
var = np.ma.array(var,mask=(var == nodata_value))

print('   Creating NetCDF ...')
nc_out = nc4.NcDataset(output_file, "w")

# create dimensions
dim_dict         = {}
dim_dict["time"] = 1
dim_dict["rlon"] = ncols
dim_dict["rlat"] = nrows
nc_out.createDimensions(dim_dict, fail=True)

# create variable
nc_out.createVariable("lat",   "single", dimensions=("rlat", "rlon"),         fail=True,zlib=True)
nc_out.variables["lat"].setncattr("units",         "degrees_north")
nc_out.variables["lat"].setncattr("long_name",     "latitude")
nc_out.variables["lat"].setncattr("standard_name", "latitude")

nc_out.createVariable("lon",   "single", dimensions=("rlat", "rlon"),         fail=True,zlib=True)
nc_out.variables["lon"].setncattr("units",         "degrees_east")
nc_out.variables["lon"].setncattr("long_name",     "longitude")
nc_out.variables["lon"].setncattr("standard_name", "longitude")

nc_out.createVariable("time",  "int32",  dimensions=("time"),                 fail=True,zlib=True)
nc_out.variables["time"].setncattr("units",         "hours since 2010-01-01 00:00:00")
nc_out.variables["time"].setncattr("long_name",     "time")
nc_out.variables["time"].setncattr("standard_name", "time")
nc_out.variables["time"].setncattr("calendar",      "gregorian")
nc_out.variables["time"].setncattr("axis",          "T")

nc_out.createVariable(varname[0], varname[1], dimensions=("time", "rlat", "rlon"), fail=True,zlib=True)
nc_out.variables[varname[0]].setncattr("units",     varname[2])
nc_out.variables[varname[0]].setncattr("long_name", varname[3])
nc_out.variables[varname[0]].setncattr("coordinates", "lon lat")
nc_out.variables[varname[0]].setncattr("missing_value", np.array([nodata_value],dtype=varname[1])[0])

if not(addinfo is None):

    # determine number of additional variables
    f = open(addinfo, 'r')
    found = False
    while not(found):
        line = f.readline().strip()
        if not(line.startswith("#")):
            found = True
            f.close()
    nvars = np.shape(line.split(';'))[0] 
    #print('line  = ',line)
    #print('nvars = ',nvars)
    
    addinfo_data = fsread(addinfo,comment='#',separator=';',snc=nvars)

    varnames    = addinfo_data.pop(0)
    vartypes    = addinfo_data.pop(0)
    units       = addinfo_data.pop(0)
    description = addinfo_data.pop(0)

    ncat = len(addinfo_data)
    dim_dict               = {}
    dim_dict["categories"] = ncat
    nc_out.createDimensions(dim_dict, fail=True)


    for ii, ivarname in enumerate(varnames):

        # datatypes: ['i8', 'f4', 'f8', 'S1', 'i2', 'i4', 'u8', 'u4', 'u1', 'u2', 'i1']
        tmp_vartypes = vartypes[ii]
        if 'char' in tmp_vartypes:
            nc_out.createVariable(varnames[ii], str, dimensions=("categories"), fail=True,zlib=True)
        else:
            nc_out.createVariable(varnames[ii], tmp_vartypes, dimensions=("categories"), fail=True,zlib=True)
        nc_out.variables[varnames[ii]].setncattr("units",     units[ii])
        nc_out.variables[varnames[ii]].setncattr("long_name", description[ii])

        addinfo_data_tmp = [ addinfo_data[jj][ii].strip() for jj in range(ncat) ]
        if 'int' in vartypes[ii]:
            addinfo_data_tmp = np.array(addinfo_data_tmp,dtype=np.int)
            nc_out.variables[varnames[ii]][:] = addinfo_data_tmp[np.newaxis,:]
        elif 'float' in vartypes[ii]:
            addinfo_data_tmp = np.array(addinfo_data_tmp,dtype=np.float64)
            nc_out.variables[varnames[ii]][:] = addinfo_data_tmp[np.newaxis,:]
        elif 'single' in vartypes[ii]:
            addinfo_data_tmp = np.array(addinfo_data_tmp,dtype=np.float32)
            nc_out.variables[varnames[ii]][:] = addinfo_data_tmp[np.newaxis,:]
        else:
            addinfo_data_tmp = np.array(addinfo_data_tmp)  # stays as string
            # string variables need to be written sequentially
            for jj, val in enumerate(addinfo_data_tmp):
                nc_out.variables[varnames[ii]][jj] = val

        

print('   Store data in NetCDF ...')
# store data in variables
nc_out.variables["time"][:]     = 0
nc_out.variables["lat"][:]      = lat[:]
nc_out.variables["lon"][:]      = lon[:]
nc_out.variables[varname[0]][:] = var[:]

# set global attributes
nc_out.setncattr('product', 'raster')
nc_out.setncattr('Conventions', 'CF-1.6')
nc_out.setncattr('License', 'These data are provided by the Canadian Surface Prediction Archive CaSPar. You should have received a copy of the License agreement with the data. Otherwise you can find them under http://caspar-data.ca/doc/caspar_license.txt or email caspar-data@uwaterloo.ca.')
nc_out.setncattr('Remarks', "This data were converted from '"+input_file+"' using 'raster2netcdf.py'")

nc_out.close()  # close output file
print('   Done.')
print('')

