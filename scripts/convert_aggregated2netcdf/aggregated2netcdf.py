#!/usr/bin/env python
from __future__ import print_function

# Copyright 2016-2020 Juliane Mai - juliane.mai(at)uwaterloo.ca
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
#    python aggregated2netcdf.py -i test-arcgis-raster.txt -g test-grid.nc -o test-arcgis-raster.nc -v test,int32,1,test variable -a addinfo.csv
#
#    python aggregated2netcdf.py -i ../../data/landcover_NALCMS_GreatLakes/landcover_NALCMS_GreatLakes_aggregated/GL_NALCMS_landcover_rdrs_v2.txt -g ../../data/meteo_forcing_RDRS-v2/grip-gl_rdrs-v2-gridonly.nc -o ../../data/landcover_NALCMS_GreatLakes/landcover_NALCMS_GreatLakes_aggregated.nc -v "perc_land_cover,float,percent,percentage of land cover type of NALCMS product for each forcing grid cell" -a ../../data/landcover_NALCMS_GreatLakes/NACLMS_legends.csv

# python aggregated2netcdf.py -i ../../data/landcover_NALCMS_GreatLakes/landcover_NALCMS_GreatLakes_aggregated/GL_NALCMS_landcover_wfdei_gem_capa.txt -g ../../data/meteo_forcings_WFDEI-GEM-CaPA/grip-gl_wfdei-gem-capa_gridonly.nc -o ../../data/landcover_NALCMS_GreatLakes/landcover_NALCMS_GreatLakes_aggregated.nc -v "perc_land_cover,float,percent,percentage of land cover type of NALCMS product for each forcing grid cell" -a ../../data/landcover_NALCMS_GreatLakes/NACLMS_legends.csv

# python aggregated2netcdf.py -i ../../data/soilclass_GSDE_GreatLakes/soilclass_GSDE_GreatLakes_aggregated/GL_GSDE_usda_soil_class_rdrs_v2.txt -g ../../data/meteo_forcing_RDRS-v2/grip-gl_rdrs-v2-gridonly.nc -o ../../data/soilclass_GSDE_GreatLakes/soilclass_GSDE_GreatLakes_aggregated_RDRS-v2.nc -l ../../data/soilclass_GSDE_GreatLakes/soilclass_GSDE_GreatLakes_aggregated/variable_info.csv

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

input_file   = 'test-arcgis-raster.dat'
grid_file    = 'test-grid.nc'
output_file  = 'test-arcgis-raster.nc'
varname      = 'test,int32,1,test variable'
addinfo      = None
varinfo      = None
nodata_value = -9999
parser      = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
              description='''Convert files from ArcGIS raster format into NetDF file usable in CaSPAr.''')
parser.add_argument('-i', '--input_file', action='store',
                    default=input_file, dest='input_file', metavar='input_file', nargs=1,
                    help='Name of input file (aggreagted txt file -  Hongren style).')
parser.add_argument('-g', '--grid_file', action='store',
                    default=grid_file, dest='grid_file', metavar='grid_file', nargs=1,
                    help='Name of file that contains the grid (NetCDF file).')
parser.add_argument('-o', '--output_file', action='store',
                    default=output_file, dest='output_file', metavar='output_file', nargs=1,
                    help='Name of output file (NetCDF file).')
parser.add_argument('-v', '--varname', action='store',
                    default=varname, dest='varname', metavar='varname', nargs=1,
                    help='Variable name, datatype, unit, longname.')
parser.add_argument('-a', '--addinfo', action='store',
                    default=addinfo, dest='addinfo', metavar='addinfo',
                    help='File containing additional information.')
parser.add_argument('-l', '--varinfo', action='store',
                    default=addinfo, dest='varinfo', metavar='varinfo',
                    help='File containing information about variables and units and longname.')
parser.add_argument('-n', '--nodata_value', action='store',
                    default=nodata_value, dest='nodata_value', metavar='nodata_value',
                    help='Specify NODATA value (default: -9999).')


args          = parser.parse_args()
input_file    = args.input_file[0]
grid_file     = args.grid_file[0]
output_file   = args.output_file[0]
varname       = np.array(args.varname[0].split(','))
addinfo       = args.addinfo
varinfo       = args.varinfo
nodata_value  = args.nodata_value

del parser, args

# ------------------------
# read grid information from NetCDF grid file
# ------------------------
print('   Reading NetCDF grid ...')
nc_grid = nc4.NcDataset(grid_file, "r")

grid_lat = nc_grid['lat'][:]
grid_lon = nc_grid['lon'][:]
ncols    = nc_grid.dimensions['rlon'].size
nrows    = nc_grid.dimensions['rlat'].size

grid_lon = np.where(grid_lon>180.,grid_lon-360.0,grid_lon) # make sure range is [-180,180]
nc_grid.close()

# ------------------------
# Create variables and dimensions
# ------------------------
print('   Creating NetCDF ...')
nc_out = nc4.NcDataset(output_file, "w")

# create dimensions
dim_dict         = {}
dim_dict["time"] = 1
dim_dict["rlon"] = ncols
dim_dict["rlat"] = nrows
nc_out.createDimensions(dim_dict, fail=True)

# ------------------------
# add additional information
# ------------------------
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
    
    addinfo_data       = fsread(addinfo,comment='#',separator=';',snc=nvars)
    column_headers_ref = fsread(addinfo,comment='#',separator=';',snc=[nvars],fill=True) # get column headers for txt file
    column_headers_ref = np.transpose(column_headers_ref[4:])[0]   # make 1d
    column_headers_ref = [ ivar.strip() for ivar in column_headers_ref ] # remove blanks

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

# ------------------------
# read variable information if available
# ------------------------
if not(varinfo is None):

    # determine number of variables
    f = open(varinfo, 'r')
    nvars = 0
    lines = f.readlines()
    for line in lines:
        if not(line.startswith("#")):
            nvars += 1
    #print('line  = ',line)
    print('nvars = ',nvars)

    varinfo_data       = np.array(fsread(varinfo,comment='#',separator=';',snc=4))
    column_headers_ref = fsread(varinfo,comment='#',separator=';',snc=[0],fill=True) # get column headers for txt file
    column_headers_ref = np.transpose(column_headers_ref[:])[0]   # make 1d
    column_headers_ref = [ ivar.strip() for ivar in column_headers_ref ] # remove blanks

    varnames    = np.array([ ivar.strip() for ivar in varinfo_data[:,0] ])
    vartypes    = np.array([ ivar.strip() for ivar in varinfo_data[:,1] ])
    units       = np.array([ ivar.strip() for ivar in varinfo_data[:,2] ])
    description = np.array([ ivar.strip() for ivar in varinfo_data[:,3] ])
                
# -------------------------------
# create variable
# -------------------------------
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

if not(addinfo is None):
    nc_out.createVariable(varname[0], varname[1], dimensions=("time", "rlat", "rlon", "categories"), fail=True,zlib=True)
    nc_out.variables[varname[0]].setncattr("units",     varname[2])
    nc_out.variables[varname[0]].setncattr("long_name", varname[3])
    nc_out.variables[varname[0]].setncattr("coordinates", "lon lat")
    nc_out.variables[varname[0]].setncattr("missing_value", np.array([nodata_value],dtype=varname[1])[0])
elif not(varinfo is None):
    for ivar in range(nvars):
        nc_out.createVariable(varnames[ivar], vartypes[ivar], dimensions=("time", "rlat", "rlon"), fail=True,zlib=True)
        nc_out.variables[varnames[ivar]].setncattr("units",     units[ivar])
        nc_out.variables[varnames[ivar]].setncattr("long_name", description[ivar])
        nc_out.variables[varnames[ivar]].setncattr("coordinates", "lon lat")
        nc_out.variables[varnames[ivar]].setncattr("missing_value", np.array([nodata_value],dtype=vartypes[ivar])[0])
else:  
    nc_out.createVariable(varname[0], varname[1], dimensions=("time", "rlat", "rlon"), fail=True,zlib=True)
    nc_out.variables[varname[0]].setncattr("units",     varname[2])
    nc_out.variables[varname[0]].setncattr("long_name", varname[3])
    nc_out.variables[varname[0]].setncattr("coordinates", "lon lat")
    nc_out.variables[varname[0]].setncattr("missing_value", np.array([nodata_value],dtype=varname[1])[0])

# find columns in aggregated data required
header   = fread(input_file,skip=1,header=True)
col_idx = header.index('Col')
row_idx = header.index('Row')
lon_idx = header.index('Gridlon')
lat_idx = header.index('Gridlat')
var_idx = [ header.index(ivar) for ivar in column_headers_ref ]
ncategories = len(var_idx) 

# read in aggregated data
col_data = fread(input_file,skip=1,nc=[col_idx])            # shape: (lines, 1)
row_data = fread(input_file,skip=1,nc=[row_idx])            # shape: (lines, 1)
lon_data = fread(input_file,skip=1,nc=[lon_idx])            # shape: (lines, 1)
lat_data = fread(input_file,skip=1,nc=[lat_idx])            # shape: (lines, 1)
agg_data = fread(input_file,skip=1,nc=var_idx)              # shape: (lines, categories)

# initialize variable
if not(varinfo is None):
    var = np.ones([nrows,ncols,nvars],dtype=float) * np.array([nodata_value],dtype=float)[0]
elif not(addinfo is None):
    var = np.ones([nrows,ncols,ncategories],dtype=varname[1]) * np.array([nodata_value],dtype=varname[1])[0]
else:
    var = np.ones([nrows,ncols],dtype=varname[1]) * np.array([nodata_value],dtype=varname[1])[0]

# set variable values
nlines = np.shape(col_data)[0]
for iline in range(nlines):

    if iline % 1000 == 999:
        print("    ",iline+1," lines of ",nlines," processed ...")

    # find right index: look for grid lat/lon that is closest to lon and lat values in aggregated file 
    res_lon = np.where( np.abs(grid_lon-lon_data[iline]) == np.min(np.abs(grid_lon-lon_data[iline])) )   # (array([ 94, 337]), array([612, 542]))
    res_lat = np.where( np.abs(grid_lat-lat_data[iline]) == np.min(np.abs(grid_lat-lat_data[iline])) )   # (array([337]), array([542]))
    res = list(set(res_lon[0]).intersection(res_lat[0])) + list(set(res_lon[1]).intersection(res_lat[1]))  # [337, 542]

    # print("res = ",res, "   --> ",row_data[iline],col_data[iline])
    # if (len(res) != 2) :

    #     print("iline   = ",iline)
    #     print("res_lon = ",res_lon)
    #     print("res_lat = ",res_lat)
    #     print("res     = ",res)
    #     print("NOT MATCHING! :(")
    #     stop
    
    grid_row_idx = np.int(row_data[iline][0])   # res[0]
    grid_col_idx = np.int(col_data[iline][0])   # res[1]

    if not(varinfo is None):
        for ivar in range(nvars):
            var[grid_row_idx,grid_col_idx,ivar] = agg_data[iline,ivar]
    elif not(addinfo is None):
        for icat in range(ncategories):
            var[grid_row_idx,grid_col_idx,icat] = agg_data[iline,icat]
    else:
        print("Not sure what to do")
        stop

if not(addinfo is None):
    var = var[np.newaxis,:,:,:]    # adding the time axis
elif not(varinfo is None):
    var = var[np.newaxis,:,:,:]    # adding the time axis
else:
    var = var[np.newaxis,:,:]      # adding the time axis
var = np.ma.array(var,mask=(var == nodata_value))


        

print('   Store data in NetCDF ...')
# store data in variables
nc_out.variables["time"][:]     = 0
nc_out.variables["lat"][:]      = grid_lat[:]
nc_out.variables["lon"][:]      = grid_lon[:]

if not(varinfo is None):
    for ivar in range(nvars):
        nc_out.variables[varnames[ivar]][:] = var[:,:,:,ivar]
elif not(addinfo is None):
    nc_out.variables[varname[0]][:] = var[:]
else:
    nc_out.variables[varname[0]][:] = var[:]

# set global attributes
nc_out.setncattr('product', 'aggregated raster')
nc_out.setncattr('Conventions', 'CF-1.6')
nc_out.setncattr('License', 'These data are provided by the Canadian Surface Prediction Archive CaSPar. You should have received a copy of the License agreement with the data. Otherwise you can find them under http://caspar-data.ca/doc/caspar_license.txt or email caspar-data@uwaterloo.ca.')
nc_out.setncattr('Remarks', "This data were converted from '"+input_file+"' using 'aggregated2netcdf.py'")

nc_out.close()  # close output file
print('   Wrote: ',output_file)
print('   Done.')
print('')

