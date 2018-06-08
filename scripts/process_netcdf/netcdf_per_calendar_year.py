#!/usr/bin/env python
from __future__ import print_function

# This function is to convert a given NetCDF file (containing multiple years of data) into
# several files thatonly contain one calendar year per file. This is for example needed for
# VIC5.0 NetCDF version.

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

# -----------------------
# add subolder scripts/lib to search path
# -----------------------
import sys
import os 
dir_path = os.path.dirname(os.path.realpath(__file__))
sys.path.append(dir_path+'/../lib')

import argparse
import numpy        as np      # to perform numerics

import netcdf4  as nc4                                         # in lib/
from   date2dec import date2dec                                # in lib/
from   dec2date import dec2date                                # in lib/

verbose          = False
input_file       = ''
output_file      = ''
year             = None
parser  = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                  description='''Extracts last time step from input file and writes it to output file. The variables and values are unchanged in the output file.''')
parser.add_argument('-i', '--input_file', action='store',
                    default=input_file, dest='input_file', metavar='input_file', nargs=1,
                    help='Name of input NetCDF file containing several years of data.')
parser.add_argument('-o', '--output_file', action='store',
                    default=output_file, dest='output_file', metavar='output_file', nargs=1,
                    help='Name of output file pattern. Files created will be <output_file>.YYYY.nc')
parser.add_argument('-y', '--year', action='store',
                    default=year, dest='year', metavar='year',
                    help='The year that needs to be extracted. If not given, all years will be extracted. ')

# example:
#     create all years: 
#         python netcdf_per_calendar_year.py -i ../../data/meteo_forcing_RDRS.nc -o ../../data/meteo_forcing_RDRS_VIC.nc
#         --> files created: [ meteo_forcing_RDRS_VIC.2010.nc, meteo_forcing_RDRS_VIC.2011.nc .... ]
#
#     create only 2012:
#         python netcdf_per_calendar_year.py -i ../../data/meteo_forcing_RDRS.nc -o ../../data/meteo_forcing_RDRS_VIC.nc -y 2012
#         --> file created: meteo_forcing_RDRS_VIC.2012

args             = parser.parse_args()
input_file       = args.input_file[0]
output_file      = args.output_file[0]
if not(args.year is None):
    year         = np.int(args.year)

# -----------------------------------------------
# Beginning of script
# -----------------------------------------------

ncin  = nc4.NcDataset(input_file, "r")

# get time
time_dim_name = 'time'
time_var_name = 'time'
ntime = len(ncin.dimensions[time_dim_name])      # number of time steps in original file

time_unit      = ncin.variables[time_var_name].units                            # e.g. u'hours since 2010-01-01 00:00:00'
time_in_file   = ncin.variables[time_var_name][:]                               # e.g. [0.0, 1.0, 2.0, ....]
ref_date       = ' '.join(time_unit.split(' ')[2:])                             # e.g. u'2010-01-01 00:00:00'
ref_date_jul   = date2dec(eng=ref_date,calendar='proleptic_gregorian')          # e.g. 733772.0

if time_unit.split(' ')[0] == 'hours':
    first_time_jul = ref_date_jul + time_in_file[ 0]/24.
    last_time_jul  = ref_date_jul + time_in_file[-1]/24.

    time_in_file_jul = np.array([ ref_date_jul + time_in_file[ii]/24. for ii in range(ntime) ])
elif time_unit.split(' ')[0] == 'days':
    first_time_jul = ref_date_jul + time_in_file[ 0]/1.
    last_time_jul  = ref_date_jul + time_in_file[-1]/1.

    time_in_file_jul = np.array([ ref_date_jul + time_in_file[ii]/1. for ii in range(ntime) ])
else:
    raise ValueError('This time unit is not implemented yet!')
    stop

time_in_file_cal = dec2date(time_in_file_jul,calendar='proleptic_gregorian')

if not(args.year is None):
    years = np.array([year])
else:
    years = np.unique(time_in_file_cal[0])   # list of all years contained in file

for iyear in years:
    # creates new filename including year 
    filename = '.'.join(output_file.split('.')[0:-1])+'.'+str(iyear)+'.nc'
    print('Producing file '+filename+' ...')

    # all the indexes for current year
    time_idx = np.where(time_in_file_cal[0]==iyear)[0]

    ncout = nc4.NcDataset(filename, "w", )

    # create dimensions in NetCDF file
    dims = ncin.dimensions.keys()
    for idim_name in dims:
        if (idim_name == time_dim_name):
            if (ncin.dimensions[idim_name].isunlimited()):
                # create unlimited time dimension
                ncout.createDimension(idim_name, None)   # 'None' means this is the unlimited dimension
            else:
                # create limited time dimension (size will be reduced)
                ncout.createDimension(idim_name, len(time_idx))
        else:
            if (ncin.dimensions[idim_name].isunlimited()):
                # create unlimited non-time dimension
                ncout.createDimension(idim_name, None)   # 'None' means this is the unlimited dimension
            else:
                # create limited non-time dimension (size stays same as before)
                ncout.createDimension(idim_name, len(ncin.dimensions[idim_name]))

    # create variable in NetCDF file
    variables = ncin.variables.keys()
    for ivar_name in variables:

        var_in   = ncin.variables[ivar_name]
        var_type = var_in.dtype
        var_dims = var_in.dimensions
        var      = ncout.createVariable(ivar_name, var_type, var_dims, zlib=True)  # zlib makes it compressed
      
        if (ivar_name == time_var_name):
            # create time variable
            var[:] = time_in_file[time_idx]
        else:
            # create non-time variable
            pos_t    = np.where(np.array(var_dims) == time_dim_name)[0]
            if (len(pos_t) == 1):
              
                # there is exactly one dimension which is time
                pos_t = pos_t[0]
              
                # check if time is last dimension
                if (pos_t != 0):
                    print("Dimensions:    ",var_dims)
                    print("Position time: ",pos_t)
                    raise ValueError("Time needs to be the first dimension!")            
              
                # copy desired time steps
                if len(var_dims) == 1:
                    var[:] = var_in[:][time_idx]
                elif len(var_dims) == 2:
                    vari = var_in[:][time_idx,:]
                    if (len(time_idx) == 1):
                        var[:] = vari[np.newaxis,:]
                    else:
                        var[:] = vari[:,:]
                elif len(var_dims) == 3:
                    vari = var_in[:][time_idx,:,:]
                    if (len(time_idx) == 1):
                        var[:] = vari[np.newaxis,:,:]
                    else:
                        var[:] = vari[:,:,:]
                else:
                    raise ValueError("This function is only implemented for variables up to 3 dimensions.")
            elif (len(pos_t) == 0):
                if np.shape(var_in) != ():  # exclude attribute-only variables like "rotated_pole":
                    # no time dimension for this variable found --> copy all
                    var[:] = var_in[:]
            else:
                raise ValueError("This is weird. Time is found in multiple dimensions.")

        # add all attributes
        attr = var_in.attributes.keys()
        attributes = {}
        for iattr in attr:
            attributes[iattr] = var_in.attributes[iattr]

        var.copyAttributes(attributes, skip=None)

    ncout.close()

ncin.close()

print('Done.')




