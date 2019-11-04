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
#    python add_leap_days.py -i 2001.nc -o 2001_leap.nc 

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
from netCDF4       import num2date, date2num
import calendar
import datetime
from random import randint    

import netcdf4     as     nc4      # in lib/
from fread         import fread    # in lib/
from fsread        import fsread   # in lib/

input_file  = '2004.nc'
output_file = '2004_leap.nc'

parser      = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
              description='''In case file conatins leap years and Feb 29 is missing, it copies Feb 28 and inserts it. ''')
parser.add_argument('-i', '--input_file', action='store',
                    default=input_file, dest='input_file', metavar='input_file', nargs=1,
                    help='Name of input file (NetCDF file).')
parser.add_argument('-o', '--output_file', action='store',
                    default=output_file, dest='output_file', metavar='output_file', nargs=1,
                    help='Name of output file (NetCDF file).')


args          = parser.parse_args()
input_file    = args.input_file[0]
output_file   = args.output_file[0]

del parser, args

# first copy file
if os.path.isfile(output_file):
    print("")
    print(os.path.isfile(output_file))
    raise ValueError("Output file exists already! Please delete first!")
shutil.copy(input_file, output_file)

print('   Creating NetCDF ...')
nc_out = nc4.NcDataset(output_file, "a")

# read dates
try :
    nc_out_t_calendar = nc_out['time'].calendar
except AttributeError : 
    nc_out_t_calendar = u"gregorian" 
nc_out_t_units    = nc_out['time'].units # get unit  "days since 1950-01-01T00:00:00Z"
nc_out_t_values   = nc_out['time'][:] # get values

datevar = num2date(nc_out_t_values,units = nc_out_t_units,calendar = nc_out_t_calendar)

# contains leap year?
leap_year = False
idx_inserted = []

years   = []
months  = []
days    = []
hours   = []
minutes = []
ntime   = len(datevar)
for idate in datevar:
    years.append(  idate.year )
    months.append( idate.month )
    days.append(   idate.day )
    hours.append(  idate.hour )
    minutes.append( idate.minute )
years   = np.array(years )
months  = np.array(months)
days    = np.array(days  )
hours   = np.array(hours )
minutes = np.array(minutes )

# contains leap days
idx_leapyears = np.where([ calendar.isleap(iyear) for iyear in years ])[0]
idx_feb_29 = np.where( (months == 2) & (days == 29) )[0]
idx_feb_28 = np.where( (months == 2) & (days == 28) )[0]

contain_leap_years = (len(idx_leapyears) > 0)
contain_leap_days  = (len(idx_feb_29)    > 0)

# data contain leap years but leap days are missing --> add day and change calendar to gregorian
if contain_leap_years and not(contain_leap_days):
    print("")
    print("   Datset contains leap years but leap days are missing!")
    print("   ToDo: Copy Feb 28 and change calendar attribute!")
    print("")

    # create list of indexes to create new variables (index at Feb 29 will point to data from feb 28)
    # insert new day in date lists
    uniq_years = np.sort(np.unique(years[idx_feb_28]))[::-1]
    idx_master = range(ntime)
    years_new   = list(years )
    months_new  = list(months)
    days_new    = list(days  )
    hours_new   = list(hours )
    minutes_new = list(minutes )
    for iyear in uniq_years:
        idx_feb_28_yr = np.where( (years == iyear) & (months == 2) & (days == 28) )[0]
        ntimetoadd = len(idx_feb_28_yr)
        for ii in idx_feb_28_yr[::-1]:
            idx_master.insert(idx_feb_28_yr[0]+ntimetoadd, ii)
            years_new.insert(  idx_feb_28_yr[0]+ntimetoadd, years[ii])
            months_new.insert( idx_feb_28_yr[0]+ntimetoadd, months[ii])
            days_new.insert(   idx_feb_28_yr[0]+ntimetoadd, days[ii]+1)  # make it 29 instead of 28
            hours_new.insert(  idx_feb_28_yr[0]+ntimetoadd, hours[ii])
            minutes_new.insert(idx_feb_28_yr[0]+ntimetoadd, minutes[ii])

        # just to check ho wmuch rain we had on Feb 28 to see if it was maybe a major event
        pr1 = np.nansum(nc_out['pr'][idx_feb_28_yr,:,:])
        print('   sum of precipitation on Feb 28: ',pr1)
        print("   average sum of preci over year: ",np.nansum(nc_out['pr'][:,:,:])/365.)
        
    years_new   = np.array(years_new )
    months_new  = np.array(months_new)
    days_new    = np.array(days_new  )
    hours_new   = np.array(hours_new )
    minutes_new = np.array(minutes_new )

    ntime_new = len(idx_master)
    print("   # of timesteps increased ",ntime,"  -->  ",ntime_new)

    # increase time variable
    nc_out_t_calendar_new = 'proleptic_gregorian'
    dates_avail = [ datetime.datetime(years_new[idate], months_new[idate], days_new[idate], hours_new[idate], minutes_new[idate]) for idate in range(ntime_new) ]
    datevar_new = date2num(dates_avail, nc_out_t_units, calendar=nc_out_t_calendar_new)
    datevar_new = np.array(datevar_new, dtype=np.int)

    # set new attributes, save new time steps
    nc_out.variables["time"].setncattr("calendar", nc_out_t_calendar_new)
    nc_out['time'][:]       = datevar_new
    
    # add new time step in each variable
    for ivar in nc_out.variables.keys():

        ivar_dims = nc_out[ivar].dimensions
        if len(ivar_dims) == 3 and 'rlat' in ivar_dims and 'rlon' in ivar_dims and 'time' in ivar_dims:

            var = nc_out[ivar][:]
            print("   var: {0:20s}  --> shape: {1}".format(ivar,np.shape(var)))  # variable has already new time steps but all are filled with NaN/nodata

            idx_time = ivar_dims.index('time')
            if idx_time == 0:
                var_new = var[idx_master,:,:]
            elif idx_time == 1:
                var_new = var[:,idx_master,:]
            elif idx_time == 2:
                var_new = var[:,:,idx_master]
            else:
                raise ValueError('Time dimension not found!')

            # check if random time steps are the same...
            rnd = [randint(0, ntime_new-1) for ii in range(10)]
            for ii in range(10):
                if ii in idx_master: # this was not a leap day index
                    
                    if idx_time == 0:
                        if np.any(np.ma.array(var_new, mask=np.isnan(var_new))[ii,:,:] != np.ma.array(var, mask=np.isnan(var))[ii,:,:]):
                            raise ValueError("DATA DONT MATCH!!    t=",ii)
                    elif idx_time == 1:
                        if np.any(np.ma.array(var_new, mask=np.isnan(var_new))[:,ii,:] != np.ma.array(var, mask=np.isnan(var))[:,ii,:]):
                            raise ValueError("DATA DONT MATCH!!    t=",ii)
                    elif idx_time == 2:
                        if np.any(np.ma.array(var_new, mask=np.isnan(var_new))[:,:,ii] != np.ma.array(var, mask=np.isnan(var))[:,:,ii]):
                            raise ValueError("DATA DONT MATCH!!    t=",ii)
                    else:
                        raise ValueError('Time dimension not found!')
                    
            # set new variable
            nc_out[ivar][:]       = var_new
    
# data contains leap years and leap days exist --> nothing to do
elif contain_leap_years and contain_leap_days:
    print("")
    print("   Datset contains leap days at leap years!")
    print("   ToDo: Nothing!")
    print("")

# data contain no leap years but calendar states "365_day" or "noleap"
elif not(contain_leap_years) and not(contain_leap_days):
    if (nc_out_t_calendar == '365_day' or nc_out_t_calendar == "noleap"):
        print("")
        print("   Datset contains no leap years but calendar attribute indicate no leap years!")
        print("   Calendar: ",nc_out_t_calendar)
        print("   ToDo: Change calendar and time values accordingly!")
        print("")

        nc_out_t_calendar_new = 'proleptic_gregorian'
        dates_avail = [ datetime.datetime(years[idate], months[idate], days[idate], hours[idate], minutes[idate]) for idate in range(ntime) ]
        datevar_new = date2num(dates_avail, nc_out_t_units, calendar=nc_out_t_calendar_new)
        datevar_new = np.array(datevar_new, dtype=np.int)

        # set calendar and new time values
        nc_out.variables["time"].setncattr("calendar", nc_out_t_calendar_new)
        nc_out['time'][:] = datevar_new
        
    else:
        print("")
        print("   Datset contains no leap years and calendar attribute seems to be fine!")
        print("   Calendar: ",nc_out_t_calendar)
        print("   ToDo: Nothing!")
        print("")

# could be that it is calendar where every year contains leap days (or something else...)
else:
    print("")
    print("   Datset contains leap days:  ",contain_leap_days)
    print("   Datset contains leap years: ",contain_leap_years)
    print("   ToDo: Don't know!")
    raise ValueError("This is something weird!")
    print("")

nc_out.close()  # close output file

