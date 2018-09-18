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

# -----------------------
# add subolder scripts/lib to search path
# -----------------------
import sys
import os 
dir_path = os.path.dirname(os.path.realpath(__file__))
sys.path.append(dir_path+'/lib')

import numpy as np
import datetime
import calendar
import copy                    # deep copy objects, arrays etc

from fread         import fread    # in lib/
from autostring    import astr          # in lib/

def convert_mm_to_m3s(time, dat, absErrorBound, lakeSurface):

    # time:          1D array of datetime.date() objects
    # dat:           1D array data in [mm]
    # absErrorBound: upper limit of returned data [m^3/s]
    # lakeSurface:   lake surface area in [m^2]

    retSet        = None
    secondsInADay = 24*60*60
    # calendar.monthrange(year, month):
    #    Returns weekday of first day of the month and
    #    number of days in month, for the specified year and month.
    dayCount = [ calendar.monthrange(time[ii].year, time[ii].month)[1] for ii in range(np.shape(dat)[0]) ]

    # magic happens here
    retSet = dat/secondsInADay/dayCount*lakeSurface/1000.

    # check of retSet values
    retSet[ np.abs(retSet) > absErrorBound ] = np.nan
    
    return(retSet)

tt  = [ [ [yy,mm] for mm in np.arange(1,13)] for yy in np.arange(1950,2016)]   # construct list of lists
ttt = np.array([item for sub in tt for item in sub])                           # flatten

# time array
times = np.array([ datetime.date(*[int(ttt[ii,0]),int(ttt[ii,1]),1]) for ii in range(np.shape(ttt)[0]) ])

# -----------------------------------------------------------
# NBS Lake Erie - MCMC samples
# -----------------------------------------------------------

# data in [mm]
#   ncols = number of months    (=792)
#   nrows = number of ensembles (=3000)
data_mm = fread("../data/monthly_posterior_distributions_Smith_Gronewold/nbs_lake_erie_samples_mm.csv")

# convert into [m^3/s]
data_m3s = copy.deepcopy(data_mm)
for ii in np.arange(np.shape(data_m3s)[0]):  # for each ensemble member
    data_m3s[ii,:] = convert_mm_to_m3s(times, data_m3s[ii,:], 10000, 27314000000)

print("")
print("Derived from samples:")
for ii in range(3):
    print('    [median,p2.5,p97.5] '+str(ii+1)+'/1950: [',np.median(data_m3s[:,ii]),',',np.percentile(data_m3s[:,ii],2.5),',',np.percentile(data_m3s[:,ii],97.5),']')

# write to file with converted units
f = open('../data/monthly_posterior_distributions_Smith_Gronewold/nbs_lake_erie_samples_m3s.csv', 'w')
for ii in range(np.shape(data_m3s)[0]):
    f.write( ','.join(astr(data_m3s[ii,:],prec=8))+'\n')
f.close()

# -----------------------------------------------------------
# NBS Lake Erie - summary statistics (derived from MCMC samples)
# -----------------------------------------------------------

# data in [mm]
#   nrows = number of months    (=792)
#   ncols = 3                   (median, 2.5 percentile, 97.5 percentile)
data_mm = fread("../data/monthly_posterior_distributions_Smith_Gronewold/nbs_lake_erie_summary-stats_mm.csv",skip=1,separator=',')[:,2:]  # only Median,2.5 Percentile,97.5 Percentile

# convert into [m^3/s]
data_m3s = copy.deepcopy(data_mm)
for ii in np.arange(np.shape(data_m3s)[1]):  # for column
    data_m3s[:,ii] = convert_mm_to_m3s(times, data_m3s[:,ii], 10000, 27314000000)

print("")
print("Read from summary stats:")
for ii in range(3):
    print('    [median,p2.5,p97.5] '+str(ii+1)+'/1950: [',data_m3s[ii,0],',',data_m3s[ii,1],',',data_m3s[ii,2],']')

# write to file with converted units
f = open('../data/monthly_posterior_distributions_Smith_Gronewold/nbs_lake_erie_summary-stats_m3s.csv', 'w')
f.write('Year,Month,Median,2.5 Percentile,97.5 Percentile\n')
for ii in range(np.shape(data_m3s)[0]):
    f.write(    astr(times[ii].year*1.0, prec=4)+','+
                astr(times[ii].month*1.0,prec=4)+','+
                astr(data_m3s[ii,0],     prec=4)+','+
                astr(data_m3s[ii,1],     prec=4)+','+
                astr(data_m3s[ii,2],     prec=4)+'\n')
f.close()



