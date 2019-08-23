#!/usr/bin/env python
from __future__ import print_function

# Copyright 2019 Juliane Mai - juliane.mai(at)uwaterloo.ca
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

import numpy as np
import glob as glob

# files = glob.glob('/Users/j6mai/Desktop/GRIP-GL_Raw_data_20190801/GRIP-GL_WSC_streamflow_raw_Objective[1-2]_old/*.txt')
files = glob.glob('/Users/j6mai/Desktop/GRIP-GL_Raw_data_20190822/all_old_format/[0-9][0-9][A-Z][A-Z][0-9][0-9][0-9].txt')

for ifile in files:
    station = ifile.split('/')[-1].split('.')[0]

    subfolder = '_'.join(ifile.split('/')[-2].split('_')[:-1])
    subfolder = 'all_new_format'
    newfile = '/'.join(ifile.split('/')[0:-2])+'/'+subfolder+'/'+station+'.csv'

    ff = open(ifile, "r")
    lines = ff.readlines()
    ff.close()

    print('')
    print('Old file: ',ifile)
    print('New file: ',newfile)

    ff = open(newfile, "w")
    ff.write("    Daily Discharge (m3/s) (PARAM = 1) \n")
    ff.write(" ID,PARAM,Date,Value,SYM \n")

    ndays_all = 0
    for ll in lines[1:]:
        # YEAR MONTH FULL_MONTH NO_DAYS MONTHLY_MEAN MONTHLY_TOTAL FIRST_DAY_MIN MIN FIRST_DAY_MAX MAX FLOW1 FLOW_SYMBOL1 FLOW2 FLOW_SYMBOL2 ...
        # 0    1     2          3       4            5             6             7   8             9   10    11           12    13           ...
        tmp = ll
        ll=ll.replace('\r\n','').replace('"','').split('\t')
        ndays = np.int(ll[3])
        vals  = ll[10::2]
        symb  = ll[11::2]

        # tmp_vals = map(float,vals)

        if len(vals) < ndays:
            print("len(vals) = ",len(vals))
            print("ndays     = ",ndays)
            raise ValueError("Somehow not enough values found!")
        if len(symb) < ndays:
            print("len(symb) = ",len(symb))
            print("ndays     = ",ndays)
            raise ValueError("Somehow not enough quality symbols found!")

        year = np.int(ll[0])
        month = np.int(ll[1])

        for iday in range(ndays):
            # 02GG002,1,1947/08/28,1.27,B
            # 02GG002,1,1947/08/29,0.481,
            # 02GG002,1,1947/08/30,0.227,
            data_str = "{:s},1,{:04d}/{:02d}/{:02d},{:s},{:s}\n".format(station,year,month,iday+1,vals[iday],symb[iday])
            ff.write(data_str)

        ndays_all += ndays

    ff.close()

    print("    ndays = ",ndays_all)
