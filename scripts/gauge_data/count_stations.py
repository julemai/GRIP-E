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
import jams

for objective in range(1,3):

    print("")
    print("----------------------------")
    print("  OBJECTIVE ",objective)
    print("----------------------------")

    # ifile = '/Users/j6mai/Desktop/GRIP-GL_Raw_data_20190801/GRIP-GL_all_gauges_info_Objective'+str(objective)+'_20190801.csv'
    # ifile = '/Users/j6mai/Desktop/GRIP-GL_Raw_data_20190812/GRIP_GL_gauge_info_20190812_objective_'+str(objective)+'.csv'
    # ifile = '/Users/j6mai/Desktop/GRIP-GL_Raw_data_20190822/GRIP_GL_gauge_info_20190822_objective_'+str(objective)+'.csv'
    ifile = '../../data/objective_'+str(objective)+'/great-lakes/gauge_info.csv'

    header = jams.fsread(ifile,separator=',',skip=1,header=True)
    ncols = len(header)

    data = np.array(jams.fsread(ifile,separator=',',skip=1,nc=0,snc=range(ncols)))

    idx_lake = header.index('Lake_info')
    idx_cal  = header.index('Calibration/Validation')
    idx_ID   = header.index('ID')
    
    lakes = np.unique(data[:,idx_lake])

    print('total # calibration stations: ',np.sum(data[:,idx_cal] == 'C'))
    print('total # validation  stations: ',np.sum(data[:,idx_cal] == 'V'))

    # print("lakes/regions: ",lakes)

    for ilake in lakes:
        n_all = np.sum(data[:,idx_lake] == ilake)
        n_cal = len(np.where( (data[:,idx_lake] == ilake) & (data[:,idx_cal] == 'C'))[0])
        n_val = len(np.where( (data[:,idx_lake] == ilake) & (data[:,idx_cal] == 'V'))[0])
        print('# stations in watershed {:>10}: {:3d} --> cal: {:3d}, val: {:3d}'.format(ilake,n_all,n_cal, n_val))

    print('')
    print("cal: ",' '.join(data[np.where( (data[:,idx_cal] == 'C'))[0],idx_ID]))
    print('')
    print("val: ",' '.join(data[np.where( (data[:,idx_cal] == 'V'))[0],idx_ID]))

    # idx_perc = header.index('Observation_percent_2000-2017')
    # perc = list(data[:,idx_perc])
    # perc = map(float,perc)
    # perc = np.array(perc)

    # print("percentage available data: ",perc)







