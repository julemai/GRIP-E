#! /bin/bash

# Copyright 2018 Juliane Mai - juliane.mai(at)uwaterloo.ca
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

set -e
#
# Converts all raw model outputs in GRIP-E project to NetCDF files.
#
#
prog=$0
pprog=$(basename ${prog})
dprog=$(dirname ${prog})
isdir=${PWD}
pid=$$

#         temporal validation:
#         bad for all models   bad for ML   bad for GR4J
stations='02HC030              04224775     04045500'

for station in ${stations} ; do
    python plot_single_gauge.py -i ../../data/objective_1/great-lakes/calibration/netcdf/all_gauges.nc -s ${station} -p ${station}.pdf -v Q -a '2001-01-01:2017-01-01' -m '../../data/objective_1/great-lakes/calibration/model/ML-LSTM/ml-lstm_phase_1_objective_1.nc ../../data/objective_1/great-lakes/calibration/model/GR4J-Raven-lp/gr4j-raven-lp_phase_1_objective_1.nc'
    pdfcrop ${station}.pdf
    mv ${station}-crop.pdf ${station}.pdf
done


exit 0
