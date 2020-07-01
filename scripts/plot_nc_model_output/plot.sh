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

datapath="../data/"

plot_models='VIC-GRU'                      # [    Lake Erie:   LBRM  HMETS-Raven-lp GR4J-Raven-lp GR4J-Raven-sd HYPE HYMOD2-DS
#                                                                GEM-Hydro WRF-Hydro MESH-SVS MESH-CLASS VIC VIC-GRU
#                                                                WATFLOOD SWAT-EPA SWAT-Guelph
#                                                                ML-ConvLSTM ML-ConvLSTM-DEM ML-ConvLSTM-LC ML-ConvLSTM-LC-DEM ML-LinReg ML-LSTM ML-XGBoost
#                                                                mHM-UFZ mHM-Waterloo
#                                                                Raven-blended
#                                            #      Great Lakes: GR4J-Raven-lp GR4J-Raven-sd LBRM-MG LBRM-ML-LSTM ML-EA-LSTM ML-LSTM ML-XGBoost]

# domain='lake-erie'                                      # [lake-erie great-lakes]
# periods='2011-01-01:2014-12-31'                         # time period(s) that should be used to derive NSE etc
# calvals='calibration'                                   # [calibration validation]  # choose ONE only

domain='lake-erie'                                      # [lake-erie great-lakes]
periods='2011-01-01:2014-12-31'                         # time period(s) that should be used to derive NSE etc
calvals='validation'                                    # [calibration validation]  # choose ONE only

# domain='great-lakes'                                      # [lake-erie great-lakes]
# periods='2001-01-01:2010-12-31 2011-01-01:2016-12-31'     # time period(s) that should be used to derive NSE etc
# calvals='calibration'                                     # [calibration validation] # choose ONE only

# domain='great-lakes'                                      # [lake-erie great-lakes]
# periods='2001-01-01:2010-12-31'                           # time period(s) that should be used to derive NSE etc
# calvals='validation'                                      # [calibration validation]  # choose ONE only

plot_obj='1 2'                                            # can be 1, 2, and/or 3
plot_phase='1'                                          # phase 0: uncalibrated, different phys. setups,
#                                                         # phase 1: calibrated,   different phys. setups,
#                                                         # phase 2: calibrated,   same phys. setups

for calval in ${calvals} ; do
    
    for imodel in ${plot_models} ; do

        imodel_lower=$( echo "$imodel" | tr '[:upper:]' '[:lower:]' )
	echo ""
	echo "-------------------------"
        echo ${imodel_lower}
	echo "-------------------------"

        for iobj in ${plot_obj} ; do

            for iphase in ${plot_phase} ; do

                for period in ${periods} ; do

		    period_str=$( echo ${period//':'/'_'} )     # "2001-01-01:2010-12-31" --> "2001-01-01_2010-12-31"

                    echo ''
                    echo 'Plot :: '${imodel}'  :: Objective #'${iobj}'  :: Phase '${iphase}'   :: Period '${period}'   :: '${calval}

                    basename="../../data/objective_${iobj}/${domain}/${calval}/model/${imodel}/${imodel_lower}_phase_${iphase}_objective_${iobj}"
                    
                    python plot_nc_model_output.py -a ${period} -i ${basename}.nc -p ${basename}_${period_str}.pdf
                    pdfcrop ${basename}_${period_str}.pdf
                    pdfsplit ${basename}_${period_str}-crop.pdf
                    mv ${basename}_${period_str}-crop1.pdf ${basename}_${period_str}_hydrographs.pdf
                    mv ${basename}_${period_str}-crop2.pdf ${basename}_${period_str}_performance.pdf
                    rm ${basename}_${period_str}-crop.pdf
                    rm ${basename}_${period_str}.pdf

                done # periods

            done # plot_phase

        done # plot_obj

    done # plot_models

done # calvals

echo ''
echo 'Done.'
echo ''

exit 0
