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

plot_models='GR4J-Raven-lp GR4J-Raven-sd' 		  # [ LBRM,  GR4J-Raven-lp GR4J-Raven-sd, HYPE, GEM-Hydro, WRF-Hydro, MESH-SVS, MESH-CLASS, VIC, VIC-GRU,
#                                         		      WATFLOOD, SWAT, ML-ConvLSTM, ML-ConvLSTM-DEM, ML-ConvLSTM-LC, ML-ConvLSTM-LC-DEM, ML-LinReg, ML-XGBoost]

# domain='lake-erie,'                     		  # [lake-erie, great-lakes]
# periods='2011-01-01:2014-12-31'         		  # time period(s) that should be used to derive NSE etc
domain='great-lakes'                      		  # [lake-erie, great-lakes]
periods='2001-01-01:2010-12-31 2011-01-01:2016-12-31'     # time period(s) that should be used to derive NSE etc

calval='calibration'                      		  # [calibration, validation]  # only for Great Lakes # choose ONE only
plot_obj='1 2'                            		  # can be 1, 2, and/or 3
plot_phase='1'                            		  # phase 0: uncalibrated, different phys. setups,
#                                         		  # phase 1: calibrated,   different phys. setups,
#                                         		  # phase 2: calibrated,   same phys. setups

for imodel in ${plot_models} ; do

    imodel_lower=$( echo "$imodel" | tr '[:upper:]' '[:lower:]' )
    echo ${imodel_lower}

    for iobj in ${plot_obj} ; do

	for iphase in ${plot_phase} ; do

	    for period in ${periods} ; do

		echo ''
		echo 'Plot :: '${imodel}'  :: Objective #'${iobj}'  :: Phase '${iphase}'   :: Period '${period}

		if [[ ( ${domain} == 'lake-erie' ) ]] ; then 
		    basename="../../data/objective_${iobj}/${domain}/model/${imodel}/${imodel_lower}_phase_${iphase}_objective_${iobj}"
		else
		    basename="../../data/objective_${iobj}/${domain}/${calval}/model/${imodel}/${imodel_lower}_phase_${iphase}_objective_${iobj}"
		fi
		python plot_nc_model_output.py -a ${period} -i ${basename}.nc -p ${basename}_${period}.pdf
		pdfcrop ${basename}_${period}.pdf
		pdfsplit ${basename}_${period}-crop.pdf
		mv ${basename}_${period}-crop1.pdf ${basename}_${period}_hydrographs.pdf
		mv ${basename}_${period}-crop2.pdf ${basename}_${period}_performance.pdf
		rm ${basename}_${period}-crop.pdf
		rm ${basename}_${period}.pdf

	    done

	done

    done

done

echo ''
echo 'Done.'
echo ''

exit 0
