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

convert_models='LBRM VIC-GRU GEM-Hydro'  # can be [LBRM, HYPE, GEM-Hydro, WRF-Hydro, MESH, VIC, VIC-GRU, WATFLOOD]
convert_obj='1 2'      # can be 1, 2, and/or 3
convert_phase='0'      # phase 0: uncalibrated, different phys. setups,
#                      # phase 1: calibrated,   different phys. setups,
#                      # phase 2: calibrated,   same phys. setups

for imodel in ${convert_models} ; do

    imodel_lower=$( echo "$imodel" | tr '[:upper:]' '[:lower:]' )
    echo ${imodel_lower}

    if [ ${imodel} == 'VIC-GRU' ] ; then
	add_inputs="-b ../../data/objective_${iobj}/model/${imodel}/subid2gauge.csv"
    else
	add_inputs=''
    fi

    for iobj in ${convert_obj} ; do

	for iphase in ${convert_phase} ; do

	    echo ''
	    echo 'Convert :: '${imodel}'  :: Objective #'${iobj}'  :: Phase '${iphase}
	    python convert_raw_to_netcdf.py -m ${imodel} -i ../../data/objective_${iobj}/model/${imodel}/${imodel_lower}_phase_${iphase}_objective_${iobj}.csv -o ../../data/objective_${iobj}/model/${imodel}/${imodel_lower}_phase_${iphase}_objective_${iobj}.nc -a ../../data/objective_${iobj}/gauge_info.csv ${add_inputs}

	done

    done

done

echo ''
echo 'Done.'
echo ''

exit 0
