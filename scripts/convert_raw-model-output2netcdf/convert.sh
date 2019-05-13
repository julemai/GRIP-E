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

convert_models='GR4J-Raven-lp GR4J-Raven-sd' # 'LBRM VIC VIC-GRU GEM-Hydro SWAT WATFLOOD GR4J-Raven-lp GR4J-Raven-sd'  # can be [LBRM, HYPE, GEM-Hydro, WRF-Hydro, MESH-SVS, MESH-CLASS, VIC, VIC-GRU, WATFLOOD]
setup_by='julie'            # Raven setup by 'julie' (outputs in separate files) or 'hongren' (outputs in one file)
convert_obj='1 2'      	    # can be 1, 2, and/or 3
convert_phase='0 1'         # phase 0: uncalibrated, different phys. setups,
#                      	    # phase 1: calibrated,   different phys. setups,
#                      	    # phase 2: calibrated,   same phys. setups

for imodel in ${convert_models} ; do

    imodel_lower=$( echo "$imodel" | tr '[:upper:]' '[:lower:]' )
    echo ${imodel_lower}

    for iobj in ${convert_obj} ; do

	for iphase in ${convert_phase} ; do

	    echo ''
	    echo 'Convert :: '${imodel}'  :: Objective #'${iobj}'  :: Phase '${iphase}

	    if [[ ( ${imodel} == 'VIC' ) || ( ${imodel} == 'VIC-GRU' ) || ( ${imodel} == 'SWAT' ) || ( ${imodel} == 'GR4J-Raven-lp' && ${setup_by} == 'julie' ) || ( ${imodel} == 'GR4J-Raven-sd' && ${setup_by} == 'julie' ) ]] ; then
		add_inputs="-b ../../data/objective_${iobj}/model/${imodel}/subid2gauge.csv"
	    else
		if [[ ( ${imodel} == 'MESH-SVS' ) || ( ${imodel} == 'MESH-CLASS' ) ]] ; then
		    add_inputs="-b ../../data/objective_${iobj}/model/${imodel}/subid2gauge.tb0"
		else
		    add_inputs=''
		fi
	    fi

	    if [[ ( ${imodel} == 'HYPE' )  || ( ${imodel} == 'GR4J-Raven-lp' && ${setup_by} == 'julie') || ( ${imodel} == 'GR4J-Raven-sd' && ${setup_by} == 'julie' ) ]] ; then
		input_csv_file=../../data/objective_${iobj}/model/${imodel}/${imodel_lower}_phase_${iphase}_objective_${iobj}_
	    else
		input_csv_file=../../data/objective_${iobj}/model/${imodel}/${imodel_lower}_phase_${iphase}_objective_${iobj}.csv
	    fi

	    if [[ ( ${imodel} == 'GR4J-Raven-lp' || ${imodel} == 'GR4J-Raven-sd' ) ]] ; then
		python convert_raw_to_netcdf.py -m ${imodel} -i ${input_csv_file} -o ../../data/objective_${iobj}/model/${imodel}/${imodel_lower}_phase_${iphase}_objective_${iobj}.nc -a ../../data/objective_${iobj}/gauge_info.csv ${add_inputs} -s ${setup_by}
	    else
		python convert_raw_to_netcdf.py -m ${imodel} -i ${input_csv_file} -o ../../data/objective_${iobj}/model/${imodel}/${imodel_lower}_phase_${iphase}_objective_${iobj}.nc -a ../../data/objective_${iobj}/gauge_info.csv ${add_inputs}
	    fi

	done

    done

done

echo ''
echo 'Done.'
echo ''

exit 0
