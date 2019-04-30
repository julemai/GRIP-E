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

plot_models='VIC' #'LBRM VIC VIC-GRU GEM-Hydro SWAT RAVEN-GR4J' # can be [LBRM, HYPE, GEM-Hydro, WRF-Hydro, MESH-SVS, MESH-CLASS, VIC, VIC-GRU, WATFLOOD]
plot_obj='1 2'                             # can be 1, 2, and/or 3
plot_phase='0 1'                           # phase 0: uncalibrated, different phys. setups,
#                                        # phase 1: calibrated,   different phys. setups,
#                                        # phase 2: calibrated,   same phys. setups

for imodel in ${plot_models} ; do

    imodel_lower=$( echo "$imodel" | tr '[:upper:]' '[:lower:]' )
    echo ${imodel_lower}

    for iobj in ${plot_obj} ; do

	for iphase in ${plot_phase} ; do

	    echo ''
	    echo 'Plot :: '${imodel}'  :: Objective #'${iobj}'  :: Phase '${iphase}
	    python plot_nc_model_output.py -a '2011-01-01:2014-12-31' -i ../../data/objective_${iobj}/model/${imodel}/${imodel_lower}_phase_${iphase}_objective_${iobj}.nc -p ../../data/objective_${iobj}/model/${imodel}/${imodel_lower}_phase_${iphase}_objective_${iobj}.pdf
	    pdfcrop ../../data/objective_${iobj}/model/${imodel}/${imodel_lower}_phase_${iphase}_objective_${iobj}.pdf
	    pdfsplit ../../data/objective_${iobj}/model/${imodel}/${imodel_lower}_phase_${iphase}_objective_${iobj}-crop.pdf
	    mv ../../data/objective_${iobj}/model/${imodel}/${imodel_lower}_phase_${iphase}_objective_${iobj}-crop1.pdf ../../data/objective_${iobj}/model/${imodel}/${imodel_lower}_phase_${iphase}_objective_${iobj}_hydrographs.pdf
	    mv ../../data/objective_${iobj}/model/${imodel}/${imodel_lower}_phase_${iphase}_objective_${iobj}-crop2.pdf ../../data/objective_${iobj}/model/${imodel}/${imodel_lower}_phase_${iphase}_objective_${iobj}_performance.pdf
	    rm ../../data/objective_${iobj}/model/${imodel}/${imodel_lower}_phase_${iphase}_objective_${iobj}-crop.pdf
	    rm ../../data/objective_${iobj}/model/${imodel}/${imodel_lower}_phase_${iphase}_objective_${iobj}.pdf

	done

    done

done

echo ''
echo 'Done.'
echo ''

exit 0
