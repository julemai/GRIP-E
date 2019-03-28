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

plot_obj='1 2'                             # can be 1, 2, and/or 3
plot_phase='1'                           # phase 0: uncalibrated, different phys. setups,
#                                        # phase 1: calibrated,   different phys. setups,
#                                        # phase 2: calibrated,   same phys. setups

for iobj in ${plot_obj} ; do

    for iphase in ${plot_phase} ; do

	echo '----------------------'
	echo "plot:  compare_models_phase_${iphase}_objective_${iobj}.pdf"
	echo '----------------------'
	files=$(\ls ../../data/objective_${iobj}/model/*/*_phase_${iphase}_objective_${iobj}.nc)
	files=$(echo ${files})
	python compare_models.py -i "${files}" -a '2011-01-01:2014-12-31' -p compare_models_phase_${iphase}_objective_${iobj}.pdf
	pdfcrop compare_models_phase_${iphase}_objective_${iobj}.pdf
	mv compare_models_phase_${iphase}_objective_${iobj}-crop.pdf compare_models_phase_${iphase}_objective_${iobj}.pdf

    done

done

echo ''
echo 'Done.'
echo ''

exit 0
