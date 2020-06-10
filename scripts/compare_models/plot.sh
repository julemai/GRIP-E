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

plot_obj='1 2'                                            # can be 1, 2, and/or 3
plot_phase='1'                                          # phase 0: uncalibrated, different phys. setups,
#                                                         # phase 1: calibrated,   different phys. setups,
#                                                         # phase 2: calibrated,   same phys. setups

# domain='lake-erie'                                        # [lake-erie, great-lakes]
# periods='2011-01-01:2014-12-31'                           # time period(s) that should be used to derive NSE etc
# calval='calibration'                                      # [calibration, validation]  # choose ONE only

domain='lake-erie'                                      # [lake-erie, great-lakes]
periods='2011-01-01:2014-12-31'                         # time period(s) that should be used to derive NSE etc
calval='validation'                                     # [calibration, validation]  # choose ONE only

# domain='great-lakes'                                    # [lake-erie, great-lakes]
# periods='2001-01-01:2010-12-31 2011-01-01:2016-12-31'   # time period(s) that should be used to derive NSE etc
# calval='calibration'                                    # [calibration, validation]  # choose ONE only

# domain='great-lakes'                                    # [lake-erie, great-lakes]
# periods='2001-01-01:2010-12-31'                         # time period(s) that should be used to derive NSE etc
# calval='validation'                                     # [calibration, validation]  # choose ONE only


for iobj in ${plot_obj} ; do

    path="../../data/objective_${iobj}/${domain}/${calval}"
    ext="_${calval}"

    for iphase in ${plot_phase} ; do

        for period in ${periods} ; do

            period_str=$( echo ${period//':'/'_'} )     # "2001-01-01:2010-12-31" --> "2001-01-01_2010-12-31"

            echo '----------------------'
            echo "plot:  compare_models_phase_${iphase}_objective_${iobj}_${domain}_${period_str}${ext}.pdf"
            echo '----------------------'
            files=$(\ls ${path}/model/*/*_phase_${iphase}_objective_${iobj}.nc)
            files=$( echo ${files} )
            
            # given -y does not sort models (y-axis)
            python compare_models.py -i "${files}" -a ${period} -p compare_models_phase_${iphase}_objective_${iobj}_${domain}_${period_str}${ext}.pdf -y
            pdfcrop compare_models_phase_${iphase}_objective_${iobj}_${domain}_${period_str}${ext}.pdf
	    pdfsplit compare_models_phase_${iphase}_objective_${iobj}_${domain}_${period_str}${ext}-crop.pdf
            mv compare_models_phase_${iphase}_objective_${iobj}_${domain}_${period_str}${ext}-crop1.pdf compare_models_phase_${iphase}_objective_${iobj}_${domain}_${period_str}${ext}_NSE.pdf
	    mv compare_models_phase_${iphase}_objective_${iobj}_${domain}_${period_str}${ext}-crop2.pdf compare_models_phase_${iphase}_objective_${iobj}_${domain}_${period_str}${ext}_PBIAS.pdf
	    rm compare_models_phase_${iphase}_objective_${iobj}_${domain}_${period_str}${ext}-crop.pdf
	    rm compare_models_phase_${iphase}_objective_${iobj}_${domain}_${period_str}${ext}.pdf
        done

    done

done

if [ -e compare_models_phase_1_objective_1_lake-erie_2011-01-01_2014-12-31_calibration_NSE.pdf ] ; then
    cp compare_models_phase_1_objective_1_lake-erie_2011-01-01_2014-12-31_calibration_NSE.pdf figure_2.pdf
    convert -density 300 -quality 95 compare_models_phase_1_objective_1_lake-erie_2011-01-01_2014-12-31_calibration_NSE.pdf compare_models_phase_1_objective_1_lake-erie_2011-01-01_2014-12-31_calibration_NSE.png
fi

if [ -e compare_models_phase_1_objective_2_lake-erie_2011-01-01_2014-12-31_calibration_NSE.pdf ] ; then
    cp compare_models_phase_1_objective_2_lake-erie_2011-01-01_2014-12-31_calibration_NSE.pdf figure_3.pdf
    convert -density 300 -quality 95 compare_models_phase_1_objective_2_lake-erie_2011-01-01_2014-12-31_calibration_NSE.pdf compare_models_phase_1_objective_2_lake-erie_2011-01-01_2014-12-31_calibration_NSE.png
fi

python barchart.py -p barchart.pdf
pdfcrop barchart.pdf
pdfsplit barchart-crop.pdf
mv barchart-crop1.pdf barchart_obj_1-2.pdf
mv barchart-crop2.pdf barchart_obj_1.pdf
mv barchart-crop3.pdf barchart_obj_2.pdf
rm barchart.pdf
rm barchart-crop.pdf

echo ''
echo 'Done.'
echo ''

exit 0
