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

convert_models='GEM-Hydro'                   # [    Lake Erie:   LBRM  HMETS-Raven-lp GR4J-Raven-lp GR4J-Raven-sd HYPE
#                                                                GEM-Hydro WRF-Hydro MESH-SVS MESH-CLASS VIC VIC-GRU
#                                                                WATFLOOD SWAT-EPA SWAT-Guelph
#                                                                ML-ConvLSTM ML-ConvLSTM-DEM ML-ConvLSTM-LC ML-ConvLSTM-LC-DEM ML-LinReg ML-XGBoost
#                                                                mHM-UFZ mHM-Waterloo
#                                                                Raven-blended
#                                            #      Great Lakes: GR4J-Raven-lp GR4J-Raven-sd LBRM-MG LBRM-ML-LSTM ML-EA-LSTM ML-LSTM ML-XGBoost]

# domain='great-lakes'                       # [lake-erie great-lakes]
# calvals='calibration validation' # validation'             # [calibration validation]  # only for Great Lakes

domain='lake-erie'                           # [lake-erie great-lakes]
calvals='calibration'                        # [calibration validation]  # only for Great Lakes

setup_by='julie'                             # Raven setup by 'julie' (outputs in separate files) or 'hongren' (outputs in one file)
convert_obj='1 2'                            # can be 1, 2, and/or 3
convert_phase='1'                          # phase 0: uncalibrated, different phys. setups,
#                                            # phase 1: calibrated,   different phys. setups,
#                                            # phase 2: calibrated,   same phys. setups

for calval in ${calvals} ; do
    
    for imodel in ${convert_models} ; do

        imodel_lower=$( echo "$imodel" | tr '[:upper:]' '[:lower:]' )
	echo ""
	echo "-------------------------"
        echo ${imodel_lower}
	echo "-------------------------"

        for iobj in ${convert_obj} ; do

            for iphase in ${convert_phase} ; do

                # ------------------
                # for great lakes split gauge-info file into only calibration and only validation stations
                # ------------------
                file_all_stations="../../data/objective_${iobj}/${domain}/gauge_info.csv"
                if [[ ${domain} == 'great-lakes' ]] ; then
                    
                    file_only_cal_or_val="../../data/objective_${iobj}/${domain}/${calval}/gauge_info.csv"
                    
                    if [[  ${calval} == 'calibration' ]] ; then
                        
                        head -1 ${file_all_stations}     > ${file_only_cal_or_val}
                        grep ',C,' ${file_all_stations} >> ${file_only_cal_or_val}
                        
                    else
                        
                        head -1 ${file_all_stations}     > ${file_only_cal_or_val}
                        grep ',V,' ${file_all_stations} >> ${file_only_cal_or_val}
                        
                    fi
                fi
                # ------------------

                echo ''
                echo 'Convert :: '${imodel}'  :: Objective #'${iobj}'  :: Phase '${iphase}

                if [[ ( ${domain} == 'lake-erie' ) ]] ; then 
                    if [[ ( ${imodel} == 'VIC' ) || ( ${imodel} == 'VIC-GRU' ) || ( ${imodel} == 'SWAT-EPA' ) || ( ${imodel} == 'SWAT-Guelph' ) || ( ${imodel} == 'GR4J-Raven-lp' && ${setup_by} == 'julie' ) || ( ${imodel} == 'GR4J-Raven-sd' && ${setup_by} == 'julie' )  || ( ${imodel} == 'HMETS-Raven-lp' && ${setup_by} == 'julie' )  || ( ${imodel} == 'Raven-blended' && ${setup_by} == 'julie' ) ]] ; then
                        add_inputs="-b ../../data/objective_${iobj}/${domain}/model/${imodel}/subid2gauge.csv"
                    else
                        if [[ ( ${imodel} == 'MESH-SVS' ) || ( ${imodel} == 'MESH-CLASS' ) ]] ; then
                            add_inputs="-b ../../data/objective_${iobj}/${domain}/model/${imodel}/subid2gauge.tb0"
                        else
                            add_inputs=''
                        fi
                    fi
                else
                    if [[ ( ${imodel} == 'VIC' ) || ( ${imodel} == 'VIC-GRU' )  || ( ${imodel} == 'SWAT-EPA' ) || ( ${imodel} == 'SWAT-Guelph' ) || ( ${imodel} == 'GR4J-Raven-lp' && ${setup_by} == 'julie' ) || ( ${imodel} == 'GR4J-Raven-sd' && ${setup_by} == 'julie' )  || ( ${imodel} == 'HMETS-Raven-lp' && ${setup_by} == 'julie' )  || ( ${imodel} == 'Raven-blended' && ${setup_by} == 'julie' ) ]] ; then
                        add_inputs="-b ../../data/objective_${iobj}/${domain}/${calval}/model/${imodel}/subid2gauge.csv"
                    else
                        if [[ ( ${imodel} == 'MESH-SVS' ) || ( ${imodel} == 'MESH-CLASS' ) ]] ; then
                            add_inputs="-b ../../data/objective_${iobj}/${domain}/${calval}/model/${imodel}/subid2gauge.tb0"
                        else
                            add_inputs=''
                        fi
                    fi
                fi

                if [[ ( ${domain} == 'lake-erie' ) ]] ; then 
                    if [[ ( ${imodel} == 'mHM-UFZ' ) || ( ${imodel} == 'mHM-Waterloo' ) || ( ${imodel} == 'HYPE' )  || ( ${imodel} == 'GR4J-Raven-lp' && ${setup_by} == 'julie') || ( ${imodel} == 'GR4J-Raven-sd' && ${setup_by} == 'julie' ) || ( ${imodel} == 'HMETS-Raven-lp' && ${setup_by} == 'julie')  || ( ${imodel} == 'Raven-blended' && ${setup_by} == 'julie' ) ]] ; then
                        input_csv_file=../../data/objective_${iobj}/${domain}/model/${imodel}/${imodel_lower}_phase_${iphase}_objective_${iobj}_
                    else
                        input_csv_file=../../data/objective_${iobj}/${domain}/model/${imodel}/${imodel_lower}_phase_${iphase}_objective_${iobj}.csv
                    fi
                else
                    if [[ ( ${imodel} == 'mHM-UFZ' ) || ( ${imodel} == 'mHM-Waterloo' ) || ( ${imodel} == 'HYPE' )  || ( ${imodel} == 'GR4J-Raven-lp' && ${setup_by} == 'julie') || ( ${imodel} == 'GR4J-Raven-sd' && ${setup_by} == 'julie' ) || ( ${imodel} == 'HMETS-Raven-lp' && ${setup_by} == 'julie')  || ( ${imodel} == 'Raven-blended' && ${setup_by} == 'julie' ) ]] ; then
                        input_csv_file=../../data/objective_${iobj}/${domain}/${calval}/model/${imodel}/${imodel_lower}_phase_${iphase}_objective_${iobj}_
                    else
                        input_csv_file=../../data/objective_${iobj}/${domain}/${calval}/model/${imodel}/${imodel_lower}_phase_${iphase}_objective_${iobj}.csv
                    fi
                fi
                
                if [[ ( ${domain} == 'lake-erie' ) ]] ; then 
                    if [[ ( ${imodel} == 'GR4J-Raven-lp' || ${imodel} == 'GR4J-Raven-sd'  || ${imodel} == 'HMETS-Raven-lp'  || ${imodel} == 'Raven-blended' ) ]] ; then
                        python convert_raw_to_netcdf.py -m ${imodel} -i ${input_csv_file} -o ../../data/objective_${iobj}/${domain}/model/${imodel}/${imodel_lower}_phase_${iphase}_objective_${iobj}.nc -a ${file_all_stations} ${add_inputs} -s ${setup_by}
                    else
                        python convert_raw_to_netcdf.py -m ${imodel} -i ${input_csv_file} -o ../../data/objective_${iobj}/${domain}/model/${imodel}/${imodel_lower}_phase_${iphase}_objective_${iobj}.nc -a ${file_all_stations} ${add_inputs}
                    fi
                else
                    if [[ ( ${imodel} == 'GR4J-Raven-lp' || ${imodel} == 'GR4J-Raven-sd'  || ${imodel} == 'HMETS-Raven-lp' || ${imodel} == 'Raven-blended'  ) ]] ; then
                        python convert_raw_to_netcdf.py -m ${imodel} -i ${input_csv_file} -o ../../data/objective_${iobj}/${domain}/${calval}/model/${imodel}/${imodel_lower}_phase_${iphase}_objective_${iobj}.nc -a ${file_only_cal_or_val} ${add_inputs} -s ${setup_by}
                    else
                        python convert_raw_to_netcdf.py -m ${imodel} -i ${input_csv_file} -o ../../data/objective_${iobj}/${domain}/${calval}/model/${imodel}/${imodel_lower}_phase_${iphase}_objective_${iobj}.nc -a ${file_only_cal_or_val} ${add_inputs}
                    fi
                fi

            done # convert_phase

        done # convert_obj

    done # convert_models

done # calvals

echo ''
echo 'Done.'
echo ''

exit 0
