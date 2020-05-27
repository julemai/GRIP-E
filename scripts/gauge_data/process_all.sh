#!/bin/bash

# Copyright 2016-2018 Juliane Mai - juliane.mai(at)uwaterloo.ca
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

# This script converts all CSV streamflow gauge data files into NetCDF format.

# Perform a cleanup if script is interupted
trap cleanup 1 2 3 6
#
prog=$0
pprog=$(basename ${prog})
dprog=$(dirname ${prog})
isdir=${PWD}
pid=$$

# tell me what to do
convert_to_single_nc=0       # converts every csv/txt into an individual NetCDF 	   --> "<station-id>.nc"
convert_to_merged_nc=1       # converts all   csv/txt's into one single NetCDF  	   --> "all_gauges.nc"
plot_merged_nc=0             # plots all stations of a NetCDF into pdf          	   --> "all_gauges.pdf"
plot_grouped_gauges=0        # plots stations grouped by watersheds (only for great-lakes) --> "all_gauges_grouped.pdf"

objectives="1 2"
domain="lake-erie"           # lake-erie or great-lakes   # choose only ONE
calvals="validation"         # calibration validation     # choose as many as you want


    

if [[ ${domain} == "lake-erie" ]] ; then
    plot_grouped_gauges=0
fi

for calval in ${calvals} ; do

    if [[ ${calval} == "None" ]] ; then
	calval=""
    fi
	
    for objective in ${objectives} ; do

	if [[ ${domain} == 'great-lakes' ]] ; then
	    gaugeinfofile="../../data/objective_${objective}/${domain}/gauge_info.csv"
	else
	    gaugeinfofile="../../data/objective_${objective}/${domain}/${calval}/gauge_info.csv"
	fi

	input_files=""
	filetype=""

	ids=$( grep -i -v 'ID,Name,Lat,Lon' "${gaugeinfofile}" | awk -F, '{ print $2 }' )
	echo ${ids}
	for ii in ${ids} ; do
	    
	    ifile="../../data/objective_${objective}/${domain}/${calval}/csv/${ii}.txt"   # this is the USGS data
	    if [ -e ${ifile} ] ; then
		# echo ""
		echo ${ifile}

		input_files+=${ifile}' '
		filetype+='USGS '
		
		# construct output file name
		tmp2=$( echo ${ifile} | rev | cut -d '.' -f 2- | cut -d '/' -f 1  | rev )
		tmp1=$( echo ${ifile} | rev | cut -d '.' -f 2- | cut -d '/' -f 3- | rev )
		ofile=$( echo "${tmp1}/netcdf/${tmp2}.nc")
		
		if [ ${convert_to_single_nc} == 1 ] ; then
		    python convert_gauge_csv_2_netcdf.py --filetype USGS --input_files ${ifile} --output_file ${ofile} --gaugeinfo_file "${gaugeinfofile}"
		fi
	    fi

	    ifile="../../data/objective_${objective}/${domain}/${calval}/csv/${ii}.csv"   # this is the WSC data
	    if [ -e ${ifile} ] ; then
		# echo ""
		echo ${ifile}

		input_files+=${ifile}' '
		filetype+='WSC '
		
		# construct output file name
		tmp2=$( echo ${ifile} | rev | cut -d '.' -f 2- | cut -d '/' -f 1  | rev )
		tmp1=$( echo ${ifile} | rev | cut -d '.' -f 2- | cut -d '/' -f 3- | rev )
		ofile=$( echo "${tmp1}/netcdf/${tmp2}.nc")
		
		if [ ${convert_to_single_nc} == 1 ] ; then
		    python convert_gauge_csv_2_netcdf.py --filetype WSC --input_files ${ifile} --output_file ${ofile} --gaugeinfo_file "${gaugeinfofile}"
		fi
	    fi
	done

	if [ ${convert_to_merged_nc} == 1 ] ; then
	    ofile=$( echo "${tmp1}/netcdf/all_gauges.nc")
	    python convert_gauge_csv_2_netcdf.py --filetype "$(echo ${filetype})" --input_files "$(echo ${input_files})" --output_file ${ofile} --gaugeinfo_file "${gaugeinfofile}"
	fi

	if [ ${plot_merged_nc} == 1 ] ; then
	    python plot_gauge_data.py -i ../../data/objective_${objective}/${domain}/${calval}/netcdf/all_gauges.nc -p ../../data/objective_${objective}/${domain}/${calval}/all_gauges.pdf
	    pdfcrop ../../data/objective_${objective}/${domain}/${calval}/all_gauges.pdf
	    mv      ../../data/objective_${objective}/${domain}/${calval}/all_gauges-crop.pdf ../../data/objective_${objective}/${domain}/${calval}/all_gauges.pdf
	fi

    done
done

if [ ${plot_grouped_gauges} == 1 ] ; then

    # time period 2010-2014
    python plot_gauge_data_grouped.py -i ../../data/objective_1/${domain}/gauge_info.csv ../../data/objective_2/${domain}/gauge_info.csv -p ../../data/all_gauges_grouped_2010-2014.pdf -v Q -y 2010:2014
    pdfcrop ../../data/all_gauges_grouped_2010-2014.pdf
    mv      ../../data/all_gauges_grouped_2010-2014-crop.pdf ../../data/all_gauges_grouped_2010-2014.pdf

    # time period 2000-2018
    python plot_gauge_data_grouped.py -i ../../data/objective_1/${domain}/gauge_info.csv ../../data/objective_2/${domain}/gauge_info.csv -p ../../data/all_gauges_grouped_2000-2018.pdf -v Q -y 2000:2018
    pdfcrop ../../data/all_gauges_grouped_2000-2018.pdf
    mv      ../../data/all_gauges_grouped_2000-2018-crop.pdf ../../data/all_gauges_grouped_2000-2018.pdf
    
fi


exit 0
