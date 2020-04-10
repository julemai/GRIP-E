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
#
# Converts ArcGIS raster format files to NetCDF.
#
set -e
#
prog=$0
pprog=$(basename ${prog})
dprog=$(dirname ${prog})
isdir=${PWD}
pid=$$

datapath="../../data/"

convert_soil_GSDE=1
convert_landcover_NACLMS=0


grid='RDRS-v2'
grid='WFDEI-GEM-CaPA'

# -------------------
# convert soil classes of GSDE
# -------------------
if [[ ${convert_soil_GSDE} -eq 1 ]] ; then
    echo "Convert aggregated soil classes and texture based on GSDE into ${grid} grid ..."

    if [[ ${grid} == 'RDRS-v2' ]] ; then
	inputfile="${datapath}soilclass_GSDE_GreatLakes/soilclass_GSDE_GreatLakes_aggregated_v1.1/GL_GSDE_usda_soil_class_rdrs_v2.txt"
	gridfile="${datapath}meteo_forcing_RDRS-v2/grip-gl_rdrs-v2-gridonly.nc"
	legendsfile=""
	varinfofile="${datapath}soilclass_GSDE_GreatLakes/soilclass_GSDE_GreatLakes_aggregated_v1.1/variable_info.csv"
    else
	if [[ ${grid} == 'WFDEI-GEM-CaPA' ]] ; then
	    inputfile="${datapath}soilclass_GSDE_GreatLakes/soilclass_GSDE_GreatLakes_aggregated_v1.1/GL_GSDE_usda_soil_class_wfdei_gem_capa.txt"
	    gridfile="${datapath}meteo_forcings_WFDEI-GEM-CaPA/grip-gl_wfdei-gem-capa_gridonly.nc"
	    legendsfile=""
	    varinfofile="${datapath}soilclass_GSDE_GreatLakes/soilclass_GSDE_GreatLakes_aggregated_v1.1/variable_info.csv"
	fi
    fi

    if [ ! -e ${inputfile} ]; then
	echo "Inputfile ${inputfile} does not exist!"
	exit
    fi

    if [ ! -e ${gridfile} ]; then
	echo "Grid defining file ${gridfile} does not exist!"
	exit
    fi

    if [[ ${legendsfile} != "" ]] ; then
	if [ ! -e ${legendsfile} ]; then
	    echo "File with legends ${legendsfile} does not exist!"
	    exit
	fi
    fi

    if [[ ${varinfofile} != "" ]] ; then
	if [ ! -e ${varinfofile} ]; then
	    echo "File with variable infos ${varinfofile} does not exist!"
	    exit
	fi
    fi

    outputfile=$( echo $( echo ${inputfile} | rev | cut -d '/' -f 2- | rev )"_${grid}.nc" )

    echo "python aggregated2netcdf.py -i ${inputfile} -g ${gridfile} -o ${outputfile} -l ${varinfofile}"
    python aggregated2netcdf.py -i ${inputfile} -g ${gridfile} -o ${outputfile} -l ${varinfofile}

fi

# -------------------
# convert aggreagted land cover classes of NACLMS
# -------------------
if [[ ${convert_landcover_NACLMS} -eq 1 ]] ; then
    echo "Convert land cover NACLMS into ${grid} grid ..."
    varname="perc_land_cover"
    vartype="float"
    unit="percent"
    description="percentage of land cover type of NALCMS product for each forcing grid cell"

    if [[ ${grid} == 'RDRS-v2' ]] ; then
	inputfile="${datapath}landcover_NALCMS_GreatLakes/landcover_NALCMS_GreatLakes_aggregated/GL_NALCMS_landcover_rdrs_v2.txt"
	gridfile="${datapath}meteo_forcing_RDRS-v2/grip-gl_rdrs-v2-gridonly.nc"
	#
	# # stores variables 4D (time, lat, lon, n_LC_types)
	# legendsfile="${datapath}landcover_NALCMS_GreatLakes/NACLMS_legends.csv"
	# varinfofile=""
	#
	# stores variables 3D (time, lat, lon) in <n_LC_types> separate variables
	legendsfile=""  
	varinfofile="${datapath}landcover_NALCMS_GreatLakes/landcover_NALCMS_GreatLakes_aggregated/variable_info.csv"
    else
	if [[ ${grid} == 'WFDEI-GEM-CaPA' ]] ; then
	    inputfile="${datapath}landcover_NALCMS_GreatLakes/landcover_NALCMS_GreatLakes_aggregated/GL_NALCMS_landcover_wfdei_gem_capa.txt"
	    gridfile="${datapath}meteo_forcings_WFDEI-GEM-CaPA/grip-gl_wfdei-gem-capa_gridonly.nc"
	    #
	    # # stores variables 4D (time, lat, lon, n_LC_types)
	    # legendsfile="${datapath}landcover_NALCMS_GreatLakes/NACLMS_legends.csv"
	    # varinfofile=""
	    #
	    # stores variables 3D (time, lat, lon) in <n_LC_types> separate variables
	    legendsfile="" 
	    varinfofile="${datapath}landcover_NALCMS_GreatLakes/landcover_NALCMS_GreatLakes_aggregated/variable_info.csv"
	fi
    fi

    if [ ! -e ${inputfile} ]; then
	echo "Inputfile ${inputfile} does not exist!"
	exit
    fi

    if [ ! -e ${gridfile} ]; then
	echo "Grid defining file ${gridfile} does not exist!"
	exit
    fi

    if [[ ${legendsfile} != "" ]] ; then
	if [ ! -e ${legendsfile} ]; then
	    echo "File with legends ${legendsfile} does not exist!"
	    exit
	fi
    fi

    if [[ ${varinfofile} != "" ]] ; then
	if [ ! -e ${varinfofile} ]; then
	    echo "File with variable infos ${varinfofile} does not exist!"
	    exit
	fi
    fi

    outputfile=$( echo $( echo ${inputfile} | rev | cut -d '/' -f 2- | rev )"_${grid}.nc" )

    # stores variables 4D (time, lat, lon, n_LC_types)
    # python aggregated2netcdf.py -i ${inputfile} -g ${gridfile} -o ${outputfile} -v "${varname},${vartype},${unit},${description}" -a ${legendsfile}

    # stores variables 3D (time, lat, lon) in <n_LC_types> separate variables
    echo "python aggregated2netcdf.py -i ${inputfile} -g ${gridfile} -o ${outputfile} -v "${varname},${vartype},${unit},${description}" -l ${varinfofile}"
    python aggregated2netcdf.py -i ${inputfile} -g ${gridfile} -o ${outputfile} -v "${varname},${vartype},${unit},${description}" -l ${varinfofile}
    

fi

exit 0



