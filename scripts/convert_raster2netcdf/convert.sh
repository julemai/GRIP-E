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

convert_dem_90m=0
convert_flowacc_90m=0
convert_flowdir_90m=0
convert_dem_500m=1
convert_flowacc_500m=1
convert_flowdir_500m=1
convert_soil=0
convert_landcover=0

# -------------------
# convert DEM (90m)
# -------------------
if [[ ${convert_dem_90m} -eq 1 ]] ; then
    echo "Convert DEM ..." 
    varname="dem"
    vartype="single"
    unit="m"
    description="digital elevation model from HydroSheds (USGS) based on conditioned, global SRTM DEM at 3 sec (90m) resolution"

    python raster2netcdf.py -i ${datapath}dem_conditioned-SRTM-90m/rect_dem_Erie.txt -o ${datapath}dem_conditioned-SRTM-90m/rect_dem_Erie.nc -v "${varname},${vartype},${unit},${description}"

fi

# -------------------
# convert flow accumulation (90m)
# -------------------
if [[ ${convert_flowacc_90m} -eq 1 ]] ; then
    echo "Convert flow accumulation ..."
    varname="flow_acc"
    vartype="int32"
    unit="1"
    description="number of cells draining into this grid cell; derived from DEM from HydroSheds by USGS based on conditioned, global SRTM DEM at 3 sec (90m) resolution"

    python raster2netcdf.py -i ${datapath}dem_conditioned-SRTM-90m/rect_flow_accumulation_Erie.txt -o ${datapath}dem_conditioned-SRTM-90m/rect_flow_accumulation_Erie.nc -v "${varname},${vartype},${unit},${description}"

fi

# -------------------
# convert flow direction (90m)
# -------------------
if [[ ${convert_flowdir_90m} -eq 1 ]] ; then
    echo "Convert flow direction ..."
    varname="flow_dir"
    vartype="int32"
    unit="1"
    description="flow direction of grid cell; 1-east, 2-southeast, 4-south, 8-southwest, 16-west, 32-northwest, 64-north, 128-northeast; derived from DEM from HydroSheds by USGS based on conditioned, global SRTM DEM at 3 sec (90m) resolution"

    python raster2netcdf.py -i ${datapath}dem_conditioned-SRTM-90m/rect_flow_direction_Erie.txt -o ${datapath}dem_conditioned-SRTM-90m/rect_flow_direction_Erie.nc -v "${varname},${vartype},${unit},${description}"

fi

# -------------------
# convert DEM (500m)
# -------------------
if [[ ${convert_dem_500m} -eq 1 ]] ; then
    echo "Convert DEM ..." 
    varname="dem"
    vartype="single"
    unit="m"
    description="digital elevation model from HydroSheds (USGS) based on conditioned, global SRTM DEM at 15 sec (500m) resolution"

    python raster2netcdf.py -i ${datapath}dem_conditioned-SRTM-500m/rect_dem_Erie.txt -o ${datapath}dem_conditioned-SRTM-500m/rect_dem_Erie.nc -v "${varname},${vartype},${unit},${description}"

fi

# -------------------
# convert flow accumulation (500m)
# -------------------
if [[ ${convert_flowacc_500m} -eq 1 ]] ; then
    echo "Convert flow accumulation ..."
    varname="flow_acc"
    vartype="int32"
    unit="1"
    description="number of cells draining into this grid cell; derived from DEM from HydroSheds by USGS based on conditioned, global SRTM DEM at 15 sec (500m) resolution"

    python raster2netcdf.py -i ${datapath}dem_conditioned-SRTM-500m/rect_flow_accumulation_Erie.txt -o ${datapath}dem_conditioned-SRTM-500m/rect_flow_accumulation_Erie.nc -v "${varname},${vartype},${unit},${description}"

fi

# -------------------
# convert flow direction (500m)
# -------------------
if [[ ${convert_flowdir_500m} -eq 1 ]] ; then
    echo "Convert flow direction ..."
    varname="flow_dir"
    vartype="int32"
    unit="1"
    description="flow direction of grid cell; 1-east, 2-southeast, 4-south, 8-southwest, 16-west, 32-northwest, 64-north, 128-northeast; derived from DEM from HydroSheds by USGS based on conditioned, global SRTM DEM at 15 sec (500m) resolution"

    python raster2netcdf.py -i ${datapath}dem_conditioned-SRTM-500m/rect_flow_direction_Erie.txt -o ${datapath}dem_conditioned-SRTM-500m/rect_flow_direction_Erie.nc -v "${varname},${vartype},${unit},${description}"

fi

# -------------------
# convert soil classes
# -------------------
if [[ ${convert_soil} -eq 1 ]] ; then
    echo "Convert soil classes ..."
    varname="soil_class"
    vartype="int32"
    unit="1"
    description="soil classes taken from FAO Harmonized World Soil Database (HWSD) v1.2 at 30 sec (1km) resolution"

    python raster2netcdf.py -i ${datapath}soilclass_HWSD/rect_HWSD_soil_class.txt -o ${datapath}soilclass_HWSD/rect_HWSD_soil_class.nc -v "${varname},${vartype},${unit},${description}" -a ${datapath}soilclass_HWSD/USDA_soil_class_legends.csv

fi

# -------------------
# convert land cover classes
# -------------------
if [[ ${convert_landcover} -eq 1 ]] ; then
    echo "Convert land cover ..."
    varname="land_cover"
    vartype="int32"
    unit="1"
    description="land cover classes taken from Global 500m MODIS MCD12Q1 product (NASA) using classification scheme 2 (UMD)"

    python raster2netcdf.py -i ${datapath}landcover_MODIS/rect_landcover_UMD_scheme_MODIS_Erie.txt -o ${datapath}landcover_MODIS/rect_landcover_UMD_scheme_MODIS_Erie.nc -v "${varname},${vartype},${unit},${description}" -a ${datapath}landcover_MODIS/UMD_scheme_MODIS_legends.csv

fi

exit 0



