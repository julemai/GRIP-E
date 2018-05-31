#!/bin/bash
#
# Produces plots for presentation.
#
set -e
#
prog=$0
pprog=$(basename ${prog})
dprog=$(dirname ${prog})
isdir=${PWD}
pid=$$

datapath="../data/"

convert_dem=0
convert_flowacc=0
convert_flowdir=0
convert_soil=1
convert_landcover=1

# -------------------
# convert DEM
# -------------------
if [[ ${convert_dem} -eq 1 ]] ; then
    echo "Convert DEM ..." 
    varname="dem"
    vartype="single"
    unit="m"
    description="digital elevation model from HydroSheds (USGS) based on conditioned, global SRTM DEM at 3 sec (90m) resolution"

    python raster2netcdf.py -i ${datapath}dem_conditioned-SRTM-90m/rect_dem_Erie.txt -o ${datapath}dem_conditioned-SRTM-90m/rect_dem_Erie.nc -v "${varname},${vartype},${unit},${description}"

fi

# -------------------
# convert flow accumulation
# -------------------
if [[ ${convert_flowacc} -eq 1 ]] ; then
    echo "Convert flow accumulation ..."
    varname="flow_acc"
    vartype="int32"
    unit="1"
    description="number of cells draining into this grid cell; derived from DEM from HydroSheds by USGS based on conditioned, global SRTM DEM at 3 sec (90m) resolution"

    python raster2netcdf.py -i ${datapath}dem_conditioned-SRTM-90m/rect_flow_accumulation_Erie.txt -o ${datapath}dem_conditioned-SRTM-90m/rect_flow_accumulation_Erie.nc -v "${varname},${vartype},${unit},${description}"

fi

# -------------------
# convert flow direction
# -------------------
if [[ ${convert_flowdir} -eq 1 ]] ; then
    echo "Convert flow direction ..."
    varname="flow_dir"
    vartype="int32"
    unit="1"
    description="flow direction of grid cell; 1-east, 2-southeast, 4-south, 8-southwest, 16-west, 32-northwest, 64-north, 128-northeast; derived from DEM from HydroSheds by USGS based on conditioned, global SRTM DEM at 3 sec (90m) resolution"

    python raster2netcdf.py -i ${datapath}dem_conditioned-SRTM-90m/rect_flow_direction_Erie.txt -o ${datapath}dem_conditioned-SRTM-90m/rect_flow_direction_Erie.nc -v "${varname},${vartype},${unit},${description}"

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



