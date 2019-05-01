## Conditioned SRTM DEM from HydroSHEDS (15 arcsec; 500m)

1.a Data (ASCII raster format)

    Conditioned SRTM DEM from HydroSHEDS, 15 second (500 m)  -->  file: rect_dem_Erie.txt
    Dervided flow accumulation from DEM using ArcGIS         -->  file: rect_flow_accumulation_Erie.txt
    Dervided flow direction    from DEM using ArcGIS         -->  file: rect_flow_direction_Erie.txt

1.b Downloaded from

    https://hydrosheds.cr.usgs.gov/index.php

1.c Prepared by

    Hongren Shen 
    Ph.D. Student
    Department of Civil and Environmental Engineering 
    University of Waterloo 
    200 University Avenue West 
    Waterloo, Ontario, Canada, N2L3G1 
    Office: E2-2359
    Phone: +1 (226)-899-3796 
    Email: hongren.shen@uwaterloo.ca, shrhongren@gmail.com

2.a Converted data (NetCDF format)

    Conditioned SRTM DEM from HydroSHEDS, 15 second (500 m)  -->  file: rect_dem_Erie.nc
    Dervided flow accumulation from DEM using ArcGIS         -->  file: rect_flow_accumulation_Erie.nc
    Dervided flow direction    from DEM using ArcGIS         -->  file: rect_flow_direction_Erie.nc

2.b Converted by

    Juliane Mai PhD
    Department Civil & Environmental Engineering
    University of Waterloo
    200 University Ave West / Waterloo ON, N2L 3G1 / Canada
    Building PHY / Room 3016
    Phone: +1 (519) 888-4567 x 30016
    Mobil: +1 (226) 505-5600
    http://www.civil.uwaterloo.ca/jmai/
    https://www.researchgate.net/profile/Juliane_Mai

2.c Data conversion

    -----------------
    DEM
    -----------------
	datapath='../../data/dem_conditioned-SRTM-500m/'
    varname="dem"
    vartype="single"
    unit="m"
    description="digital elevation model from HydroSheds (USGS) based on conditioned, global SRTM DEM at 15 sec (500m) resolution"

    python raster2netcdf.py -i ${datapath}dem_conditioned-SRTM-500m/rect_dem_Erie.txt \
		     	        -o ${datapath}dem_conditioned-SRTM-500m/rect_dem_Erie.nc \
		     	        -v "${varname},${vartype},${unit},${description}"

    -----------------
    Flow accumulation
    -----------------
	datapath='../../data/dem_conditioned-SRTM-500m/'
    varname="flow_acc"
    vartype="int32"
    unit="1"
    description="number of cells draining into this grid cell; \
                 derived from DEM from HydroSheds by USGS based on \
                 conditioned, global SRTM DEM at 15 sec (500m) resolution"

    python raster2netcdf.py -i ${datapath}dem_conditioned-SRTM-500m/rect_flow_accumulation_Erie.txt
		     	        -o ${datapath}dem_conditioned-SRTM-500m/rect_flow_accumulation_Erie.nc
		     	        -v "${varname},${vartype},${unit},${description}"

    -----------------
    Flow direction
    -----------------
	datapath='../../data/dem_conditioned-SRTM-500m/'
    varname="flow_dir"
    vartype="int32"
    unit="1"
    description="flow direction of grid cell; 1-east, 2-southeast, 4-south, 8-southwest, \
                 16-west, 32-northwest, 64-north, 128-northeast; \
                 derived from DEM from HydroSheds by USGS based on  \
                 conditioned, global SRTM DEM at 15 sec (500m) resolution"

    python raster2netcdf.py -i ${datapath}dem_conditioned-SRTM-500m/rect_flow_direction_Erie.txt
		     	        -o ${datapath}dem_conditioned-SRTM-500m/rect_flow_direction_Erie.nc
		     	        -v "${varname},${vartype},${unit},${description}"
