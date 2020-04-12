## Conditioned SRTM DEM from HydroSHEDS (3 arcsec; 90m)

External download from:

* DEM, flow direction and accumulation (3'', ~90m grid) in ASCII raster format<br>
  http://www.civil.uwaterloo.ca/jmai/GRIP-GL/dem_HydroSHEDS-90m_GreatLakes.zip (577.5 MB)
* Slope derived from DEM (3'', ~90m grid) in TIF format<br>
  http://www.civil.uwaterloo.ca/jmai/GRIP-GL/slope_HydroSHEDS-90m_GreatLakes.tif_v1.1.zip (322.1 MB)
* Aggregated slope to forcing grids in tabulated TXT format<br>
  https://github.com/julemai/GRIP-E/blob/master/data/dem_HydroSHEDS-90m_GreatLakes/slope_HydroSHEDS-90m_GreatLakes_aggregated_v1.1.zip
* Aggregated slope to forcing grids in NetCDF format<br>
  https://github.com/julemai/GRIP-E/blob/master/data/dem_HydroSHEDS-90m_GreatLakes/slope_HydroSHEDS-90m_GreatLakes_aggregated_v1.1_RDRS-v2.nc <br>
  https://github.com/julemai/GRIP-E/blob/master/data/dem_HydroSHEDS-90m_GreatLakes/slope_HydroSHEDS-90m_GreatLakes_aggregated_v1.1_WFDEI-GEM-CaPA.nc
  

1.a Data (ASCII raster format)

    Conditioned SRTM DEM from HydroSHEDS, 3 second (90 m)    -->  file: rect_dem_GreatLakes.txt
    Dervided flow accumulation from DEM using ArcGIS         -->  file: rect_flow_accumulation_GreatLakes.txt
    Dervided flow direction    from DEM using ArcGIS         -->  file: rect_flow_direction_GreatLakes.txt

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

    Conditioned SRTM DEM from HydroSHEDS, 3 second (90 m)    -->  file: rect_dem_GreatLakes.nc
    Dervided flow accumulation from DEM using ArcGIS         -->  file: rect_flow_accumulation_GreatLakes.nc
    Dervided flow direction    from DEM using ArcGIS         -->  file: rect_flow_direction_GreatLakes.nc

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
	datapath='../../data/'
    varname="dem"
    vartype="single"
    unit="m"
    description="digital elevation model from HydroSheds (USGS) based on conditioned, global SRTM DEM at 3 sec (90m) resolution"

    python raster2netcdf.py -i ${datapath}dem_HydroSHEDS-90m_GreatLakes/rect_dem_GreatLakes.txt \
		     	        -o ${datapath}dem_HydroSHEDS-90m_GreatLakes/rect_dem_GreatLakes.nc \
		     	        -v "${varname},${vartype},${unit},${description}"

    -----------------
    Flow accumulation
    -----------------
	datapath='../../data/'
    varname="flow_acc"
    vartype="int32"
    unit="1"
    description="number of cells draining into this grid cell; derived from DEM from HydroSheds by USGS based on conditioned, global SRTM DEM at 3 sec (90m) resolution"

    python raster2netcdf.py -i ${datapath}dem_HydroSHEDS-90m_GreatLakes/rect_flow_accumulation_GreatLakes.txt
		     	        -o ${datapath}dem_HydroSHEDS-90m_GreatLakes/rect_flow_accumulation_GreatLakes.nc
		     	        -v "${varname},${vartype},${unit},${description}"

    -----------------
    Flow direction
    -----------------
	datapath='../../data/'
    varname="flow_dir"
    vartype="int32"
    unit="1"
    description="flow direction of grid cell; 1-east, 2-southeast, 4-south, 8-southwest, 16-west, 32-northwest, 64-north, 128-northeast; derived from DEM from HydroSheds by USGS based on conditioned, global SRTM DEM at 3 sec (90m) resolution"

    python raster2netcdf.py -i ${datapath}dem_HydroSHEDS-90m_GreatLakes/rect_flow_direction_GreatLakes.txt
		     	        -o ${datapath}dem_HydroSHEDS-90m_GreatLakes/rect_flow_direction_GreatLakes.nc
		     	        -v "${varname},${vartype},${unit},${description}"
