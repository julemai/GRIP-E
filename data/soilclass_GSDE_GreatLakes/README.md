## Harmonized World soil data (FAO), 30 second (1 km) 

### 1.a Data (ASCII raster format; classification information given in separate file)

    Global Soil Dataset for Earth System Modeling, 30 second (1 km)   --> file: rect_GSDE_soil_class.txt
    Classification information (formatting adjusted)
    --> file: USDA_soil_class_legends.csv

	The soil classes were identified from the GSDE-based soil texture
	data (percentages of silt, sand and clay) in each grid cell, and
	the classification was based on the USDA soil texture triangle. 

	The resolution of the processed soil class data is the same to the GSDE data set.

	The procedures to generate soil classes for grid cells are briefly
    described as follows:
	1. The GSDE data set provides soil texture data for eight soil
		layers, defined as
		0- 0.045, 0.045- 0.091, 0.091- 0.166, 0.166- 0.289, 0.289- 0.493, 0.493- 0.829, 0.829- 1.383
		and 1.383- 2.296 m. 
		Each layer of the data are clipped by a rectangular polygon
		that fully covered the Great Lakes domain.
	2.  Soil texture data of the eight layers are aggregated as one
        value at each grid cell for silt, sand, and clay,
        respectively. The averaging is based on weights of soil
        depths.
	3.  The USDA-based soil classes are identified based on the soil texture triangle. 
	Each grid cell corresponds to only one soil class.

	Dataset got updated March 2, 2020:
	Updates:

	(1) The extent of the domain has been enlarged to contain both the wfdei-gem-capa and the rdrs_version2 forcing grids.
	(2) An important bug has been fixed when aggregate the 8-layer
	    soil texture into one layer.
	    Soil class checked more consistent with the HWSD dataset after this modification.


### 1.b Downloaded from

    http://globalchange.bnu.edu.cn/research/soilw

### 1.c Prepared by

    Hongren Shen 
    Ph.D. Student
    Department of Civil and Environmental Engineering 
    University of Waterloo 
    200 University Avenue West 
    Waterloo, Ontario, Canada, N2L3G1 
    Office: E2-2359
    Phone: +1 (226)-899-3796 
    Email: hongren.shen@uwaterloo.ca, shrhongren@gmail.com

### 2.a Converted data (NetCDF format including classification information)

     Global Soil Dataset for Earth System Modeling, 30 second (1 km)   --> file: rect_GSDE_soil_class.nc

### 2.b Converted by

    Juliane Mai PhD
    Department Civil & Environmental Engineering
    University of Waterloo
    200 University Ave West / Waterloo ON, N2L 3G1 / Canada
    Building PHY / Room 3016
    Phone: +1 (519) 888-4567 x 30016
    Mobil: +1 (226) 505-5600
    http://www.civil.uwaterloo.ca/jmai/
    https://www.researchgate.net/profile/Juliane_Mai

### 2.c Data conversion

    -----------------
    Soil class
    -----------------
	datapath='../../data/'
    varname="soil_class"
    vartype="int32"
    unit="1"
    description="soil classes derived from Global Soil Dataset for Earth System Modeling (GSDE) at 30 sec (1km) resolution"

    python raster2netcdf.py -i ${datapath}soilclass_GSDE/rect_GSDE_soil_class.txt
			    -o ${datapath}soilclass_GSDE/rect_GSDE_soil_class.nc
			    -v "${varname},${vartype},${unit},${description}"
			    -a ${datapath}soilclass_GSDE/USDA_soil_class_legends.csv

### 3 Aggregated products to RDRS-v2 and WFDEI-GEM-CaPA grids

    file: soilclass_GSDE_GreatLakes_aggregated.zip
	
    Mapping the soil texture data to forcing grids:

	(1) Two forcing grids are applied: 
		a. The wfdei_gem_capa forcing grids, 0.125 degree (about 12 km).
		b. The RDRS forcing grids, 10 km.
		
	(2) Soil texture (i.e. silt, sand and clay in percent) data are aggregated from the original 1 km grid cell to the forcing grids with respect to 8 soil layers. 

	(3) Variable definitions:
		FID: Order number of grid cells.
		FGID: Index of grid cells in the forcing netCDF file (in the two-dimensional matrix).
		Row: Row number of grid cells in the forcing netCDF file.
		Col: Column number of the grid cells in the forcing netCDF file.
		Gridlat: Lattitude of the grid central point.
		Gridlon: Longitude of the grid central point.
		Area: Area of the grid cells (km2).
		SILT_L1 to SILT_L8: Percentage of silt in soil layer 1 to 8. -9999 denotes NaN values from the original GSDE data set.
		SAND_L1 to SAND_L8: Percentage of sand in soil layer 1 to 8. -9999 denotes NaN values from the original GSDE data set.
		CLAY_L1 to CLAY_L8: Percentage of clay in soil layer 1 to 8. -9999 denotes NaN values from the original GSDE data set.
		BD_L1 to BD_L8: Bulk density (0.01 g/cm3) of soil layer 1 to 8. -9999 denotes NaN values from the original GSDE data set.
		USDA_class_L1 to USDA_class_L8: The USDA-based soil class identified from soil layer 1 to 8.
		USDA_class_8layers: The USDA-based soil class identified from the full soil column.
			
	(4) Note
		When checking the aggregated soil texture results, it was found that there were a large part of lake-grid cells having 
		normal soil texture values instead of -9999, which means these lake-grid cells were identified as soil grids. This is 
		caused by the nature of GSDE data set.

		There are some grid cells with very small area. These grids are situated at the four edges of the netCDF domain.
		Since there are no more central points beyond these edges, we cannot get FULL grid polygon for those grid cells. 
		However, these grid cells are chosen as buffer, which will not be used in the hydrological modeling. 
		Thus, these small grids may look weird when checking the grid area info, but will not affect our modeling work.
		
		
	If you have any questions about the soil texture mapping results, please contact Hongren Shen: hongren.shen@uwaterloo.ca
		
	2020-03-01

	
