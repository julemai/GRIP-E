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
