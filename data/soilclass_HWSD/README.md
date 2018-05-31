--------------------------------------------------------------------------------------------------------
data (ASCII raster format; classification information given in separate file)

    Harmonized World soil data (FAO), 30 second (1 km)                --> file: rect_landcover_UMD_scheme_MODIS_Erie.txt
    Classification information (formatting adjusted)                  --> file: USDA_soil_class_legends.csv

downloaded from

    http://www.fao.org/soils-portal/soil-survey/soil-maps-and-databases/harmonized-world-soil-database-v12/en/

prepared by

    Hongren Shen 
    Ph.D. Student
    Department of Civil and Environmental Engineering 
    University of Waterloo 
    200 University Avenue West 
    Waterloo, Ontario, Canada, N2L3G1 
    Office: E2-2359
    Phone: +1 (226)-899-3796 
    Email: hongren.shen@uwaterloo.ca, shrhongren@gmail.com

--------------------------------------------------------------------------------------------------------
converted data (NetCDF format including classification information)

    Harmonized World soil data (FAO), 30 second (1 km)                --> file: rect_landcover_UMD_scheme_MODIS_Erie.nc

converted by

    Juliane Mai PhD
    Department Civil & Environmental Engineering
    University of Waterloo
    200 University Ave West / Waterloo ON, N2L 3G1 / Canada
    Building PHY / Room 3016
    Phone: +1 (519) 888-4567 x 30016
    Mobil: +1 (226) 505-5600
    http://www.civil.uwaterloo.ca/jmai/
    https://www.researchgate.net/profile/Juliane_Mai

data conversion

    -----------------
    Soil class
    -----------------
    varname="soil_class"
    vartype="int32"
    unit="1"
    description="soil classes taken from FAO Harmonized World Soil Database (HWSD) v1.2 at 30 sec (1km) resolution"

    python raster2netcdf.py -i ${datapath}soilclass_HWSD/rect_HWSD_soil_class.txt
                            -o ${datapath}soilclass_HWSD/rect_HWSD_soil_class.nc
			    -v "${varname},${vartype},${unit},${description}"
			    -a ${datapath}soilclass_HWSD/USDA_soil_class_legends.csv