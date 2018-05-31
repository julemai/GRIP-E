--------------------------------------------------------------------------------------------------------
data (ASCII raster format; classification information given in separate file)

    MODIS/Terra+Aqua Land Cover Type Yearly L3 Global 500m SIN Grid V006 (MCD12Q1_006), 500 m   --> file: rect_landcover_UMD_scheme_MODIS_Erie.txt
    Classification information (formatting adjusted)                                            --> file: UMD_scheme_MODIS_legends.csv

downloaded from

    https://lpdaac.usgs.gov/dataset_discovery/modis/modis_products_table/mcd12q1_v006

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

    MODIS/Terra+Aqua Land Cover Type Yearly L3 Global 500m SIN Grid V006 (MCD12Q1_006), 500 m   --> file: rect_landcover_UMD_scheme_MODIS_Erie.nc

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
    Land cover
    -----------------
    varname="land_cover"
    vartype="int32"
    unit="1"
    description="land cover classes taken from Global 500m MODIS MCD12Q1 \
                 product (NASA) using classification scheme 2 (UMD)"

    python raster2netcdf.py -i ${datapath}landcover_MODIS/rect_landcover_UMD_scheme_MODIS_Erie.txt
                            -o ${datapath}landcover_MODIS/rect_landcover_UMD_scheme_MODIS_Erie.nc
			    -v "${varname},${vartype},${unit},${description}"
			    -a ${datapath}landcover_MODIS/UMD_scheme_MODIS_legends.csv