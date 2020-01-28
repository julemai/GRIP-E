## NALCMS, 19 land cover classes for North America (30m, Landsat, 2010 from Mexico and Canada, 2011 for U.S.)

External download from:
http://www.civil.uwaterloo.ca/jmai/GRIP-GL/landcover_NALCMS_GreatLakes.txt.zip (ASCII raster format; 345.5 MB)
http://www.civil.uwaterloo.ca/jmai/GRIP-GL/landcover_NALCMS_GreatLakes.tif.zip (TIF format incl projection information; 370.9 MB)

1.a Data (ASCII raster format; classification information given in separate file)

    NALCMS, 19 land cover classes for North America (30m, Landsat, 2010 from Mexico and Canada, 2011 for U.S.)   --> file: rect_landcover_NACLMS_GreatLakes.txt
    Classification information (formatting adjusted)                                            --> file: NACLMS_legends.csv

1.b Downloaded from

	https://www.mrlc.gov/data/north-american-land-change-monitoring-system

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

2.a Converted data (NetCDF format including classification information)

    NOT converted since lat/lon are in [m] and not [degree]

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

2.c Data conversion into NetCDF

	NOT CONVERTED!!!!

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

2.d Data conversion into TIFF
