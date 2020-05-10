## Global Soil Dataset for Earth System Modeling (GSDE), 30 second (1 km) 

External Download links:
* Soil classes (1km) in ASCII raster format <br>
  https://github.com/julemai/GRIP-E/blob/master/data/soilclass_GSDE_GreatLakes/rect_GSDE_soil_class.txt.zip
* Soil classes (1km) in NetCDF format <br>
  https://github.com/julemai/GRIP-E/blob/master/data/soilclass_GSDE_GreatLakes/rect_GSDE_soil_class.nc
* Soil textures (sand, silt, clay, bulk density; 1km) in ASCII raster format <br>
  https://github.com/julemai/GRIP-E/blob/master/data/soilclass_GSDE_GreatLakes/rect_GSDE_soil_texture.zip
* Soil classses and textures aggregated to forcing grids in tabulated TXT format <br>
  https://github.com/julemai/GRIP-E/blob/master/data/soilclass_GSDE_GreatLakes/soilclass_GSDE_GreatLakes_aggregated_v1.1.zip
* Soil classses and textures aggregated to forcing grids in NetCDF format <br>
  https://github.com/julemai/GRIP-E/blob/master/data/soilclass_GSDE_GreatLakes/soilclass_GSDE_GreatLakes_aggregated_v1.1_RDRS-v2.nc <br>
  https://github.com/julemai/GRIP-E/blob/master/data/soilclass_GSDE_GreatLakes/soilclass_GSDE_GreatLakes_aggregated_v1.1_WFDEI-GEM-CaPA.nc
  

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

    python raster2netcdf.py -i ${datapath}soilclass_GSDE_GreatLakes/rect_GSDE_soil_class.txt
                -o ${datapath}soilclass_GSDE_GreatLakes/rect_GSDE_soil_class.nc
                -v "${varname},${vartype},${unit},${description}"
                -a ${datapath}soilclass_GSDE_GreatLakes/USDA_soil_class_legends.csv

    -----------------
    Soil texture (sand)
    -----------------
    datapath='../../data/'
    varname="sand"
    vartype="float"
    unit="percent"
    description="Sand content as percent of weight "

    python raster2netcdf.py -i ${datapath}soilclass_GSDE_GreatLakes/rect_GSDE_soil_texture/rect_GSDE_soil_texture_sand_8_layers_avg.txt
                -o ${datapath}soilclass_GSDE_GreatLakes/rect_GSDE_soil_texture/rect_GSDE_soil_texture_sand_8_layers_avg.nc
                -v "${varname},${vartype},${unit},${description}"

    -----------------
    Soil texture (silt)
    -----------------
    datapath='../../data/'
    varname="silt"
    vartype="float"
    unit="percent"
    description="Silt content as percent of weight "

    python raster2netcdf.py -i ${datapath}soilclass_GSDE_GreatLakes/rect_GSDE_soil_texture/rect_GSDE_soil_texture_silt_8_layers_avg.txt
                -o ${datapath}soilclass_GSDE_GreatLakes/rect_GSDE_soil_texture/rect_GSDE_soil_texture_silt_8_layers_avg.nc
                -v "${varname},${vartype},${unit},${description}"

    -----------------
    Soil texture (clay)
    -----------------
    datapath='../../data/'
    varname="clay"
    vartype="float"
    unit="percent"
    description="Clay content as percent of weight "

    python raster2netcdf.py -i ${datapath}soilclass_GSDE_GreatLakes/rect_GSDE_soil_texture/rect_GSDE_soil_texture_clay_8_layers_avg.txt
                -o ${datapath}soilclass_GSDE_GreatLakes/rect_GSDE_soil_texture/rect_GSDE_soil_texture_clay_8_layers_avg.nc
                -v "${varname},${vartype},${unit},${description}"

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


### 4 History of versions

    Mapping the soil texture data to forcing grids:

    #### 2020-05-07 Updates:

    In this version, the GSDE 8-layer soil texture aggregation results at the 30
    sec (1km, GSDE original grids) are updated.

    (1) A bug is fixed related to the soil class identification:

        This bug may occur in lake part where the soil textures are
        -9999 and thus, soil class are set as -9999 as well.

    Files revised/created:
    - rect_GSDE_soil_class_v1.1.txt/nc
        
    (2) In this version, soil textures (silt, clay and sand) aggregation
    results at the 30 sec resolution are added.

    Files revised/created:
    - rect_GSDE_soil_texture.zip --> rect_GSDE_soil_texture_clay_8_layers_avg.txt/nc
    - rect_GSDE_soil_texture.zip --> rect_GSDE_soil_texture_sand_8_layers_avg.txt/nc
    - rect_GSDE_soil_texture.zip --> rect_GSDE_soil_texture_silt_8_layers_avg.txt/nc

    (3) We would like to clarify the 8-layer aggregation weights we used
    in the processing. These weights are computed based on the GSDE soil
    depths of each layer: 
        
        layer 1: 0-4.5 cm
        layer 2: 4.5-9.1 cm
        layer 3: 9.1-16.6 cm
        layer 4: 16.6-28.9 cm
        layer 5: 28.9-49.3 cm
        layer 6: 49.3-82.9 cm
        layer 7: 82.9-138.3 cm
        layer 8: 138.3-229.6 cm

        The soil depths for layer1 to layer8 are 4.5, 4.6, 7.5, 12.3,
        20.4, 33.6, 55.4, and 91.3 cm, respectively. 

        Thus, we can compute the weights of each soil layer by
        dividing each layer depth by the total depth, yielding ratios
        0.0195993 , 0.02003484, 0.03266551, 0.05357143, 0.08885017,
        0.14634146, 0.2412892 , and 0.39764808 respectively.
        
    #### 2020-04-21 Updates:

    In v1.2 products, we optimized the output soil texture and class format:

    In grid-cells whose original GSDE soil textures are missing, the
    aggregated soil texture data are set as -9999. 

    In soil layers where at least one of the soil textures is missing, the
    coorresponding soil class is set as -9999. 

    Please note that in the GSDE metadata, there are grid-cells that the
    sum of their soil texture percentages (silt, sand and clay) may not be
    equal to 100%. Therefore, a compulsory re-scale is adopted in this
    product to ensure the sum is 100%. 

    #### 2020-04-10 Updates:

    To distinguish different versions of the soil mapping products, we
    started to name the latest (current) product as v1.1. The previous
    (untill Mar 12, 2020) product is named as v1.0. 

    In v1.1 products, two bugs has been fixed:

    (1) We noticed that in the RDRS forcing file, lattitude and longitude
    values are defined by a float number with 6 digits. However, in v1.0
    product the lat/lon values of some grid  
    cells were exported with 5-digit precision. This issue may cause
    problem when users want to locate those grid cells in the forcing
    domain. 6-digit values have been used for lat/lon 
    since v1.1. 

    (2) A bug related to bulk density and soil texture spatial-join has
    been fixed. 

    In addition, the exported ASCII file format was optimized.

    After fixing the bugs, both the [wfdei_gem_capa] and [RDRS] forcing
    grids based soil aggregation were re-run to keep everything else
    consistent.  

    Please use the latest version in your modelling. 

    #### 2020-03-12 Updates:

    We noticed that using the strict soil texture aggregation algrithom,
    i.e. considering bulk density when upscaling soil texture data, may
    cause an issue in some  
    grids whose bulk density values are missing ([NaN]) but soil texture
    values are normal. This algrithom, however, identifies the aggregated
    soil texture data in those grids  
    as [NaN] too, and this issue may impact the parameterization in
    distributed hydrological models. Thus, we optimized the aggregation in
    those grids that we use the  
    weighted texture data without bulk density and then scale the three
    texture values to ensure their sum is equal to 1.0. 

    #### 2020-03-01:

    (1) Two forcing grids are applied: 
        a. The wfdei_gem_capa forcing grids, 0.125 degree (about 12 km).
        b. The RDRS forcing grids, 10 km.
        
    (2) Soil texture (i.e. silt, sand and clay in percent) data are
    aggregated from the original 1 km grid cell to the forcing grids with
    respect to 8 soil layers.  

    (3) Variable definitions:
        FID:             Order number of grid cells.
        FGID:            Index of grid cells in the forcing netCDF
                             file (in the two-dimensional matrix). 
        Row:             Row number of grid cells in the forcing
                             netCDF file. 
        Col:             Column number of the grid cells in the
                     forcing netCDF file. 
        Gridlat:         Lattitude of the grid central point.
        Gridlon:         Longitude of the grid central point.
        Area:            Area of the grid cells (km2).
        SILT_L1 to SILT_L8:  Percentage of silt in soil layer 1 to
                             8. -9999 denotes NaN values from the
                             original GSDE data set.  
        SAND_L1 to SAND_L8:  Percentage of sand in soil layer 1 to
                             8. -9999 denotes NaN values from the
                             original GSDE data set.  
        CLAY_L1 to CLAY_L8:  Percentage of clay in soil layer 1 to
                             8. -9999 denotes NaN values from the
                             original GSDE data set.  
        BD_L1 to BD_L8:      Bulk density (0.01 g/cm3) of soil layer 1 to
                             8. -9999 denotes NaN values from the
                             original GSDE data set.  
        USDA_class_L1 to USDA_class_L8: The USDA-based soil class
                             identified from soil layer 1 to 8. 
        USDA_class_8layers:  The USDA-based soil class identified from
                             the full soil column. 
            
    (4) Note
        When checking the aggregated soil texture results, it was
        found that there were a large part of lake-grid cells having  
        normal soil texture values instead of -9999, which means these
        lake-grid cells were identified as soil grids. This is  
        caused by the nature of GSDE data set.

        There are some grid cells with very small area. These grids
        are situated at the four edges of the netCDF domain. 
        Since there are no more central points beyond these edges, we
        cannot get FULL grid polygon for those grid cells.  
        However, these grid cells are chosen as buffer, which will not
        be used in the hydrological modeling.  
        Thus, these small grids may look weird when checking the grid
        area info, but will not affect our modeling work. 
        
        
    If you have any questions about the soil texture mapping results,
    please contact Hongren Shen: hongren.shen@uwaterloo.ca 
    

