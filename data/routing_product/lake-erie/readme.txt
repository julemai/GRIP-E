GRIP-E spatial validation gauges are selected from the GRIP-GL data sets. Details of the gauge selection procedures can be referred to the April 2020 teleconference slides.

#####################################################################################################
2020-05-25 Updates:

GRIP-E spatial validation gauge list is updated based on the May 2020 monthly tele-conference.

(1) Six gauges from the GRIP-GL data set are selected. Routing network of these six gauges are extracted from the GRIP-GL routing product v1.2.

(2) Two new gauges from the USGS are added into the validation list. Routing network of these two gauges are newly delineated based on the GRIP-GL datas sets. This new routing network is independent to the GRIP-GL network, but they differ in the catchment numbers.

Since these eight catchments and their catchments are archived seperately, it will not affect hydrological modelling because of the different network system (catchment numbers). 

#####################################################################################################
2020-04-18 Updates:

GRIP-E spatial calidation catchment data include:

(1) [GRIP_GL_for_Lake_Erie_spatial_validation_gauge_info.csv] This is a gauge info file, containing the basic information for the selected gauges. This file comes from GRIP-GL gauge info file.

(2) [GRIP-E_spatial_validation_routing_network] This folder contains the three versions of routing network cooresponding to the GRIP-GL routing network v1.1.1.
In this routing network, the catchment polygons and the catchment-forcing-overlaying polygons for each gauge are seperated.
The RDRS forcing used in GRIP-E is used to generated this routing network. To distinguish this RDRS with the upcoming new RDRS, we hereafter use "RDRS_v1" to indicate the older version.
Variable names and their meanings in the attribute tables can be referred to GRIP-GL routing network.

(3) [GRIP-E_spatial_validation_RDRS_v1_24h_LocalTime_2010-2014] This folder contains the aggregated RDRS_v1 data at catchment level.
Data are archived in ASCII format. It can be directly used in Raven.
Data time step: 24-hr. Data range: 2000-01-01 to 2014-12-31, 1826 steps in total.
Variables aggregated are: 
    Cumulated daily precipitaion (PRECIP), mm/day.
    Daily minimum air temperature (TEMP_MIN), Celsius degree.
    Daily maximum air temperature (TEMP_MAX), Celsius degree.
    Daily mean air temperature (TEMP_DAILY_AVE), Celsius degree.

Please feel free to let us know if there is any question.

Hongren Shen 

hongren.shen@uwaterloo.calidation