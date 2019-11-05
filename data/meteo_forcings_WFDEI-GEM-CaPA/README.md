## WFDEI-GEM-CaPA dataset

External download from:
http://www.civil.uwaterloo.ca/jmai/GRIP-GL/grip-gl_wfdei-gem-capa_gridonly.nc                 (236 KB)
http://www.civil.uwaterloo.ca/jmai/GRIP-GL/grip-gl_wfdei-gem-capa_2000-2016_leap.nc           (8.4 GB)
http://www.civil.uwaterloo.ca/jmai/GRIP-GL/grip-gl_wfdei-gem-capa_2000-2016_leap/2000_leap.nc (470 MB)
http://www.civil.uwaterloo.ca/jmai/GRIP-GL/grip-gl_wfdei-gem-capa_2000-2016_leap/2001_leap.nc (470 MB)
...
http://www.civil.uwaterloo.ca/jmai/GRIP-GL/grip-gl_wfdei-gem-capa_2000-2016_leap/2016_leap.nc (470 MB)



1.a Data

	Data provided by Mohamed Elshamy.
	3 hourly, 10~km regular grid, 1979-2016, North America
	UTC time $\curvearrowright$ need to be shifted by 6h to get "local" time
	Variables available:
	- hus: 				 Specific Humidity at Lowest Model Level (approx 40 m)
	- pr:  				 Precipitation (surface)
	- ps:  				 Surface Pressure (surface)
	- rlds:				 Surface Downwelling Longwave Radiation
	- rsds:				 Surface Downwelling Shortwave Radiation
	- rsds_thresholded:  Surface Downwelling Shortwave Radiation (rsds) variable post-processed setting values to zeros if they were zero in WFDEI
	- wind_speed:        Zonal (Eastward) Wind (approx 40 m)
	- ta:                Air Temperature (approx 40 m)

1.b Downloaded from

	Cuizinart

1.c Prepared by

    Juliane Mai PhD
    Department Civil & Environmental Engineering
    University of Waterloo
    200 University Ave West / Waterloo ON, N2L 3G1 / Canada
    Building PHY / Room 3016
    Phone: +1 (519) 888-4567 x 30016
    Mobil: +1 (226) 505-5600
    http://www.civil.uwaterloo.ca/jmai/
    https://www.researchgate.net/profile/Juliane_Mai

1.d Leap year data 

    python add_leap_days.py -i 2001.nc -o 2001_leap.nc
    ...

	--> scripts/add_leap/days/process.sh
