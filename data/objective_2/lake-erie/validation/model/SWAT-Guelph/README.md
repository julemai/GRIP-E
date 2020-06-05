## Validation procedure

1. validation gauges that are located at sub-basin outlet:
- area-ratio method
- will be applied to stations: 04185000, 04195500, 04208000
 
2. validation gauges that are NOT located at sub-basin outlet:
- area-ratio method
- will be applied to stations: 04201500, 02GE003, 04167000, 04168000


Narayan email
-------------------
We don't need to re-run the model. You can extract the data from the
OUTPUT.RCH files that I sent for Objectives 1 &2 (I hope you still
have them!). To extract, you need the corresponding reach number (or
sub-basin number) and a factor (for area-weightage
calculation). Please find the attached excel sheet which contains
information regarding this. For the gauging station (04201500, ROCKY R
NR BEREA OH), as it is exactly located at the confluence of two
reaches, you have to add FLOW_OUT for both these reaches. I have
explicitly mentioned this in the attached excel sheet (information
below).

--> done in script: "extract_validation_results.py"



FID	ID	Name	Lat	Lon	lat_snap_90m	lon_snap_90m	Area	Total delineated area at associated sub-basin outlet	Delineated area at gauging location	Factor	Associated Subbasin	Comment
0	4185000	Tiffin River at Stryker OH	41.50449567	-84.4296719	41.4920998	-84.441803	1061.9	1011.533494	969.3743551	0.9583	442	Take FLOWOUT of reach 442 and multiply by the factor
1	4195500	PORTAGE R AT WOODVILLE OH	41.449496	-83.3613156	41.4522018	-83.3509979	1108.52	1127.804036	1121.264384	1.0058	471	Take FLOWOUT of reach 471 and multiply by the factor
2	4201500	ROCKY R NR BEREA OH	41.4075483	-81.88263759	41.4188004	-81.8684998	691.53	683.1983687	683.1983688	1.0000	558	This station is right at the confluence of 558 and 521. Please extract the FLOW_OUT of these two reaches and add them up.
3	4208000	CUYAHOGA R AT INDEPENDENCE OH	41.39533087	-81.6298478	41.3653984	-81.6092987	1831.13	1719.614609	1709.470481	1.0059	505	Take FLOWOUT of reach 505 and multiply by the factor
4	02GE003	THAMES RIVER AT THAMESVILLE	42.54486084	-81.9672699	42.5429001	-81.985199	4370	4330.817975	4330.39966	1.0001	235	Take FLOWOUT of reach 235 and multiply by the factor
5	4168000	LOWER RIVER ROUGE AT INKSTER MI	42.300593	-83.300208	42.28289576	-83.38604706	219.2	293.329166	274.4004306	1.0690	309	Take FLOWOUT of reach 309 and multiply by the factor
6	4167000	MIDDLE RIVER ROUGE NEAR GARDEN CITY MI	42.348093	-83.311597	42.33136272	-83.24857193	229.5	195.2274214	176.78672	1.1043	295	Take FLOWOUT of reach 295 and multiply by the factor
