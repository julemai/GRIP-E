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

--> python extract_validation_results.py -m SWAT-Guelph -i ../../data/objective_1/lake-erie/validation/model/SWAT-Guelph/swat-guelph_phase_1_objective_1.csv -o ../../data/objective_1/lake-erie/validation/model/SWAT-Guelph/swat-guelph_phase_1_objective_1.nc -a ../../data/objective_1/lake-erie/validation/gauge_info.csv -b ../../data/objective_1/lake-erie/validation/model/SWAT-Guelph/subid2gauge.csv

--> python extract_validation_results.py -m SWAT-Guelph -i ../../data/objective_2/lake-erie/validation/model/SWAT-Guelph/swat-guelph_phase_1_objective_2.csv -o ../../data/objective_2/lake-erie/validation/model/SWAT-Guelph/swat-guelph_phase_1_objective_2.nc -a ../../data/objective_2/lake-erie/validation/gauge_info.csv -b ../../data/objective_2/lake-erie/validation/model/SWAT-Guelph/subid2gauge.csv


Use SWAT reach 442 with weight 0.9583 to derive streamflow at gauge 04185000
Use SWAT reach 471 with weight 1.0058 to derive streamflow at gauge 04195500
Use SWAT reach 558 with weight 1.0000 to derive streamflow at gauge 04201500
Use SWAT reach 521 with weight 1.0000 to derive streamflow at gauge 04201500
Use SWAT reach 505 with weight 1.0059 to derive streamflow at gauge 04208000
Use SWAT reach 235 with weight 1.0001 to derive streamflow at gauge 02GE003
Use SWAT reach 309 with weight 1.0690 to derive streamflow at gauge 04168000
Use SWAT reach 295 with weight 1.1043 to derive streamflow at gauge 04167000
