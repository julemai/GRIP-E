## Monthly Posterior Distributions 

### 1.a Data delivery

Data were received by email on May 30, 2018 from 

    Joeseph Paul Smith
    B. Sci (2012) Informatics at University of Michigan - Data Mining and Information Analysis
    General Programmer/Analyst
    Cooperative Institute for Great Lakes Research (CIGLR)
    at NOAA Great Lakes Environmental Research Laboratory (NOAA-GLERL)
    4840 S. State Rd. Ann Arbor, MI 48108
    +1 (734) 741-2252
    joeseph@umich.edu or joeseph.smith@noaa.gov

### 1.b Content of Files

#### Summary statistics of monthly posterior distributions

In units of cubic meters per second (m^3/s)
- Lake St. Clair inflow (= miHuronOutflow)         <br> --> file: inflow_lake_st-clair_summary-stats.csv
- Lake St. Clair net basin supply                  <br> --> file: nbs_lake_st-clair_summary-stats.csv
- Lake Erie inflow (= clairOutflow)                <br> --> file: inflow_lake_erie_summary-stats.csv
- Lake Erie net basin supply                       <br> --> file: nbs_lake_erie_summary-stats.csv

These are CSV files, with the columns from left to right being:
- Year
- Month
- Median of posterior distribution
-  2.5 percentile
- 97.5 percentile

#### All ensemble members of monthly posterior distributions

In units of cubic meters per second (m^3/s)
- Lake St. Clair inflow (= miHuronOutflow)         <br> --> file: inflow_lake_st-clair_samples.csv
- Lake St. Clair net basin supply                  <br> --> file: nbs_lake_st-clair_samples.csv
- Lake Erie inflow (= clairOutflow)                <br> --> file: inflow_lake_erie_samples.csv
- Lake Erie net basin supply                       <br> --> file: nbs_lake_erie_samples.csv

These are CSV files with the T=792 columns and N=3000 rows. Each row
is one ensemble member, each column is one month starting at Jan 1950.

### 1.c Period of data
- Jan 1950 to Dec 2015
- monthly data


### 1.d Notes

The net basin supply for Lake Erie was initially provided in [mm] instead of
[m^3/s]. The conversion was performed by Juliane Mai in collaboration
with Joeseph Smith.

`` new = old * multiplier ``

where

`` multiplier = 1./(seconds_this_month)*lake_area_m2/1000``

with the lake area (lake_area_m2) was set to be 27314*1000^2
[m^2]. The seconds per month (seconds_this_month) was determined for
every month considering leap years and different lengths of
months. The multiplier for a 30-day month (e.g. April) would be around 10.5378.

The original files are named "nbs_lake_erie_summary-stats_mm.csv" and "nbs_lake_erie_samples_mm.csv". The conversion was done using this script [here](https://github.com/julemai/GRIP-E/blob/master/scripts/monthly_nbs_inflows/convert_nbs_lake_erie_to_m3s.py). 
