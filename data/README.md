## Datasets used for project

### Geophysical Datasets 
The geophysical datasets (DEM, land cover, soil) shared here were used by some models in the project. At the end of the project we agreed on a common geophysical dataset indicated under "Common dataset" in the Wiki [here](https://github.com/julemai/GRIP-E/wiki/Datasets).

## Meteorologic forcings 
The meteorologic forcings shared here are used by all models in the GRIP-E project. They can be downloaded using the external download links. The fully gridded dataset is [here](http://www.civil.uwaterloo.ca/jmai/GRIP-E/RDRS_CaPA24hr_forcings_final.zip) (1.1GB). The forcings aggregated to basin level are [here](http://www.civil.uwaterloo.ca/jmai/GRIP-E/RDRS_CaPA24hr_forcings_final_aggregated.zip) (1.4 MB).

## Observed streamflow
The data for observed streamflow are available in raw csv/txt format (US stations are in ft<sup>3</sup>/s; CA stations in m<sup>3</sup>/s) and in harmonized NetCDF format (all stations in m<sup>3</sup>/s). The observations are split into stations ob objective 1 (low-human impact) and objective 2 (most downstream). There are stations that appear in both objcetives (low-human impact and most downstream). 

* Observed streamflow of objective 1 stations in calibration: [here](https://github.com/julemai/GRIP-E/tree/master/data/objective_1/lake-erie/calibration) (list of stations [here](https://github.com/julemai/GRIP-E/blob/master/data/objective_1/lake-erie/calibration/gauge_info.csv))
* Observed streamflow of objective 1 stations in validation: [here](https://github.com/julemai/GRIP-E/tree/master/data/objective_1/lake-erie/validation) (list of stations [here](https://github.com/julemai/GRIP-E/blob/master/data/objective_1/lake-erie/validation/gauge_info.csv))
* Observed streamflow of objective 2 stations in calibration: [here](https://github.com/julemai/GRIP-E/tree/master/data/objective_2/lake-erie/calibration) (list of stations [here](https://github.com/julemai/GRIP-E/blob/master/data/objective_2/lake-erie/calibration/gauge_info.csv))
* Observed streamflow of objective 2 stations in validation: [here](https://github.com/julemai/GRIP-E/tree/master/data/objective_2/lake-erie/validation) (list of stations [here](https://github.com/julemai/GRIP-E/blob/master/data/objective_2/lake-erie/validation/gauge_info.csv))

## Model results
The model results are available in native model output format and harmonized NetCDF format. The latter is the same for each file and certainly preferable for post-processing. The results are in folders `data/objective_XXX/lake-erie/YYY/model/ZZZ` where `XXX in {1,2}`, `YYY in {calibration, validation}`, and `ZZZ` is the model name. The netCDF file names follow the pattern `ZZZ_phase_WWW_objective_XXX.nc` where `ZZZ` is the model name (lowercase), `WWW in {0,1}` (phase of project; 0=pre-calibration/ default run; 1=after calibration), and `XXX in {1,2}` (objective; 1=low-human impact; 2=most-downstream gauges).

Example:
The results of the HYPE model for objective 1 for the calibration stations after calibration are in the file
`data/objective_1/lake-erie/calibration/model/HYPE/hype_phase_1_objective_1.nc` ([here](https://github.com/julemai/GRIP-E/blob/master/data/objective_1/lake-erie/calibration/model/HYPE/hype_phase_1_objective_1.nc)).
