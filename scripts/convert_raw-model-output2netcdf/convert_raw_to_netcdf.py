#!/usr/bin/env python
from __future__ import print_function

# Copyright 2016-2018 Juliane Mai - juliane.mai(at)uwaterloo.ca
#
# License
# This file is part of Juliane Mai's personal code library.
#
# Juliane Mai's personal code library is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Juliane Mai's personal code library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.

# You should have received a copy of the GNU Lesser General Public License
# along with Juliane Mai's personal code library.  If not, see <http://www.gnu.org/licenses/>.
#

# run:
#
#    ------------
#    LBRM
#    ------------
#    python convert_raw_to_netcdf.py -m LBRM -i ../../data/objective_1/lake-erie/calibration/model/LBRM/lbrm_phase_0_objective_1.csv -o ../../data/objective_1/lake-erie/calibration/model/LBRM/lbrm_phase_0_objective_1.nc -a ../../data/objective_1/lake-erie/calibration/gauge_info.csv
#
#    ------------
#    VIC
#    ------------
#    python convert_raw_to_netcdf.py -m VIC -i ../../data/objective_1/lake-erie/calibration/model/VIC/vic_phase_0_objective_1.csv -o ../../data/objective_1/lake-erie/calibration/model/VIC/vic_phase_0_objective_1.nc -a ../../data/objective_1/lake-erie/calibration/gauge_info.csv -b ../../data/objective_1/lake-erie/calibration/model/VIC/subid2gauge.csv

#    ------------
#    VIC-GRU
#    ------------
#    python convert_raw_to_netcdf.py -m VIC-GRU -i ../../data/objective_1/lake-erie/calibration/model/VIC-GRU/vic-gru_phase_0_objective_1.csv -o ../../data/objective_1/lake-erie/calibration/model/VIC-GRU/vic-gru_phase_0_objective_1.nc -a ../../data/objective_1/lake-erie/calibration/gauge_info.csv -b ../../data/objective_1/lake-erie/calibration/model/VIC-GRU/subid2gauge.csv

#    ------------
#    GEM-Hydro
#    ------------
#    python convert_raw_to_netcdf.py -m GEM-Hydro -i ../../data/objective_1/lake-erie/calibration/model/GEM-Hydro/gem-hydro_phase_0_objective_1.csv -o ../../data/objective_1/lake-erie/calibration/model/GEM-Hydro/gem-hydro_phase_0_objective_1.nc -a ../../data/objective_1/lake-erie/calibration/gauge_info.csv

#    ------------
#    -- Lake Erie   :: 'ML-ConvLSTM' or 'ML-ConvLSTM-DEM' or 'ML-ConvLSTM-LC' or 'ML-ConvLSTM-LC-DEM' or 'ML-LinReg' or 'ML-XGBoost'    
#    -- Great Lakes :: 'ML-EA-LSTM' or 'ML-LSTM'
#    ------------
#    python convert_raw_to_netcdf.py -m ML-ConvLSTM-w-LC -i ../../data/objective_1/lake-erie/calibration/model/ML-ConvLSTM-w-LC/ml-convlstm-w-lc_phase_1_objective_1.csv -o ../../data/objective_1/lake-erie/calibration/model/ML-ConvLSTM-w-LC/ml-convlstm-w-lc_phase_1_objective_1.nc -a ../../data/objective_1/lake-erie/calibration/gauge_info.csv

#    ------------
#    HYPE
#    ------------
#    python convert_raw_to_netcdf.py -m HYPE -i ../../data/objective_1/lake-erie/calibration/model/HYPE/hype_phase_0_objective_1_ -o ../../data/objective_1/lake-erie/calibration/model/HYPE/hype_phase_0_objective_1.nc -a ../../data/objective_1/lake-erie/calibration/gauge_info.csv

#    ------------
#    HYMOD2-DS
#    ------------
#    python convert_raw_to_netcdf.py -m HYMOD2-DS -i ../../data/objective_1/lake-erie/calibration/model/HYMOD2-DS/hymod_phase_0_objective_1_ -o ../../data/objective_1/lake-erie/calibration/model/HYMOD2-DS/hymod_phase_0_objective_1.nc -a ../../data/objective_1/lake-erie/calibration/gauge_info.csv

#    ------------
#    HMETS-Raven-lp
#    ------------
#    python convert_raw_to_netcdf.py -m HMETS-Raven-lp -i ../../data/objective_1/lake-erie/calibration/model/HMETS-Raven-lp/raven-hmets-lp_phase_1_objective_1_ -o ../../data/objective_1/lake-erie/calibration/model/HMETS-Raven-lp/raven-hmets-lp_phase_0_objective_1.nc -a ../../data/objective_1/lake-erie/calibration/gauge_info.csv -b ../../data/objective_1/lake-erie/calibration/model/HMETS-Raven-lp/subid2gauge.csv -s julie

#    ------------
#    Raven-blended
#    ------------
#    python convert_raw_to_netcdf.py -m Raven-blended -i ../../data/objective_1/lake-erie/calibration/model/Raven-blended/raven-blended_phase_1_objective_1_ -o ../../data/objective_1/lake-erie/calibration/model/Raven-blended/raven-blended_phase_0_objective_1.nc -a ../../data/objective_1/lake-erie/calibration/gauge_info.csv -b ../../data/objective_1/lake-erie/calibration/model/Raven-blended/subid2gauge.csv -s julie

#    ------------
#    GR4J-Raven-lp
#    ------------
#    python convert_raw_to_netcdf.py -m GR4J-Raven-lp -i ../../data/objective_1/lake-erie/calibration/model/GR4J-Raven-lp/raven-gr4j-lp_phase_0_objective_1_ -o ../../data/objective_1/lake-erie/calibration/model/GR4J-Raven-lp/raven-gr4j-lp_phase_0_objective_1.nc -a ../../data/objective_1/lake-erie/calibration/gauge_info.csv -b ../../data/objective_1/lake-erie/calibration/model/GR4J-Raven-lp/subid2gauge.csv -s julie
#    python convert_raw_to_netcdf.py -m GR4J-Raven-lp -i ../../data/objective_1/lake-erie/calibration/model/GR4J-Raven/raven-gr4j-lp_phase_0_objective_1.csv -o ../../data/objective_1/lake-erie/calibration/model/GR4J-Raven-lp/raven-gr4j-lp_phase_0_objective_1.nc -a ../../data/objective_1/lake-erie/calibration/gauge_info.csv -s hongren

#    ------------
#    GR4J-Raven-sd
#    ------------
#    python convert_raw_to_netcdf.py -m GR4J-Raven-sd -i ../../data/objective_1/lake-erie/calibration/model/GR4J-Raven-sd/raven-gr4j-sd_phase_0_objective_1_ -o ../../data/objective_1/lake-erie/calibration/model/GR4J-Raven-lp/raven-gr4j-sd_phase_0_objective_1.nc -a ../../data/objective_1/lake-erie/calibration/gauge_info.csv -b ../../data/objective_1/lake-erie/calibration/model/GR4J-Raven-sd/subid2gauge.csv -s julie

#    ------------
#    SWAT-EPA
#    ------------
#    python convert_raw_to_netcdf.py -m SWAT-EPA -i ../../data/objective_1/lake-erie/calibration/model/SWAT-EPA/swat-epa_phase_0_objective_1.csv -o ../../data/objective_1/lake-erie/calibration/model/SWAT-EPA/swat-epa_phase_0_objective_1.nc -a ../../data/objective_1/lake-erie/calibration/gauge_info.csv -b ../../data/objective_1/lake-erie/calibration/model/SWAT-EPA/subid2gauge.csv

#    ------------
#    SWAT-Guelph
#    ------------
#    python convert_raw_to_netcdf.py -m SWAT-Guelph -i ../../data/objective_1/lake-erie/calibration/model/SWAT-Guelph/swat-guelph_phase_0_objective_1.csv -o ../../data/objective_1/lake-erie/calibration/model/SWAT-Guelph/swat-guelph_phase_0_objective_1.nc -a ../../data/objective_1/lake-erie/calibration/gauge_info.csv -b ../../data/objective_1/lake-erie/calibration/model/SWAT-Guelph/subid2gauge.csv

#    ------------
#    WATFLOOD
#    ------------
#    python convert_raw_to_netcdf.py -m WATFLOOD -i ../../data/objective_1/lake-erie/calibration/model/WATFLOOD/watflood_phase_0_objective_1.csv -o ../../data/objective_1/lake-erie/calibration/model/WATFLOOD/watflood_phase_0_objective_1.nc -a ../../data/objective_1/lake-erie/calibration/gauge_info.csv

#    ------------
#    MESH-SVS
#    ------------
#    python convert_raw_to_netcdf.py -m MESH-SVS -i ../../data/objective_1/lake-erie/calibration/model/MESH-SVS/mesh-svs_phase_0_objective_1.csv -o ../../data/objective_1/lake-erie/calibration/model/MESH-SVS/mesh-svs_phase_0_objective_1.nc -a ../../data/objective_1/lake-erie/calibration/gauge_info.csv -b ../../data/objective_1/lake-erie/calibration/model/MESH-SVS/subid2gauge.tb0

#    ------------
#    MESH-CLASS
#    ------------
#    python convert_raw_to_netcdf.py -m MESH-CLASS -i ../../data/objective_1/lake-erie/calibration/model/MESH-CLASS/mesh-class_phase_0_objective_1.csv -o ../../data/objective_1/lake-erie/calibration/model/MESH-CLASS/mesh-class_phase_0_objective_1.nc -a ../../data/objective_1/lake-erie/calibration/gauge_info.csv -b ../../data/objective_1/lake-erie/calibration/model/MESH-CLASS/subid2gauge.tb0

# -----------------------
# add subolder scripts/lib to search path
# -----------------------
import sys
import os 
dir_path = os.path.dirname(os.path.realpath(__file__))
sys.path.append(dir_path+'/../lib')

import argparse
import textwrap                # nicer formatting of help text in parser
import numpy as np             # to perform numerics
import shutil                  # file operations
import copy                    # deep copy objects, arrays etc
import datetime                # converting dates
import netCDF4 as nc           # NetCDF writing
import glob                    # for file listing
import xarray as xr
import scipy.io                # to read *.mat files

#import netcdf4     as     nc4         # in lib/
from fread         import fread        # in lib/
from fsread        import fsread       # in lib/
from writenetcdf   import writenetcdf  # in lib/
from dec2date      import dec2date     # in lib/
from date2dec      import date2dec     # in lib/

model                      = ['LBRM']
input_file                 = ['lbrm_phase_0_objective_1.csv']
output_file                = ['lbrm_phase_0_objective_1.nc']
gaugeinfo_file             = ['gauge_info.csv']
mapping_subbasinID_gaugeID = ['']
setup_by                   = None
parser      = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
              description='''Convert LBRM raw streamflow model outputs into NetDF format (consistent across all models in GRIP-E).''')
parser.add_argument('-i', '--input_file', action='store',
                    default=input_file, dest='input_file', metavar='input_file', nargs=1,
                    help='Name of input file (raw model output file; for HYPE, HYMOD2-DS, and Raven basename (filename to be expected to be basename<gaugeId>.txt).')
parser.add_argument('-o', '--output_file', action='store',
                    default=output_file, dest='output_file', metavar='output_file', nargs=1,
                    help='Name of output file (NetCDF file).')
parser.add_argument('-m', '--model', action='store',
                    default=model, dest='model', metavar='model', nargs=1,
                    help='Model (e.g., LBRM, VIC, VIC-GRU, HYPE, GEM-Hydro).')
parser.add_argument('-a', '--gaugeinfo_file', action='store',
                    default=gaugeinfo_file, dest='gaugeinfo_file', metavar='gaugeinfo_file', nargs=1,
                    help='File containing additional information about gauges (e.g., darinage area, lat/lon).')
parser.add_argument('-b', '--mapping_subbasinID_gaugeID', action='store',
                    default=mapping_subbasinID_gaugeID, dest='mapping_subbasinID_gaugeID', metavar='mapping_subbasinID_gaugeID',nargs=1,
                    help='File containing mapping of subbasin ID (col 1) to gauge ID (col 2). All other columns are ignored. One header line. Only required for VIC and VIC-GRU.')
parser.add_argument('-s', '--setup_by', action='store',
                    default=setup_by, dest='setup_by', metavar='setup_by',
                    help='Model was setup by this person. Outputs might vary. E.g. Raven setup by Julie dumps one file for each gauge. Raven setup by Hongren dumps all in same file.')

args                       = parser.parse_args()
model                      = args.model[0]
input_file                 = args.input_file[0]
output_file                = args.output_file[0]
gaugeinfo_file             = args.gaugeinfo_file[0]
mapping_subbasinID_gaugeID = args.mapping_subbasinID_gaugeID[0]
setup_by                   = args.setup_by

del parser, args

# nodata
nodata = -9999.0

if ( (model != 'LBRM')                 and
     (model != 'VIC')                  and
     (model != 'VIC-GRU')              and
     (model != 'GEM-Hydro')            and
     (model != 'HYPE')                 and
     (model != 'HYMOD2-DS')            and  
     (model != 'ML-ConvLSTM')          and          # Lake Erie
     (model != 'ML-ConvLSTM-DEM')      and          # Lake Erie
     (model != 'ML-ConvLSTM-LC')       and          # Lake Erie
     (model != 'ML-ConvLSTM-LC-DEM')   and          # Lake Erie
     (model != 'ML-LinReg')            and          # Lake Erie
     (model != 'ML-XGBoost')           and          # Lake Erie
     (model != 'ML-EA-LSTM')           and          # Great Lakes
     (model != 'ML-LSTM')              and          # Great Lakes
     (model != 'LBRM-ML-LSTM')         and          # Great Lakes
     (model != 'LBRM-MG')              and          # Great Lakes
     (model != 'HMETS-Raven-lp')       and          # Great Lakes
     (model != 'GR4J-Raven-lp')        and
     (model != 'GR4J-Raven-sd')        and
     (model != 'Raven-blended')        and          # Lake Erie
     (model != 'SWAT-EPA')             and
     (model != 'SWAT-Guelph')          and
     (model != 'WATFLOOD')             and
     (model != 'MESH-SVS')             and    
     (model != 'MESH-CLASS')           and    
     (model != 'mHM-Waterloo')         and    
     (model != 'mHM-UFZ') ):
    raise ValueError('This model is not supported yet!')

if ( ((model == 'VIC-GRU')                                and (mapping_subbasinID_gaugeID == '')) or
     ((model == 'VIC')                                    and (mapping_subbasinID_gaugeID == '')) or
     ((model == 'HMETS-Raven-lp' and setup_by == 'julie') and (mapping_subbasinID_gaugeID == '')) or
     ((model == 'Raven-blended'  and setup_by == 'julie') and (mapping_subbasinID_gaugeID == '')) or
     ((model == 'GR4J-Raven-lp'  and setup_by == 'julie') and (mapping_subbasinID_gaugeID == '')) or
     ((model == 'GR4J-Raven-sd'  and setup_by == 'julie') and (mapping_subbasinID_gaugeID == '')) or
     ((model == 'SWAT-EPA')                               and (mapping_subbasinID_gaugeID == '')) or
     ((model == 'SWAT-Guelph')                            and (mapping_subbasinID_gaugeID == '')) or
     ((model == 'MESH-SVS')                               and (mapping_subbasinID_gaugeID == '')) or
     ((model == 'MESH-CLASS')                             and (mapping_subbasinID_gaugeID == '')) ):
    raise ValueError('For VIC, SWAT-EPA, SWAT-Guelph, and RAVEN model CSV file containing the mapping of subbasin ID (col 1) to gauge ID (col 2) needs to be provided. All other columns in that file will be ignored. Exactly one header line needs to be provided.\n For MESH-SVS and MESH-CLASS the file is assumed to be a model setup tb0 file where only the line with :ColumnName is read. It should contain the gauge names. The order of the gauges in :ColumnName is assumed to be the order of the columns in the MESH csv output files.')

if ((model == 'GR4J-Raven-lp' or model == 'GR4J-Raven-sd' or model == 'HMETS-Raven-lp' or model == 'Raven-blended') and (setup_by is None)):
    raise ValueError('For GR4J-Raven-lp and GR4J-Raven-sd and HMETS-Raven-lp and Raven-blended the person who has setup the model needs to be named.')
if ( ((model == 'GR4J-Raven-lp' or model == 'GR4J-Raven-sd' or model == 'HMETS-Raven-lp' or model == 'Raven-blended') and not(setup_by == 'julie' or setup_by == 'hongren')) ):
    raise ValueError('Person who has setup GR4J-Raven-lp or GR4J-Raven-sd and HMETS-Raven-lp and Raven-blended must be "julie" or "hongren".')

if (model == 'mHM-UFZ' or model == 'mHM-Waterloo'):
    # ---------------
    # read model outputs
    # - every gauge is in a separate file (NetCDF)
    # ---------------
    input_files    = glob.glob(input_file+"*.nc")
    model_stations = [ ii.split(input_file)[1].split('.')[0] for ii in input_files ]
    model_data     = [ [] for ii in input_files ]
    model_dates    = None
    for ii,iinput_file in enumerate(input_files):

        tmp = xr.open_dataset(iinput_file)

        # find variable that contains simulation results
        var_sim=[s for s in tmp.variables.keys() if 'Qsim' in s][0]

        model_data[ii] = tmp[var_sim].data

        # make sure all model dates are same in all files
        if model_dates is None:
            # save first file's dates
            model_dates = tmp[var_sim].time.data
            model_dates = [ itt.astype('M8[D]').astype('O') for itt in model_dates ]
        else:
            # check if dates are same as already saved
            tmp_dates = tmp[var_sim].time.data
            tmp_dates = [ itt.astype('M8[D]').astype('O') for itt in tmp_dates ]
            if not np.all(tmp_dates == model_dates):
                print('Time steps first file: ',input_files[0])
                print('     ',model_dates)
                print('Time steps current file: ',input_files[ii])
                print('     ',tmp_dates)
                raise ValueError('Time step in files must be all the same!')

    model_data  = np.transpose(np.array(model_data))
    model_dates = np.transpose(np.array(model_dates))
        

# read model output file
if (model == 'HYPE'):
    # ---------------
    # read model outputs
    # - every gauge is in a separate file (ASCII)
    # ---------------
    input_files    = glob.glob(input_file+"*.txt")
    model_stations = [ ii.split(input_file)[1].split('.')[0] for ii in input_files ]
    model_data     = [ [] for ii in input_files ]
    model_dates    = None
    for ii,iinput_file in enumerate(input_files):
        
        # find column containing discharge "cout"
        head = fread(iinput_file,skip=2,cskip=1,header=True)
        idx = head[0].index('cout')
        
        model_data[ii]  = fread(iinput_file,skip=2,cskip=1,header=False)[:,idx]

        # make sure all model dates are same in all files
        if model_dates is None:
            # save first file's dates
            model_dates = fsread(iinput_file,skip=2,snc=1)
            model_dates = [ datetime.datetime( int(str(iii[0])[0:4]),int(str(iii[0])[5:7]),int(str(iii[0])[8:10]),0,0 ) for iii in model_dates ]
        else:
            # check if dates are same as already saved
            tmp_dates = fsread(iinput_file,skip=2,snc=1)
            tmp_dates = [ datetime.datetime( int(str(iii[0])[0:4]),int(str(iii[0])[5:7]),int(str(iii[0])[8:10]),0,0 ) for iii in tmp_dates ]
            if not np.all(tmp_dates == model_dates):
                print('Time steps first file: ',input_files[0])
                print('     ',model_dates)
                print('Time steps current file: ',input_files[ii])
                print('     ',tmp_dates)
                raise ValueError('Time step in files must be all the same!')
            

    model_data  = np.transpose(np.array(model_data))
    model_dates = np.transpose(np.array(model_dates))


# read model output file
if (model == 'HYMOD2-DS'):
    # ---------------
    # read model outputs
    # - every gauge is in a separate file (MAT)
    # ---------------
    input_files    = glob.glob(input_file+"*.mat")
    model_stations = [ ii.split(input_file)[1].split('.')[0] for ii in input_files ]
    model_data     = [ [] for ii in input_files ]
    model_dates    = None
    for ii,iinput_file in enumerate(input_files):

        # read matlab file
        mat = scipy.io.loadmat(iinput_file)

        model_data[ii]  = mat['hymod_result']['simulated_flow'][0][0][:,0] 

        # make sure all model dates are same in all files
        if model_dates is None:
            # save first file's dates
            model_dates = mat['hymod_result']['simulation_period'][0][0]  # formatted as, e.g. 01-Jan-2010
            model_dates = [ datetime.datetime.strptime(iii, '%d-%b-%Y') for iii in model_dates ]
        else:
            # check if dates are same as already saved
            tmp_dates = mat['hymod_result']['simulation_period'][0][0]  # formatted as, e.g. 01-Jan-2010
            tmp_dates = [ datetime.datetime.strptime(iii, '%d-%b-%Y') for iii in tmp_dates ]
            if not np.all(tmp_dates == model_dates):
                print('Time steps first file: ',input_files[0])
                print('     ',model_dates)
                print('Time steps current file: ',input_files[ii])
                print('     ',tmp_dates)
                raise ValueError('Time step in files must be all the same!')
            

    model_data  = np.transpose(np.array(model_data))
    model_dates = np.transpose(np.array(model_dates))
    

# read model output file
if (model == 'GR4J-Raven-lp' or model == 'GR4J-Raven-sd' or model == 'HMETS-Raven-lp' or model == 'Raven-blended'):

    if (setup_by == 'julie'):
        input_files    = glob.glob(input_file+"*.csv")
        model_stations = [ ii.split(input_file)[1].split('.')[0] for ii in input_files ]
        model_data     = [ [] for ii in input_files ]
        model_dates    = None

        # ---------------
        # read mapping info subbasin ID --> gauge station ID
        # ---------------
        mapping = fsread(mapping_subbasinID_gaugeID,skip=0,snc=2)
        mapping = np.array(mapping)

        # ---------------------------------------------------------------------------
        # outputs are distributed in separte files (RAVEN format) == JULIE
        # ---------------------------------------------------------------------------
        for ii,iinput_file in enumerate(input_files):
          
            # find column containing discharge "cout"
            head = fread(iinput_file,skip=1,cskip=4,header=True)

            # find column with subbasin ID matching the gauge ID in file name (saved in 'model_stations')
            subID = mapping[np.where(mapping[:,1]==model_stations[ii])[0][0]][0]
            desired_column_header = 'sub'+subID+' [m3/s]'

            if desired_column_header in head:
                idx = head.index(desired_column_header)
            else:
                desired_column_header = 'raven-weighted [m3/s]'
                if desired_column_header in head:
                    idx = head.index(desired_column_header)
                else:
                    raise ValueError('Column header not found in '+iinput_file)
          
            model_data[ii]  = fread(iinput_file,skip=1,cskip=4,header=False,fill=True,fill_value=nodata)[:,idx]

            # make sure all model dates are same in all files
            if model_dates is None:
                # save first file's dates
                model_dates    = fsread(iinput_file,skip=1,cskip=1,snc=2)
                model_dates    = [ datetime.datetime( int(str(mm[0])[0:4]),int(str(mm[0])[5:7]),int(str(mm[0])[8:10]),int(str(mm[1])[0:2]),int(str(mm[1])[3:5]) ) - datetime.timedelta(days=1) for mm in model_dates ]
            else:
                # check if dates are same as already saved
                tmp_dates = fsread(iinput_file,skip=1,cskip=1,snc=2)
                tmp_dates = [ datetime.datetime( int(str(mm[0])[0:4]),int(str(mm[0])[5:7]),int(str(mm[0])[8:10]),int(str(mm[1])[0:2]),int(str(mm[1])[3:5]) ) - datetime.timedelta(days=1) for mm in tmp_dates ]
                if not np.all(tmp_dates == model_dates):
                    print('Time steps first file: ',input_files[0])
                    print('     ',model_dates)
                    print('Time steps current file: ',input_files[ii])
                    print('     ',tmp_dates)
                    raise ValueError('Time step in files must be all the same!')
                
        model_data  = np.transpose(np.array(model_data))
        model_dates = np.transpose(np.array(model_dates))
        
    elif (setup_by == 'hongren'):

        # ---------------
        # read model outputs
        # ---------------
        # header names look like "'ObsQ_02GA018 [m3/s]', 'SimQ_02GA038 [m3/s]', ..."   --> THIS IS NOT RAVEN OUTPUT FORMAT!!!!
        head = fread(input_file,skip=1,cskip=1,header=True)

        # find column with 'SimQ_*'
        idx = [ii for ii, ss in enumerate(head) if 'SimQ' in ss]

        # save station names
        model_stations = np.array([ii.split('_')[1].split(' ')[0] for ii in head ])[idx]
        model_stations = list(model_stations)

        # save simulated data
        model_data  = fread(input_file,skip=1,cskip=1,header=False,fill=True,fill_value=nodata)[:,idx]
        model_data  = np.array(model_data,dtype=np.float32)

        # model dates --> THIS IS NOT RAVEN OUTPUT FORMAT!!!!
        # shifting time because RAven reports period beginning (but observations are period ending)
        model_dates    = fsread(input_file,skip=1,cskip=1,snc=[0])
        model_dates    = [ datetime.datetime( int(str(mm[0])[0:4]),int(str(mm[0])[5:7]),int(str(mm[0])[8:10]),int(str(mm[0])[11:13]),int(str(mm[0])[14:16]) ) - datetime.timedelta(days=1) for mm in model_dates ]

    else:

        raise ValueError('Person who has setup GR4J-Raven must be "julie" or "hongren".')

    # ---------------------------------------------------------------------------
    # HONGREN'S outputs
    # ---------------------------------------------------------------------------
    # # ---------------
    # # read mapping info subbasin ID --> gauge station ID
    # # ---------------
    # mapping = fsread(mapping_subbasinID_gaugeID,skip=0,snc=2)
    # mapping = np.array(mapping)

    # # ---------------
    # # read model outputs
    # # - every gauge is in a separate file
    # # ---------------
    # for ii,iinput_file in enumerate(input_files):
        
    #     # find column containing discharge "cout"
    #     head = fread(iinput_file,skip=1,cskip=4,header=True)

    #     # find column with subbasin ID matching the gauge ID in file name (saved in 'model_stations')
    #     desired_column_header = 'sub'+mapping[np.where(mapping[:,1] == model_stations[ii])[0][0]][0]+' [m3/s]'
        
    #     idx = head.index(desired_column_header)
        
    #     model_data[ii]  = fread(iinput_file,skip=1,cskip=4,header=False,fill=True,fill_value=nodata)[:,idx]

    #     # make sure all model dates are same in all files
    #     if model_dates is None:
    #         # save first file's dates
    #         model_dates    = fsread(iinput_file,skip=1,cskip=1,snc=2)
    #         model_dates    = [ datetime.datetime( int(str(mm[0])[0:4]),int(str(mm[0])[5:7]),int(str(mm[0])[8:10]),int(str(mm[1])[0:2]),int(str(mm[1])[3:5]) ) for mm in model_dates ]
    #     else:
    #         # check if dates are same as already saved
    #         tmp_dates = fsread(iinput_file,skip=1,cskip=1,snc=2)
    #         tmp_dates = [ datetime.datetime( int(str(mm[0])[0:4]),int(str(mm[0])[5:7]),int(str(mm[0])[8:10]),int(str(mm[1])[0:2]),int(str(mm[1])[3:5]) ) for mm in tmp_dates ]
    #         if not np.all(tmp_dates == model_dates):
    #             print('Time steps first file: ',input_files[0])
    #             print('     ',model_dates)
    #             print('Time steps current file: ',input_files[ii])
    #             print('     ',tmp_dates)
    #             raise ValueError('Time step in files must be all the same!')
            
    # model_data  = np.transpose(np.array(model_data))
    # model_dates = np.transpose(np.array(model_dates))
    
if (model == 'LBRM'):
    # ---------------
    # read model outputs
    # ---------------
    model_stations = fread(input_file,skip=1,cskip=1,header=True)
    model_data     = fread(input_file,skip=1,cskip=1,header=False)
    model_data     = np.array(model_data,dtype=np.float32)
    model_dates    = fsread(input_file,skip=1,snc=1)
    model_dates    = [ datetime.datetime( int(str(ii[0])[0:4]),int(str(ii[0])[5:7]),int(str(ii[0])[8:10]),0,0 ) for ii in model_dates ]

if (model == 'WATFLOOD'):
    # ---------------
    # read model outputs
    # ---------------

    input_f = open(input_file, "r")
    dump = input_f.readlines()
    input_f.close()
    
    dump = [ dd.strip() for dd in dump ]    # remove leading and trailing blanks and '\n'

    # indexes of lines that start with "location"
    idx_location_line = np.where( [ line.startswith('location') for line in dump ] )[0]

    # col_1 = # of gauges in following block
    # col_2 = max hour
    # col_3 = max hour
    # col_4 = delta hour
    # next block will have <col_1> lines and some columns; first column contains gauge ID
    # block after has <col_1*2 + 1> columns; first is time step in hours; then obs_1, mod_1, obs_2, mod_2, ..., obs_ngauge, mod_ngauge
    #                 <col_3> / <col_4> rows
    ref_dates  = [ ' '.join(dump[ii-2].split()).split('\\')[-1].split('_')[0] for ii in idx_location_line ]
    ref_dates  = [ datetime.datetime( int(str(ii)[0:4]),int(str(ii)[4:6]),int(str(ii)[6:8]),0,0 ) for ii in ref_dates ]
    block_info = np.array( [ ' '.join(dump[ii-1].split()).split(' ') for ii in idx_location_line],dtype=np.int )

    # stations per block
    model_stations_block = [ [ ' '.join(dump[jj].split()).split(' ')[0] for jj in range(ii+1, ii+block_info[iii,0]+1)] for iii,ii in enumerate(idx_location_line) ]

    # unique list of stations (first flatten list of lists)
    model_stations_uniq = [item for sublist in model_stations_block for item in sublist]
    model_stations_uniq  = np.unique( model_stations_uniq )

    # model dates are ref_dates plus hours (first column in second block)
    hours_block = [ [ int(' '.join(dump[jj].split()).split(' ')[0]) for jj in range(ii+block_info[iii,0]+2,
                                                                      ii+block_info[iii,0]+2+block_info[iii,2]/block_info[iii,3])]
                                                      for iii,ii in enumerate(idx_location_line) ]
    # ASK FRANK IF FIRST DAY IS JAN 2 (-0) OR JAN 1 (-24) 
    model_dates_block = [ [ iref + datetime.timedelta(hours=jj-24) for jj in hours_block[ii] ] for ii,iref in enumerate(ref_dates) ]

    # unique list of dates (first flatten list of lists)
    model_dates_uniq = [item for sublist in model_dates_block for item in sublist]
    model_dates_uniq  = np.unique( model_dates_uniq )

    # data of each block; not sorted yet following model_stations_uniq
    model_data_block = [ [ list(np.array(' '.join(dump[jj].split()).split(' ')[2::2],dtype=np.float32)) for jj in range(ii+block_info[iii,0]+2,
                                                                      ii+block_info[iii,0]+2+block_info[iii,2]/block_info[iii,3])] for iii,ii in enumerate(idx_location_line) ]

    # finally put data
    model_data = np.ones([np.shape(model_dates_uniq)[0], np.shape(model_stations_uniq)[0]]) * -9999.9

    for iblock in range(len(model_data_block)):
        row_idx = np.array([ list(model_dates_uniq).index(ii) for ii in model_dates_block[iblock] ])
        col_idx = np.array([ list(model_stations_uniq).index(ii) for ii in model_stations_block[iblock] ])

        model_data[row_idx[:, np.newaxis],col_idx] = np.array(model_data_block[iblock])

    model_stations = list(model_stations_uniq)
    model_dates    = list(model_dates_uniq)

if (model == 'ML-ConvLSTM' or model == 'ML-ConvLSTM-DEM' or model == 'ML-ConvLSTM-LC' or model == 'ML-ConvLSTM-LC-DEM' or model == 'ML-LinReg' or model == 'ML-XGBoost' or model == 'ML-EA-LSTM' or model == 'ML-LSTM' or model == 'LBRM-MG' or model == 'LBRM-ML-LSTM' ):    #'ANN-LinReg'
    # ---------------
    # read model outputs
    # - model outputs in pickle exported to CSV
    # - stations are blockwise stored one after each other
    # - ,date,prediction,actual,station
    # ---------------
    model_stations = np.unique(fsread(input_file,skip=1,header=False,snc=[4]))
    nstations      = np.shape(model_stations)[0]

    # nodata lines are just not written
    # --> create date range from first to last date
    first_date = np.sort(np.unique(fsread(input_file,skip=1,snc=[1])))[0]
    first_date = datetime.datetime(int(str(first_date)[0:4]),int(str(first_date)[5:7]),int(str(first_date)[8:10]),0,0)
    last_date  = np.sort(np.unique(fsread(input_file,skip=1,snc=[1])))[-1]
    last_date  = datetime.datetime(int(str(last_date)[0:4]),int(str(last_date)[5:7]),int(str(last_date)[8:10]),0,0)
    ndays = (last_date-first_date).days+1
    model_dates = np.array([ first_date + datetime.timedelta(ii) for ii in range(ndays) ])

    # read just everything (f=floats, s=strings)
    all_data_f, all_data_s = fsread(input_file,skip=1,snc=[1,4],nc=[0,2,3])
    all_data_s = np.array(all_data_s) # make array
    # all_data_s[:,0] = np.array([ datetime.datetime(int(str(ii)[0:4]),int(str(ii)[5:7]),int(str(ii)[8:10]),0,0) for ii in all_data_s[:,0] ])

    # allocate final array
    model_data = np.ones([ndays,nstations]) * nodata

    # find which column contains 'prediction'
    head_input = fsread(input_file,cskip=1,skip=1,header=True)
    idx_prediction = head_input.index('prediction')

    # map data onto data array
    for istation in model_stations:

        idx_station = np.where(all_data_s[:,1] == istation)[0]   # all lines of that station
        iistation   = np.where(model_stations == istation)[0][0]
        iitime = [ np.where(model_dates == idate)[0][0]  for idate in [ datetime.datetime(int(str(ii)[0:4]),int(str(ii)[5:7]),int(str(ii)[8:10]),0,0) for ii in all_data_s[idx_station,0] ] ]
        iitime = np.array(iitime)
        dats = all_data_f[idx_station,idx_prediction] # columns are [running ID, prediction, actual] or [running ID, actual, prediction]
        dats = np.where(dats<0.0,0.0,dats) 
        model_data[iitime,iistation] = dats 

    model_stations = list(model_stations)

if (model == 'GEM-Hydro'):
    # ---------------
    # read model outputs
    # - model outputs already pre-processes using "strf_graphs_scores.py"
    # ---------------
    model_stations = fread(input_file,skip=1,cskip=2,header=True)
    model_data     = fread(input_file,skip=1,cskip=2,header=False)
    model_data     = np.array(model_data,dtype=np.float32)
    model_dates    = fsread(input_file,skip=1,snc=2)
    model_dates    = [ datetime.datetime( int(str(ii[0])[0:4]),int(str(ii[0])[5:7]),int(str(ii[0])[8:10]),int(str(ii[1])[0:2]),int(str(ii[1])[3:5]) ) for ii in model_dates ]

if (model == 'VIC-GRU' or model == 'VIC'):
    # ---------------
    # read model outputs
    # - model outputs contain subbasin ID and not gauge ID --> need to be remapped
    # ---------------
    model_stations = fread(input_file,skip=1,cskip=4,header=True)   # this is subbasin IDs    
    model_data     = fread(input_file,skip=1,cskip=4,header=False,fill=True,fill_value=nodata)
    model_data     = np.array(model_data,dtype=np.float32)
    model_dates    = fsread(input_file,skip=1,cskip=1,snc=2)
    model_dates    = [ datetime.datetime( int(str(ii[0])[0:4]),int(str(ii[0])[5:7]),int(str(ii[0])[8:10]),int(str(ii[1])[0:2]),int(str(ii[1])[3:5]) ) - datetime.timedelta(days=1) for ii in model_dates ]

    # ---------------
    # read mapping info subbasin ID --> gauge station ID
    # ---------------
    mapping = fsread(mapping_subbasinID_gaugeID,skip=0,snc=2)
    mapping = np.array(mapping)
    
    # ---------------
    # map subbasin ID's to gauging stations IDs
    # ---------------
    for ii,isubbasin in enumerate(model_stations):  # they look like "sub676 [m3/s]" --> "676"

        if not('observed' in isubbasin):
            subbasin_ID = isubbasin.split(' ')[0].split('sub')[1]
            idx = np.where(mapping[:,0]==subbasin_ID)[0][0]
            gauge_id = mapping[idx,1]

            model_stations[ii] = gauge_id

if (model == 'SWAT-EPA'):

    # ---------------
    # read mapping info subbasin ID --> gauge station ID
    # ---------------
    mapping = fsread(mapping_subbasinID_gaugeID,skip=1,snc=[0,1,2])
    mapping = np.array(mapping)
    mapping[:,0] = np.array( [ii.strip() for ii in mapping[:,0]] )   # remove trailing blanks: subbasin ID
    mapping[:,1] = np.array( [ii.strip() for ii in mapping[:,1]] )   # remove trailing blanks: USGS/WSC gauge ID
    mapping[:,2] = np.array( [ii.strip() for ii in mapping[:,2]] )   # remove trailing blanks: variable name in output file
    nstations = np.shape(mapping)[0]
    
    # ---------------
    # read model outputs
    # - model outputs contain subbasin ID and not gauge ID --> need to be remapped
    # ---------------
    #
    # col 1 :: tmp[idx,0] :: subbasin ID: 1-199 (basins that are important are in mapping)
    # col 2 :: tmp[idx,1] :: Year
    # col 3 :: tmp[idx,2] :: day of the year (1-366)
    # col ? :: tmp[idx,3] :: FLOW_OUTcms = Q[m^3/s]

    header = fread(input_file, skip=1,cskip=0,header=True)
    
    required_data_vars = np.unique(mapping[:,2])   # ['FLOW_OUTcms','FLOW_INcms']
    required_data_idx  = [ header.index(rr) for rr in required_data_vars ]  # [5,6]
    
    tmp = fread(input_file,nc=[0,1,2]+required_data_idx, skip=1,cskip=0,header=False)

    model_stations_all = np.array(np.array(tmp[:,0],np.int),dtype=str)   # this is subbasin IDs
    ntime = np.shape(np.where(model_stations_all==model_stations_all[0])[0])[0]

    model_data = np.ones([ntime,nstations],dtype=np.float32) * -9999.9
    model_stations = list(np.array(['None']*nstations,dtype=str))
    for istation,station in enumerate([ str(int(ii)) for ii in mapping[:,0] ]):

        idx = np.where(model_stations_all==station)[0]   # lines containing data for current station
        idx_var = np.where( required_data_vars == mapping[istation,2] )[0]   # find which column the current variable is in
        # print('var = ',mapping[istation,2], ' in column #',3+idx_var, ' of tmp')
        imodel_data = tmp[idx,3+idx_var]                 # FLOW_OUTcms for current station
        ndat = np.shape(imodel_data)[0]                  # number of data points
        if ndat != ntime:
            print('ERROR :: station: ',station, '   # of time steps = ',ndat, '         expected # of time steps = ',ntime)
            stop
        
        model_data[:,istation] = imodel_data
        model_stations[istation] = station   # this is subbasin IDs    

        year = tmp[idx,1]  # YEAR
        doy  = tmp[idx,2]  # MON = day of the year (1...366)

    model_dates = [ dec2date(date2dec(yr=year[ii],mo=1,dy=1)-1.0 + doy[ii]) for ii in range(len(year)) ]
    model_dates = [ datetime.datetime( ii[0],ii[1],ii[2],0,0 ) for ii in model_dates ]
    
    # ---------------
    # map subbasin ID's to gauging stations IDs
    # ---------------
    for ii,isubbasin in enumerate(model_stations):  # they look like "sub676 [m3/s]" --> "676"

        subbasin_ID = str(isubbasin) # isubbasin.split(' ')[0].split('sub')[1]
        idx = np.where(mapping[:,0]==subbasin_ID)[0][0]
        gauge_id = mapping[idx,1]

        model_stations[ii] = gauge_id

if (model == 'SWAT-Guelph'):

    # ---------------
    # read mapping info subbasin ID --> gauge station ID
    # ---------------
    mapping = fsread(mapping_subbasinID_gaugeID,skip=1,snc=[0,1,2])
    mapping = np.array(mapping)
    mapping[:,0] = np.array( [ii.strip() for ii in mapping[:,0]] )   # remove trailing blanks: subbasin ID
    mapping[:,1] = np.array( [ii.strip() for ii in mapping[:,1]] )   # remove trailing blanks: USGS/WSC gauge ID
    mapping[:,2] = np.array( [ii.strip() for ii in mapping[:,2]] )   # remove trailing blanks: variable name in output file
    nstations = np.shape(mapping)[0]
    
    # ---------------
    # read model outputs
    # - model outputs contain subbasin ID and not gauge ID --> need to be remapped
    # ---------------
    #
    # col 1 :: REACH  (only string)
    # col 2 :: reach ID: 1-699 (something)
    # col 3 :: GIS    (whole time 0)
    # col 4 :: MON day of the year (1-366)
    # col ? :: FLOW_INcms   = Q[m^3/s]
    # col ? :: FLOW_OUTcms  = Q[m^3/s]

    header = fread(input_file, skip=1,cskip=0,full_header=True,header=True,separator=' ',strip=True,squeeze=True)[0]
    header = ' '.join(header.strip().split())   # merge several blanks to one blank
    header = header.split()

    required_data_vars = np.unique(mapping[:,2])   # ['FLOW_OUTcms','FLOW_INcms']
    required_data_idx  = [ header.index(rr) for rr in required_data_vars ]  # [5,6]
    
    tmp = fread(input_file,nc=[1,3]+required_data_idx, skip=1,cskip=0,header=False)          # No YEAR!!!

    model_stations_all = np.array(np.array(tmp[:,0],np.int),dtype=str)   # this is subbasin IDs
    ntime = np.shape(np.where(model_stations_all==model_stations_all[0])[0])[0]

    model_data = np.ones([ntime,nstations],dtype=np.float32) * -9999.9
    model_stations = list(np.array(['None']*nstations,dtype=str))
    for istation,station in enumerate([ str(int(ii)) for ii in mapping[:,0] ]):

        idx = np.where(model_stations_all==station)[0]   # lines containing data for current station
        idx_var = np.where( required_data_vars == mapping[istation,2] )[0]   # find which column the current variable is in
        # print('var = ',mapping[istation,2], ' in column #',3+idx_var, ' of tmp')
        imodel_data = tmp[idx,2+idx_var]                 # FLOW_OUTcms for current station
        ndat = np.shape(imodel_data)[0]                  # number of data points
        if ndat != ntime:
            print('ERROR :: station: ',station, '   # of time steps = ',ndat, '         expected # of time steps = ',ntime)
            stop
        
        model_data[:,istation] = imodel_data
        model_stations[istation] = station   # this is subbasin IDs

        doy  = tmp[idx,1]  # MON = day of the year (1...366)

        # create YEAR array since it is missing in file output
        year = np.ones(ntime,dtype=np.int) * 2011                        # <<<<<<<<<<<< ASSUMED start year 2011
        for itime in range(1,ntime):
            if doy[itime] != 1:
                year[itime] = year[itime-1]
            else:
                year[itime] = year[itime-1] + 1        

    model_dates = [ dec2date(date2dec(yr=year[ii],mo=1,dy=1)-1.0 + doy[ii]) for ii in range(len(year)) ]
    model_dates = [ datetime.datetime( ii[0],ii[1],ii[2],0,0 ) for ii in model_dates ]
    
    # ---------------
    # map subbasin ID's to gauging stations IDs
    # ---------------
    for ii,isubbasin in enumerate(model_stations):  # they look like "sub676 [m3/s]" --> "676"

        subbasin_ID = str(isubbasin) # isubbasin.split(' ')[0].split('sub')[1]
        idx = np.where(mapping[:,0]==subbasin_ID)[0][0]
        gauge_id = mapping[idx,1]

        model_stations[ii] = gauge_id

if (model == 'MESH-SVS' or model == 'MESH-CLASS'):
    # ---------------
    # read model outputs
    # - model outputs are sorted in a way we dont know yet (is in tb0 setup file) --> need remapping
    # ---------------
    model_stations = fread(input_file,skip=1,cskip=2,header=True)   # this is just generic column names [QOMEAS1, QOSIM1, QOMEAS2, QOSIM2, ....]   
    model_data     = fread(input_file,skip=1,cskip=2,header=False,fill=True,fill_value=nodata)
    model_data     = np.array(model_data,dtype=np.float32)[:,1::2] # only every second column, others are obeservations
    model_dates    = fread(input_file,skip=1,cskip=0,header=False,fill=True,fill_value=nodata)
    model_dates    = np.array(model_dates,dtype=np.float32)[:,0:2] # only first two columns: col1=year, col2=doy
    model_dates    = [ datetime.datetime( int(ii[0]-1),12,31,0,0 ) + datetime.timedelta(int(ii[1])) for ii in model_dates ]

    # ---------------
    # read names of stations
    # ---------------
    input_f = open(mapping_subbasinID_gaugeID, "r")
    dump = input_f.readlines()
    input_f.close()

    dump = [ dd.strip() for dd in dump ]    # remove leading and trailing blanks and '\n'

    found = False
    for dd in dump:
        # print("dd: ",dd)
        first = dd.split()[0]
        if first.upper() == ':COLUMNNAME':
            mapping = dd.split()[1:]
            found = True

    if not(found):
        raise ValueError('Line starting with :ColumnName not found in file '+mapping_subbasinID_gaugeID)
    else:
        model_stations = mapping
    

# get gauge station file
gaugeinfo_header = fsread(gaugeinfo_file,comment='#',separator=',',skip=1,header=True)
ncols = len(gaugeinfo_header)
gaugeinfo_data   = fsread(gaugeinfo_file,comment='#',separator=',',skip=1,snc=ncols,header=False)

# read gauge information from CSV file
# [NO,ID,Name,Lat,Lon,Country,Drainage_area,...]
gaugeinfo_all  = []
station_id     = model_stations                                   # station IDs
with open(gaugeinfo_file) as f:
    content = f.readlines()

for cc in content[1:]:
    ccc = ' '.join(cc.strip().split())

    # split at "," but not when in quotes
    import re
    PATTERN = re.compile(r'''((?:[^,"']|"[^"]*"|'[^']*')+)''')
    # print( PATTERN.split(ccc)[1::2] )
    ccc = PATTERN.split(ccc)[1::2]
    
    gaugeinfo_all.append(ccc[1:7])

# slim down gauge information to only gauges present (in right order)
station_id_all = [ ii[0] for ii in gaugeinfo_all ]

# ------------------------
# only dump stations to NetCDF that are in gauge_info.csv
# ------------------------
# check if data read from model are only stations required (see "station_id_all"); if not, delete data from:
#      model_stations
#      model_data
nn = len(model_stations)
ii = 0
while ii < nn:
    if not(model_stations[ii] in station_id_all):
        # print('not in: ',model_stations[ii])
        model_stations.pop(ii)
        model_data = np.delete(model_data,ii,axis=1)
        nn -= 1
    else:
        # print('in:     ',model_stations[ii])
        ii += 1

# ------------------------
# do not dump stations multiple times
# ------------------------
unique_stations = list(np.unique(model_stations))
for uu in unique_stations:
    
    idx = np.sort(np.where(np.array(model_stations) == uu)[0])    
    for ii in idx[:-1]: #idx[:0:-1]:
        # delete that column from
        #      model_stations
        #      model_data
        model_stations.pop(ii)
        model_data = np.delete(model_data,ii,axis=1)

idx_stations   = [ station_id_all.index(ss) for ss in station_id ]
gaugeinfo      = [ gaugeinfo_all[idx] for idx in idx_stations ]   # gauge information of only the requested gauges: [ID,Name,Lat,Lon,Country,Drainage_area]
station_info   = [[] for istation in range(len(gaugeinfo))]
for istation in range(len(gaugeinfo)):
    station_info[istation] = gaugeinfo[istation][0]+' : '+gaugeinfo[istation][1]+' ('+gaugeinfo[istation][4]+')'

# times
start_day            = model_dates[0]
ref_date             = 'days since '+str(start_day)
times_since_in_days = [ ((tt - start_day).total_seconds())/(60.*60.*24.) for tt in model_dates ]

# ----------------------------------------------
# write NetCDF
# ----------------------------------------------

# open netcdf file and add some general information
print("Write '"+output_file+"' ...")
fh = nc.Dataset( output_file, 'w', 'NETCDF4' )
# File attributes
FiAtt   = ([['description', 'Gauging station file created from '+input_file],
            ['history'    , 'Created by Juliane Mai'],
            ['Conventions', 'CF-1.6'],
            ['featureType', 'timeSeries']])
handle  = writenetcdf(fh, fileattributes=FiAtt)

# dummy arrays for dimensions and dimension lengths
dim_name     = np.array([ 'dimmmmmmmm_'+str(ii) for ii in range(100) ])  # ['dim_0', 'dim_1', 'dim_2', ...., 'dim_99']
dims         = np.zeros(np.shape(dim_name),dtype=int)
dd = 0

# ------------------
# station IDs
# ------------------
nstations = len(station_id)

varName      = 'station_id'
varAtt       = ([['long_name', 'station ID'],
                 ['units',     '1'],
                 ['cf_role',   'timeseries_id']])
dims[dd]     = nstations
dim_name[dd] = 'nstations'

# create variable
dh = fh.createDimension(dim_name[dd], dims[dd])
vh = fh.createVariable(varName, str, tuple(['nstations']), zlib=True)
for ii in range(len(varAtt)):
    vh.setncattr(varAtt[ii][0], varAtt[ii][1])
    
# set variable values
vh[:] = np.array(station_id)

dd += 1

# ------------------
# station info
# ------------------
nstations = len(station_id)

varName      = 'station_info'
varAtt       = ([['long_name', 'station long information'],
                 ['units',     '1']])
dims[dd]     = nstations
dim_name[dd] = 'nstations'

# create variable
vh = fh.createVariable(varName, str, tuple(['nstations']), zlib=True)
for ii in range(len(varAtt)):
    vh.setncattr(varAtt[ii][0], varAtt[ii][1])
    
# set variable values
vh[:] = np.array(station_info)

dd += 1

# only if additional info about gauges is given
if (not(gaugeinfo_file is None)):
    
    # ------------------
    # latitudes
    # ------------------
    varName    = 'lat'
    attributes = {"long_name":     "latitude",
                  "standard_name": "latitude",
                  "units":         "degrees_north"}
        
    # set variable values
    # gaugeinfo = [ID,Name,Lat,Lon,Country,Drainage_area]
    lat   = np.array([ np.float(gg[2]) for gg in gaugeinfo ], dtype=np.float32)

    arr_shape = list([len(lat)])
    idx = [ np.where(dims==arr_shape[ddd])[0][0] for ddd in np.arange(len(arr_shape)) ]

    # write values to dimension
    vh   = writenetcdf(fh, name=varName, dims=dim_name[idx],
                       comp=True,
                       attributes=attributes, 
                       vartype=lat.dtype)
    writenetcdf( fh, vh, var = lat)

    # ------------------
    # longitudes
    # ------------------
    varName    = "lon"
    attributes = {"long_name":     "longitude",
                  "standard_name": "longitude",
                  "units":         "degrees_east"}
        
    # set variable values
    # gaugeinfo = [ID,Name,Lat,Lon,Country,Drainage_area]
    lon   = np.array([ np.float(gg[3]) for gg in gaugeinfo ], dtype=np.float32)

    arr_shape = list([len(lon)])
    idx = [ np.where(dims==arr_shape[ddd])[0][0] for ddd in np.arange(len(arr_shape)) ]

    # write values to dimension
    vh   = writenetcdf(fh, name=varName, dims=dim_name[idx],
                       comp=True,
                       attributes=attributes, 
                       vartype=lon.dtype)
    writenetcdf( fh, vh, var = lon)
    
    # ------------------
    # drainage area
    # ------------------
    varName    = "drainage_area"
    attributes = {"long_name":     "drainage area",
                  "units":         "km**2"}
        
    # set variable values
    # gaugeinfo = [ID,Name,Lat,Lon,Country,Drainage_area]
    area   = np.array([ np.float(gg[5]) for gg in gaugeinfo ], dtype=np.float32)

    arr_shape = list([len(area)])
    idx = [ np.where(dims==arr_shape[ddd])[0][0] for ddd in np.arange(len(arr_shape)) ]

    # write values to dimension
    vh   = writenetcdf(fh, name=varName, dims=dim_name[idx],
                       comp=True,
                       attributes=attributes, 
                       vartype=area.dtype)
    writenetcdf( fh, vh, var = area)

# ------------------
# time
# ------------------
# set attributes
varName     = 'time'
varAtt  = ([['units',         ref_date],
            ['calendar',      'gregorian'],
            ['standard_name', 'time']])
ntime        = np.shape(times_since_in_days)[0]
dims[dd]     = ntime
dim_name[dd] = 'time'

# create dimension; dims=None makes it unlimited
th = writenetcdf(fh, name=dim_name[dd], dims=None, attributes=varAtt, isdim=True)
# write values to dimension
writenetcdf(fh, th, time=list(range(ntime)), var=times_since_in_days)

dd += 1

# ------------------
# the discharge data itself
# ------------------

# set attributes
varName     = 'Q'
attributes = {'long_name': "discharge",
              'units':     'm**3 s**-1',
              '_FillValue': np.float32(nodata)}

model_data = np.array(model_data,np.float32)
model_data = np.transpose(model_data) # to match order of dimensions of observations

arr_shape = np.shape(model_data)
idx = [ np.where(dims==arr_shape[ddd])[0][0] for ddd in np.arange(len(arr_shape)) ]

# write values to dimension
vh   = writenetcdf(fh, name=varName, dims=dim_name[idx],
                       comp=True,
                       attributes=attributes, 
                       vartype=model_data.dtype)
writenetcdf( fh, vh, var = model_data)

# set global attributes
fh.setncattr('product', 'raster')
fh.setncattr('Conventions', 'CF-1.6')
fh.setncattr('License', 'These data are provided by the GWF funded IMPC project A5 - Hydrologic model inter-comparison and multi-model analysis for improved prediction. The data are under a GWF licence.')
fh.setncattr('Remarks', "This data were converted from "+input_file+" using convert_raw_to_netcdf.py")

# close netcdf
fh.close() 

