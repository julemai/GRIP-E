#! /usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from netcdf4 import NcDataset
import netCDF4


def concatDimension(fnames, dim, dimvar=None, outfname=None, reference_date=None, **kwargs):
    """
    Arguments:
    ----------
    datasets (List[str]):          datasets to concatenate
    dim (String):                  name of the dimension to concatenate
    dimvar (Optional[String]):     name of the variable defining values
                                   for the concatenation dimension, 
                                   default: dimvar = dim
    outfname (Optional[str]):      file name of an output dataset
    kwargs (Optional[Any]):        paramaters to pass to the createVariable
                                   Method of NcDataset
    reference_date                 datetime object with reference date for time
                                   if not given first timestamp of first file is taken as reference
                                

    Return
    ------
    Optional[netcdf4.NcDataset]

    Purpose
    -------
    Concatenate the given datasets into a new one along the given dimensions
    'dim'. Values of the dimensions need to be given in 'dimvar' (default
    variable name is the same as the concatenation dimension).
    If 'outfname' is given the concatenation is file based, otherwise an in-memory
    dataset will be returned.
    """

    if reference_date is None:
        # if reference date not given take reference first time step of first file as reference
        with NcDataset(fnames[0], "r") as nc:
            time_var = nc.variables['time']
            reference_date = netCDF4.num2date(time_var[0],time_var.units)

    def getDimVarValues(fnames, vardim, reference_date=None):
        vals = []
        types = []
        for fname in fnames:
            with NcDataset(fname, "r") as nc:
                # TODO:
                # check if values are given, otherwise make up
                # some fake data simpling starting at 0
                if vardim == 'time':
                    time_var = nc.variables['time']
                    # if not reference date chosen take first date in first file
                    if reference_date is None:
                        reference_date = netCDF4.num2date(time_var[0],time_var.units)
                    # hours since reference date
                    data = np.array([ np.int(tt.days*24 + tt.seconds/60./60.) for tt in netCDF4.num2date(time_var[:],time_var.units) - reference_date ], dtype=time_var.dtype)
                else:
                    data = nc.variables[vardim][:]
                types.append(data.dtype)
                vals.extend(data.tolist())

        nvals_all    = len(vals)
        nvals_unique = len(np.unique(vals))

        if (nvals_all != nvals_unique):

            print("-----------------------------------------------------------")
            print("                    WARNING                                ")
            print("-----------------------------------------------------------")
            print("Be arware that there are timesteps that appear in multiple ")
            print("files (e.g., overlapping lead times of multiple forecasts).")
            print("These time steps will be overwritten by each other.        ")
            print("The value in the last file will be kept.                   ")
            print("-----------------------------------------------------------")
                
        return np.unique(vals).astype(dtype=np.find_common_type(types, []))
    
    def setConcatDimension(nc, dim, dimvar, values):
        nc.createDimension(dim, None)
        var = nc.createVariable(dimvar, type(values[0]), (dim,))
        var[:] = values

    def getCommonVarTypes(fnames):
        vars = {}
        for fname in fnames:
            with NcDataset(fname, "r") as nc:
                for vname, var in nc.variables.items():
                    try:
                        vars[vname].append(var.dtype)
                    except KeyError:
                        vars[vname] = [var.dtype]
        return {
            vname: np.find_common_type(vtypes, [])
            for vname, vtypes in vars.items()}
            
    if not dimvar:
        dimvar = dim

    dimvals = getDimVarValues(fnames, dimvar, reference_date=reference_date)
    dtypes = getCommonVarTypes(fnames)

    # the new concat dimension
    out = NcDataset(outfname, "w")
    out.createDimension(dim, None, fail=False)
    var = out.createVariable(dimvar, dimvals.dtype, (dim,), fail=False, **kwargs)
    var[:] = dimvals

    for i, fname in enumerate(fnames):

        nc = NcDataset(fname, "r")
        print ("file: " + fname)
        
        out.copyDimensions(nc.dimensions, skip=out.dimensions, fix=True)
        out.copyAttributes(nc.attributes)

        for vname, var in nc.variables.items():

            outvar = out.copyVariable(
                    var, dtype=dtypes[vname], data=False, fail=False, **kwargs
                )

            if dimvar == 'time' and vname == 'time':
                out.variables[vname].setncattr("units", 'hours since '+reference_date.strftime('%Y-%m-%d %H:%M:%S'))

            # build up the indices
            if dim in var.dimensions:
                if dimvar == 'time':
                    time_var = nc.variables['time']
                    vv = np.array([ np.int(tt.days*24 + tt.seconds/60./60.) for tt in netCDF4.num2date(time_var[:],time_var.units) - reference_date ], dtype=time_var.dtype)
                else:
                    vv = nc.variables[dimvar]
                idx = [np.where(dimvals == v)[0][0] for v in vv]
                slices = [slice(None,)] * len(var.dimensions)
                slices[var.dimensions.index(dim)] = idx
            else:
                # writing them again is stupid...
                slices = slice(None)

            if vname == 'time' and dimvar == 'time':
                data = dimvals[idx]
            else:
                data = var[:]

            # this seems to be necessary, because netCDF4 sometimes comes
            # up with a masked without a explicit user definition
            # (i.e. no fill_value, valid_range, whatever is set)
            if isinstance(data, np.ma.MaskedArray):
                data = data.data

            outvar[slices] = data

    if outfname is None:
        return out
    out.close()
