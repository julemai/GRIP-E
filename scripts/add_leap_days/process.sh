#! /bin/bash

# Copyright 2019 Juliane Mai - juliane.mai(at)uwaterloo.ca
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
#
# Adds leap days to all files (or just changes calendar if not leap year).
#
set -e
#
prog=$0
pprog=$(basename ${prog})
dprog=$(dirname ${prog})
isdir=${PWD}
pid=$$

add_leap=0
merge_files=1

datapath="/Users/j6mai/Documents/GWF/data/grip-gl_wfdei-gem-capa/"

# -------------------
# add leap years
# -------------------
if [[ ${add_leap} -eq 1 ]] ; then
    files=$( \ls ${datapath}20[0-9][0-9].nc )
    for oldfile in ${files} ; do
	
	newfile=$( echo ${oldfile} | rev | cut -d '.' -f 2- | rev )'_leap.nc'
	
	echo ''
	echo '------------------------------------------'
	echo ${newfile}
	echo '------------------------------------------'
	
	if [ -e ${newfile} ] ; then
	    rm ${newfile}
	fi
	
	python add_leap_days.py -i ${oldfile} -o ${newfile}
	
    done
fi

# -------------------
# merging all files together
# -------------------
if [[ ${merge_files} -eq 1 ]] ; then

    files=$( \ls ${datapath}20[0-9][0-9]'_leap.nc' )
    ncrcat ${files} ${datapath}/ grip-gl_wfdei-gem-capa_2000-2016_leap.nc

fi

exit 0



