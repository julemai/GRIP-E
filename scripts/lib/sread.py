#!/usr/bin/env python
from __future__ import print_function
import numpy as np

__all__ = ['sread']

def sread(infile, nc=0, cname=None, skip=0, cskip=0, hskip=0, hstrip=True, separator=None,
          squeeze=False, reform=False, skip_blank=False, comment=None,
          fill=False, fill_value='', strip=None,
          header=False, full_header=False,
          transpose=False, strarr=False):
    """
        Read strings into string array from a file.
        Lines or columns can be skipped.
        Columns can also be picked specifically.
        Blank (only whitespace) and comment lines can be excluded.
        The header of the file can be read separately.

        This routines is exactly the same as fread but reads
        everything as strings except of floats.


        Definition
        ----------
        def sread(infile, nc=0, cname=None, skip=0, cskip=0, hskip=0, hstrip=True, separator=None,
                  squeeze=False, reform=False, skip_blank=False, comment=None,
                  fill=False, fill_value='', strip=None,
                  header=False, full_header=False,
                  transpose=False, strarr=False):


        Input
        -----
        infile         source file name


        Optional Input Parameters
        -------------------------
        nc           number of columns to be read (default: all (nc<=0))
                     nc can be a vector of column indeces,
                     starting with 0; cskip will be ignored then.
        cname        columns can alternatively be chosen by the values in the first header line;
                     must be iterable with strings.
        skip         number of lines to skip at the beginning of file (default: 0)
        cskip        number of columns to skip at the beginning of each line (default: 0)
        hskip        number of lines in skip that do not belong to header (default: 0)
        hstrip       If true strip header cells to match with cname (default: True)
        separator    column separator
                     If not given, columns separator are (in order):
                     comma, semi-colon, whitespace
        comment      line gets excluded if first character of line is in comment sequence
                     sequence can be e.g. string, list or tuple
        fill_value   value to fill in if not enough columns in line
                     and fill=True (default: '')
        strip        Strip strings with str.strip(strip) (default: None)


        Options
        -------
        squeeze      True:  2-dim array will be cleaned of degenerated
                            dimension, i.e. results in vector
                     False: array will be two-dimensional as read (default)
        reform       Same as squeeze.
        skip_blank   True:  continues reading after blank line
                     False: stops reading at first blank line (default)
        fill         True:  fills in fill_value if not enough columns in line
                     False: stops execution and returns None if not enough
                            columns in line (default)
        header       True:  header strings will be returned
                     False  numbers in file will be returned (default)
        full_header  True:  header is a string vector of the skipped rows
                     False: header will be split in columns,
                            exactly as the data, and will hold only the
                            selected columns (default)
        transpose    True:  column-major format output(0:ncolumns,0:nlines)
                     False: row-major format output(0:nlines,0:ncolumns) (default)
        strarr       True:  return as numpy array of strings
                     False: return as list


        Output
        ------
        Depending on options:
            list of strings if header=False
            2D-array of strings if strarr=True and header=False
            String array of file header if header=True
            String vector of file header if header=True and full_header=True


        Restrictions
        ------------
        If header=True then skip is counterintuitive because it is
          actally the number of header rows to be read. This is to
          be able to have the exact same call of the function, once
          with header=False and once with header=True.
        If fill=True, blank lines are not filled but are expected
          end of file.
        transpose=True has no effect on 1D output such as 1 header line


        Examples
        --------
        >>> # Create some data
        >>> filename = 'test.dat'
        >>> ff = open(filename, 'w')
        >>> ff.writelines('head1 head2 head3 head4\\n')
        >>> ff.writelines('1.1 1.2 1.3 1.4\\n')
        >>> ff.writelines('2.1 2.2 2.3 2.4\\n')
        >>> ff.close()

        >>> # Read sample file in different ways
        >>> # data
        >>> print(sread(filename, skip=1))
        [['1.1', '1.2', '1.3', '1.4'], ['2.1', '2.2', '2.3', '2.4']]
        >>> print(sread(filename, skip=2))
        [['2.1', '2.2', '2.3', '2.4']]
        >>> print(sread(filename, skip=1, cskip=1))
        [['1.2', '1.3', '1.4'], ['2.2', '2.3', '2.4']]
        >>> print(sread(filename, nc=2, skip=1, cskip=1))
        [['1.2', '1.3'], ['2.2', '2.3']]
        >>> print(sread(filename, nc=[1,3], skip=1))
        [['1.2', '1.4'], ['2.2', '2.4']]
        >>> print(sread(filename, nc=1, skip=1))
        [['1.1'], ['2.1']]
        >>> print(sread(filename, nc=1, skip=1, reform=True))
        ['1.1', '2.1']

        >>> # header
        >>> print(sread(filename, nc=2, skip=1, header=True))
        ['head1', 'head2']
        >>> print(sread(filename, nc=2, skip=1, header=True, full_header=True))
        ['head1 head2 head3 head4']
        >>> print(sread(filename, nc=1, skip=2, header=True))
        [['head1'], ['1.1']]
        >>> print(sread(filename, nc=1, skip=2, header=True, squeeze=True))
        ['head1', '1.1']
        >>> print(sread(filename, nc=1, skip=2, header=True, squeeze=True, strarr=True))
        ['head1' '1.1']
        >>> print(sread(filename, nc=1, skip=2, header=True, squeeze=True, transpose=True))
        ['head1', '1.1']

        >>> # skip blank lines
        >>> ff = open(filename, 'a')
        >>> ff.writelines('\\n')
        >>> ff.writelines('3.1 3.2 3.3 3.4\\n')
        >>> ff.close()
        >>> print(sread(filename, skip=1))
        [['1.1', '1.2', '1.3', '1.4'], ['2.1', '2.2', '2.3', '2.4']]
        >>> print(sread(filename, skip=1, skip_blank=True))
        [['1.1', '1.2', '1.3', '1.4'], ['2.1', '2.2', '2.3', '2.4'], ['3.1', '3.2', '3.3', '3.4']]
        >>> print(sread(filename, skip=1, strarr=True))
        [['1.1' '1.2' '1.3' '1.4']
         ['2.1' '2.2' '2.3' '2.4']]

        >>> print(sread(filename, skip=1, strarr=True, transpose=True))
        [['1.1' '2.1']
         ['1.2' '2.2']
         ['1.3' '2.3']
         ['1.4' '2.4']]
        >>> print(sread(filename, skip=1, transpose=True))
        [['1.1', '2.1'], ['1.2', '2.2'], ['1.3', '2.3'], ['1.4', '2.4']]

        >>> # skip comment lines
        >>> ff = open(filename, 'a')
        >>> ff.writelines('# First\\n')
        >>> ff.writelines('! Second second comment\\n')
        >>> ff.writelines('4.1 4.2 4.3 4.4\\n')
        >>> ff.close()
        >>> print(sread(filename, skip=1))
        [['1.1', '1.2', '1.3', '1.4'], ['2.1', '2.2', '2.3', '2.4']]
        >>> print(sread(filename, skip=1, skip_blank=True, comment='#'))
        [['1.1', '1.2', '1.3', '1.4'], ['2.1', '2.2', '2.3', '2.4'], ['3.1', '3.2', '3.3', '3.4'], ['!', 'Second', 'second', 'comment'], ['4.1', '4.2', '4.3', '4.4']]
        >>> print(sread(filename, skip=1, skip_blank=True, comment='#!'))
        [['1.1', '1.2', '1.3', '1.4'], ['2.1', '2.2', '2.3', '2.4'], ['3.1', '3.2', '3.3', '3.4'], ['4.1', '4.2', '4.3', '4.4']]
        >>> print(sread(filename, skip=1, skip_blank=True, comment=('#','!')))
        [['1.1', '1.2', '1.3', '1.4'], ['2.1', '2.2', '2.3', '2.4'], ['3.1', '3.2', '3.3', '3.4'], ['4.1', '4.2', '4.3', '4.4']]
        >>> print(sread(filename, skip=1, skip_blank=True, comment=['#','!']))
        [['1.1', '1.2', '1.3', '1.4'], ['2.1', '2.2', '2.3', '2.4'], ['3.1', '3.2', '3.3', '3.4'], ['4.1', '4.2', '4.3', '4.4']]

        >>> # Create some more data with escaped numbers
        >>> filename2 = 'test2.dat'
        >>> ff = open(filename2, 'w')
        >>> ff.writelines('"head1" "head2" "head3" "head4"\\n')
        >>> ff.writelines('"1.1" "1.2" "1.3" "1.4"\\n')
        >>> ff.writelines('2.1 nan Inf "NaN"\\n')
        >>> ff.close()
        >>> print(sread(filename2, skip=1, strarr=True, transpose=True, strip='"'))
        [['1.1' '2.1']
         ['1.2' 'nan']
         ['1.3' 'Inf']
         ['1.4' 'NaN']]

        >>> # Create some more data with an extra (shorter) header line
        >>> filename3 = 'test3.dat'
        >>> ff = open(filename3, 'w')
        >>> ff.writelines('Extra header\\n')
        >>> ff.writelines('head1 head2 head3 head4\\n')
        >>> ff.writelines('1.1 1.2 1.3 1.4\\n')
        >>> ff.writelines('2.1 2.2 2.3 2.4\\n')
        >>> ff.close()

        >>> print(sread(filename3, skip=2))
        [['1.1', '1.2', '1.3', '1.4'], ['2.1', '2.2', '2.3', '2.4']]
        >>> print(sread(filename3, skip=2, hskip=1))
        [['1.1', '1.2', '1.3', '1.4'], ['2.1', '2.2', '2.3', '2.4']]
        >>> print(sread(filename3, nc=2, skip=2, hskip=1, header=True))
        ['head1', 'head2']

        >>> # Create some more data with missing values
        >>> filename4 = 'test4.dat'
        >>> ff = open(filename4, 'w')
        >>> ff.writelines('Extra header\\n')
        >>> ff.writelines('head1,head2,head3,head4\\n')
        >>> ff.writelines('1.1,1.2,1.3,1.4\\n')
        >>> ff.writelines('2.1,,2.3,2.4\\n')
        >>> ff.close()

        >>> print(sread(filename4, skip=2))
        [['1.1', '1.2', '1.3', '1.4'], ['2.1', '', '2.3', '2.4']]
        >>> print(sread(filename4, skip=2, fill=True, fill_value='-1'))
        [['1.1', '1.2', '1.3', '1.4'], ['2.1', '-1', '2.3', '2.4']]

        >>> # cname
        >>> print(sread(filename, cname='head2', skip=1, skip_blank=True, comment='#!', squeeze=True))
        ['1.2', '2.2', '3.2', '4.2']
        >>> print(sread(filename, cname=['head1','head2'], skip=1, skip_blank=True, comment='#!'))
        [['1.1', '1.2'], ['2.1', '2.2'], ['3.1', '3.2'], ['4.1', '4.2']]
        >>> print(sread(filename, cname=['head1','head2'], skip=1, skip_blank=True, comment='#!', header=True))
        ['head1', 'head2']
        >>> print(sread(filename, cname=['head1','head2'], skip=1, skip_blank=True, comment='#!', header=True, full_header=True))
        ['head1 head2 head3 head4']
        >>> print(sread(filename, cname=['  head1','head2'], skip=1, skip_blank=True, comment='#!', hstrip=False))
        [['1.2'], ['2.2'], ['3.2'], ['4.2']]

        >>> # Clean up doctest
        >>> import os
        >>> os.remove(filename)
        >>> os.remove(filename2)
        >>> os.remove(filename3)
        >>> os.remove(filename4)


        License
        -------
        This file is part of the JAMS Python package.

        The JAMS Python package is free software: you can redistribute it and/or modify
        it under the terms of the GNU Lesser General Public License as published by
        the Free Software Foundation, either version 3 of the License, or
        (at your option) any later version.

        The JAMS Python package is distributed in the hope that it will be useful,
        but WITHOUT ANY WARRANTY; without even the implied warranty of
        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
        GNU Lesser General Public License for more details.

        You should have received a copy of the GNU Lesser General Public License
        along with the JAMS Python package (cf. gpl.txt and lgpl.txt).
        If not, see <http://www.gnu.org/licenses/>.

        Copyright 2009-2015 Matthias Cuntz


        History
        -------
        Written,  MC, Jul 2009
        Modified, MC, Feb 2012 - transpose
                  MC, Feb 2013 - ported to Python 3
                  MC, Nov 2014 - bug when nc is list and contains 0
                  MC, Nov 2014 - hskip
                  MC, Feb 2015 - no lif, nc can be tuple
                               - large re-write
                  MC, Nov 2017 - use range instead of np.arange for producing indexes
                  MC, Nov 2017 - cname, sname, file->infile, hstrip
    """
    #
    # Open file
    try:
        f = open(infile, 'r')
    except IOError:
        raise IOError('Cannot open file '+infile)
    #
    # Read header and Skip lines
    if hskip > 0:
        ihskip = 0
        while ihskip < hskip:
            tmp = f.readline().rstrip()
            ihskip += 1
    if skip > 0:
        head = ['']*(skip-hskip)
        iskip = 0
        while iskip < (skip-hskip):
            head[iskip] = f.readline().rstrip()
            iskip += 1
    #
    # read first line to determine nc and separator (if not set)
    split = -1
    while True:
        s = f.readline().rstrip()
        if len(s) == 0:
            if skip_blank:
                continue
            else:
                break
        if comment is not None:
            if (s[0] in comment): continue
        break
    if separator is None:
        sep = ','
        res = s.split(sep)
        nres = len(res)
        if nres == 1:
            sep = ';'
            res = s.split(sep)
            nres = len(res)
            if nres == 1:
                sep = None
                res = s.split(sep)
                nres = len(res)
    else:
        sep = separator
        res = s.split(sep)
        nres = len(res)
    #
    # Determine indices
    if nc != 0 and cname is not None:
        f.close()
        raise ValueError('nc and cname are mutually exclusive.')
    if cname is not None:
        # from first header line
        if (skip-hskip) <= 0:
            f.close()
            raise IOError('No header line left for choosing columns by name.')
        if not isinstance(cname, (list, tuple, np.ndarray)): cname = [cname]
        if hstrip: cname = [ h.strip() for h in cname ]
        hres = head[0].split(sep)
        if hstrip: hres = [ h.strip() for h in hres ]
        iinc = []
        for k in range(len(hres)):
            if hres[k] in cname: iinc.append(k)
        nnc = len(iinc)
    else:
        # from nc keyword
        if isinstance(nc, (list, tuple, np.ndarray)):
            nnc  = len(nc)
            iinc = tuple(nc)
        else:
            if nc <= 0:
                iinc = range(cskip,nres)
                nnc  = nres-cskip
            else:
                iinc = range(cskip,cskip+nc)
                nnc = nc
    miinc = max(iinc)
    #
    # Header
    if header:
        # Split header
        var = None
        if (skip-hskip) > 0:
            if full_header:
                var = head
            else:
                var = list()
                k = 0
                while k < (skip-hskip):
                    hres = head[k].split(sep)
                    nhres = len(hres)
                    if (miinc >= nhres) and (not fill):
                        f.close()
                        raise ValueError('Line has not enough columns to index: '+head[k])
                    null = line2var(hres, var, iinc, strip)
                    k += 1
                if (skip-hskip) == 1: var = var[0]
        f.close()
        if strarr:
            var = np.array(var, dtype=np.str)
            if transpose: var = var.T
            if squeeze or reform: var = var.squeeze()
            if fill: var = np.where(var=='', fill_value, var)
        else:
            if fill:
                var = [ [ fill_value if i=='' else i for i in row ] for row in var ]
            if squeeze or reform:
                maxi = max([ len(i) for i in var])
                if maxi==1: var = [ i[0] for i in var ]
            if transpose and isinstance(var[0], list):
                var = [list(i) for i in zip(*var)] # transpose
        return var
    #
    # Values - first line
    if (miinc >= nres) and (not fill):
        f.close()
        raise ValueError('Line has not enough columns to index: '+s)
    var = list()
    null = line2var(res, var, iinc, strip)
    #
    # Values - rest of file
    for line in f:
        s = line.rstrip()
        if len(s) == 0:
            if skip_blank:
                continue
            else:
                break
        if comment is not None:
            if (s[0] in comment): continue
        res = s.split(sep)
        nres = len(res)
        if (miinc >= nres) and (not fill):
            f.close()
            raise ValueError('Line has not enough columns to index: '+s)
        null = line2var(res, var, iinc, strip)
    f.close()
    if strarr:
        var = np.array(var, dtype=np.str)
        if transpose: var = var.T
        if squeeze or reform: var = var.squeeze()
        if fill: var = np.where(var=='', fill_value, var)
    else:
        if fill:
            var = [ [ fill_value if i=='' else i for i in row ] for row in var ]
        if squeeze or reform:
            maxi = max([ len(i) for i in var])
            if maxi==1: var = [ i[0] for i in var ]
        if transpose and isinstance(var[0], list):
            var = [list(i) for i in zip(*var)] # transpose

    return var


# Helper for append var with current line already splitted into list
def line2var(res, var, iinc, strip):
    nres = len(res)
    if strip is None:
        tmp = [res[i].strip('"').strip("'") for i in iinc if i < nres]
    elif not strip:
        tmp = [res[i] for i in iinc if i < nres]
    else:
        tmp = [res[i].strip(strip) for i in iinc if i < nres]
    rest = len([ i for i in iinc if i >= nres ])
    if rest > 0:
        tmp.extend(['']*rest)
    var.append(tmp)
    return


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)
