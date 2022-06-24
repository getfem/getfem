#!/bin/sh
#
# Copyright (C) 2001-2020 Yves Renard
#
# This file is a part of GetFEM
#
# GetFEM  is  free software;  you  can  redistribute  it  and/or modify it
# under  the  terms  of the  GNU  Lesser General Public License as published
# by  the  Free Software Foundation;  either version 2.1 of the License,  or
# (at your option) any later version.
# This program  is  distributed  in  the  hope  that it will be useful,  but
# WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
# or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
# License for more details.
# You  should  have received a copy of the GNU Lesser General Public License
# along  with  this program;  if not, write to the Free Software Foundation,
# Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.




echo "running checks..."
echo "  srcdir='$srcdir'"

#export MATLABPATH=$(pwd)/../../src/matlab:$(pwd)/../../tests/matlab:${srcdir}:${srcdir}/../../src/matlab

#echo "  export MATLAB_ROOT=$MATLAB_ROOT"
#echo "  export MATLABPATH=$MATLABPATH"
#echo "  setenv LD_LIBRARY_PATH $LD_LIBRARY_PATH"
#echo "  setenv PATH $PATH"
#export PATH;

# s=$(echo "s=getenv('MATLABPATH'); while (length(s)), [a,s]=strtok(s,':'); addpath(a); end; disp(pwd); check_all; pause(1)" | ${MLAB} 2>&1);

s=$(echo "addpath('../../src/octave'); addpath('${srcdir}/../../src/octave'); disp(pwd); pause(1); check_all; pause(1);" | octave 2>&1);

# echo $s

k=`echo "$s" | grep "All tests succeeded"`;
if test x"$k" = x""; then
  echo "failure: "
  echo "$s" | grep FAIL
  exit 1;
else
  echo "$s" | grep -i "succe";
fi;

