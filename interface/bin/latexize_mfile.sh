#!/bin/sh
#
# Copyright (C) 2004-2009 Yves Renard, Julien Pommier.
#                                                       
# This file is a part of GETFEM++                                         
#                                                                         
# GetFEM++  is  free software;  you  can  redistribute  it  and/or modify it
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
#


IN=$1
OUT=$2
bf=$(basename $IN)
TMP=$(mktemp /tmp/${bf}.XXXXXX) || exit 1
#TMP=./demo_laplacian.tmp
echo '\\begin{mcode}' > $TMP
sed -e 's/{/\\{/g' -e 's/}/\\}/g' < $IN >> $TMP
echo '\\end{mcode}' >> $TMP

perl $(dirname $0)/latexize_mcode.pl < $TMP > $OUT || ( rm $OUT; exit 1)
