# Copyright (C) 2015-2015 Yves Renard
#
# This file is a part of GETFEM++
#
# Getfem++  is  free software;  you  can  redistribute  it  and/or modify it
# under  the  terms  of the  GNU  Lesser General Public License as published
# by  the  Free Software Foundation;  either version 3 of the License,  or
# (at your option) any later version along with the GCC Runtime Library
# Exception either version 3.1 or (at your option) any later version.
# This program  is  distributed  in  the  hope  that it will be useful,  but
# WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
# or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
# License and GCC Runtime Library Exception for more details.
# You  should  have received a copy of the GNU Lesser General Public License
# along  with  this program;  if not, write to the Free Software Foundation,
# Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.

$bin_dir = "$ENV{srcdir}/../bin";
$tmp = `$bin_dir/createmp cont.param`;

sub catch { `rm -f $tmp`; print "error caught..\n"; exit(1); }
$SIG{INT} = 'catch';

open(TMPF, ">$tmp") or die "Open file impossible : $!\n";

print TMPF <<""
NX = 10;
FEM_TYPE = 'FEM_PK(1,1)';
INTEGRATION = 'IM_GAUSS1D(3)';
DATAPATH = 'data/';
IND_BRANCH = 2;
DIRECTION = 1.;
LAMBDA0 = 0.;
NBSTEP = 80;
SINGULARITIES = 2;
H_INIT = 2E-2;
H_MAX = 2E-1;
H_MIN = 2E-5;
H_INC = 1.3;
H_DEC = 0.5;
MAXITER = 5;
THR_ITER = 4;
RESIDUAL = 1E-6;
DIFFERENCE = 1E-6;
COS = 0.997;
RESIDUAL_SOLVE = 1E-8
NOISY = 1;

;
close(TMPF);

$er = 0;
open F, "./test_continuation $tmp 2>&1 |" or die "could not open $tmp\n";
while (<F>) {
  #print $_; #uncomment this line in case of problem..
  if ($_ =~ /smooth bifurcation point/) {
    ($a, $b) = split(',', $_);
    if (abs($b - 2) > 0) { print "\nWrong number of bifurcation points\n"; $er = 1; }
  }
}
close(F); if ($?) { `rm -f $tmp`; exit(1); }
if ($er == 1) { `rm -f $tmp`; exit(1); }
`rm -f $tmp`;


