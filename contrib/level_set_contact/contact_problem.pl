# Copyright (C) 2012-2015 Yves Renard
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

$bin_dir = "$ENV{srcdir}/../../bin";
$tmp = `$bin_dir/createmp elas.param`;

sub catch { `rm -f $tmp`; exit(1); }
$SIG{INT} = 'catch';

open(TMPF, ">$tmp") or die "Open file impossible : $!\n";
print TMPF <<""
RESIDUAL = 1e-8;
N = 2;
NSTEP = 2;
APPLIED_DISP = -1;
LS_OFFSET = -0.12;
DIVxM = 4;
DIVyM = 8;
DIVzM = 1;
xM = 0.3;
yM = 0.92;
zM = 0.1;
LxM = 0.4;
LyM = 0.5;
LzM = 0.1;
APPROX_ORDER_MASTER = 1;
LM_INT_TYPE = 'IM_STRUCTURED_COMPOSITE(IM_GAUSS1D(2),4)';
INT_ORDER_MASTER = 2;
MESH_TYPE_MASTER = 'QK';
LAMBDA_MASTER = 110.0;
MU_MASTER = 70.0;
DIVxS = 20;
DIVyS = 
20;
DIVzS = 4;
xS = 0;
yS = 0;
zS = 0;
LxS = 1;
LyS = 1;
LzS = 0.3;
APPROX_ORDER_SLAVE = 1;
INT_ORDER_SLAVE = 2;
MESH_TYPE_SLAVE = 'QK';
LAMBDA_SLAVE = 0;
MU_SLAVE = 7;

;
close(TMPF);

$er = 0;
open F, "./test_contact $tmp 2>&1 |" or die;
while (<F>) {
  #print $_;
  if ($_ =~ /error has been detected/)
  {
    $er = 1;
    print "============================================\n";
    print $_, <F>;
  }
}
close(F); if ($?) { `rm -f $tmp`; exit(1); }
if ($er == 1) { `rm -f $tmp`; exit(1); }
`rm -f $tmp`;


