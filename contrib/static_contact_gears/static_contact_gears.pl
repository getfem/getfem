# Copyright (C) 2012-2015 Yves Renard
#
# This file is a part of GetFEM++
#
# GetFEM++  is  free software;  you  can  redistribute  it  and/or modify it
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
$tmp = `$bin_dir/createmp test_static_contact_gears.param`;

sub catch {
`rm -f $tmp`;
exit(1);
}
$SIG{INT} = 'catch';

open(TMPF, ">$tmp") or die "Open file impossible : $!\n";
print TMPF <<""
MU = 0.83E+5;
LAMBDA = 1.18E+5;
ROT_ANGLE = -1.5E-2;
RESIDUAL = 1E-6;
FRICTION_COEFFICIENT = 0.0E+0;
MESHNAME_GEAR1 = 'gmsh:$ENV{srcdir}/gear1.msh';
MESHNAME_GEAR2 = 'gmsh:$ENV{srcdir}/gear2.msh';
CONTACT_FACES_1 = [113];
CONTACT_FACES_2 = [113];
DIRICHLET_FACES_1 = [133,142,143,173,182,183];
DIRICHLET_FACES_2 = [133,142,143,173,182,183];
FEM_TYPE = 'FEM_QK(3, 1)';
INTEGRATION = 'IM_HEXAHEDRON(5)';
ROOTFILENAME = 'static_contact_gears';
CONTACT_ALGO = 5;
MULT_FEM_TYPE = 'FEM_QK(3, 1)';

;
close(TMPF);

$er = 0;
open F, "./static_contact_gears $tmp 2>&1 |" or die;
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


