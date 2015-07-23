# Copyright (C) 2001-2015 Yves Renard
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

$bin_dir = "$ENV{srcdir}/../bin";
$tmp = `$bin_dir/createmp elas.param`;

sub catch { `rm -f $tmp`; exit(1); }
$SIG{INT} = 'catch';

open(TMPF, ">$tmp") or die "Open file impossible : $!\n";
print TMPF <<""
MU = 1.0;	        % Lamé coefficient.
LAMBDA = 0.0;   	% Lamé coefficient.
EPSILON = 0.01;          % thickness of the plate
PRESSURE = 0.01;         % pressure on the top surface of the plate.
MESH_TYPE = 'GT_QK(2,1)'; % linear rectangles
NX = 10;            	          % space step.
MESH_NOISED = 0; % Set to one if you want to "shake" the mesh
FEM_TYPE_UT = 'FEM_QK(2,1)';
FEM_TYPE_U3 = 'FEM_QK(2,2)';
FEM_TYPE_THETA= 'FEM_QK(2,1)';
DATA_FEM_TYPE = 'FEM_QK(2,1)';
INTEGRATION = 'IM_GAUSS_PARALLELEPIPED(2,4)';
INTEGRATION_CT = 'IM_GAUSS_PARALLELEPIPED(2,4)';
RESIDUAL = 1E-9;     	% residu for conjugate gradient.
ROOTFILENAME = 'plate';     % Root of data files.
VTK_EXPORT = 0 % export solution to a .vtk file ?
MIXED = 0;
MITC = 0;
SYMMETRIZED = 0;
SOL_REF = 0;
ETA = 0.;
DX_EXPORT = 0;

;
close(TMPF);

$er = 0;
open F, "./plate $tmp 2>&1 |" or die;
while (<F>) {
  # print $_;
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


