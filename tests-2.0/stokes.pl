# Copyright (C) 2001-2009 Yves Renard
#
# This file is a part of GETFEM++
#
# Getfem++  is  free software;  you  can  redistribute  it  and/or modify it
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

$bin_dir = "$ENV{srcdir}/../bin";
$tmp = `$bin_dir/createmp elas.param`;

sub catch { `rm -f $tmp`; exit(1); }
$SIG{INT} = 'catch';

open(TMPF, ">$tmp") or die "Open file impossible : $!\n";
print TMPF <<
LX = 1.0;		% size in X.
LY = 1.0;	        % size in Y.
LZ = 1.0;		% size in Z.
NU = 1.0;	        % Lamé coefficient.
MESH_TYPE = 'GT_QK(2,1)'; % linear rectangles
NX = 10;            	          % space step.
MESH_NOISED = 0; % Set to one if you want to "shake" the mesh
NOISY=0;
FEM_TYPE = 'FEM_QK(2,2)';
FEM_TYPE_P = 'FEM_QK(2,1)';
DATA_FEM_TYPE = 'FEM_QK(2,1)';
INTEGRATION = 'IM_GAUSS_PARALLELEPIPED(2,4)';
RESIDUAL = 1E-9;     	% residu for conjugate gradient.
ROOTFILENAME = 'stokes';     % Root of data files.
VTK_EXPORT = 0 % export solution to a .vtk file ?

;
close(TMPF);

$er = 0;
open F, "./stokes $tmp 2>&1 |" or die;
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


