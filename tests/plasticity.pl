# Copyright (C) 2001-2012 Yves Renard
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
print TMPF <<""
LX =100;   % size in X in mm.        %2.0; %1.0;		
LY =20;   % size in Y in mm.        %0.5; %1.0;
LZ =20;   % size in Z in mm.        %0.5;
MU = 80769.; % Lamé coefficient in N/mm^2.	      % 385 Lamé coefficient.
LAMBDA = 121150.;  % Lamé coefficient in N/mm^2.      % 330 pour plane_stress,
INCLINE = 0;            % Incline of the mesh.
MESH_TYPE = 'GT_PK(2,1)';         % linear triangles
NX =10 ; %5            	          % space step.
NY =10 ;
NZ =3 ;
MESH_NOISED = 0; % Set to one if you want to "shake" the mesh
FEM_TYPE = 'FEM_PK(2,1)';  % P1 for triangles
FEM_TYPE_SIGMA = 'FEM_PK(2,0)'; 
INTEGRATION = 'IM_TRIANGLE(6)'; % quadrature rule for polynomials up
GENERIC_DIRICHLET = 0;  % Generic Dirichlet condition for non-lagrangian elts.
ROOTFILENAME = 'plasticity';     % Root of data files.
STRESS_THRESHOLD =8000.;  % plasticity stress_threshold
RESIDUAL=1E-6;                      % RESIDUAL for iterative solvers
FLAG_HYP=0;     % option for the calculation hypothesis : 1 for stress plane
FORCE=330;

;
close(TMPF);



$er = 0;
open F, "./plasticity $tmp 2>&1 |" or die;
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


