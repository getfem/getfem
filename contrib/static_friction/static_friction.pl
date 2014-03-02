# Copyright (C) 2001-2012 Yves Renard
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
LX = 1.0;		% size in X.
LY = 1.0;	        % size in Y.
LZ = 1.0;		% size in Z.
%MU = 7700;	        % Lam%/1€Œiso8859-15é coefficient.
%LAMBDA = 11500;   	% Lam%/1€Œiso8859-15é coefficient.
MU = 5;
LAMBDA = 10;
FRICTION_COEF = 0.5;    % Friction coefficient.
% PG = 9810; 		% gravitation constante (on earth) (mm/s^2).
% PG = 1000000; 	% gravitation constante (on jupiter !) (mm/s^2).
PG=12000
RHO = 6e-6;     	% "realistic" density for steel
%%%%%   discretisation parameters  :                     	      %%%%%
MESH_TYPE = 'GT_PK(2,1)';         % linear triangles
% MESH_TYPE = 'GT_QK(3,1)';       % 
% MESH_TYPE = 'GT_PRISM(3,1)';    % 3D prisms
NX = 20;            	          % space step.
MESH_NOISE = 0;         % Set to one if you want to "shake" the mesh
RESIDUAL = 1E-9;     	% residual for Newton.
METHOD = 0;             % 0 = Newton.
			% 1 = genetic for 2D problem only.
			% 2 = Additive Schwarz Newton
NOISY = 3;
POPULATION = 100;       % Parameter for genetic algorithm
R = 100.0;              % Augmentation parameter
DIRICHLET = 1;          % 0 = no Dirichlet boundary
			% 1 = Dirichlet boundary on the top
			% 2 = Dirichlet boundary on the left
NEUMANN = 0;            % 0 = no non homogeneous Neumann Boundary
			% 1 = Non homogeneous Neumann Boudary on the top
NEUMANN_INTENSITY = -0.0;
DIRICHLET_RATIO = -0.1;  % parametre pour la condition de Dirichlet
CONTACT_CONDITION = 0;  % 0 = Condition almost conformal in u
			% 1 = Condition almost conformal in forces on contact
			%     boundary with FEM_TYPE_L for the multipliers
FEM_TYPE = 'FEM_PK(2, 1)';      % Main FEM
FEM_TYPE_L = 'FEM_PK(2, 1)';    % FEM fo the multipliers
%DATA_FEM_TYPE = 'FEM_PK(2,1)'; % must be defined for non-Lagrangian main FEM
INTEGRATION = 'IM_TRIANGLE(6)'; % Quadrature rule
% INTEGRATION = 'IM_GAUSS_PARALLELEPIPED(3,6)'; % Quadrature rule
MESHNAME='splx:';
% MESHNAME='meshes/donut_regulier_8_elements_288ddl.mesh';
% MESHNAME='donut_regulier_64_elements_1920ddl.mesh';
% MESHNAME='donut_regulier_512_elements_13824ddl.mesh';
% MESHNAME='donut_regulier_32_elements.mesh';
% MESHNAME='donut_regulier_72_elements.mesh';
% MESHNAME='donut_regulier_128_elements.mesh';
% MESHNAME='donut_regulier_200_elements.mesh';
% MESHNAME='donut_regulier_288_elements.mesh';
% MESHNAME='donut_regulier_392_elements.mesh';
% MESHNAME='donut_regulier_512_elements.mesh';
% MESHNAME='donut_regulier_648_elements.mesh';
% MESHNAME='donut_regulier_800_elements.mesh';
ROOTFILENAME = 'dynamic_friction';     % Root of data files.
DX_EXPORT = 0;
ML_EXPORT = 0;

;
close(TMPF);

$er = 0;
open F, "./static_friction $tmp 2>&1 |" or die;
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


