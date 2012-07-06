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



$srcdir = "$ENV{srcdir}";
$bin_dir = "$srcdir/../bin";
$tmp = `$bin_dir/createmp elas.param`;

sub catch { `rm -f $tmp`; exit(1); }
$SIG{INT} = 'catch';

open(TMPF, ">$tmp") or die "Open file impossible : $!\n";
print TMPF <<""
N = 2;
MU = 7700;	        % Lamé coefficient.
LAMBDA = 11500;   	% Lamé coefficient.
R = 1;                  % Augmentation parameter
FRICTION_COEF = 0.1;    % Friction coefficient.
NOISY = 1;
if (N == 2)
  MESHNAME1='meshes/disc_with_a_hole.mesh';
  MESHNAME2='structured:GT="GT_PK(2,1)";ORG=[-0.5,0];SIZES=[1,0.1];NSUBDIV=[20,2]';
  FEM_TYPE = 'FEM_PK(2, 2)';      % Main FEM
  MULT_FEM_TYPE = 'FEM_PK(2, 1)'; % FEM for multipiers
  DATA_FEM_TYPE = 'FEM_PK(2, 2)'; % must be defined for non-Lagrangian main FEM
  INTEGRATION = 'IM_TRIANGLE(4)'; % Quadrature rule
end;
ROOTFILENAME = 'large_sliding_contact';     % Root of data files.

;

close(TMPF);





$er = 0;
open F, "./test_large_sliding_contact $tmp 2>&1 |" or die;
while (<F>) {
  # print $_;
  if ($_ =~ /error has been detected/)
  {
    $er = 1;
    print " =============================================================\n";
    print $_, <F>;
  }
}
close(F); if ($?) { exit(1); }
if ($er == 1) { exit(1); }


