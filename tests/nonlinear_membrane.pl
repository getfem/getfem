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
$tmp = `$bin_dir/createmp memb.param`;

sub catch { `rm -f $tmp`; exit(1); }
$SIG{INT} = 'catch';

open(TMPF, ">$tmp") or die "Open file impossible : $!\n";
print TMPF <<""
LX =50;
LY =25;
NX =10;
NY =6;
P0 =1500;
P1 =0.01;
P2 =1500;
P3 =0;
bdyUp=0;	
bdyLow=0;	
bdyLeft=0;	
bdyRight=0;
bdy_type=-1;
bdy_element1= 9;%431 ;
bdy_face1= 1;%0;
bdy_element2= 29;%432 ;
bdy_face2= 1;%1;
bdy_element3=235;% 435 ;
bdy_face3=-1;%0 ;
bdy_element4=218;% 436 ;
bdy_face4= -1;%1 ;
bdy_element5=219;%439 ;
bdy_face5= -1;%0;
bdy_element6= -1;
bdy_face6= -1;
bdy_element7= -1 ;
bdy_face7= -1 ;
bdy_element8= -1 ;
bdy_face8= -1 ;
dx=0;
dy=0;
dz=0;
opposite_bdy_reversed=0;
INITIAL_DISP=1;		
initialDispAmplitude=1e-1;
PRINT_CONVEXES=0;
PRINT_STRESSES=0;
PRETENSION_X=0;
PRETENSION_Y=0;
src_type=0;
FORCEX = 0;
FORCEY = 0;
FORCEZ = 10;
PUNCTUAL_DOF1=0;
PUNCTUAL_FORCE1=0;
PUNCTUAL_DOF2=0;
PUNCTUAL_FORCE2=0;
MESH_TYPE = 'GT_PK(2,1)';
MESH_NOISED = 0;
FEM_TYPE = 'FEM_PK(2,1)';
FEM_TYPE_P = 'FEM_PK(2,1)';  % P1 for triangles
VONMISES_FEM_TYPE = 'FEM_PK_DISCONTINUOUS(2,1)';
INTEGRATION = 'IM_TRIANGLE(3)';
RESIDUAL = 1E-8;
MAXITER = 15;
DIRICHLET_VERSION=0;
NBSTEP =10;
ROOTFILENAME = 'nonlinear_membrane';
VTK_EXPORT = 1;
NOISY=1;

;
close(TMPF);



$er = 0;
open F, "./nonlinear_membrane $tmp 2>&1 |" or die;
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


