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

$bin_dir = "$ENV{srcdir}/../bin";
$tmp = `$bin_dir/createmp elas.param`;

sub catch { `rm -f $tmp`; print "error caught..\n"; exit(1); }
$SIG{INT} = 'catch';

open(TMPF, ">$tmp") or die "Open file impossible : $!\n";

print TMPF <<""
MU = 1.0;
LAMBDA = 1.0;
MESH_NOISED = 0;
MESH_TYPE = 'GT_PK(2,1)';
NX = 16;
MIXED_PRESSURE=0;
VECTORIAL_ENRICHMENT = 0;
INTEGRATION = 'IM_TRIANGLE(6)';
SIMPLEX_INTEGRATION = 'IM_TRIANGLE(6)';
SINGULAR_INTEGRATION = 'IM_STRUCTURED_COMPOSITE(IM_GAUSS_PARALLELEPIPED(2, 6), 9)';
FEM_TYPE = 'FEM_PK(2, 1)';
MORTAR_FEM_TYPE = FEM_TYPE;
DATA_FEM_TYPE = 'FEM_PK(2,1)';
INTEGRATION = 'IM_TRIANGLE(6)';
MODE = 2;
RESIDUAL = 1E-9;
CUTOFF_FUNC = 0;
CUTOFF=0.3;
CUTOFF1 = 0.2;
CUTOFF0 = 0.45;
ADDITIONAL_CRACK = 0;
ENRICHMENT_OPTION = 2;
FEM_TYPE = 'FEM_PK(2,1)';
FEM_TYPE_P = 'FEM_PK_DISCONTINUOUS(2,0)';
DATA_FEM_TYPE = 'FEM_PK(2,1)';
RADIUS_ENR_AREA = 0.2;
SPIDER_RADIUS =  0.2;
SPIDER_NR =  3;
SPIDER_NTHETA = 5;
SPIDER_K=1;
ROOTFILENAME = 'crack';
DIRICHLET_VERSION = 2;
VTK_EXPORT = 0;

;
close(TMPF);

$er = 0;
open F, "./crack $tmp 2>&1 |" or die "could not open $tmp\n";
while (<F>) {
  #print $_; #uncomment this line in case of problem..
  if ($_ =~ /H1 ERROR/) {
    ($a, $b) = split(':', $_);
    if ($b > 0.08) { print "\nError too large\n"; $er = 1; }
  }
  if ($_ =~ /L2 ERROR/) {
    ($a, $b) = split(':', $_);
    if ($b > 0.0015) { print "\nError too large\n"; $er = 1; }
  }
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


