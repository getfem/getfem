# Copyright (C) 2001-2015 Yves Renard
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
MU = 1.0;
LAMBDA = 1.0;
FT = 0.1;
SOL_SING = 0;
REFINE = 1;
N = 2;
MESH_FILE='structured:GT="GT_PK(2,1)";SIZES=[1,1];NOISED=1';
MIXED_PRESSURE=0;
INTEGRATION = 'IM_TRIANGLE(13)';
NX = 30;
RESIDUAL = 1E-9;
FEM_TYPE = 'FEM_PK(2,1)';
ROOTFILENAME = 'elasto_static';
DIRICHLET_VERSION = 0;

;
close(TMPF);


sub start_program {
  my $def   = $_[0];

 # print "def = $def\n";

  my $h1err = "null";
  open F, "./elastostatic $tmp $def 2>&1 |" or die("elastostatic not found");
  while (<F>) {
    if ($_ =~ /H1 error/) {
      ($a, $h1err) = split('=', $_); 
      $h1err =~ s/\n//;
      #print "La norme en question :", $h1err;
    }
  }
  close(F);
  if ($?) {
    #`rm -f $tmp`; 
    print "./elastostatic $tmp $def 2>&1 failed\n";
    exit(1);
  }
  return $h1err;
}

$err1 = start_program("");
if ($err1 > 0.027) {
  print "error too large\n"; exit(1);
}
print ".";
$err1 = start_program(" -d NX=4 -d 'FEM_TYPE=\"FEM_PK(2,1)\"'");
$err2 = start_program(" -d NX=8 -d 'FEM_TYPE=\"FEM_PK(2,1)\"'");

if ($err2 > $err1 / 1.6) {
  print "Convergence error: P1: $err1 $err2\n";
  exit(1);
}
print ".";

$err1 = start_program(" -d NX=4 -d 'FEM_TYPE=\"FEM_PK(2,2)\"'");
$err2 = start_program(" -d NX=8 -d 'FEM_TYPE=\"FEM_PK(2,2)\"'");

if ($err2 > $err1 / 2.5) {
  #maybe a problem with a neumann corner?
  print "Convergence error: P2: $err1 $err2\n";
  exit(1);
}
print ".";

$err1 = start_program(" -d NX=4 -d 'FEM_TYPE=\"FEM_PK(2,3)\"'");
$err2 = start_program(" -d NX=8 -d 'FEM_TYPE=\"FEM_PK(2,3)\"'");

if ($err2 > $err1 / 5.9) {
  print "Convergence error: P3: $err1 $err2\n";
  exit(1);
}

print ".";

$err3 = start_program(" -d REFINE=1 -d NX=8 -d 'FEM_TYPE=\"FEM_PK(2,3)\"'");

if ($err3 > $err1) {
  print "Refinement error";
  exit(1);
}

print ".";

$err1 = start_program(" -d REFINE=0 -d 'FEM_TYPE=\"FEM_PK(2,2)\"' -d 'MESH_FILE=\"$srcdir/../tests/meshes/disc_P2_h4.mesh\"'");
$err2 = start_program(" -d REFINE=0 -d 'FEM_TYPE=\"FEM_PK(2,2)\"' -d 'MESH_FILE=\"$srcdir/../tests/meshes/disc_P2_h2.mesh\"'");

if ($err2 > $err1 / 4) {
  print "Convergence error: P2 isoparametric: $err1 $err2\n";
  exit(1);
}

`rm -f $tmp`;

print ".\n";
