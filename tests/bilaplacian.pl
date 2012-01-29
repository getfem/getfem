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

$srcdir = "$ENV{srcdir}";
$bin_dir = "$srcdir/../bin";
$tmp = `$bin_dir/createmp elas.param`;

sub catch { `rm -f $tmp`; exit(1); }
$SIG{INT} = 'catch';

open(TMPF, ">$tmp") or die "Open file impossible : $!\n";
print TMPF <<""
LX = 1.0; LY = LX; LZ = LX;
MESH_NOISED = 0;
MESH_TYPE = 'GT_PK(2,1)';
NEWBRICKS=1;
INTEGRATION = 'IM_TRIANGLE(13)';
NX = 10;
RESIDUAL = 1E-9;
FEM_TYPE = 'FEM_ARGYRIS';
DIRICHLET_FEM_TYPE = 'FEM_PK(2,3)';
DATA_FEM_TYPE = 'FEM_PK(2, 5)';
ROOTFILENAME = 'bilaplacian';
DIRICHLET_VERSION = 0;
D = 1.0;
KL = 0;

;

close(TMPF);

sub start_program {
  my $def   = $_[0];

 # print "def = $def\n";

  my $h2err = "null";
  open F, "./bilaplacian $tmp $def 2>&1 |" or die("bilaplacian not found");
  while (<F>) {
    if ($_ =~ /H2 error/) {
      ($a, $h2err) = split('=', $_); 
      $h2err =~ s/\n//;
      #print "La norme en question :", $h2err;
    }
  }
  close(F);
  if ($?) {
    #`rm -f $tmp`; 
    print "./bilaplacian $tmp $def 2>&1 failed\n";
    exit(1);
  }
  return $h2err;
}

$err1 = start_program("");
if ($err1 > 0.027) {
  print "error too large: $err1\n"; exit(1);
}
print ".";
$err1 = start_program(" -d NX=4");
$err2 = start_program(" -d NX=8");

if ($err2 > $err1 / 1.6) {
  print "Convergence error: P1: $err1 $err2\n";
  exit(1);
}
print ".\n";

