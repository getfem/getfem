$srcdir = "$ENV{srcdir}";
$bin_dir = "$srcdir/../bin";
$tmp = `$bin_dir/createmp elas.param`;

sub catch { `rm -f $tmp`; exit(1); }
$SIG{INT} = 'catch';

open(TMPF, ">$tmp") or die "Open file impossible : $!\n";
print TMPF <<
MU = 1.0;
LAMBDA = 1.0;
LX = 1.0;
LY = 1.0;
LZ = 1.0;
INCLINE = 0.0;
MESH_NOISED = 1;
FT = 0.1;
MESH_TYPE = 'GT_PK(2,1)';
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

if ($err2 > $err1 / 3.6) {
  #maybe a problem with a neumann corner?
  print "Convergence error: P2: $err1 $err2\n";
  exit(1);
}
print ".";

$err1 = start_program(" -d NX=4 -d 'FEM_TYPE=\"FEM_PK(2,3)\"'");
$err2 = start_program(" -d NX=8 -d 'FEM_TYPE=\"FEM_PK(2,3)\"'");

if ($err2 > $err1 / 7.6) {
  print "Convergence error: P3: $err1 $err2\n";
  exit(1);
}

print ".";

$err1 = start_program(" -d 'FEM_TYPE=\"FEM_PK(2,2)\"' -d 'MESH_FILE=\"$srcdir/../tests/meshes/disc_P2_h8.mesh\"'");
$err2 = start_program(" -d 'FEM_TYPE=\"FEM_PK(2,2)\"' -d 'MESH_FILE=\"$srcdir/../tests/meshes/disc_P2_h4.mesh\"'");

if ($err2 > $err1 / 4) {
  print "Convergence error: P2 isoparametric: $err1 $err2\n";
  exit(1);
}

print ".\n";
