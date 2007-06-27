$srcdir = "$ENV{srcdir}";
$bin_dir = "$srcdir/../bin";
$tmp = `$bin_dir/createmp xfem.param`;

sub catch { `rm -f $tmp`; exit(1); }
$SIG{INT} = 'catch';



sub start_program {
  my $def = $_[0];
  my $nx = $_[1];

  open(TMPF, ">$tmp") or die "Open file impossible : $!\n";
  print TMPF <<
  RADIUS = .167;
  LEVEL_SET_DEGREE = 2;
  MASS_ELIMINATION = 1;
  ELIMINATION_THRESHOLD = 1e-8
  IM = 'IM_STRUCTURED_COMPOSITE(IM_TRIANGLE(6),1)';
  IM_SIMPLEX = 'IM_STRUCTURED_COMPOSITE(IM_TRIANGLE(6),1)';
  FEM = 'FEM_PK(2,2)';
  FEM_RHS = FEM;
  FEM_MULT = 'FEM_PK(2,1)';

    ;

  print TMPF "MESH_FILE=\'structured:GT=\"GT_PK(2,1)\";SIZES=[1,1];NOISED=0;NSUBDIV=[$nx,$nx]\';\n";

  close(TMPF);


 # print "def = $def\n";

  my $h1err = "null";
  open F, "./xfem_dirichlet_demo $tmp $def 2>&1 |" or die("xfem_dirichlet_demo not found");
  while (<F>) {
      # print $_;
      if ($_ =~ /H1 error/) {
      ($a, $h1err) = split(':', $_);
      $h1err += 0.0;
      # print "La norme en question :", $h1err;
    }
  }
  close(F);
  if ($?) {
    #`rm -f $tmp`; 
    print "./xfem_dirichlet_demo $tmp $def 2>&1 failed\n";
    exit(1);
  }
  return $h1err;
}

$err1 = start_program("", 15);
if ($err1 > 7e-4) {
  print "error too large\n"; exit(1);
}

$err2 = start_program("", 30);

if ($err2 > $err1 / 3.4) {
  print "Convergence error: P1: $err1 $err2\n";
  exit(1);
}


