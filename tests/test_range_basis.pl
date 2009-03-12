$bin_dir = "$ENV{srcdir}/../bin";
$tmp = `$bin_dir/createmp test_range_basis.param`;
#$tmp=toto;
sub catch { `rm -f $tmp`; exit(1); }
$SIG{INT} = 'catch';

open(TMPF, ">$tmp") or die "Open file $tmp impossible : $!\n";
print TMPF <<
LX = 1.0;		          % size in X.
LY = 1.0;	                  % size in Y.
LZ = 1.0;	                  % size in Z.
INCLINE = 0;                      % Incline of the mesh.
FT = 0.1;                         % parameter for the exact solution.
MESH_TYPE = 'GT_PK(2,1)';         % linear triangles
NX = 10;            	          % space step.
MESH_NOISED = 0;                  % Set to one if you want to "shake" the mesh
FEM_TYPE = 'FEM_PK(2,1)';         % P1 for triangles
MULT_FEM_TYPE = 'FEM_PK(2,2)';    % P2 for triangles
INTEGRATION = 'IM_TRIANGLE(6)';   % quadrature rule for polynomials up
                                  % to degree 6 on triangles
RESIDUAL = 1E-9;     	          % residu for conjugate gradient.
ROOTFILENAME = 'test_range_basis';       % Root of data files.

;
close(TMPF);


$er = 0;

sub start_program { # (N, K, NX, OPTION, SOLVER)

  my $def   = $_[0];

 # print "def = $def\n";

  open F, "./test_range_basis $tmp $def 2>&1 |" or die("test_range_basis not found");
  while (<F>) {
    if ($_ =~ /L2 error/) {
  #    print $_;
      ($a, $b) = split('=', $_);
      # print "La norme en question :", $b;
      if ($b > 0.01) { 
	print "\nError too large: $b\n"; 
	print "./test_range_basis $tmp $def 2>&1 failed\n";
	$er = 1; 
      }
    }
    if ($_ =~ /error has been detected/) {
      $er = 1;
      print "============================================\n";
      print $_, <F>;
    }
 # print $_;
  }
  close(F);
  if ($?) {
    #`rm -f $tmp`;
    print "./test_range_basis $tmp $def 2>&1 failed\n";
    exit(1);
  }
}

start_program("");

#`rm -f $tmp`;
if ($er == 1) { exit(1); }


