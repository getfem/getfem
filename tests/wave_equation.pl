$bin_dir = "$ENV{srcdir}/../bin";
$tmp = `$bin_dir/createmp wave.param`;
#$tmp=toto;
sub catch { `rm -f $tmp`; exit(1); }
$SIG{INT} = 'catch';

open(TMPF, ">$tmp") or die "Open file $tmp impossible : $!\n";
print TMPF <<
LX = 1.0;		          % size in X.
LY = 1.0;	                  % size in Y.
LZ = 1.0;	                  % size in Z.
INCLINE = 0;                      % Incline of the mesh.
MESH_TYPE = 'GT_PK(2,1)';         % linear triangles
MESH_NOISED = 1;                  % Set to one if you want to "shake" the mesh
FEM_TYPE = 'FEM_PK(2,2)';         % P1 for triangles
INTEGRATION = 'IM_TRIANGLE(6)';   % quadrature rule for polynomials up
                                  % to degree 6 on triangles
RESIDUAL = 1E-9;     	          % residu for conjugate gradient.
ROOTFILENAME = 'wave';            % Root of data files.
DIRICHLET_VERSION = 1;      	  % 0 = With Lagrange multipliers
			    	  % 1 = penalization.
DIRICHLET_COEFFICIENT = 1E10;	  % Penalization coefficient.
FT = 2;                           % parameter for the exact solution.
C = 1.5;
EXPORT_SOLUTION = 0;
SQRTN = 1.4142135623730951;       % N = 2
% SQRTN = 1.7320508075688772;     % N = 3
T = 5/(2*SQRTN*FT*C);             % Final time (1.25 period)
DT = T/50;                        % Time step 
THETA = 0.5;
NX = 20;            	          % space step.
SCHEME = 1;                       % Time integration scheme
BETA = 0.25;   GAMMA = 0.5;       % Newmark scheme parameters

;
close(TMPF);


$er = 0;

sub start_program { # (N, K, NX, OPTION, SOLVER)

  my $def   = $_[0];

 # print "def = $def\n";

  open F, "./wave_equation $tmp $def 2>&1 |" or die("wave_equation not found");
  while (<F>) {
    if ($_ =~ /L2 error/) {
  #    print $_;
      ($a, $b) = split('=', $_);
      # print "La norme en question :", $b;
      if ($b > 0.004) {
	print "\nError too large: $b\n";
	print "./wave_equation $tmp $def 2>&1 failed\n";
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
    print "./wave_equation $tmp $def 2>&1 failed\n";
    exit(1);
  }
}

start_program("");
#print ".";
#start_program("-d 'MESH_TYPE=\"GT_PK(1,1)\"' -d 'FEM_TYPE=\"FEM_PK(1,2)\"' -d 'INTEGRATION=\"IM_GAUSS1D(3)\"' -d NX=10 -d FT=1.0");
#print ".";
#start_program("-d 'MESH_TYPE=\"GT_PK(3,1)\"' -d 'FEM_TYPE=\"FEM_PK(3,2)\"' -d 'INTEGRATION=\"IM_TETRAHEDRON(5)\"' -d NX=3 -d FT=0.01");
#print ".";
#start_program("-d 'MESH_TYPE=\"GT_PK(2,1)\"' -d 'FEM_TYPE=\"FEM_PK(2,2)\"' -d 'INTEGRATION=\"IM_TRIANGLE(4)\"' -d NX=5 -d GENERIC_DIRICHLET=0");
#print ".";
#start_program("-d 'INTEGRATION=\"IM_TRIANGLE(2)\"'");
#print ".";
#start_program("-d 'INTEGRATION=\"IM_TRIANGLE(19)\"'");
#print ".";
#start_program("-d 'MESH_TYPE=\"GT_QK(2,1)\"' -d 'FEM_TYPE=\"FEM_QK(2,1)\"' -d 'INTEGRATION=\"IM_NC_PARALLELEPIPED(2,2)\"'");
##start_program("-d INTEGRATION=1  -d MESH_TYPE=1");
#print ".";
#start_program("-d 'MESH_TYPE=\"GT_QK(2,1)\"' -d 'FEM_TYPE=\"FEM_QK(2,1)\"' -d 'INTEGRATION=\"IM_QUAD(3)\"'");
##start_program("-d INTEGRATION=33 -d MESH_TYPE=1");
#print ".";
#start_program("-d 'MESH_TYPE=\"GT_QK(2,1)\"' -d 'FEM_TYPE=\"FEM_QK(2,1)\"' -d 'INTEGRATION=\"IM_QUAD(17)\"'");
##start_program("-d INTEGRATION=35 -d MESH_TYPE=1");
#print ".";
#start_program("-d 'MESH_TYPE=\"GT_PRISM(3,1)\"' -d 'FEM_TYPE=\"FEM_PK_PRISM(3,1)\"' -d 'INTEGRATION=\"IM_NC_PRISM(3,2)\"' -d NX=3 -d FT=0.01 -d GENERIC_DIRICHLET=0");
##start_program("-d N=3 -d INTEGRATION=1 -d MESH_TYPE=2 -d NX=3 -d FT=0.01");
#print ".";

#`rm -f $tmp`;
if ($er == 1) { exit(1); }


