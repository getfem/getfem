$bin_dir = "$ENV{srcdir}/../bin";
$tmp = `$bin_dir/createmp elas.param`;

sub catch { `rm -f $tmp`; exit(1); }
$SIG{INT} = 'catch';

open(TMPF, ">$tmp") or die "Open file impossible : $!\n";
print TMPF <<
LX = .5;		% size in X.
LY = 1.0;	        % size in Y.
LZ = 2.0;		% size in Z.
P1 = 1.;	        % First elastic coefficient.
P2 = 1.;   	        % Second elastic coefficient.
P3 = -1.4;   	        % Third elastic coefficient.
LAW = 2;
FORCEX = 0          % Amplitude of the external force
FORCEY = 0;
FORCEZ = 0;
TORSION = 2.65; %3.14;
EXTENSION = 0.0;
MESH_TYPE = 'GT_PK(3,1)';         % linear triangles
NX = 4;            	          % space step.
NZ = 10;
MESH_NOISE = 0; % Set to one if you want to "shake" the mesh
FEM_TYPE = 'FEM_PK(3,2)';  % P1 for triangles
FEM_TYPE_P = 'FEM_PK(3,1)';  % P1 for triangles
DATA_FEM_TYPE = 'FEM_PK(3,1)'
INTEGRATION = 'IM_TETRAHEDRON(6)'
RESIDU = 1E-6;     	% residu for iterative solvers.
MAXITER = 15;
USE_MULTIPLIERS=0
NBSTEP = 50;
ROOTFILENAME = 'nonlinear_elastostatic';     % Root of data files.
VTK_EXPORT = 1 % export solution to a .vtk file ?
NOISY=2

;
close(TMPF);



$er = 0;
open F, "./nonlinear_elastostatic $tmp 2>&1 |" or die;
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


