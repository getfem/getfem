$bin_dir = "$ENV{srcdir}/../bin";
$tmp = `$bin_dir/createmp elas.param`;

sub catch { `rm -f $tmp`; exit(1); }
$SIG{INT} = 'catch';

open(TMPF, ">$tmp") or die "Open file impossible : $!\n";
print TMPF <<
LX = 1.0;		% size in X.
LY = 1.0;	        % size in Y.
LZ = 1.0;		% size in Z.
NU = 1.0;	        % Lamé coefficient.
MESH_TYPE = 'GT_QK(2,1)'; % linear rectangles
NX = 10;            	          % space step.
MESH_NOISED = 0; % Set to one if you want to "shake" the mesh
NOISY=0;
FEM_TYPE = 'FEM_QK(2,2)';
FEM_TYPE_P = 'FEM_QK(2,1)';
DATA_FEM_TYPE = 'FEM_QK(2,1)';
INTEGRATION = 'IM_GAUSS_PARALLELEPIPED(2,4)';
RESIDUAL = 1E-9;     	% residu for conjugate gradient.
ROOTFILENAME = 'stokes';     % Root of data files.
VTK_EXPORT = 0 % export solution to a .vtk file ?

;
close(TMPF);

$er = 0;
open F, "./stokes $tmp 2>&1 |" or die;
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


