$bin_dir = "$ENV{srcdir}/../bin";
$tmp = `$bin_dir/createmp schw.param`;

sub catch { `rm -f $tmp`; exit(1); }
$SIG{INT} = 'catch';

open(TMPF, ">$tmp") or die "Open file impossible : $!\n";
print TMPF <<
N = 3;                  % dimension.
PG = 9.81;		% gravity constant.
RHO = 0.1;     	        % mass density
MU = 1.0;	        % shear elastic stiffness.
LAMBDA = 1.0;   	% elastic stiffness.
LX = 1.0;		% size in X.
LY = 1.0;	        % size in Y.
LZ = 1.0;		% size in Z.
D = 0.00;		% Dirichlet condition.
%%%%%   discretisation parameters :        			      %%%%%
K = 1;         % Degree of the finite element method.
NX = 10;       % space step.
NXCOARSE = 4;  % space step for the coarse mesh.
USECOARSE = 1; % use a coarse mesh or not.
RESIDU = 1E-7;  %
SOLVER = 1;     % 0 = C.G.
		% 1 = additive Schwarz with global and local CG
		% 2 = additive Schwarz with global et local Gmres
CGRADIENT = 20;
NSDMX = 4;     % Nomber of sub-domains in x direction
NSDMY = 4;     % Nomber of sub-domains in y direction
NSDMZ = 4;     % Nomber of sub-domains in z direction
OVERLAP = 0.0; % overlap between sub-domains in %
MESHNAME = '';

;
close(TMPF);


$er = 0;
open F, "./schwarz_additive $tmp 2>&1 |" or die;
while (<F>) {
  # print $_;
  if ($_ =~ /error has been detected/)
  {
    $er = 1;
    print " =============================================================\n";
    print $_, <F>;
  }
}
close(F); if ($?) { exit(1); }
if ($er == 1) { exit(1); }



