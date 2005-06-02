$bin_dir = "$ENV{srcdir}/../../bin";
$tmp = `$bin_dir/createmp elas.param`;

sub catch { `rm -f $tmp`; exit(1); }
$SIG{INT} = 'catch';

open(TMPF, ">$tmp") or die "Open file impossible : $!\n";
print TMPF <<
LX = 1.0;
LY = 1.0;
LZ = 1.0;
NU = 1.0;
MESH_TYPE = 'GT_QK(2,1)';
NX = 6;
DT = 0.001;
T = 0.01;
MESH_NOISE = 0;
NOISY=0;
OPTION = 3;
FEM_TYPE = 'FEM_QK(2,2)';
FEM_TYPE_P = 'FEM_QK(2,1)';
DATA_FEM_TYPE = 'FEM_QK(2,1)';
INTEGRATION = 'IM_GAUSS_PARALLELEPIPED(2,4)';
RESIDU = 1E-9;
ROOTFILENAME = 'navier_stokes';
DX_EXPORT = 0;
DT_EXPORT = 0.01;

;
close(TMPF);

$er = 0;
open F, "./navier_stokes $tmp 2>&1 |" or die;
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


