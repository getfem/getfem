$bin_dir = "$ENV{srcdir}/../bin";
$tmp = `$bin_dir/createmp helm.param`;

sub catch { `rm -f $tmp`; exit(1); }
$SIG{INT} = 'catch';

open(TMPF, ">$tmp") or die "Open file impossible : $!\n";
print TMPF <<
WAVENUM_R = 3;
WAVENUM_I = 0;
R0 = 2.;
R1 = 10.;
GTDEGREE = 3;
NTHETA = 16;
NR = 10;
DIRICHLET_WITH_MULTIPLIERS = 0;
FEM_TYPE = 'FEM_QK(2,2)';
DATA_FEM_TYPE = 'FEM_QK(2,4)';
INTEGRATION = 'IM_GAUSS_PARALLELEPIPED(2,12)';
RESIDU = 1E-6;
ROOTFILENAME = 'helmholtz';
VTK_EXPORT = 2

;
close(TMPF);



$er = 0;
open F, "./helmholtz $tmp 2>&1 |" or die;
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


