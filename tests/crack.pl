$bin_dir = "$ENV{srcdir}/../bin";
$tmp = `$bin_dir/createmp elas.param`;

sub catch { `rm -f $tmp`; exit(1); }
$SIG{INT} = 'catch';

open(TMPF, ">$tmp") or die "Open file impossible : $!\n";
print TMPF <<
MU = 1.0;
LAMBDA = 1.0;
MESH_NOISE = 0;
MESH_TYPE = 'GT_PK(2,1)';
MIXED_PRESSURE=0;
INTEGRATION = 'IM_TRIANGLE(6)';
SIMPLEX_INTEGRATION = 'IM_TRIANGLE(6)';
NX = 16;
RESIDUE = 1E-9;
CUTOFF=0.3;
ADDITIONAL_CRACK = 0;
ENRICHMENT_OPTION = 2;
FEM_TYPE = 'FEM_PK(2,1)';
FEM_TYPE_P = 'FEM_PK_DISCONTINUOUS(2,0)';
DATA_FEM_TYPE = 'FEM_PK(2,1)';
RADIUS_ENR_AREA = 0.2;
SPIDER_RADIUS =  0.2;
SPIDER_NR =  3;
SPIDER_NTHETA = 5;
SPIDER_K=1;
ROOTFILENAME = 'crack';
DIRICHLET_VERSION = 2;
VTK_EXPORT = 0;

;
close(TMPF);



$er = 0;
open F, "./crack $tmp 2>&1 |" or die;
while (<F>) {
  # print $_;
  if ($_ =~ /H1 ERROR/) {
    ($a, $b) = split(':', $_);
    if ($b > 0.12) { print "\nError too large\n"; $er = 1; }
  }
  if ($_ =~ /L2 ERROR/) {
    ($a, $b) = split(':', $_);
    if ($b > 0.0025) { print "\nError too large\n"; $er = 1; }
  }
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


