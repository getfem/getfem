$bin_dir = "$ENV{srcdir}/../bin";
$tmp = `$bin_dir/createmp elas.param`;

sub catch { `rm -f $tmp`; exit(1); }
$SIG{INT} = 'catch';

open(TMPF, ">$tmp") or die "Open file impossible : $!\n";
print TMPF "MU = 1.0;\n";
print TMPF "LAMBDA = 1.0;\n";
print TMPF "MESH_NOISE = 0.0;\n";
print TMPF "MESH_TYPE = 'GT_PK(2,1)';\n";
print TMPF "MIXED_PRESSURE=0;\n";
print TMPF "INTEGRATION = 'IM_TRIANGLE(4)';\n";
print TMPF "SIMPLEX_INTEGRATION = 'IM_TRIANGLE(4)';\n";
print TMPF "NX = 5;\n";
print TMPF "RESIDU = 1E-9;\n";
print TMPF "FEM_TYPE = 'FEM_PK(2,1)';\n";
print TMPF "ROOTFILENAME = 'crack';\n";
print TMPF "\n\n";
close(TMPF);



$er = 0;
open F, "./crack $tmp 2>&1 |" or die;
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


