$bin_dir = "$ENV{srcdir}/../bin";
$tmp = `$bin_dir/createmp geo.param`;

sub catch { `rm -f $tmp`; exit(1); }
$SIG{INT} = 'catch';

open(TMPF, ">$tmp") or die "Open file impossible : $!\n";
print TMPF "N = 2;\n";
print TMPF "LX = 1.0;\n";
print TMPF "LY = 1.0;\n";
print TMPF "LZ = 1.0;\n";
print TMPF "MESH_TYPE = 0;\n";
print TMPF "NX = 10;\n";
print TMPF "NB_POINTS = 1000;\n";
print TMPF "BASE = 10;\n";
print TMPF "\n\n";
close(TMPF);

$er = 0;
open F, "geo_trans_inv $tmp 2>&1 |" or die;
while (<F>) {
  # print $_;
  if ($_ =~ /error has been detected/)
  {
    $er = 1;
    print "=============================================================\n";
    print $_, <F>;
  }
}
if ($er == 1) { exit(1); }
`geo_trans_inv $tmp`;
if ($?) { exit(1); }
