$bin_dir = "$ENV{srcdir}/../bin";


$tmp = `$bin_dir/createmp test.param`;
# print "TMP = $tmp\n";
sub catch { `rm -f $tmp`; exit(1); }
$SIG{INT} = 'catch';

open(TMPF, ">$tmp") or die "Open file impossible : $!\n";
print TMPF "N = 2;\n";
print TMPF "LX = 1.0\n";
print TMPF "LY = 1.0\n";
print TMPF "LZ = 1.0\n";
print TMPF "K = 1;\n";
print TMPF "KI = 2;\n";
print TMPF "INTEGRATION = 15;\n";
print TMPF "NX1 = 7;\n";
print TMPF "NX2 = 6;\n";
print TMPF "\n\n";
close(TMPF);


$er = 0;
open F, "./test_link_fem $tmp 2>&1 |" or die;
while (<F>) {
  # print $_;
    if ($_ =~ /error has been detected/) {
    $er = 1;
    print "=============================================================\n";
    print $_, <F>;
  }

}
if ($er == 1) { `rm -f $tmp`; exit(1); }
`./test_link_fem $tmp`;
`rm -f $tmp`;
if ($?) { exit(1); }
