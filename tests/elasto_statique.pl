$bin_dir = "$ENV{srcdir}/../bin";
$tmp = `$bin_dir/createmp elas.param`;

sub catch { `rm -f $tmp`; exit(1); }
$SIG{INT} = 'catch';

open(TMPF, ">$tmp") or die "Open file impossible : $!\n";
print TMPF "N = 2;\n";
print TMPF "PG = 9.81;\n";
print TMPF "G = 1.0;\n";
print TMPF "LAMBDA = 1.0;\n";
print TMPF "LX = 1.0;\n";
print TMPF "LY = 1.0;\n";
print TMPF "LZ = 1.0;\n";
print TMPF "FT = 0.1;\n";
print TMPF "K = 2;\n";
print TMPF "MESH_TYPE = 0;\n";
print TMPF "KI = 1;\n";
print TMPF "INTEGRATION = 0;\n";
print TMPF "NX = 5;\n";
print TMPF "RESIDU = 1E-9;\n";
print TMPF "FEM_TYPE = 0;\n"; 
print TMPF "ROOTFILENAME = 'elasto_statique';\n";
print TMPF "\n\n";
close(TMPF);



$er = 0;
open F, "./elasto_statique $tmp 2>&1 |" or die;
while (<F>) {
  # print $_;
  if ($_ =~ /error has been detected/)
  {
    $er = 1;
    print "============================================\n";
    print $_, <F>;
  }
}
if ($er == 1) { exit(1); }
`./elasto_statique $tmp`;
if ($?) { exit(1); }

