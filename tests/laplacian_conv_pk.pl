
# -*- perl -*-
eval 'exec perl -S $0 "$@"'
  if 0;

# Effectue un test de convergence pour des elements PK
# mettre bin_dir = ../bin ou ../../bin selon l'usage
$bin_dir = "../../bin";
$tmp = `$bin_dir/createmp laplacian.param`;
$tmp_res = `$bin_dir/createmp laplacian.res`;
$tmp_gnuplot = `$bin_dir/createmp laplacian.gnuplot`;

sub catch { `rm -f $tmp $tmp_res $tmp_gnuplot`; exit(1); }
$SIG{INT} = 'catch';

open(TMPF, ">$tmp") or die "Open file impossible : $!\n";
print TMPF "N = 2;\n";
print TMPF "LX = 1.0\n";
print TMPF "LY = 1.0\n";
print TMPF "LZ = 1.0\n";
print TMPF "INCLINE = 0.0\n";
print TMPF "FT = 30.0\n";
print TMPF "MESH_TYPE = 0;\n";
print TMPF "K = 1;\n";
print TMPF "KI = 1;\n";
print TMPF "INTEGRATION = 0;\n";
print TMPF "NX = 7;\n";
print TMPF "RESIDU = 1E-17;\n";
print TMPF "FEM_TYPE = 0;\n"; 
print TMPF "ROOTFILENAME = 'laplacian';\n";
print TMPF "\n\n";
close(TMPF);

$linferror = 0.0;

sub start_program # (N, K, NX, OPTION, SOLVER)
{
  my $def   = $_[0];

  # print ("def = $def\n");

  open F, "laplacian $tmp $def 2>&1 |" or die;
  while (<F>) {
    if ($_ =~ /Linfty error/) {
      ($a, $b) = split('=', $_);
      chomp $b;
      $linferror = $b;
    }
    if ($_ =~ /error has been detected/) {
      $er = 1;
      print "============================================\n";
      print $_, <F>;
    }
    # print $_;
  }
}

$INTE = 0;

$NMAX = 20400;

while ($INTE < 3) {

open(RES, ">$tmp_res");
$K = 1; $N = 1; $FT = 1.0; $NX = 1;
while ($NX <= $NMAX) {

  $K = 1;
  print "Test for NX = $NX \t";
  print RES 1.0/$NX;
  while ($K * $NX * $K * $K <= 4*$NMAX && $K <= 12) {

    start_program("-d N=$N -d NX=$NX -d K=$K -d FT=$FT -d INTEGRATION=$INTE");
    print RES "$linferror ";
    print ".";

    ++$K;
  }
  print RES "\n"; print "\n";
  if ($NX >= 5) { $NX = int($NX * 2); } else { ++$NX; }
}
close(RES);

open(GNF, ">$tmp_gnuplot");
print GNF "set data style line\n";
print GNF "set logscale\n";
print GNF "plot '$tmp_res' using (\$1):2,";
print GNF "     '$tmp_res' using ((\$1)/2):3,";
print GNF "     '$tmp_res' using ((\$1)/4):4,";
print GNF "     '$tmp_res' using ((\$1)/5):5,";
print GNF "     '$tmp_res' using ((\$1)/6):6,";
print GNF "     '$tmp_res' using ((\$1)/7):7,";
print GNF "     '$tmp_res' using ((\$1)/8):8,";
print GNF "     '$tmp_res' using ((\$1)/9):9,";
print GNF "     '$tmp_res' using ((\$1)/10):10,";
print GNF "     '$tmp_res' using ((\$1)/11):11,";
print GNF "     '$tmp_res' using ((\$1)/12):12\n";
print GNF "pause -1;\n";
print GNF "set output 'laplacian_$INTE.ps'\n";
print GNF "set term postscript\n";
print GNF "replot\n";

close(GNF);
`gnuplot $tmp_gnuplot`;

$INTE += 2;
}

`rm -f $tmp $tmp_res $tmp_gnuplot`;

# `rm -f $tmp $tmp_gnuplot`;
