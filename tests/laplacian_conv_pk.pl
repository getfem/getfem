
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
print TMPF "FT = 3.0\n";
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


$NDDLMAX = 10400;
$PAUSE = 0;
$SKIP = 0;
$FT = 10.0;

##########################################################################
print "   TESTS EN DIMENSION 1, ET ELEMENTS PK                         \n";
##########################################################################
$FEM_TYPE = 0;
$INTE = 0;
while ($INTE < 3 && $SKIP < 1) {
open(RES, ">$tmp_res");
$K = 1; $N = 1;  $NX = 1;
while ($NX**$N <= $NDDLMAX) {
  $K = 1;
  print "Test for NX = $NX \t"; print RES $NX**$N;
  while ((($K * $NX)**$N) * $K * $K <= 4*$NDDLMAX && $K <= 12) {
    start_program("-d N=$N -d NX=$NX -d K=$K -d FT=$FT -d INTEGRATION=$INTE -d FEM_TYPE=$FEM_TYPE");
    print RES "$linferror "; print ".";
    ++$K;
  }
  print RES "\n"; print "\n";
  if ($NX >= 5) { $NX = int($NX * 2); } else { ++$NX; }
}
close(RES);

open(GNF, ">$tmp_gnuplot");
print GNF "set data style line\n";
print GNF "set logscale\n";
print GNF "set xlabel 'number of dof'\n";
print GNF "set ylabel 'L-infinity error'\n";
print GNF "plot '$tmp_res' using (\$1):2 title 'PK(1,1)',";
print GNF "     '$tmp_res' using ((\$1)*2):3 title 'PK(1,2)',";
print GNF "     '$tmp_res' using ((\$1)*3):4 title 'PK(1,3)',";
print GNF "     '$tmp_res' using ((\$1)*4):5 title 'PK(1,4)',";
print GNF "     '$tmp_res' using ((\$1)*5):6 title 'PK(1,5)',";
print GNF "     '$tmp_res' using ((\$1)*6):7 title 'PK(1,6)',";
print GNF "     '$tmp_res' using ((\$1)*7):8 title 'PK(1,7)',";
print GNF "     '$tmp_res' using ((\$1)*8):9 title 'PK(1,8)',";
print GNF "     '$tmp_res' using ((\$1)*9):10 title 'PK(1,9)',";
print GNF "     '$tmp_res' using ((\$1)*10):11 title 'PK(1,10)',";
print GNF "     '$tmp_res' using ((\$1)*11):12 title 'PK(1,11)',";
print GNF "     '$tmp_res' using ((\$1)*12):13 title 'PK(1,12)'\n";
if ($PAUSE) { print GNF "pause -1;\n"; }
print GNF "set output 'laplacian_1D_$INTE.ps'\n";
print GNF "set term postscript color\n";
print GNF "replot\n";

close(GNF);
`gnuplot $tmp_gnuplot`;

$INTE += 1;
}


##########################################################################
print "   TESTS EN DIMENSION 1, ET ELEMENTS PK HIERARCHIQUES           \n";
##########################################################################
$FEM_TYPE = 2;
$INTE = 0;
while ($INTE < 3 && $SKIP < 2) {
open(RES, ">$tmp_res");
$K = 1; $N = 1; $NX = 1;
while ($NX**$N <= $NDDLMAX) {
  $K = 1;
  print "Test for NX = $NX \t"; print RES $NX**$N;
  while ((($K * $NX)**$N) * $K * $K <= 4*$NDDLMAX && $K <= 12) {
    start_program("-d N=$N -d NX=$NX -d K=$K -d FT=$FT -d INTEGRATION=$INTE -d FEM_TYPE=$FEM_TYPE");
    print RES "$linferror "; print ".";
    ++$K;
  }
  print RES "\n"; print "\n";
  if ($NX >= 5) { $NX = int($NX * 2); } else { ++$NX; }
}
close(RES);

open(GNF, ">$tmp_gnuplot");
print GNF "set data style line\n";
print GNF "set logscale\n";
print GNF "set xlabel 'number of dof'\n";
print GNF "set ylabel 'L-infinity error'\n";
print GNF "plot '$tmp_res' using (\$1):2 title 'PK(1,1)',";
print GNF "     '$tmp_res' using ((\$1)*2):3 title 'PK(1,2)',";
print GNF "     '$tmp_res' using ((\$1)*3):4 title 'PK(1,3)',";
print GNF "     '$tmp_res' using ((\$1)*4):5 title 'PK(1,4)',";
print GNF "     '$tmp_res' using ((\$1)*5):6 title 'PK(1,5)',";
print GNF "     '$tmp_res' using ((\$1)*6):7 title 'PK(1,6)',";
print GNF "     '$tmp_res' using ((\$1)*7):8 title 'PK(1,7)',";
print GNF "     '$tmp_res' using ((\$1)*8):9 title 'PK(1,8)',";
print GNF "     '$tmp_res' using ((\$1)*9):10 title 'PK(1,9)',";
print GNF "     '$tmp_res' using ((\$1)*10):11 title 'PK(1,10)',";
print GNF "     '$tmp_res' using ((\$1)*11):12 title 'PK(1,11)',";
print GNF "     '$tmp_res' using ((\$1)*12):13 title 'PK(1,12)'\n";
if ($PAUSE) { print GNF "pause -1;\n"; }
print GNF "set output 'laplacian_1D_hier_$INTE.ps'\n";
print GNF "set term postscript color\n";
print GNF "replot\n";

close(GNF);
`gnuplot $tmp_gnuplot`;

$INTE += 1;
}

##########################################################################
print "   TESTS EN DIMENSION 2, ET ELEMENTS PK                         \n";
##########################################################################
$NDDLMAX = 41600; $FT = 0.1;
$FEM_TYPE = 0;
$INTE = 0;
while ($INTE < 2 && $SKIP < 3) {
open(RES, ">$tmp_res");
$K = 1; $N = 2; $NX = 1;
while ($NX**$N <= $NDDLMAX) {
  $K = 1;
  print "Test for NX = $NX \t"; print RES $NX**$N;
  while ((($K * $NX)**$N) * $K <= 2*$NDDLMAX && $K <= 9) {
    start_program("-d N=$N -d NX=$NX -d K=$K -d FT=$FT -d INTEGRATION=$INTE -d FEM_TYPE=$FEM_TYPE");
    print RES "$linferror "; print ".";
    ++$K;
  }
  print RES "\n"; print "\n";
  if ($NX >= 5) { $NX = int($NX * 1.5); } else { ++$NX; }
}
close(RES);

open(GNF, ">$tmp_gnuplot");
print GNF "set data style line\n";
print GNF "set logscale\n";
print GNF "set xlabel 'number of dof'\n";
print GNF "set ylabel 'L-infinity error'\n";
print GNF "plot '$tmp_res' using (\$1):2 title 'PK(2,1)',";
print GNF "     '$tmp_res' using ((\$1)*4):3 title 'PK(2,2)',";
print GNF "     '$tmp_res' using ((\$1)*9):4 title 'PK(2,3)',";
print GNF "     '$tmp_res' using ((\$1)*16):5 title 'PK(2,4)',";
print GNF "     '$tmp_res' using ((\$1)*25):6 title 'PK(2,5)',";
print GNF "     '$tmp_res' using ((\$1)*36):7 title 'PK(2,6)',";
print GNF "     '$tmp_res' using ((\$1)*49):8 title 'PK(2,7)',";
print GNF "     '$tmp_res' using ((\$1)*64):9 title 'PK(2,8)',";
print GNF "     '$tmp_res' using ((\$1)*81):10 title 'PK(2,9)',";
print GNF "     '$tmp_res' using ((\$1)*100):11 title 'PK(2,10)',";
print GNF "     '$tmp_res' using ((\$1)*121):12 title 'PK(2,11)',";
print GNF "     '$tmp_res' using ((\$1)*144):13 title 'PK(2,12)'\n";
if ($PAUSE) { print GNF "pause -1;\n"; }
print GNF "set output 'laplacian_2D_$INTE.ps'\n";
print GNF "set term postscript color\n";
print GNF "replot\n";

close(GNF);
`gnuplot $tmp_gnuplot`;

$INTE += 1;
}


##########################################################################
print "   TESTS EN DIMENSION 2, ET ELEMENTS PK HIERARCHIQUES           \n";
##########################################################################
$FEM_TYPE = 2;
$INTE = 0;
while ($INTE < 2 && $SKIP < 4) {
open(RES, ">$tmp_res");
$K = 1; $N = 2; $FT = 1.0; $NX = 1;
while ($NX**$N <= $NDDLMAX) {
  $K = 1;
  print "Test for NX = $NX \t"; print RES $NX**$N;
  while ((($K * $NX)**$N) * $K <= 2*$NDDLMAX && $K <= 9) {
    start_program("-d N=$N -d NX=$NX -d K=$K -d FT=$FT -d INTEGRATION=$INTE -d FEM_TYPE=$FEM_TYPE");
    print RES "$linferror "; print ".";
    ++$K;
  }
  print RES "\n"; print "\n";
  if ($NX >= 5) { $NX = int($NX * 1.5); } else { ++$NX; }
}
close(RES);

open(GNF, ">$tmp_gnuplot");
print GNF "set data style line\n";
print GNF "set logscale\n";
print GNF "set xlabel 'number of dof'\n";
print GNF "set ylabel 'L-infinity error'\n";
print GNF "plot '$tmp_res' using (\$1):2 title 'PK(2,1)',";
print GNF "     '$tmp_res' using ((\$1)*4):3 title 'PK(2,2)',";
print GNF "     '$tmp_res' using ((\$1)*9):4 title 'PK(2,3)',";
print GNF "     '$tmp_res' using ((\$1)*16):5 title 'PK(2,4)',";
print GNF "     '$tmp_res' using ((\$1)*25):6 title 'PK(2,5)',";
print GNF "     '$tmp_res' using ((\$1)*36):7 title 'PK(2,6)',";
print GNF "     '$tmp_res' using ((\$1)*49):8 title 'PK(2,7)',";
print GNF "     '$tmp_res' using ((\$1)*64):9 title 'PK(2,8)',";
print GNF "     '$tmp_res' using ((\$1)*81):10 title 'PK(2,9)',";
print GNF "     '$tmp_res' using ((\$1)*100):11 title 'PK(2,10)',";
print GNF "     '$tmp_res' using ((\$1)*121):12 title 'PK(2,11)',";
print GNF "     '$tmp_res' using ((\$1)*144):13 title 'PK(2,12)'\n";
if ($PAUSE) { print GNF "pause -1;\n"; }
print GNF "set output 'laplacian_2D_hier_$INTE.ps'\n";
print GNF "set term postscript color\n";
print GNF "replot\n";

close(GNF);
`gnuplot $tmp_gnuplot`;

$INTE += 1;
}


















`rm -f $tmp $tmp_res $tmp_gnuplot`;

# `rm -f $tmp $tmp_gnuplot`;
