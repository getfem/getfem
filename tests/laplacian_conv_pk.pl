
# -*- perl -*-
eval 'exec perl -S $0 "$@"'
  if 0;

# Effectue un test de convergence pour des elements PK
# mettre bin_dir = ../bin ou ../../bin selon l'usage
$bin_dir = "../../bin";
$tmp = `$bin_dir/createmp laplacian.param`;
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
print TMPF "GENERIC_DIRICHLET = 0;\n";
print TMPF "\n\n";
close(TMPF);

sub start_program # (N, K, NX, OPTION, SOLVER)
{
  my $def   = $_[0];

  $linferror = 100.0;

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
    if ($AFFICH) { print $_; }
  }
}


# $NDDLMAX = 20800;
$NDDLMAX = 4800;
$PAUSE = 0;
$SKIP = 2;
$GRAPHONLY=1;
$FT = 20.0;
$GENDIR = 0;
$AFFICH = 0;

@Ks=(1, 2, 3, 4, 6, 9, 12, 15, 18, 24);

##########################################################################
print "   TESTS EN DIMENSION 1, ET ELEMENTS PK                         \n";
##########################################################################
$FEM_TYPE = 0;
$INTE = 0;
while ($INTE < 3 && $SKIP < 1) {
  if (!($GRAPHONLY)) {
    open(RES, ">laplacian_1D_$INTE.res");
    $N = 1;  $NX = 1;
    while ($NX**$N <= $NDDLMAX) {
      print "Test for NX = $NX \t"; print RES $NX**$N;
      foreach $K (@Ks) {
	if ((($K * $NX)**$N) * $K <= 2*$NDDLMAX) {
	  start_program("-d N=$N -d NX=$NX -d K=$K -d FT=$FT -d INTEGRATION=$INTE -d FEM_TYPE=$FEM_TYPE");
	  print RES "$linferror "; print ".";
	}
      }
      print RES "\n"; print "\n";
      if ($NX >= 5) { $NX = int($NX * 2); } else { ++$NX; }
    }
    close(RES);
  }

open(GNF, ">$tmp_gnuplot");
print GNF "set data style line\n";
print GNF "set logscale\n";
print GNF "set xlabel 'number of dof'\n";
print GNF "set ylabel 'L-infinity error'\n";
print GNF "plot ";
$first = 0; $rank = 2;
foreach $K (@Ks) {
  if ($first) { print GNF ", "; }
  print GNF " 'laplacian_1D_$INTE.res' using ((\$1)*$K):$rank title 'PK(1,$K)'";
  $first = 1; ++$rank;
}
print GNF "\n";
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
  if (!($GRAPHONLY)) {
open(RES, ">laplacian_1D_hier_$INTE.res");
$K = 1; $N = 1; $NX = 1;
while ($NX**$N <= $NDDLMAX) {
  print "Test for NX = $NX \t"; print RES $NX**$N;
  foreach $K (@Ks) {
    if ((($K * $NX)**$N) * $K <= 2*$NDDLMAX) {
      start_program("-d N=$N -d NX=$NX -d K=$K -d FT=$FT -d INTEGRATION=$INTE -d FEM_TYPE=$FEM_TYPE");
      print RES "$linferror "; print ".";
    }
  }
  print RES "\n"; print "\n";
  if ($NX >= 5) { $NX = int($NX * 2); } else { ++$NX; }
}
close(RES);
}

open(GNF, ">$tmp_gnuplot");
print GNF "set data style line\n";
print GNF "set logscale\n";
print GNF "set xlabel 'number of dof'\n";
print GNF "set ylabel 'L-infinity error'\n";
print GNF "plot ";
$first = 0; $rank = 2;
foreach $K (@Ks) {
  if ($first) { print GNF ", "; }
  print GNF " 'laplacian_1D_hier_$INTE.res' using ((\$1)*$K):$rank title 'HIERARCHICAL_PK(1,$K)'";
  $first = 1; ++$rank;
}
print GNF "\n";
if ($PAUSE) { print GNF "pause -1;\n"; }
print GNF "set output 'laplacian_1D_hier_$INTE.ps'\n";
print GNF "set term postscript color\n";
print GNF "replot\n";

close(GNF);
`gnuplot $tmp_gnuplot`;

$INTE += 1;
}

##########################################################################
print "   TESTS EN DIMENSION 2, ET ELEMENTS PK                        \n";
##########################################################################
@Ks=(1, 2, 3, 4, 6, 9, 12, 15);
$NDDLMAX = 100000; $FT = 10.0;
$FEM_TYPE = 0;
$INTE = 0;
while ($INTE < 2 && $SKIP < 3) {
if (!($GRAPHONLY)) {
open(RES, ">laplacian_2D_$INTE.res");
$K = 1; $N = 2; $NX = 1;
while ($NX**$N <= $NDDLMAX) {
  print "Test for NX = $NX \t"; print RES $NX**$N;
  foreach $K (@Ks) {
    if ((($K * $NX)**$N) * $K <= 2*$NDDLMAX) {
      start_program("-d N=$N -d NX=$NX -d K=$K -d FT=$FT -d INTEGRATION=$INTE -d FEM_TYPE=$FEM_TYPE");
      print RES "$linferror "; print ".";
    }
  }
  print RES "\n"; print "\n";
  $NX = int($NX * 2.001);
}
close(RES);
}

open(GNF, ">$tmp_gnuplot");
print GNF "set data style line\n";
print GNF "set logscale\n";
print GNF "set xlabel 'number of dof'\n";
print GNF "set ylabel 'L-infinity error'\n";
print GNF "plot ";
$first = 0; $rank = 2;
foreach $K (@Ks) {
  if ($first) { print GNF ", "; }
  $KK = $K * $K;
  print GNF " 'laplacian_2D_$INTE.res' using ((\$1)*$KK):$rank title 'PK(2,$K)'";
  $first = 1; ++$rank;
}
print GNF "\n";
if ($PAUSE) { print GNF "pause -1;\n"; }
print GNF "set output 'laplacian_2D_$INTE.ps'\n";
print GNF "set term postscript color\n";
print GNF "replot\n";

close(GNF);
`gnuplot $tmp_gnuplot`;

$INTE += 1;
}


##########################################################################
print "   TESTS EN DIMENSION 2, ET ELEMENTS PK HIERARCHIQUES          \n";
##########################################################################
$NDDLMAX = 100000; $FT = 10.0;
$FEM_TYPE = 2;
$INTE = 0;
$GENDIR = 1;

while ($INTE < 2 && $SKIP < 3) {
if (!($GRAPHONLY)) {
open(RES, ">laplacian_2D_hier_$INTE.res");
$K = 1; $N = 2; $NX = 1;
while ($NX**$N <= $NDDLMAX) {
  print "Test for NX = $NX \t"; print RES $NX**$N;
  foreach $K (@Ks) {
    if ((($K * $NX)**$N) * $K <= 2*$NDDLMAX) {
      start_program("-d N=$N -d NX=$NX -d K=$K -d FT=$FT -d INTEGRATION=$INTE -d FEM_TYPE=$FEM_TYPE -d GENERIC_DIRICHLET=$GENDIR");
      print RES "$linferror "; print ".";
    }
  }
  print RES "\n"; print "\n";
  $NX = int($NX * 2.001);
}
close(RES);
}

open(GNF, ">$tmp_gnuplot");
print GNF "set data style line\n";
print GNF "set logscale\n";
print GNF "set xlabel 'number of dof'\n";
print GNF "set ylabel 'L-infinity error'\n";
print GNF "plot ";
$first = 0; $rank = 2;
foreach $K (@Ks) {
  if ($first) { print GNF ", "; }
  $KK = $K * $K;
  print GNF " 'laplacian_2D_hier_$INTE.res' using ((\$1)*$KK):$rank title 'HIERARCHICAL_PK(2,$K)'";
  $first = 1; ++$rank;
}
print GNF "\n";
if ($PAUSE) { print GNF "pause -1;\n"; }
print GNF "set output 'laplacian_2D_hier_$INTE.ps'\n";
print GNF "set term postscript color\n";
print GNF "replot\n";

close(GNF);
`gnuplot $tmp_gnuplot`;

$INTE += 1;
}


##########################################################################
print "   TESTS EN DIMENSION 3, ET ELEMENTS PK                        \n";
##########################################################################
@Ks=(1, 2, 3, 4, 6, 9);
$NDDLMAX = 100000; $FT = 2.0;
$FEM_TYPE = 0;
$INTE = 1;
while ($INTE < 2 && $SKIP < 4) {
if (!($GRAPHONLY)) {
open(RES, ">laplacian_3D_$INTE.res");
$K = 1; $N = 3; $NX = 1;
while ($NX**$N <= $NDDLMAX) {
  print "Test for NX = $NX \t"; print RES $NX**$N;
  foreach $K (@Ks) {
    if ((($K * $NX)**$N) * $K <= 2*$NDDLMAX) {
      start_program("-d N=$N -d NX=$NX -d K=$K -d FT=$FT -d INTEGRATION=$INTE -d FEM_TYPE=$FEM_TYPE");
      print RES "$linferror "; print ".";
    }
  }
  print RES "\n"; print "\n";
  $NX = int($NX * 2.001);
}
close(RES);
}

open(GNF, ">$tmp_gnuplot");
print GNF "set data style line\n";
print GNF "set logscale\n";
print GNF "set xlabel 'number of dof'\n";
print GNF "set ylabel 'L-infinity error'\n";
print GNF "plot ";
$first = 0; $rank = 2;
foreach $K (@Ks) {
  if ($first) { print GNF ", "; }
  $KK = $K * $K * $K;
  print GNF " 'laplacian_3D_$INTE.res' using ((\$1)*$KK):$rank title 'PK(3,$K)'";
  $first = 1; ++$rank;
}
print GNF "\n";
if ($PAUSE) { print GNF "pause -1;\n"; }
print GNF "set output 'laplacian_3D_$INTE.ps'\n";
print GNF "set term postscript color\n";
print GNF "replot\n";

close(GNF);
`gnuplot $tmp_gnuplot`;

$INTE += 1;
}

##########################################################################
print "   TESTS EN DIMENSION 4, ET ELEMENTS PK                        \n";
##########################################################################
@Ks=(1, 2, 3, 4, 6);
$FEM_TYPE = 0;
$INTE = 1;
while ($INTE < 2 && $SKIP < 4) {
if (!($GRAPHONLY)) {
open(RES, ">laplacian_4D_$INTE.res");
$K = 1; $N = 4; $NX = 1;
while ($NX**$N <= $NDDLMAX) {
  print "Test for NX = $NX \t"; print RES $NX**$N;
  foreach $K (@Ks) {
    if ((($K * $NX)**$N) * $K <= 2*$NDDLMAX) {
      start_program("-d N=$N -d NX=$NX -d K=$K -d FT=$FT -d INTEGRATION=$INTE -d FEM_TYPE=$FEM_TYPE");
      print RES "$linferror "; print ".";
    }
  }
  print RES "\n"; print "\n";
  $NX = int($NX * 2.001);
}
close(RES);
}

open(GNF, ">$tmp_gnuplot");
print GNF "set data style line\n";
print GNF "set logscale\n";
print GNF "set xlabel 'number of dof'\n";
print GNF "set ylabel 'L-infinity error'\n";
print GNF "plot ";
$first = 0; $rank = 2;
foreach $K (@Ks) {
  if ($first) { print GNF ", "; }
  $KK = $K * $K * $K * $K;
  print GNF " 'laplacian_4D_$INTE.res' using ((\$1)*$KK):$rank title 'PK(4,$K)'";
  $first = 1; ++$rank;
}
print GNF "\n";
if ($PAUSE) { print GNF "pause -1;\n"; }
print GNF "set output 'laplacian_4D_$INTE.ps'\n";
print GNF "set term postscript color\n";
print GNF "replot\n";

close(GNF);
`gnuplot $tmp_gnuplot`;

$INTE += 1;
}



`rm -f $tmp $tmp_gnuplot`;


