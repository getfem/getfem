
$tmp = `../bin/createmp laplacian.param`;
# print "TMP = $tmp\n";
sub catch { `rm -f $tmp`; }
$SIG{INT} = 'catch';

open(TMPF, ">$tmp") or die "Open file impossible : $!\n";
print TMPF "N = 2;\n";
print TMPF "LX = 1.0\n";
print TMPF "LY = 1.0\n";
print TMPF "LZ = 1.0\n";
print TMPF "INCLINE = 0.0\n";
print TMPF "FT = 0.1\n";
print TMPF "MESH_TYPE = 0;\n";
print TMPF "K = 1;\n";
print TMPF "INTEGRATION = 0;\n";
print TMPF "NX = 10;\n";
print TMPF "RESIDU = 1E-9;\n";
print TMPF "FEM_TYPE = 0;\n"; 
print TMPF "ROOTFILENAME = 'laplacian';\n";
print TMPF "\n";
print TMPF "\n";
close(TMPF);


$er = 0;

sub start_program # (N, K, NX, OPTION, SOLVER)
{
  my $def   = $_[0];

  # print ("def = $def\n");

  open F, "laplacian $tmp $def 2>&1 |" or die;
  while (<F>) {
    if ($_ =~ /L2 error/) {
      ($a, $b) = split('=', $_);
      # print "La norme en question :", $b;
      if ($b > 0.01) { print "Error too large\n"; $er = 1; }
    }
    if ($_ =~ /error has been detected/) {
      $er = 1;
      print "============================================\n";
      print $_, <F>;
    }
 # 
 #   print $_;
  }
  `laplacian $tmp $def`;
  if ($?) { `rm -f $tmp`; exit(1); }
}

start_program("");
print ".";
start_program("-d N=1 -d NX=10 -d FT=1.0");
print ".";
start_program("-d N=3 -d NX=5");
print ".";
start_program("-d N=3 -d INTEGRATION=25 -d NX=5");
print ".";
start_program("-d K=2 -d NX=5");
print ".";
start_program("-d K=2 -d NX=5");
print ".";
start_program("-d INTEGRATION=12");
print ".";
start_program("-d INTEGRATION=15");
print ".";
start_program("-d INTEGRATION=17");
print ".";
start_program("-d INTEGRATION=1  -d MESH_TYPE=1");
print ".";
start_program("-d INTEGRATION=33 -d MESH_TYPE=1");
print ".";
start_program("-d INTEGRATION=35 -d MESH_TYPE=1");
print ".";
start_program("-d N=3 -d INTEGRATION=1 -d MESH_TYPE=2 -d NX=5 -d FT=0.01");
print ".";
start_program("-d INTEGRATION=2 -d MESH_TYPE=1 -d INCLINE=0.5");
print ".\n";

`rm -f $tmp`;
if ($er == 1) { exit(1); }


