$bin_dir = "$ENV{srcdir}/../bin";

$er = 0;

sub start_program 
{
  my $def   = $_[0];

  # print ("def = $def\n");

  open F, "test_assembly $def 2>&1 |" or die;
  while (<F>) {
    if ($_ =~ /FAILED/) {
      $er = 1;
      print "============================================\n";
      print $_, <F>;
    }
 # 
 #   print $_;
  }
}

start_program("-d NX=8 -d NDIM=2");
print ".";
start_program("-d NX=3 -d NDIM=3 -d K=1");
print ".";
start_program("-d NX=6 -d NDIM=2 -d Kdata=3");
print ".\n";

if ($er == 1) { exit(1); }


