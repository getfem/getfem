$er = 0;
open F, "laplacien laplacien.param 2>&1 |" or die;
while (<F>) {
  if ($_ =~ /L2 error/) {
    ($a, $b) = split('=', $_);
    # print "La norme en question :", $b;
    if ($b > 0.1) { $er = 1; }
  }
  if ($_ =~ /error has been detected/) {
    $er = 1;
    print "=============================================================\n";
    print $_, <F>;
  }
  
  # print $_;
}
if ($er == 1) { exit(1); }
