$er = 0;
open F, "integration 2>&1 |" or die;
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
