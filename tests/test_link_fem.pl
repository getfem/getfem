$er = 0;
open F, "test_link_fem test_link_fem.param 2>&1 |" or die;
while (<F>) {
  # print $_;
    if ($_ =~ /error has been detected/) {
    $er = 1;
    print "=============================================================\n";
    print $_, <F>;
  }

}
if ($er == 1) { exit(1); }
`test_link_fem test_link_fem.param`;
if ($?) { exit(1); }
