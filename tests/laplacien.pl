open F, "laplacien laplacien.param |" or die;
while (<F>) {
  if ($_ =~ /L2 error/) {
    ($a, $b) = split('=', $_);
    print "La norme en question :", $b;
    if ($b > 0.1) { exit(1); }
  }
  print $_;
}

