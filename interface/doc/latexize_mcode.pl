# -*- perl -*-
eval 'exec perl -S $0 "$@"'
  if 0;

# recognize strings
sub mark_strings {
  local($s) = $_[0];
  while ($s =~ /([[ ,(={]|^)'[^']*'/) {
    $s =~ s/'([^']*)'/\\str{\1}/;
  }
  # escape % in strings
  while ($s =~ /\\str{(.*[^\\])%(.*)}/) {
    $s =~ s/\\str{(.*[^\\])%(.*)}/\\str{\1\\%\2}/;
  }
  return $s;
}

sub mark_comments {
  local($s) = $_[0];
  if ($s =~ /[^\\]%/) {
    $s =~ s/([^\\])%(.*$)/\1\\mlabcomment{\\%\2}/;
  }
  return $s;
}

sub transforme {
  local($s, $mode) = @_;
  if ($mode != 0) {
    $s = &mark_strings($s);
    $s =~ s/@\{/\{/g;
    $s =~ s/@\}/\}/g;
  }
  if ($mode == 2 || $mode == 3) {
    if ($mode == 3) {
      if ($s =~ />>/) {
	$s =~ s/>>/\\mlabprompt/;
      } elsif (length($s)) {
	$s = "\\mlaboutput{".$s."}";
      }
    }
    if ($mode == 2) {
      ($s,$c,@b) = split('(%)',$s);
      while ($s =~ /[ \t,](function|error|try|catch|for|end)/) {
	$k = $1;
	$s =~ s/([ \t,])$k/\1\\mlabkeyword{$k}/;
      }
      $s = $s.$c;
      while (@b) { $s .= shift(@b); }
    }
    $s = "  ".$s;
    $s = mark_comments($s);
  }
  #recongnize getfem functions
  if ($mode != 0) {
    if ($s =~ /gf/) {
      local($tmp) = "";
      while ($s =~ /([^\'{\\]|^)gf/) {
	($bef, $rest) = split('gf', $s, 2);
	$rest =~ s/([a-zA-Z\\\_#]+)//;
	local($a) = "gf".$1;
	local($b) = $a; $b =~ s/\_//g; $b =~ s/\\//g;
	$tmp .= $bef."\\kwl{".$b."}{$a}";
	$s = $rest;
      }
      $s=$tmp.$s;
    }
  }
  return($s);
}

print "% file generated automatically -- do not edit this one\n";
my $cmode = 0;
while (my $ligne = <STDIN>) {
  chop($ligne);
  @t = split(/(##|@@|\\begin{mcode}|\\end{mcode}|\\begin{matlab}|\\end{matlab})/,$ligne);
  while (@t) {
    $tag = shift(@t);
#    print "\ntag = $tag";
#    print "TRANS[".&transforme($debut, $cmode)."]";
    if ($tag eq "##") {
      $ligne = shift(@t);
      $ligne =~ s/([a-zA-Z\\\_#]+)//;
      local($a) = $1;
      if ($a =~ /\\[^_]/) { # for case such as ##gf~mesh\\
	$a =~ s/(\\[^_].*)//; $ligne = $1.$ligne;
      }
      local($b) = $a; $b =~ s/\\_//g; $b =~ s/_//g;
      $ligne = "\\kwl{".$b."}{$a}".$ligne;
      print $ligne;
    } elsif ($tag eq "@@") {
      if ($cmode == 0) {
	$cmode = 1;
	print "\\inlinematlab{"
      } elsif ($cmode == 1) { 
	$cmode = 0;
	print "}";
      } else { die "imbricated @@ tags"; }
    } elsif ($tag eq "\\begin{mcode}") {
      print $tag;
      if ($cmode == 0) {
	$cmode = 2;
      } else { die "imbricated mcode tags"; }
    } elsif ($tag eq "\\end{mcode}") {
      print $tag;
      if ($cmode == 2) {
	$cmode = 0;
      } else { die "imbricated mcode tags"; }
    } elsif ($tag eq "\\begin{matlab}") {
      print $tag;
      if ($cmode == 0) {
	$cmode = 3;
      } else { die "imbricated matlab tags"; }
    } elsif ($tag eq "\\end{matlab}") {
      print $tag;
      if ($cmode == 3) {
	$cmode = 0;
      } else { die "imbricated mcode tags"; }
    } else {
      print &transforme($tag, $cmode);
    }
  }
  print "\n";
}

