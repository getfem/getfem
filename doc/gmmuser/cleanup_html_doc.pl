# -*- perl -*-
# Copyright (C) 2001-2009 Yves Renard
#
# This file is a part of GETFEM++
#
# Getfem++  is  free software;  you  can  redistribute  it  and/or modify it
# under  the  terms  of the  GNU  Lesser General Public License as published
# by  the  Free Software Foundation;  either version 2.1 of the License,  or
# (at your option) any later version.
# This program  is  distributed  in  the  hope  that it will be useful,  but
# WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
# or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
# License for more details.
# You  should  have received a copy of the GNU Lesser General Public License
# along  with  this program;  if not, write to the Free Software Foundation,
# Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.

eval 'exec perl -S $0 "$@"'
  if 0;

open(CONTENTF, "gmmuser_2.html") or die "Open input file impossible : $!\n";

my $content = "";
my %hrefs=();
my @flist;
my $in_li=0;
while ($li = <CONTENTF>) {
  chomp($li);
#  if ($li=~/<li>/ || $li =~ /<\/ul>/) {
#    if ($in_li) { $li = "</li>\n".$li; } # close tags for hyperlatex..
#    $in_li = 1;
#  } elsif ($li =~ /<ul>/) { $in_li = 0; }
  if ($li=~/<ul>.*/ || $li=~/<li>.*/ || $li=~/<\/ul>.*/ || $li=~/<\/li>.*/) {
    $_ = $li;
    if (/href="(.*)"/) {
      my $fname = $1;
      if ($1 =~ /#/) {
      } else {
	push(@flist, "$fname");
      }
    }
    $_ = $li;
    if (/Contents/) {
    } else {
      $_ = $li;
      if (/<a/) {
	if (/\#/) {
	} else {
	  $href = $li; $href =~ s/.*href=\"([^"]+)\".*/href=\"\1\"/;
	  $title = $li; $title =~ s/<a(.*)>(.*)<\/a>/\2/;
	  $hrefs{$href} = $title;
	}
	$li =~ s/<a(.*)>(.*)<\/a>/<a title="\2"\1>\2<\/a>/;
      }
      #if (/<li>/) { $li .= "</li>"; } #.. hyperlatex claims to produce valid xhtml..
      $content .= "$li\n";
    }
  }
}
print $content;

sub transform_line {
  local($li) = $_[0];
  local($nextli) = $_[1];
  $_ = $li;

  if ($li =~ /using Hyperlatex v 2.6/) {
    $li.="modified with a perl script, cleaned up with tidy for xhtml conformance..\n";
  }

#  $li =~ s/rel=stylesheet/rel=\"stylesheet\"/g;
#  $li =~ s/(<a name=\"[^\"]*\")>/\1 \/>/g; # fix missing slash for <a name="..">
#  $li =~ s/<\/A>//g; # remove all </A> don't know where they come from ... brain dead hyperlatex ...
#  $li =~ s/<p>/<p \/>/g;
  
  # replace <font color="#dfd"> (not xhtml valid) with <span style="color:#dfd">
#  $li =~ s/<font color=\"/<span style=\"color:/g; $li =~ s/<\/font>/<\/span>/g;
  # do the same for <font size="+x">
#  $li =~ s/<font size=\"/<span style=\"font-size:/g;

#  $li =~ s/.css\" type=\"text\/css\">/.css\" type=\"text\/css\" \/>/; # fix missing slash for <link rel=stylesheet..>
  if (/<pre>/) { $inpre=1; }
  if (/<\/pre>/) { $inpre=1; }
  if ( $inpre == 1 && /^  / ) { $li = substr($li,2); } #hyperlatex insert 2 whitespaces in pre blocks
  if (/<\/head>/) {
    if ($prevfile) { print FOUT "<link rel=\"prev\" href=\"$prevfile\" />\n"; }
    if ($nextfile) { print FOUT "<link rel=\"next\" href=\"$nextfile\" />\n"; }
  }
  if (/<body>/) {
    print FOUT "<body>\n<div id=\"menu\">\n";
    print FOUT "<p><a href=\"http://home.gna.org/getfem/gmm_intro\"><img src=\"gmmlogo_small.png\" title=\"getfem documentation index\" alt=\"getfem documentation index\"></img></a></p>\n";
    print FOUT "<h1>Gmm++ User Documentation</h1>\n";
    print FOUT $content;
    print FOUT "</div><div id=\"content\">\n";
  } elsif (/<\/body>/) {
    print FOUT "</div>\n";
    print FOUT "<div id=\"navbar\">";
    if ($prevfile) { print FOUT "<a title=\"Prev\" href=\"$prevfile\">&lsaquo;</a>"; }
    if ($nextfile) { print FOUT "<a title=\"Next\" href=\"$nextfile\">&rsaquo;</a>"; }
    print FOUT "</div>\n";
    print FOUT "$li";
  } else {
    $_ = $nextli;
    if (/<\/pre>/) { #hyperlatex inserts a bad carriage return before its </pre>
      chomp($li);
    }
    print FOUT $li;
  }
}


#foreach $fname (@flist) {
for ($i=0; $i<@flist; $i=$i+1) {
  if ($i > 0) { $prevfile = $flist[$i-1]; }
  $fname = $flist[$i];
  if ($i < @flist-1) { $nextfile = $flist[$i+1]; }
  my $fnameout = "m-".$fname;
  print "doing file $fname\n";
  open(FIN, $fname) or die "Open input file impossible : $!\n";
  open(FOUT, ">$fnameout") or die "Open output file impossible : $!\n";
  $pli=<FIN>;
  $inpre = 0;
  while ($li = <FIN>) {
    transform_line($pli,$li);
    $pli = $li;
  }
  transform_line($pli,"");
  close(FIN); close(FOUT);
  system("tidy -q -clean < $fnameout > $fname; rm '$fnameout'");
  #rename ("$fnameout", "$fname") || die "Cannot rename --> $fnameout $fname $!\n";
}
