# Copyright (C) 2001-2012 Yves Renard
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

$bin_dir = "$ENV{srcdir}/../bin";

$er = 0;

sub start_program 
{
  my $def   = $_[0];

  # print ("def = $def\n");

  open F, "./test_assembly $def 2>&1 |" or die;
  while (<F>) {
    if ($_ =~ /FAILED/) {
      $er = 1;
      print "============================================\n";
      print $_, <F>;
    }
 # 
 #   print $_;
  }
  close(F); if ($?) { exit(1); }
}

start_program("-d NX=8 -d NDIM=2");
print ".";
start_program("-d NX=3 -d NDIM=3 -d K=1");
print ".";
start_program("-d NX=6 -d NDIM=2 -d Kdata=3");
print ".\n";

if ($er == 1) { exit(1); }


