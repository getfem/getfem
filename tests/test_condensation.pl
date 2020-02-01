# Copyright (C) 2019-2020 Yves Renard
#
# This file is a part of GetFEM++
#
# GetFEM++  is  free software;  you  can  redistribute  it  and/or modify it
# under  the  terms  of the  GNU  Lesser General Public License as published
# by  the  Free Software Foundation;  either version 3 of the License,  or
# (at your option) any later version along with the GCC Runtime Library
# Exception either version 3.1 or (at your option) any later version.
# This program  is  distributed  in  the  hope  that it will be useful,  but
# WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
# or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
# License and GCC Runtime Library Exception for more details.
# You  should  have received a copy of the GNU Lesser General Public License
# along  with  this program;  if not, write to the Free Software Foundation,
# Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.

$er = 0;

sub start_program
{
  my $def   = $_[0];

  # print ("def = $def\n");

  open F, "./test_condensation $def 2>&1 |" or die;
  while (<F>) {
    if ($_ =~ /FAILED/) {
      $er = 1;
      print "============================================\n";
      print $_, <F>;
    }
    print $_;
  }
  close(F); if ($?) { exit(1); }
}

start_program("-d NX=2 -d NX=3 -d DIFFICULTY=0");
print ".";
start_program("-d NX=2 -d NX=3 -d DIFFICULTY=1");
print ".";
start_program("-d NX=2 -d NX=3 -d DIFFICULTY=10");
print ".";
start_program("-d NX=2 -d NX=3 -d DIFFICULTY=11");
print ".";
start_program("-d NX=2 -d NX=3 -d DIFFICULTY=20");
print ".";
start_program("-d NX=2 -d NX=3 -d DIFFICULTY=21");
print ".";
start_program("-d NX=2 -d NX=3 -d DIFFICULTY=30");
print ".";
start_program("-d NX=2 -d NX=3 -d DIFFICULTY=31");
print ".";
start_program("-d NX=2 -d NX=3 -d DIFFICULTY=100");
print ".";
start_program("-d NX=2 -d NX=3 -d DIFFICULTY=101");
print ".";
start_program("-d NX=2 -d NX=3 -d DIFFICULTY=110");
print ".";
start_program("-d NX=2 -d NX=3 -d DIFFICULTY=111");
print ".";
start_program("-d NX=2 -d NX=3 -d DIFFICULTY=120");
print ".";
start_program("-d NX=2 -d NX=3 -d DIFFICULTY=121");
print ".";
start_program("-d NX=2 -d NX=3 -d DIFFICULTY=130");
print ".";
start_program("-d NX=2 -d NX=3 -d DIFFICULTY=131");
print ".";
start_program("-d NX=2 -d NX=3 -d DIFFICULTY=1000");
print ".";
start_program("-d NX=2 -d NX=3 -d DIFFICULTY=1001");
print ".";
start_program("-d NX=2 -d NX=3 -d DIFFICULTY=1010");
print ".";
start_program("-d NX=2 -d NX=3 -d DIFFICULTY=1011");
print ".";
start_program("-d NX=2 -d NX=3 -d DIFFICULTY=1020");
print ".";
start_program("-d NX=2 -d NX=3 -d DIFFICULTY=1021");
print ".";
start_program("-d NX=2 -d NX=3 -d DIFFICULTY=1030");
print ".";
start_program("-d NX=2 -d NX=3 -d DIFFICULTY=1031");
print ".";
start_program("-d NX=2 -d NX=3 -d DIFFICULTY=1100");
print ".";
start_program("-d NX=2 -d NX=3 -d DIFFICULTY=1101");
print ".";
start_program("-d NX=2 -d NX=3 -d DIFFICULTY=1110");
print ".";
start_program("-d NX=2 -d NX=3 -d DIFFICULTY=1111");
print ".";
start_program("-d NX=2 -d NX=3 -d DIFFICULTY=1120");
print ".";
start_program("-d NX=2 -d NX=3 -d DIFFICULTY=1121");
print ".";
start_program("-d NX=2 -d NX=3 -d DIFFICULTY=1130");
print ".";
start_program("-d NX=2 -d NX=3 -d DIFFICULTY=1131");
print ".\n";

if ($er == 1) { exit(1); }


