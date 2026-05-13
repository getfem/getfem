# Copyright (C) 2001-2026 Yves Renard
#
# This file is a part of GetFEM
#
# GetFEM  is  free software;  you  can  redistribute  it  and/or modify it
# under  the  terms  of the  GNU  Lesser General Public License as published
# by  the  Free Software Foundation;  either version 3 of the License,  or
# (at your option) any later version along with the GCC Runtime Library
# Exception either version 3.1 or (at your option) any later version.
# This program  is  distributed  in  the  hope  that it will be useful,  but
# WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
# or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
# License and GCC Runtime Library Exception for more details.
# You  should  have received a copy of the GNU Lesser General Public License
# along  with  this program.  If not, see https://www.gnu.org/licenses/.

$er = 0;

sub start_program
{
  my $def   = $_[0];

  # print ("def = $def\n");

  open F, "./test_assembly_assignment $def 2>&1 |" or die;
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

start_program("-d N=10 -d RESULT=0.04478308527");
print ".\n";
start_program("-d N=15 -d RESULT=0.02998533866");
print ".\n";

if ($er == 1) { exit(1); }


