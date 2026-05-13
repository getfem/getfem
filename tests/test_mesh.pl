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
open F, "./test_mesh 2>&1 |" or die;
while (<F>) {
  # print $_;
    if ($_ =~ /error has been detected/) {
    $er = 1;
    print "=============================================================\n";
    print $_, <F>;
  }
}
close(F); if ($?) { exit(1); }
if ($er == 1) { exit(1); }

