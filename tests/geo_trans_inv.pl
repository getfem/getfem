# Copyright (C) 2001-2015 Yves Renard
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

$bin_dir = "$ENV{srcdir}/../bin";
$tmp = `$bin_dir/createmp geo.param`;

sub catch { `rm -f $tmp`; exit(1); }
$SIG{INT} = 'catch';

open(TMPF, ">$tmp") or die "Open file impossible : $!\n";
print TMPF "N = 2;\n";
print TMPF "LX = 1.0;\n";
print TMPF "LY = 1.0;\n";
print TMPF "LZ = 1.0;\n";
print TMPF "MESH_TYPE = 0;\n";
print TMPF "NX = 10;\n";
print TMPF "NB_POINTS = 100;\n";
print TMPF "BASE = 10;\n";
print TMPF "\n\n";
close(TMPF);

$er = 0;
open F, "./geo_trans_inv $tmp 2>&1 |" or die;
while (<F>) {
  # print $_;
  if ($_ =~ /error has been detected/)
  {
    $er = 1;
    print "=============================================================\n";
    print $_, <F>;
  }
}
close(F);
if ($?) { `rm -f $tmp`; exit(1); }
if ($er == 1) { `rm -f $tmp`; exit(1); }
`rm -f $tmp`;
