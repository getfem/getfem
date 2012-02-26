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


$tmp = `$bin_dir/createmp test.param`;
# print "TMP = $tmp\n";
sub catch { `rm -f $tmp`; exit(1); }
$SIG{INT} = 'catch';

open(TMPF, ">$tmp") or die "Open file impossible : $!\n";
print TMPF "N = 2;\n";
print TMPF "LX = 1.0\n";
print TMPF "LY = 1.0\n";
print TMPF "LZ = 1.0\n";
print TMPF "K = 1;\n";
print TMPF "KI = 2;\n";
print TMPF "INTEGRATION = 15;\n";
print TMPF "NX1 = 7;\n";
print TMPF "NX2 = 6;\n";
print TMPF "\n\n";
close(TMPF);


$er = 0;
open F, "./test_interpolated_fem $tmp 2>&1 |" or die;
while (<F>) {
  # print $_;
    if ($_ =~ /error has been detected/) {
    $er = 1;
    print "=============================================================\n";
    print $_, <F>;
  }

}
close(F); if ($?) { `rm -f $tmp`; exit(1); }
if ($er == 1) { `rm -f $tmp`; exit(1); }
`rm -f $tmp`;
