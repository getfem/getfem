# Copyright (C) 2015-2015 Yves Renard
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
$tmp = `$bin_dir/createmp elas.param`;

sub catch { `rm -f $tmp`; print "error caught..\n"; exit(1); }
$SIG{INT} = 'catch';

open(TMPF, ">$tmp") or die "Open file impossible : $!\n";

print TMPF <<""
epsilon = 1.;
E = 21E6;
nu = 0.3;
F = 100E2;
kappa = 4.;
D = 10.;
air_temp = 20.;
alpha_th = 16.6E-6;
T0 = 20.;
rho_0 = 1.754E-8;
alpha = 0.0039;
h = 4.;
elements_degree = 2;
export_mesh = 1;
solve_in_two_steps = 1;

;
close(TMPF);

$er = 0;
open F, "./thermo_elasticity_electrical_coupling $tmp 2>&1 |" or die "could not open $tmp\n";
while (<F>) {
  #print $_; #uncomment this line in case of problem..
  if ($_ =~ /L2 norm of temperature/) {
    ($a, $b) = split('=', $_);
    if (abs($b - 964) > 2) { print "\nWrong temperature norm\n"; $er = 1; }
  }
}
close(F); if ($?) { `rm -f $tmp`; exit(1); }
if ($er == 1) { `rm -f $tmp`; exit(1); }
`rm -f $tmp`;


