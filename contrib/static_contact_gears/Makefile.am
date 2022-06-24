#  Copyright (C) 1999-2020 Yves Renard
#
#  This file is a part of GetFEM
#
#  GetFEM  is  free software;  you  can  redistribute  it  and/or modify it
#  under  the  terms  of the  GNU  Lesser General Public License as published
#  by  the  Free Software Foundation;  either version 3 of the License,  or
#  (at your option) any later version along with the GCC Runtime Library
#  Exception either version 3.1 or (at your option) any later version.
#  This program  is  distributed  in  the  hope  that it will be useful,  but
#  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
#  or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
#  License and GCC Runtime Library Exception for more details.
#  You  should  have received a copy of the GNU Lesser General Public License
#  along  with  this program;  if not, write to the Free Software Foundation,
#  Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.

# SUBDIRS = 

check_PROGRAMS = static_contact_gears static_contact_gears_u1_u2

CLEANFILES = 						 \
	static_contact_gears.U static_contact_gears.vtk  \
	static_contact_gears.RHS static_contact_gears.VM \
	static_contact_gears.mfd static_contact_gears.mf \
	static_contact_gears.mesh


static_contact_gears_SOURCES = static_contact_gears.cc
static_contact_gears_u1_u2_SOURCES = static_contact_gears_u1_u2.cc

AM_CPPFLAGS = -I$(top_srcdir)/src -I../../src
LDADD    = ../../src/libgetfem.la -lm @SUPLDFLAGS@

TESTS = static_contact_gears.pl

EXTRA_DIST = \
	gear1.msh \
	gear1_2D.msh \
	gear1_2teeth.msh \
	gear2.msh \
	gear2_2D.msh \
	gear2_2teeth.msh \
	static_contact_planetary_1.msh \
	static_contact_planetary_2.msh \
	static_contact_planetary_3.msh \
	static_contact_planetary_4.msh \
	static_contact_planetary_5.msh \
	static_contact_planetary.py \
	static_contact_gears.param \
	static_contact_gears.pl \
	static_contact_gears_2D.param \
	static_contact_gears_2teeth.param

LOG_COMPILER = perl
