#  Copyright (C) 1999-2017 Yves Renard
#
#  This file is a part of GetFEM++
#
#  GetFEM++  is  free software;  you  can  redistribute  it  and/or modify it
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

check_PROGRAMS =                 \
        test_contact

CLEANFILES = \
	slave_0.vtk slave_1.vtk master_0.vtk master_2.vtk master_1.vtk \
        slave_2.vtk


test_contact_SOURCES = contact_problem.cpp test_contact.cpp contact_problem.h

AM_CPPFLAGS = -I$(top_srcdir)/src -I../../src
LDADD    = ../../src/libgetfem.la -lm @SUPLDFLAGS@ @BOOST_LIBS@

TESTS = contact_problem.pl

EXTRA_DIST = \
	contact_problem.pl              \
	contact_2D.param		\
	contact_3D.param

LOG_COMPILER = perl
