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

if METIS
optprogs = xfem_dirichlet
else
optprogs = 
endif

check_PROGRAMS = xfem_contact xfem_stokes $(optprogs)

CLEANFILES = 

xfem_contact_SOURCES = xfem_contact.cc
if METIS
xfem_dirichlet_SOURCES = xfem_dirichlet.cc
endif
xfem_stokes_SOURCES = xfem_stokes.cc

AM_CPPFLAGS = -I$(top_srcdir)/src -I../../src
LDADD    = ../../src/libgetfem.la -lm @SUPLDFLAGS@

TESTS = xfem_contact.pl

EXTRA_DIST = \
	xfem_contact.pl                  \
	xfem_contact.param               \
	xfem_dirichlet.param             \
	xfem_stokes.param

LOG_COMPILER = perl
