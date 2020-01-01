#  Copyright (C) 1999-2020 Yves Renard
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

if QHULL
check_PROGRAMS = crack_mindlin crack_bilaplacian
else
check_PROGRAMS =
endif

CLEANFILES = 

if QHULL
crack_mindlin_SOURCES = crack_mindlin.cc
crack_bilaplacian_SOURCES = crack_bilaplacian.cc crack_bilaplacian_singularities.cc crack_bilaplacian_problem.cc crack_bilaplacian.h crack_bilaplacian_moment.cc crack_bilaplacian_tools.cc crack_bilaplacian_sif.cc
endif

AM_CPPFLAGS = -I$(top_srcdir)/src -I../../src
LDADD    = ../../src/libgetfem.la -lm @SUPLDFLAGS@

if QHULL
TESTS = crack_mindlin.pl
else
TESTS =
endif

EXTRA_DIST =		 		\
	crack_mindlin.pl		\
	crack_mindlin.param 		\
	crack_bilaplacian.param 	\
	crack_bilaplacian_mixed.param

#	serie.param \
#	serie_bis.param \
#

LOG_COMPILER = perl
