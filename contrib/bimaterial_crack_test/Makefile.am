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
	bimaterial_crack_test    \
	crack

CLEANFILES = 

bimaterial_crack_test_SOURCES = bimaterial_crack_test.cc getfem_Xfem.cc
crack_SOURCES = crack.cc crack_exact_solution.cc crack_exact_solution.h getfem_Xfem.cc getfem_spider_fem.h getfem_Xfem.h


AM_CPPFLAGS = -I$(top_srcdir)/src -I../../src 
LDADD    = ../../src/libgetfem.la -lm @SUPLDFLAGS@

TESTS = bimaterial_crack_test.pl \
        crack.pl

EXTRA_DIST = \
	bimaterial_crack_test.pl                  \
	bimaterial_crack_test.param          \
	crack.pl                  \
	crack.param 

LOG_COMPILER = perl
