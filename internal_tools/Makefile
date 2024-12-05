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






all : simplexification_refelt


morley_base : morley_base.cc
	g++ -O3 -I ../src/getfem -I ../src/ -I ../src/gmm morley_base.cc -o morley_base ../src/.libs/libgetfem.a -lblas -lm

simplexification_refelt.o : simplexification_refelt.cc
	g++ -O3 -I ../src -I ../gcc/src simplexification_refelt.cc -c

simplexification_refelt : simplexification_refelt.o
	g++ simplexification_refelt.o -o simplexification_refelt ../gcc/src/.libs/libgetfem.a -lqhull -lblas -lg2c

