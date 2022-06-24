#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Python GetFEM interface
#
# Copyright (C) 2009-2020 Luis Saavedra.
#
# This file is a part of GetFEM
#
# GetFEM  is  free software;  you  can  redistribute  it  and/or modify it
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
#
############################################################################
"""  test user global functions.

  This program is used to check that python-getfem is working. This is
  also a good example of use of python-getfem..

  $Id$
"""
import os

import getfem as gf

# mesh fem to export
m = gf.Mesh('triangles grid', [-1, -0.5, 0, 0.5, 1], [-1, -0.5, 0, 0.5, 1])
# m = gf.Mesh('import','gmsh','quad.msh')
mf = gf.MeshFem(m)
mf.set_fem(gf.Fem('FEM_PK(2,1)'))
PTs = mf.basic_dof_nodes()

# crack:
ck0  = gf.GlobalFunction('crack',0)
ck1  = gf.GlobalFunction('crack',1)
ck2  = gf.GlobalFunction('crack',2)
ck3  = gf.GlobalFunction('crack',3)
ck4  = gf.GlobalFunction('crack',4)
ck5  = gf.GlobalFunction('crack',5)
ck6  = gf.GlobalFunction('crack',6)
ck7  = gf.GlobalFunction('crack',7)
ck8  = gf.GlobalFunction('crack',8)
ck9  = gf.GlobalFunction('crack',9)
ck10 = gf.GlobalFunction('crack',10)
ck11 = gf.GlobalFunction('crack',11)
mf.export_to_pos( 'check_global_functions0.pos', ck0(PTs), 'ck0')
mf.export_to_pos( 'check_global_functions1.pos', ck1(PTs), 'ck1')
mf.export_to_pos( 'check_global_functions2.pos', ck2(PTs), 'ck2')
mf.export_to_pos( 'check_global_functions3.pos', ck3(PTs), 'ck3')
mf.export_to_pos( 'check_global_functions4.pos', ck4(PTs), 'ck4')
mf.export_to_pos( 'check_global_functions5.pos', ck5(PTs), 'ck5')
mf.export_to_pos( 'check_global_functions6.pos', ck6(PTs), 'ck6')
mf.export_to_pos( 'check_global_functions7.pos', ck7(PTs), 'ck7')
mf.export_to_pos( 'check_global_functions8.pos', ck8(PTs), 'ck8')
mf.export_to_pos( 'check_global_functions9.pos', ck9(PTs), 'ck9')
mf.export_to_pos('check_global_functions10.pos',ck10(PTs),'ck10')
mf.export_to_pos('check_global_functions11.pos',ck11(PTs),'ck11')

# cutoff:
co0 = gf.GlobalFunction('cutoff',-1,0.4,0.01,0.4)
co1 = gf.GlobalFunction('cutoff', 0,0.4,0.01,0.4)
co2 = gf.GlobalFunction('cutoff', 1,0.4,0.01,0.4)
co3 = gf.GlobalFunction('cutoff', 2,0.4,0.01,0.4)
mf.export_to_pos('check_global_functions12.pos',co0(PTs),'cutoff -1')
mf.export_to_pos('check_global_functions13.pos',co1(PTs),'cutoff  0')
mf.export_to_pos('check_global_functions14.pos',co2(PTs),'cutoff  1')
mf.export_to_pos('check_global_functions15.pos',co3(PTs),'cutoff  2')


# parser:
p0 = gf.GlobalFunction('parser','0')
p1 = gf.GlobalFunction('parser','1')
p2 = gf.GlobalFunction('parser','2')
p3 = gf.GlobalFunction('parser','3')
p00 = gf.GlobalFunction('parser','x','[1;0]')
p11 = gf.GlobalFunction('parser','y','[0;1]')
p22 = gf.GlobalFunction('parser','r','[x/r;y/r]','[y*y/(r*r*r);-x*y/(r*r*r);-x*y/(r*r*r);x*x/(r*r*r)]')
p33 = gf.GlobalFunction('parser','theta','[-y/(r*r);x/(r*r)]','[2*x*y/(pow(r,4));(y*y-x*x)/(pow(r,4));(y*y-x*x)/(pow(r,4));-2*x*y/(pow(r,4))]')

mf.export_to_pos('check_global_functions16.pos',p0(PTs),'0')
mf.export_to_pos('check_global_functions17.pos',p1(PTs),'1')
mf.export_to_pos('check_global_functions18.pos',p2(PTs),'2')
mf.export_to_pos('check_global_functions19.pos',p3(PTs),'3')
mf.export_to_pos('check_global_functions20.pos',p00(PTs),'x')
mf.export_to_pos('check_global_functions21.pos',p11(PTs),'y')
mf.export_to_pos('check_global_functions22.pos',p22(PTs),'r')
mf.export_to_pos('check_global_functions23.pos',p33(PTs),'theta')
mf.export_to_pos('check_global_functions24.pos',p00.grad(PTs),'grad(x)')
mf.export_to_pos('check_global_functions25.pos',p11.grad(PTs),'grad(y)')
mf.export_to_pos('check_global_functions26.pos',p22.grad(PTs),'grad(r)')
mf.export_to_pos('check_global_functions27.pos',p33.grad(PTs),'grad(theta)')

# product:
p0 = ck0*ck1
p1 = ck1*ck2
p2 = ck2*ck3
mf.export_to_pos('check_global_functions28.pos',p0(PTs),'ck0*ck1')
mf.export_to_pos('check_global_functions29.pos',p1(PTs),'ck1*ck2')
mf.export_to_pos('check_global_functions30.pos',p2(PTs),'ck2*ck3')

# add:
ad0 = ck0+ck1
ad1 = ck1+ck2
ad2 = ck2+ck3
mf.export_to_pos('check_global_functions31.pos',ad0(PTs),'ck0+ck1')
mf.export_to_pos('check_global_functions32.pos',ad1(PTs),'ck1+ck2')
mf.export_to_pos('check_global_functions33.pos',ad2(PTs),'ck2+ck3')

for i in range(34):
  os.remove("check_global_functions%i.pos" % i);
