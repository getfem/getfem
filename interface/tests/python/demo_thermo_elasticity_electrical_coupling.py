#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Python GetFEM interface
#
# Copyright (C) 2015-2026 Yves Renard.
#
# This file is a part of GetFEM
#
# GetFEM  is  free software;  you  can  redistribute  it  and/or modify it
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
#
############################################################################

import numpy as np
import getfem as gf

# Deformation of a plate under the coupling of thermal, elasticity, and
# electric effects.
#
#
#     ______________________________________
#   /|         __       __       __         |->
#   /|        /  \     /  \     /  \        |->
#   /|       |    |   |    |   |    |       |-> F
#   /|        \__/     \__/     \__/        |->
#   /|______________________________________|->
#     
#
# Elastic problem: The plate is clamped at rhe left boundary and a
#   traction density of force F is prescribed at the right boundary.
#   Plane stress conditions are assumed.
# Electric problem: The potential is prescribed to be 0V at the right
#   boundary and 0.1V at the left boundary.
# Thermal problem: A thermal insulation condition is prescribed at the
#   left, right, and hole boudnaries. The remaining boundary and the
#   plate front and back surfaces are supposed to transfer heat by
#   convection with respect to the surrounding air at 20 deg C.
# Coupling terms:
#   - Joule heating: source term  1/rho ||Grad_V||^2
#   - Dependance of the thermal resistivity on temperature :
#     rho = rho_0(1+alpha(T-T0))
#     with T0 = 20 deg C, rho_0 the resistivity at T0
#     and alpha the resistivity-temperature coefficient.
#   - Thermal expansion:
#     stress_tensor = E/(1+nu) ( nu/(1-nu) (div(u) - 2 alpha_th DT) I
#                               + (epsilon(u) - alpha_th DT I) )
#     with alpha_th being the thermal expansion coefficient.
# The first two coupling terms are nonlinear ones.


#
# Physical parameters
#

t = 1.             # Thickness of the plate (cm)
E = 21E6           # Young Modulus (N/cm^2)
nu = 0.3           # Poisson ratio
F = 100E2          # Force density at the right boundary (N/cm^2)
kappa = 4.         # Thermal conductivity (W/(cm K))
D = 10.            # Heat transfer coefficient (W/(K cm^2))
air_temp = 20.     # Temperature of the air in deg C
alpha_th = 16.6E-6 # Thermal expansion coefficient (1/K)
T0 = 20.           # Reference temperature in deg C
rho_0 = 1.754E-8   # Resistivity at T0
alpha = 0.0039     # Resistivity-temperature coefficient


#
# Numerical parameters
#

h = 2.                    # Approximate mesh size
elements_degree = 2       # Degree of the finite element methods
export_mesh = False       # Export the mesh after mesh generation or not
solve_in_two_steps = 2    # Solve the elasticity problem separately (1)
                          # or in a coupled way (0) or both and compare (2)


#
# Mesh generation. Meshes can also be imported in various formats.
#

mo1 = gf.MesherObject('rectangle', [0., 0.], [100., 25.])
mo2 = gf.MesherObject('ball', [25., 12.5], 8.)
mo3 = gf.MesherObject('ball', [50., 12.5], 8.)
mo4 = gf.MesherObject('ball', [75., 12.5], 8.)
mo5 = gf.MesherObject('union', mo2, mo3, mo4)
mo  = gf.MesherObject('set minus', mo1, mo5)

print('Mesh generation')
gf.util('trace level', 2)   # No trace for mesh generation
mesh = gf.Mesh('generate', mo, h, 2)


#
# Boundary selection
#

fb1 = mesh.outer_faces_in_box([1., 1.], [99., 24.])    # Boundary of the holes
fb2 = mesh.outer_faces_with_direction([ 1., 0.], 0.01) # Right boundary
fb3 = mesh.outer_faces_with_direction([-1., 0.], 0.01) # Left boundary
fb4 = mesh.outer_faces_with_direction([0.,  1.], 0.01) # Top boundary
fb5 = mesh.outer_faces_with_direction([0., -1.], 0.01) # Bottom boundary
fb6 = mesh.outer_faces_in_ball([25., 12.5], 8.+0.01*h) # Left hole boundary
fb7 = mesh.outer_faces_in_ball([50., 12.5], 8.+0.01*h) # Center hole boundary
fb8 = mesh.outer_faces_in_ball([75., 12.5], 8.+0.01*h) # Right hole boundary

RIGHT_BOUND=1; LEFT_BOUND=2; TOP_BOUND=3; BOTTOM_BOUND=4; HOLE_BOUND=5;
HOLE1_BOUND = 6; HOLE2_BOUND = 7; HOLE3_BOUND = 8;
mesh.set_region( RIGHT_BOUND, fb2)
mesh.set_region(  LEFT_BOUND, fb3)
mesh.set_region(   TOP_BOUND, fb4)
mesh.set_region(BOTTOM_BOUND, fb5)
mesh.set_region(  HOLE_BOUND, fb1)
mesh.set_region( HOLE1_BOUND, fb6)
mesh.set_region( HOLE2_BOUND, fb7)
mesh.set_region( HOLE3_BOUND, fb8)
mesh.region_subtract( RIGHT_BOUND, HOLE_BOUND)
mesh.region_subtract(  LEFT_BOUND, HOLE_BOUND)
mesh.region_subtract(   TOP_BOUND, HOLE_BOUND)
mesh.region_subtract(BOTTOM_BOUND, HOLE_BOUND)
mesh.region_merge(HOLE1_BOUND, HOLE2_BOUND)
mesh.region_merge(HOLE1_BOUND, HOLE3_BOUND)
np.testing.assert_array_equal(mesh.region(HOLE_BOUND), mesh.region(HOLE1_BOUND))

if (export_mesh):
    mesh.export_to_vtu('mesh.vtu');
    print('\nYou can view the mesh for instance with');
    print('mayavi2 -d mesh.vtu -f ExtractEdges -m Surface \n');


#
# Definition of finite element methods and integration method
#

mfu = gf.MeshFem(mesh, 2)  # Finite element for the elastic displacement
mfu.set_classical_fem(elements_degree)
mft = gf.MeshFem(mesh, 1)  # Finite element for temperature and electrical field
mft.set_classical_fem(elements_degree)
mfvm = gf.MeshFem(mesh, 1) # Finite element for Von Mises stress interpolation
mfvm.set_classical_discontinuous_fem(elements_degree)
mim = gf.MeshIm(mesh, elements_degree*2)   # Integration method


#
# Model definition
#

md=gf.Model('real');
md.add_fem_variable('u', mfu)   # Displacement of the structure
md.add_fem_variable('T', mft)   # Temperature
md.add_fem_variable('V', mft)   # Electric potential
md.add_initialized_data('t', t)

# Membrane elastic deformation and thermal expansion
md.add_initialized_data('E', E)
md.add_initialized_data('nu', nu)
md.add_initialized_data('alpha_th', alpha_th)
md.add_initialized_data('T0', T0)
md.add_macro('sigma',
             'E/(1+nu)*( nu/(1-nu)*(Div(u)-2*alpha_th*T)*Id(2)'
                       '+(Sym(Grad(u))-alpha_th*T*Id(2)) )')
md.add_linear_term(mim, 't*sigma:Grad(Test_u)')
md.add_Dirichlet_condition_with_multipliers(mim, 'u', elements_degree-1, LEFT_BOUND)
md.add_initialized_data('F', F)
md.add_linear_term(mim, '-t*F*Test_u(1)', RIGHT_BOUND)

# Electric field
md.add_initialized_data('rho_0', rho_0)
md.add_initialized_data('alpha', alpha)
md.add_macro('rho', 'rho_0*(1+alpha*(T-T0))')
md.add_nonlinear_term(mim, 't/rho * Grad(V).Grad(Test_V)')
md.add_Dirichlet_condition_with_multipliers(mim, 'V', elements_degree-1, RIGHT_BOUND)
md.add_initialized_data('DdataV', 0.1)
md.add_Dirichlet_condition_with_multipliers(mim, 'V', elements_degree-1, LEFT_BOUND, 'DdataV')

# Thermal problem
md.add_initialized_data('kappaT', kappa)
md.add_initialized_data('D', D)
md.add_initialized_data('T_air', air_temp)
md.add_linear_term(mim, 't*kappaT*Grad(T).Grad(Test_T) + 2*D*(T-T_air)*Test_T')
md.add_linear_term(mim, 't*D*(T-T_air).Test_T', TOP_BOUND)
md.add_linear_term(mim, 't*D*(T-T_air).Test_T', BOTTOM_BOUND)
# Joule heating term
md.add_nonlinear_term(mim, '-t/rho * Norm_sqr(Grad(V))*Test_T')


#
# Model solve
#

if (solve_in_two_steps >= 1):
  md.disable_variable('u')
  print('First problem with', md.nbdof(), ' dofs')
  md.solve('max_res', 1E-9, 'max_iter', 100, 'noisy')
  md.enable_variable('u')
  md.disable_variable('T')
  md.disable_variable('V')
  print('Second problem with ', md.nbdof(), ' dofs')
  md.solve('max_res', 1E-9, 'max_iter', 100, 'noisy')
  md.enable_variable('T')
  md.enable_variable('V')
  U1 = md.variable('u')
  
if (solve_in_two_steps == 0):
  print('Global problem with ', md.nbdof(), ' dofs')
  md.solve('max_res', 1E-9, 'max_iter', 100, 'noisy')

if (solve_in_two_steps == 2):
  print('Global problem with ', md.nbdof(), ' dofs')
  md.set_variable('u', md.variable('u')*0.);
  md.solve('max_res', 1E-9, 'max_iter', 100, 'noisy')
  U2 = md.variable('u')
  print (np.linalg.norm(U2-U1));
  if (np.linalg.norm(U2-U1) > 1E-10):
      print("Too big difference between solve in one and two steps")
      exit(1)

#
# Solution export
#

print("L2 norm of temperature =", np.sqrt(gf.asm_generic(mim, 0, "sqr(T)", -1, md)))
VM = md.local_projection(mim, "sqrt(Norm_sqr(sigma)+sqr(sigma(1,2))-sigma(1,1)*sigma(2,2))", mfvm)
#VM = md.local_projection(mim, "sqrt(1.5)*Norm_sqr([1,0;0,1;0,0]*sigma*[1,0,0;0,1,0]-1/3*Trace(sigma)*Id(3))", mfvm)
CO = md.interpolation('-t/rho * Grad(V)', mfvm).reshape(2, mfvm.nbdof(), order='F')
mfvm.export_to_vtu('displacement_with_von_mises.vtu', mfvm,  VM, 'Von Mises Stresses',
                                                      mfu, md.variable('u'), 'Displacements')
mft.export_to_vtu('temperature.vtu', mft, md.variable('T'), 'Temperature')
mft.export_to_vtu('electric_potential.vtu', mft, md.variable('V'), 'Electric potential')
print('You can view solutions with for instance:\nmayavi2 -d displacement_with_von_mises.vtu -f WarpVector -m Surface')
print('mayavi2 -d temperature.vtu -f WarpScalar -m Surface')
print('mayavi2 -d electric_potential.vtu -f WarpScalar -m Surface')
