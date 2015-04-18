#!/usr/bin/env python
# -*- coding: UTF8 -*-
# Python GetFEM++ interface
#
# Copyright (C)  2015-2015 Yves Renard.
#
# This file is a part of GetFEM++
#
# Getfem++  is  free software;  you  can  redistribute  it  and/or modify it
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

import getfem as gf
import numpy as np

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
# Electric problem: The potential is prescribed to be 0V at the right
#   boundary and 0.1V at the left boundary.
# Thermal problem: A thermal insulation condition is prescribed at the
#   left and hole boudnaries. The remaining boundary and the plate itself
#   is supposed to be submitted to an heat transfert with respect to the
#   air at 20oC.
# Coupling terms:
#   - Joule heating: source term  sigma|Grad_V|^2
#   - Dependance of the thermal conductivity in temperature :
#     sigma = 1/(rho_0(1+alpha(theta-T0)))
#     with T0 = 20oC, rho_0 the resistance temperature coefficient at T0
#     and alpha the second resistance temperature coefficient.
#   - Thermal expansion:
#     stress_tensor = clambdastar div(u) I + 2 cmu epsilon(u) - beta theta I
#     with beta = alpha_th E/(1-2nu), alpha_th being the thermal
#     expansion coefficient.
# The first two coupling terms are nonlinear ones.


#
# Physical parameters
#
epsilon = 1.       # Thickness of the plate (cm)
E = 21E6           # Young Modulus (N/cm^2)
nu = 0.3           # Poisson ratio
clambda = E*nu/((1+nu)*(1-2*nu)) # First Lame coefficient (N/cm^2)
cmu = E/(2*(1+nu))               # Second Lame coefficient (N/cm^2)
clambdastar = 2*clambda*cmu/(clambda+2*cmu) # Lame coefficient for Plane stress (N/cm^2)
F = 100E2          # Force density at the right boundary (N/cm^2)
kappa = 4.         # Thermal conductivity (W/(cm K))
D = 10.            # Heat transfert coefficient (W/(K cm^2))
air_temp = 20.     # Temperature of the air in oC.
alpha_th = 16.6E-6 # Thermal expansion coefficient (/K).
T0 = 20.           # Reference temperature in oC.
rho_0 = 1.754E-8   # Resistance temperature coefficient at T0 = 20oC
alpha = 0.0039     # Second resistance temperature coefficient.

#
# Numerical parameters
#
h = 2.                     # Approximate mesh size
elements_degree = 2        # Degree of the finite element methods
export_mesh = True         # Draw the mesh after mesh generation or not
solve_in_two_steps = True  # Solve the elasticity problem separately or not

#
# Mesh generation. Meshes can also been imported from several formats.
#
mo1 = gf.MesherObject('rectangle', [0., 0.], [100., 25.])
mo2 = gf.MesherObject('ball', [25., 12.5], 8.)
mo3 = gf.MesherObject('ball', [50., 12.5], 8.)
mo4 = gf.MesherObject('ball', [75., 12.5], 8.)
mo5 = gf.MesherObject('union', mo2, mo3, mo4)
mo  = gf.MesherObject('set minus', mo1, mo5)

print 'Mesh generation';
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

RIGHT_BOUND=1; LEFT_BOUND=2; TOP_BOUND=3; BOTTOM_BOUND=4; HOLE_BOUND=5;

mesh.set_region( RIGHT_BOUND, fb2)
mesh.set_region(  LEFT_BOUND, fb3)
mesh.set_region(   TOP_BOUND, fb4)
mesh.set_region(BOTTOM_BOUND, fb5)
mesh.set_region(  HOLE_BOUND, fb1)
mesh.region_subtract( RIGHT_BOUND, HOLE_BOUND)
mesh.region_subtract(  LEFT_BOUND, HOLE_BOUND)
mesh.region_subtract(   TOP_BOUND, HOLE_BOUND)
mesh.region_subtract(BOTTOM_BOUND, HOLE_BOUND)

if (export_mesh):
    mesh.export_to_vtk('mesh.vtk');
    print ('\nYou can view the mesh for instance with');
    print ('mayavi2 -d mesh.vtk -f ExtractEdges -m Surface \n');

#
# Definition of finite elements methods and integration method
#

mfu = gf.MeshFem(mesh, 2)  # Finite element for the elastic displacement
mfu.set_classical_fem(elements_degree)
mft = gf.MeshFem(mesh, 1)  # Finite element for temperature and electrical field
mft.set_classical_fem(elements_degree)
mfvm = gf.MeshFem(mesh, 1) # Finite element for Von Mises stress interpolation
mfvm.set_classical_discontinuous_fem(elements_degree)
mim = gf.MeshIm(mesh, pow(elements_degree,2))   # Integration method


#
# Model definition
#

md=gf.Model('real');
md.add_fem_variable('u', mfu)       # Displacement of the structure
md.add_fem_variable('theta', mft)   # Temperature
md.add_fem_variable('V', mft)       # Electric potential

# Membrane elastic deformation
md.add_initialized_data('cmu', [cmu])
md.add_initialized_data('clambdastar', [clambdastar])
md.add_isotropic_linearized_elasticity_brick(mim, 'u', 'clambdastar', 'cmu')

md.add_Dirichlet_condition_with_multipliers(mim, 'u', elements_degree-1, LEFT_BOUND)
md.add_initialized_data('Fdata', [F*epsilon, 0])
md.add_source_term_brick(mim, 'u', 'Fdata', RIGHT_BOUND)

# Electrical field
sigmaeps = '(eps/(rho_0*(1+alpha*(theta-T0))))'
md.add_initialized_data('eps', [epsilon])
md.add_initialized_data('rho_0', [rho_0])
md.add_initialized_data('alpha', [alpha])
md.add_initialized_data('T0', [T0])
md.add_nonlinear_generic_assembly_brick(mim, sigmaeps+'*(Grad_V.Grad_Test_V)')
md.add_Dirichlet_condition_with_multipliers(mim, 'V', elements_degree-1, RIGHT_BOUND)
md.add_initialized_data('DdataV', [0.1])
md.add_Dirichlet_condition_with_multipliers(mim, 'V', elements_degree-1, LEFT_BOUND, 'DdataV')

# Thermal problem
md.add_initialized_data('kaeps', [kappa*epsilon])
md.add_generic_elliptic_brick(mim, 'theta', 'kaeps')
md.add_initialized_data('D2', [D*2])
md.add_initialized_data('D2airt', [air_temp*D*2])
md.add_mass_brick(mim, 'theta', 'D2')
md.add_source_term_brick(mim, 'theta', 'D2airt')
md.add_initialized_data('Deps', [D/epsilon])
md.add_initialized_data('Depsairt', [air_temp*D/epsilon])
md.add_Fourier_Robin_brick(mim, 'theta', 'Deps', TOP_BOUND)
md.add_source_term_brick(mim, 'theta', 'Depsairt', TOP_BOUND)
md.add_Fourier_Robin_brick(mim, 'theta', 'Deps', BOTTOM_BOUND)
md.add_source_term_brick(mim, 'theta', 'Depsairt', BOTTOM_BOUND)

# Joule heating term
md.add_nonlinear_generic_assembly_brick(mim, '-'+sigmaeps+'*Norm_sqr(Grad_V)*Test_theta')

# Thermal expansion term
md.add_initialized_data('beta', [alpha_th*E/(1-2*nu)])
md.add_linear_generic_assembly_brick(mim, 'beta*(T0-theta)*Trace(Grad_Test_u)')


#
# Model solve
#
if (solve_in_two_steps):
  md.disable_variable('u')
  print 'First problem with', md.nbdof(), ' dofs'
  md.solve('max_res', 1E-9, 'max_iter', 100, 'noisy')
  md.enable_variable('u')
  md.disable_variable('theta')
  md.disable_variable('V')
  print 'Second problem with ', md.nbdof(), ' dofs'
  md.solve('max_res', 1E-9, 'max_iter', 100, 'noisy')
else:
  print 'Global problem with ', md.nbdof(), ' dofs'
  md.solve('max_res', 1E-9, 'max_iter', 100, 'noisy')


#
# Solution export
#  
U = md.variable('u')
V = md.variable('V')
THETA = md.variable('theta')
VM = md.compute_isotropic_linearized_Von_Mises_or_Tresca('u', 'clambdastar', 'cmu', mfvm)
CO = np.reshape(md.interpolation('-'+sigmaeps+'*Grad_V', mfvm), (2, mfvm.nbdof()), 'F')

mfvm.export_to_vtk('displacement_with_von_mises.vtk', mfvm,  VM, 'Von Mises Stresses', mfu, U, 'Displacements')
print ('You can view solutions with for instance:\nmayavi2 -d displacement_with_von_mises.vtk -f WarpVector -m Surface')
mft.export_to_vtk('temperature.vtk', mft, THETA, 'Temperature')
print ('mayavi2 -d temperature.vtk -f WarpScalar -m Surface')
mft.export_to_vtk('electric_potential.vtk', mft, V, 'Electric potential')
print ('mayavi2 -d electric_potential.vtk -f WarpScalar -m Surface')










