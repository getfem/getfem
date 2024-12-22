// Copyright (C) 2010-2020 Yann COLLETTE.
//
// This file is a part of GetFEM
//
// GetFEM  is  free software;  you  can  redistribute  it  and/or modify it
// under  the  terms  of the  GNU  Lesser General Public License as published
// by  the  Free Software Foundation;  either version 3 of the License,  or
// (at your option) any later version along with the GCC Runtime Library
// Exception either version 3.1 or (at your option) any later version.
// This program  is  distributed  in  the  hope  that it will be useful,  but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
// License and GCC Runtime Library Exception for more details.
// You  should  have received a copy of the GNU Lesser General Public License
// along  with  this program;  if not, write to the Free Software Foundation,
// Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.

lines(0);

path = get_absolute_file_path('demo_plate.sce');

printf('demo plate started\n');

if getos()=='Windows' then
  // Under Windows, all the trace messages are available in the dos console
  // Under Linuxs, all the trace messages are redirected to the Scilab console
  consolebox('on');
end
gf_util('trace level',3);
gf_util('warning level',3);


// Simple supported Mindlin-Reissner plate


Emodulus = 1;          // Young Modulus
nu       = 0.5;        // Poisson Coefficient
epsilon  = 0.001;      // Plate thickness
kappa     = 5/6;       // Shear correction factor
f = -5*epsilon^3;      // Prescribed force on the top of the plate

variant = 2;           // 0 : not reduced, 1 : with reduced integration, 2 : MITC reduction
quadrangles = %t;      // Locking free only on quadrangle for the moment
K = 1;                 // Degree of the finite element method
with_Mindlin_brick = %t; // Uses the Reissner-Mindlin predefined brick or not
dirichlet_version = 1; // 0 = simplification, 1 = with multipliers, 2 = penalization

plot_mesh = %f;
draw_solution = %t;

// trace on;
gf_workspace('clear all');
NX = 80;
if (quadrangles)
  m = gf_mesh('cartesian',[0:1/NX:1],[0:1/NX:1]);
else
  m=gf_mesh('import','structured',sprintf('GT=""GT_PK(2,1)"";SIZES=[1,1];NOISED=0;NSUBDIV=[%d,%d];', NX, NX));
end

// Create a mesh_fem of for a 2 dimension vector field
mftheta = gf_mesh_fem(m,2);
mfu = gf_mesh_fem(m,1);
// Assign the QK or PK fem to all convexes of the mesh_fem, and define an
// integration method
if (quadrangles)
  gf_mesh_fem_set(mftheta,'fem',gf_fem(sprintf('FEM_QK(2,%d)', K)));
  gf_mesh_fem_set(mfu,'fem',gf_fem(sprintf('FEM_QK(2,%d)', K)));
  mim = gf_mesh_im(m, gf_integ('IM_GAUSS_PARALLELEPIPED(2,6)'));
  mim_reduced = gf_mesh_im(m, gf_integ('IM_GAUSS_PARALLELEPIPED(2,1)'));
else
  gf_mesh_fem_set(mftheta,'fem',gf_fem(sprintf('FEM_PK(2,%d)', K)));
  gf_mesh_fem_set(mfu,'fem',gf_fem(sprintf('FEM_PK(2,%d)', K)));
  mim = gf_mesh_im(m, gf_integ('IM_TRIANGLE(6)'));
  mim_reduced = gf_mesh_im(m, gf_integ('IM_TRIANGLE(1)'));
end

// detect the border of the mesh
border = gf_mesh_get(m,'outer faces');
// mark it as boundary #1
gf_mesh_set(m, 'boundary', 1, border);
if (plot_mesh)
  gf_plot_mesh(m, 'regions', [1]); // the boundary edges appears in red
  sleep(1000);
end

md=gf_model('real');
gf_model_set(md, 'add fem variable', 'u', mfu);
gf_model_set(md, 'add fem variable', 'theta', mftheta);
gf_model_set(md, 'add initialized data', 'E', Emodulus);
gf_model_set(md, 'add initialized data', 'nu', nu);
gf_model_set(md, 'add initialized data', 'epsilon', epsilon);
gf_model_set(md, 'add initialized data', 'kappa', kappa);


if (with_Mindlin_brick)
  gf_model_set(md, 'add Mindlin Reissner plate brick', mim, mim_reduced, 'u', 'theta', 'E', 'nu', 'epsilon', 'kappa', variant);
else
  gf_model_set(md, 'add elementary rotated RT0 projection', 'RT0_projection');
  gf_model_set(md, 'add linear generic assembly brick', mim, '(E*epsilon*epsilon*epsilon*(1-nu)/(48 * (1 - nu*nu))) * ((Grad_theta+Grad_theta''):(Grad_Test_theta+Grad_Test_theta''))');
  gf_model_set(md, 'add linear generic assembly brick', mim, '(E*epsilon*epsilon*epsilon*nu/(12 * (1 - nu*nu))) * (Trace(Grad_theta)*Trace(Grad_Test_theta))');
  if (variant == 0)
    gf_model_set(md, 'add linear generic assembly brick', mim, '(E*kappa*epsilon/(1 + nu)) * ((Grad_u + theta).Grad_Test_u) + (E*kappa*epsilon/(1 + nu)) * ((Grad_u + theta).Test_theta)');
  elseif (variant == 1)
    gf_model_set(md, 'add linear generic assembly brick', mim_reduced, '(E*kappa*epsilon/(1 + nu)) * ((Grad_u + theta).Grad_Test_u) + (E*kappa*epsilon/(1 + nu)) * ((Grad_u + theta).Test_theta)');
  else
    gf_model_set(md, 'add linear generic assembly brick', mim, '(E*kappa*epsilon/(1 + nu)) * ((Grad_u + Elementary_transformation(theta,RT0_projection)).Grad_Test_u) + (E*kappa*epsilon/(1 + nu)) * ((Grad_u + Elementary_transformation(theta, RT0_projection)).(Elementary_transformation(Test_theta, RT0_projection)))');  
  end
end

gf_model_set(md, 'add initialized data', 'VolumicData', f);

gf_model_set(md, 'add source term brick', mim, 'u', 'VolumicData');
gf_model_set(md, 'add initialized data', 'DirichletData', 0);
select (dirichlet_version)
  case 0,
    gf_model_set(md, 'add Dirichlet condition with simplification', 'u', 1, 'DirichletData');   
  case 1, 
    gf_model_set(md, 'add Dirichlet condition with multipliers', mim, 'u', mfu, 1, 'DirichletData');
  case 2,
    gf_model_set(md, 'add Dirichlet condition with penalization', mim, 'u', r, 1, 'DirichletData');
end
gf_model_get(md, 'solve');
U = gf_model_get(md, 'variable', 'u');

if (draw_solution)
  hh = scf();
  hh.color_map = jetcolormap(255);
  gf_plot(mfu,U,'mesh','off', 'zplot', 'on'); 
  colorbar(min(U),max(U)); title('computed solution');
end























