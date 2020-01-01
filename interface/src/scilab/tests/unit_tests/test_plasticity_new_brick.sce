// Copyright (C) 2009-2020 Yann Colette
// 
//  This file is a part of GetFEM++
// 
//  GetFEM++  is  free software;  you  can  redistribute  it  and/or modify it
//  under  the  terms  of the  GNU  Lesser General Public License as published
//  by  the  Free Software Foundation;  either version 3 of the License,  or
//  (at your option) any later version along with the GCC Runtime Library
//  Exception either version 3.1 or (at your option) any later version.
//  This program  is  distributed  in  the  hope  that it will be useful,  but
//  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
//  or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
//  License and GCC Runtime Library Exception for more details.
//  You  should  have received a copy of the GNU Lesser General Public License
//  along  with  this program;  if not, write to the Free Software Foundation,
//  Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.

gf_workspace('clear all');

//

// We try to compute a plasticity problem with a Von Mises crierion
// For convenience we consider an homogenous Dirichlet condition on the left
// of the domain and an easy computed Neumann Condition on the right

//

// Initialize used data
L      = 100;
H      = 20;
lambda = 121150;
mu     = 80769;
von_mises_threshold = 8000;
f = [0 -330]';
t = [0 0.9032 1 1.1 1.3 1.5 1.7 1.74 1.7 1.5 1.3 1.1 1 0.9032 0.7 0.5 0.3 0.1 0];

// Create the mesh
m = gf_mesh('triangles grid', [0:4:L], [0:2:H]);

// Plotting
h_graph = scf();
h_graph.color_map = jetcolormap(256);
drawlater;
gf_plot_mesh(m, 'vertices', 'on', 'convexes', 'on');
drawnow;

// Define used MeshIm
mim = gf_mesh_im(m);  gf_mesh_im_set(mim, 'integ', gf_integ('IM_TRIANGLE(6)')); // Gauss methods on triangles

// Define used MeshFem
mf_u     = gf_mesh_fem(m,2); gf_mesh_fem_set(mf_u,     'fem', gf_fem('FEM_PK(2,2)'));
mf_data  = gf_mesh_fem(m);   gf_mesh_fem_set(mf_data,  'fem', gf_fem('FEM_PK_DISCONTINUOUS(2,0)'));
mf_sigma = gf_mesh_fem(m,4); gf_mesh_fem_set(mf_sigma, 'fem', gf_fem('FEM_PK_DISCONTINUOUS(2,1)'));
mf_err   = gf_mesh_fem(m);   gf_mesh_fem_set(mf_err,   'fem', gf_fem('FEM_PK(2,0)'));
mf_vm    = gf_mesh_fem(m);   gf_mesh_fem_set(mf_vm,    'fem', gf_fem('FEM_PK_DISCONTINUOUS(2,1)'));
mf_pl    = gf_mesh_fem(m);   gf_mesh_fem_set(mf_pl,    'fem', gf_fem('FEM_PK_DISCONTINUOUS(2,1)'));

// Find the border of the domain
P = gf_mesh_get(m, 'pts');
pidleft  = find(abs(P(1,:))<1e-6); // Retrieve index of points which x near to 0
pidright = find(abs(P(1,:) - L)<1e-6); // Retrieve index of points which x near to L
fleft    = gf_mesh_get(m, 'faces from pid', pidleft); 
fright   = gf_mesh_get(m, 'faces from pid', pidright);

// Decomposed the mesh into 2 regions with different values of Lamé coeff
CV       = gf_mesh_get(m, 'cvid');
CVbottom = find(CV <= 250); // Retrieve index of convex located at the bottom
CVtop    = find(CV > 250);  // Retrieve index of convex located at the top

// Definition of Lame coeff
lambda(CVbottom) = 121150; // Stell
lambda(CVtop)    = 84605;  // Iron
//lambda(CV) = 84605;
mu(CVbottom) = 80769; // Stell
mu(CVtop)    = 77839; // Iron
//mu(CV) = 77839;
von_mises_threshold(CVbottom) = 7000;
von_mises_threshold(CVtop)    = 8000;

// Assign boundary numbers
gf_mesh_set(m,'boundary',1,fleft);  // for Dirichlet condition
gf_mesh_set(m,'boundary',2,fright); // for Neumann condition

// Create the model
md = gf_model('real');

// Declare that u is the unknown of the system on mf_u
// 2 is the number of version of the data stored, for the time integration scheme 
gf_model_set(md, 'add fem variable', 'u', mf_u, 2); 

// Declare that lambda is a data of the system on mf_data
gf_model_set(md, 'add initialized fem data', 'lambda', mf_data, lambda);

// Declare that mu is a data of the system on mf_data
gf_model_set(md, 'add initialized fem data', 'mu', mf_data, mu);

// Declare that von_mises_threshold is a data of the system on mf_data
gf_model_set(md, 'add initialized fem data', 'von_mises_threshold', mf_data, von_mises_threshold);

// Declare that sigma is a data of the system on mf_sigma
// 2 is the number of version of the data stored, for the time integration scheme
gf_model_set(md, 'add fem data', 'sigma', mf_sigma);

// Add plasticity brick on u
gf_model_set(md, 'add elastoplasticity brick', mim, 'VM', 'u', 'lambda', 'mu', 'von_mises_threshold', 'sigma');

// Add homogeneous Dirichlet condition to u on the left hand side of the domain
gf_model_set(md, 'add Dirichlet condition with multipliers', mim, 'u', mf_u, 1);

// Add a source term to the system
gf_model_set(md,'add initialized fem data', 'VolumicData', mf_data, gf_mesh_fem_get_eval(mf_data, list(list(f(1,1)),list(f(2,1)*t(1)))));
gf_model_set(md, 'add source term brick', mim, 'u', 'VolumicData', 2);

VM = zeros(1,gf_mesh_fem_get(mf_vm, 'nbdof'));

nbstep = size(t,2);

dd = gf_mesh_fem_get(mf_err, 'basic dof from cvid');

h_graph_2 = scf();
h_graph_2.color_map = jetcolormap(256);

h_graph_3 = scf();
h_graph_3.color_map = jetcolormap(256);

for step=1:nbstep,
  if step > 1 then
    gf_model_set(md, 'variable', 'VolumicData', gf_mesh_fem_get_eval(mf_data, list(list(f(1,1)),list(f(2,1)*t(step)))));
  end

  // Solve the system
  gf_model_get(md, 'solve','lsolver', 'superlu', 'lsearch', 'simplest',  'alpha min', 0.8, 'very noisy', 'max_iter', 100, 'max_res', 1e-6);

  // Retrieve the solution U
  U = gf_model_get(md, 'variable', 'u', 0);
    
  // Compute new plasticity constraints used to compute 
  // the Von Mises or Tresca stress
  gf_model_get(md, 'elastoplasticity next iter', mim, 'u', 'VM', 'lambda', 'mu', 'von_mises_threshold', 'sigma');
  plast = gf_model_get(md, 'compute plastic part', mim, mf_pl, 'u', 'VM', 'lambda', 'mu', 'von_mises_threshold', 'sigma');
      
  // Compute Von Mises or Tresca stress
  VM = gf_model_get(md, 'compute elastoplasticity Von Mises or Tresca', 'sigma', mf_vm, 'Von Mises');
  
  scf(h_graph_2);
  drawlater;
  clf(h_graph_2);
  subplot(2,1,1);
  gf_plot(mf_vm,VM,'deformed_mesh', 'on', 'deformation', U, 'deformation_mf', mf_u, 'refine', 4, 'deformation_scale',1); 
  colorbar(min(VM),max(VM));

  n = t(step);
  title(['Von Mises criterion for t = ', string(n)]);
  
  ERR = gf_compute(mf_u, U, 'error estimate', mim);
  E = ERR; E(dd) = ERR;

  subplot(2,1,2);
  gf_plot(mf_err, E, 'mesh','on', 'refine', 1); 
  colorbar(min(E),max(E));
  title('Error estimate');
  drawnow;

  scf(h_graph_3);
  drawlater;
  clf(h_graph_3);
  gf_plot(mf_pl,plast,'deformed_mesh','on', 'deformation',U,'deformation_mf',mf_u,'refine', 4, 'deformation_scale',1); 
  colorbar(min(plast),max(plast));
  n = t(step);
  title(['Plastification for t = ', string(n)]);
  drawnow;

  sleep(1000);
end
