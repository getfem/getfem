// Matlab GetFEM interface
//
// Copyright (C) 2009-2020 Yves Renard.
//
// This file is a part of GetFEM
//
// GetFEM  is  free software;  you  can  redistribute  it  and/or modify it
// under  the  terms  of the  GNU  Lesser General Public License as published
// by  the  Free Software Foundation;  either version 2.1 of the License,  or
// (at your option) any later version.
// This program  is  distributed  in  the  hope  that it will be useful,  but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
// License for more details.
// You  should  have received a copy of the GNU Lesser General Public License
// along  with  this program;  if not, write to the Free Software Foundation,
// Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.
//
// Static equilibrium of an elastic solid in contact with a rigid foundation.
// Tests the different contact/friction formulations of Getfem.
//
// This program is used to check that matlab-getfem is working. This is also
// a good example of use of GetFEM.
//

gf_workspace('clear all');
path = get_absolute_file_path('demo_static_contact.sce');

if getos()=='Windows' then
  // Under Windows, all the trace messages are available in the dos console
  // Under Linuxs, all the trace messages are redirected to the Scilab console
  consolebox('on');
end
gf_util('trace level',3);
gf_util('warning level',3);

// Import the mesh : disc
// m = gf_mesh('load', path + 'data/disc_P2_h2.mesh');
// m = gf_mesh('load', path + 'data/disc_P2_h1.mesh');
// m = gf_mesh('load', path + 'data/disc_P2_h0.5.mesh');
// m = gf_mesh('load', path + 'data/disc_P2_h0.25.mesh');
// m = gf_mesh('load', path + 'data/disc_P2_h0.15.mesh');

// Import the mesh : sphere
// m = gf_mesh('load', path + 'data/sphere_with_quadratic_tetra_8_elts.mesh');
m = gf_mesh('load', path + 'data/sphere_with_quadratic_tetra_80_elts.mesh');
// m = gf_mesh('load', path + 'data/sphere_with_quadratic_tetra_400_elts.mesh');
// m = gf_mesh('load', path + 'data/sphere_with_quadratic_tetra_2000_elts.mesh');
// m = gf_mesh('load', path + 'data/sphere_with_quadratic_tetra_16000_elts.mesh');

d = gf_mesh_get(m, 'dim'); // Mesh dimension

// Parameters of the model
clambda = 1;           // Lame coefficient
cmu     = 1;           // Lame coefficient
friction_coeff = 0.4;  // coefficient of friction
vertical_force = 0.05; // Volumic load in the vertical direction
r = 10;                // Augmentation parameter
condition_type = 0;    // 0 = Explicitely kill horizontal rigid displacements
                       // 1 = Kill rigid displacements using a global penalization
                       // 2 = Add a Dirichlet condition on the top of the structure
penalty_parameter = 1E-6; // Penalization coefficient for the global penalization

if (d == 2) then
    cpoints = [0, 0];   // constraigned points for 2d
    cunitv  = [1, 0];   // corresponding constraigned directions for 2d
else
    cpoints = [0, 0, 0,   0, 0, 0,   5, 0, 5];  // constraigned points for 3d
    cunitv  = [1, 0, 0,   0, 1, 0,   0, 1, 0];  // corresponding constraigned directions for 3d
end

niter     = 100;   // Maximum number of iterations for Newton's algorithm.
plot_mesh = %t;
version = 13; // 1 : frictionless contact and the basic contact brick
              // 2 : contact with 'static' Coulomb friction and basic contact brick
              // 3 : frictionless contact and the contact with a
              //     rigid obstacle brick
              // 4 : contact with 'static' Coulomb friction and the contact with a
              //     rigid obstacle brick
              // 5 : frictionless contact and the integral brick
              //     Newton and Alart-Curnier augmented lagrangian,
              //     unsymmetric version
              // 6 : frictionless contact and the integral brick
              //     Newton and Alart-Curnier augmented lagrangian, symmetric
              //     version.
              // 7 : frictionless contact and the integral brick
              //     Newton and Alart-Curnier augmented lagrangian,
              //     unsymmetric version with an additional augmentation.
              // 8 : frictionless contact and the integral brick
              //     New unsymmetric method.
              // 9 : frictionless contact and the integral brick : Uzawa
              //     on the Lagrangian augmented by the penalization term.
              // 10 : contact with 'static' Coulomb friction and the integral brick
              //     Newton and Alart-Curnier augmented lagrangian,
              //     unsymmetric version.
              // 11 : contact with 'static' Coulomb friction and the integral brick
              //     Newton and Alart-Curnier augmented lagrangian,
              //     nearly symmetric version.
              // 12 : contact with 'static' Coulomb friction and the integral brick
              //     Newton and Alart-Curnier augmented lagrangian,
              //     unsymmetric version with an additional augmentation.
              // 13 : contact with 'static' Coulomb friction and the integral brick
              //     New unsymmetric method.
              // 14 : contact with 'static' Coulomb friction and the integral brick : Uzawa
              //     on the Lagrangian augmented by the penalization term.
              // 15 : penalized contact with 'static' Coulomb friction (r is the penalization
              //     coefficient).
 // Signed distance representing the obstacle
if (d == 2) then obstacle = 'y'; else obstacle = 'z'; end;

// Selection of the contact and Dirichlet boundaries
GAMMAC = 1; GAMMAD = 2;

border  = gf_mesh_get(m,'outer faces');
normals = gf_mesh_get(m, 'normal of faces', border);
contact_boundary = border(:, find(normals(d, :) < -0.01));
gf_mesh_set(m, 'region', GAMMAC, contact_boundary);
contact_boundary = border(:, find(normals(d, :) > 0.01));
gf_mesh_set(m, 'region', GAMMAD, contact_boundary);

// Finite element methods
u_degree      = 2;
lambda_degree = 2;

mfu = gf_mesh_fem(m, d);
gf_mesh_fem_set(mfu, 'classical fem', u_degree);
mfd = gf_mesh_fem(m, 1);
gf_mesh_fem_set(mfd, 'classical fem', u_degree);
mflambda = gf_mesh_fem(m, 1); // used only by versions 5 to 13
gf_mesh_fem_set(mflambda, 'classical fem', lambda_degree);
mfvm = gf_mesh_fem(m, 1);
gf_mesh_fem_set(mfvm, 'classical discontinuous fem', u_degree-1);

// Integration method
mim = gf_mesh_im(m, 4);
if (d == 2) then
  mim_friction = gf_mesh_im(m, ...
      gf_integ('IM_STRUCTURED_COMPOSITE(IM_TRIANGLE(4),4)'));
else
   mim_friction = gf_mesh_im(m, ...
      gf_integ('IM_STRUCTURED_COMPOSITE(IM_TETRAHEDRON(5),4)')); 
end

// Plot the mesh
if (plot_mesh) then
  scf(1); clf();
  gf_plot_mesh(m, 'regions', [GAMMAC]);
  title('Mesh and contact boundary (in red)');
  sleep(100);
end

// Volumic density of force
nbdofd = gf_mesh_fem_get(mfd, 'nbdof');
nbdofu = gf_mesh_fem_get(mfu, 'nbdof');
F = zeros(nbdofd*d, 1);
F(d:d:nbdofd*d) = -vertical_force;

// Elasticity model
md = gf_model('real');
gf_model_set(md, 'add fem variable', 'u', mfu);
gf_model_set(md, 'add initialized data', 'cmu', [cmu]);
gf_model_set(md, 'add initialized data', 'clambda', [clambda]);
gf_model_set(md, 'add isotropic linearized elasticity brick', mim, 'u', ...
                 'clambda', 'cmu');
gf_model_set(md, 'add initialized fem data', 'volumicload', mfd, F);
gf_model_set(md, 'add source term brick', mim, 'u', 'volumicload');

if (condition_type == 2) then
  Ddata = zeros(1, d); Ddata(d) = -5;
  gf_model_set(md, 'add initialized data', 'Ddata', Ddata);
  gf_model_set(md, 'add Dirichlet condition with multipliers', mim, 'u', u_degree, GAMMAD, 'Ddata');
elseif (condition_type == 0) then
  gf_model_set(md, 'add initialized data', 'cpoints', cpoints);
  gf_model_set(md, 'add initialized data', 'cunitv', cunitv);
  gf_model_set(md, 'add pointwise constraints with multipliers', 'u', 'cpoints', 'cunitv');
elseif (condition_type == 1) then
  // Small penalty term to avoid rigid motion (should be replaced by an
  // explicit treatment of the rigid motion with a constraint matrix)
  gf_model_set(md, 'add initialized data', 'penalty_param', ...
              [penalty_parameter]);          
  gf_model_set(md, 'add mass brick', mim, 'u', 'penalty_param');
end

// The contact condition

cdof = gf_mesh_fem_get(mfu, 'dof on region', GAMMAC);
nbc  = size(cdof, 2) / d;

if (nbc <= 0) then
    disp('No contact zone');
    return;
end

solved = %f;
if (version == 1 | version == 2) then // defining the matrices BN and BT by hand
  contact_dof   = cdof(d:d:nbc*d);
  contact_nodes = gf_mesh_fem_get(mfu, 'basic dof nodes', contact_dof);
  BN   = spzeros(nbc, nbdofu);
  ngap = zeros(nbc, 1);
  for i = 1:nbc
    BN(i, contact_dof(i)) = -1.0;
    ngap(i) = contact_nodes(d, i);
  end

  if (version == 2) then
    BT = spzeros(nbc*(d-1), nbdofu);
    for i = 1:nbc
      for j = 1:d-1
        BT(j+(i-1)*(d-1), contact_dof(i)-d+j) = 1.0;
      end
    end
  end

  gf_model_set(md, 'add variable', 'lambda_n', nbc);
  gf_model_set(md, 'add initialized data', 'r', [r]);
  if (version == 2) then
    gf_model_set(md, 'add variable', 'lambda_t', nbc * (d-1));
    gf_model_set(md, 'add initialized data', 'friction_coeff', ...
                 [friction_coeff]);
  end
  gf_model_set(md, 'add initialized data', 'ngap', ngap);
  gf_model_set(md, 'add initialized data', 'alpha', ones(nbc, 1));
  if (version == 1) then
    gf_model_set(md, 'add basic contact brick', 'u', 'lambda_n', 'r', ...
        BN, 'ngap', 'alpha', 1);
  else
    gf_model_set(md, 'add basic contact brick', 'u', 'lambda_n', ...
		 'lambda_t', 'r', BN, BT, 'friction_coeff', 'ngap', 'alpha', 1);
  end
elseif (version == 3 | version == 4) then // BN and BT defined by contact brick

  gf_model_set(md, 'add variable', 'lambda_n', nbc);
  gf_model_set(md, 'add initialized data', 'r', [r]);
  if (version == 3) then
    gf_model_set(md, 'add nodal contact with rigid obstacle brick', mim, 'u', ...
	         'lambda_n', 'r', GAMMAC, obstacle, 1);
  else
    gf_model_set(md, 'add variable', 'lambda_t', nbc * (d-1));
    gf_model_set(md, 'add initialized data', 'friction_coeff', ...
		 [friction_coeff]);
    gf_model_set(md, 'add nodal contact with rigid obstacle brick', mim, 'u', ...
	         'lambda_n', 'lambda_t', 'r', 'friction_coeff', GAMMAC, ...
		 obstacle, 1);
  end

elseif (version >= 5 & version <= 8) then // The integral version, Newton
 
  ldof = gf_mesh_fem_get(mflambda, 'dof on region', GAMMAC);
  mflambda_partial = gf_mesh_fem('partial', mflambda, ldof);
  gf_model_set(md, 'add fem variable', 'lambda_n', mflambda_partial);
  gf_model_set(md, 'add initialized data', 'r', [r]);
  OBS = gf_mesh_fem_get_eval(mfd, list(obstacle));
  gf_model_set(md, 'add initialized fem data', 'obstacle', mfd, OBS);
  gf_model_set(md, 'add integral contact with rigid obstacle brick', ...
      mim_friction, 'u', 'lambda_n', 'obstacle', 'r', GAMMAC, version-4);
          
elseif (version == 9) then // The integral version, Uzawa on the augmented Lagrangian
    
  ldof = gf_mesh_fem_get(mflambda, 'dof on region', GAMMAC);
  mflambda_partial = gf_mesh_fem('partial', mflambda, ldof);
  nbc = gf_mesh_fem_get(mflambda_partial, 'nbdof');
  OBS = gf_mesh_fem_get_eval(mfd, list(obstacle));
  M = gf_asm('mass matrix', mim, mflambda_partial, mflambda_partial, GAMMAC);
  lambda_n = zeros(1, nbc);
  gf_model_set(md, 'add initialized fem data', 'lambda_n', mflambda_partial, lambda_n);
  gf_model_set(md, 'add initialized data', 'r', [r]);
  OBS = gf_mesh_fem_get_eval(mfd, list(obstacle));
  gf_model_set(md, 'add initialized fem data', 'obstacle', mfd, OBS);
  gf_model_set(md, 'add penalized contact with rigid obstacle brick', mim_friction, 'u', ...
	         'obstacle', 'r', GAMMAC, 'lambda_n', 2);
  
  for ii=1:100
      printf('iteration %d\n', ii);
      gf_model_get(md, 'solve', 'max_res', 1E-9, 'max_iter', niter); // , 'very noisy');
      U = gf_model_get(md, 'variable', 'u');
      lambda_n_old = lambda_n;
      lambda_n = (M\ gf_asm('integral contact Uzawa projection', GAMMAC, mim_friction, mfu, U, mflambda_partial, lambda_n, mfd, OBS, r))';
      gf_model_set(md, 'variable', 'lambda_n', lambda_n);
      difff = max(abs(lambda_n-lambda_n_old));
      printf('diff : %g\n', difff/max(abs(lambda_n)));
      // pause;
      if (difff/max(abs(lambda_n)) < penalty_parameter) then break; end;
  end
  
  solved = %t;
  
elseif (version >= 10 & version <= 13) then // The integral version with friction, Newton
 
  gf_mesh_fem_set(mflambda, 'qdim', d);
  ldof = gf_mesh_fem_get(mflambda, 'dof on region', GAMMAC);
  mflambda_partial = gf_mesh_fem('partial', mflambda, ldof);
  gf_model_set(md, 'add fem variable', 'lambda', mflambda_partial);
  gf_model_set(md, 'add initialized data', 'r', [r]);
  gf_model_set(md, 'add initialized data', 'friction_coeff', [friction_coeff]);
  OBS = gf_mesh_fem_get_eval(mfd, list(obstacle));
  gf_model_set(md, 'add initialized fem data', 'obstacle', mfd, OBS);
  gf_model_set(md, 'add integral contact with rigid obstacle brick', mim_friction, 'u', ...
	         'lambda', 'obstacle', 'r', 'friction_coeff', GAMMAC, version-9);

elseif (version == 14) then // The integral version, Uzawa on the augmented Lagrangian with friction
  
  gf_mesh_fem_set(mflambda, 'qdim', d);
  ldof = gf_mesh_fem_get(mflambda, 'dof on region', GAMMAC);
  mflambda_partial = gf_mesh_fem('partial', mflambda, ldof);
  nbc = gf_mesh_fem_get(mflambda_partial, 'nbdof');
  gf_model_set(md, 'add initialized data', 'friction_coeff', [friction_coeff]);
  OBS = gf_mesh_fem_get_eval(mfd, list(obstacle));
  M = gf_asm('mass matrix', mim, mflambda_partial, mflambda_partial, GAMMAC);
  lambda = zeros(1, nbc);
  gf_model_set(md, 'add initialized fem data', 'lambda', mflambda_partial, lambda);
  gf_model_set(md, 'add initialized data', 'r', [r]);
  OBS = gf_mesh_fem_get_eval(mfd, list(obstacle));
  gf_model_set(md, 'add initialized fem data', 'obstacle', mfd, OBS);
  gf_model_set(md, 'add penalized contact with rigid obstacle brick', mim_friction, 'u', ...
	         'obstacle', 'r', 'friction_coeff', GAMMAC, 2, 'lambda');
  
  for ii=1:100
      printf('iteration %d\n', ii);
      gf_model_get(md, 'solve', 'max_res', 1E-9, 'max_iter', niter); // , 'very noisy');
      U = gf_model_get(md, 'variable', 'u');
      lambda_old = lambda;
      lambda = (M\ gf_asm('integral contact Uzawa projection', GAMMAC, mim_friction, mfu, U, mflambda_partial, lambda, mfd, OBS, r, friction_coeff))';
      gf_model_set(md, 'variable', 'lambda', lambda);
      difff = max(abs(lambda-lambda_old));
      printf('diff : %g\n', difff/max(abs(lambda)));
      // pause;
      if (difff/max(abs(lambda)) < penalty_parameter) then break; end;
  end
  
  solved = %t;

elseif (version == 15) then
 
  gf_model_set(md, 'add initialized data', 'r', [r]);
  gf_model_set(md, 'add initialized data', 'friction_coeff', [friction_coeff]);
  OBS = gf_mesh_fem_get_eval(mfd, list(obstacle));
  gf_model_set(md, 'add initialized fem data', 'obstacle', mfd, OBS);
  gf_model_set(md, 'add penalized contact with rigid obstacle brick', mim_friction, 'u', ...
	         'obstacle', 'r', 'friction_coeff', GAMMAC);
    
else
  error('Inexistent version');
end

// Solve the problem
if (~solved) then
  gf_model_get(md, 'solve', 'max_res', 1E-9, 'very noisy', 'max_iter', niter); // ,  'lsearch', 'simplest'); // , 'with pseudo potential');
end

U = gf_model_get(md, 'variable', 'u');
// lambda_n = gf_model_get(md, 'variable', 'lambda_n');
VM = gf_model_get(md, 'compute_isotropic_linearized_Von_Mises_or_Tresca', ...
		  'u', 'clambda', 'cmu', mfvm);
    

// set a custom colormap
// r=[0.7 .7 .7]; l = r($,:); s=63; s1=20; s2=25; s3=48;s4=55; for i=1:s, c1 = max(min((i-s1)/(s2-s1),1),0);c2 = max(min((i-s3)/(s4-s3),1),0); r($+1,:)=(1-c2)*((1-c1)*l + c1*[1 0 0]) + c2*[1 .8 .2]; end; colormap(r);

h = scf(2); clf();
h.color_map = gf_colormap('chouette');

if (d == 3) then
  c = [0.1;0;20]; 
  x = [1;0;0]; 
  y = [0;1;0];
  z = [0;0;1];
  // Whole boundary
  // sl2 = gf_slice(list('boundary',list('none')), m, 5);
  // Slice, 3 planes
  // sl2 = gf_slice(list('boundary',list('union',list('planar',+1,c,x),list('planar',+1,c,y),list('planar',+1,c,z))),m,5);
  // Slice, 2 planes
  // sl2 = gf_slice(list('boundary',list('union',list('planar',+1,c,x),list('planar',+1,c,y))),m,5);
  // Slice, 1 plane
  sl2 = gf_slice(list('boundary',list('planar',+1,c,x)), m, 5);

  P = gf_slice_get(sl2,'pts'); dP = gf_compute(mfu,U,'interpolate on',sl2);
  gf_slice_set(sl2, 'pts', P+dP);
  VMsl = gf_compute(mfvm,VM,'interpolate on',sl2);
  //set(gcf,'renderer','zbuffer');
  gf_plot_slice(sl2,'mesh','on','mesh_slice_edges','off','data',VMsl);
  //view(-80,-15); axis on; camlight;
  // map=[1:-1/10:0]'*[1 1 1]; colormap(map); // for NB
    
  // gf_plot(mfvm, VM, 'mesh', 'off', 'cvlst', ...
  //        gf_mesh_get(mfu,'outer faces'), 'deformation', U, ...
  //        'deformation_mf', mfu, 'deformation_scale', 1, 'refine', 8);
  // view(-5,-10); camlight; colormap(map);
  xlabel('x'); ylabel('y'); zlabel('z');
  title('Sliced deformed configuration (not really a small deformation of course ...)');
else
  gf_plot(mfvm, VM, 'deformed_mesh', 'on', 'deformation', U, ...
          'deformation_mf', mfu, 'deformation_scale', 1, 'refine', 8);
  xlabel('x'); ylabel('y');
  title('Deformed configuration (not really a small deformation of course ...)');
end

//colorbar;
sleep(100);
