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
stacksize('max');

gf_workspace('clear all');


path = get_absolute_file_path('demo_plasticity.sce');

if getos()=='Windows' then
  // Under Windows, all the trace messages are available in the dos console
  // Under Linuxs, all the trace messages are redirected to the Scilab console
  consolebox('on');
end
gf_util('trace level',3);
gf_util('warning level',3);



// We compute a plasticity problem with a Von Mises criterion with or
// without kinematic hardening
// For convenience we consider an homogenous Dirichlet condition on the left
// of the domain and an easy computed Neumann Condition on the right


with_hardening = 1;
bi_material = %f;
test_tangent_matrix = %f;
do_plot = %t;



// Initialize used data
LX = 100;
LY = 20;
NX = 50;
NY = 20;

// alpha is parameter of the generalized integration algorithms.
// The choice alpha = 1/2 yields the mid point method and alpha = 1 leads to
// backward Euler integration
alpha = 1.0;





f = [0 -600]';
t = [0 0.5 0.6 0.7 0.8 0.9 1 0.9 0.8 0.7 0.6 0.5 0.4 0.3 0.2 0.1 0];
if (with_hardening == 1)
  f = [15000 0]';
  t = [0 0.5 0.6 0.7 0.8 0.9 1 0.9 0.8 0.7 0.6 0.5 0.4 0.3 0.2 0.1 0 -0.1 -0.2 -0.4 -0.6 -0.8 -0.6 -0.4 -0.2 0];
end

// Create the mesh
// m = gf_mesh('triangles grid', [0:(LX/NX):LX], [0:(LY/NY):LY]);
m = gf_mesh('import','structured',sprintf('GT=""GT_PK(2,1)"";SIZES=[%d,%d];NOISED=0;NSUBDIV=[%d,%d];', LX, LY, NX, NY));
N = gf_mesh_get(m, 'dim');
  
// Plotting
// gf_plot_mesh(m, 'vertices', 'on', 'convexes', 'on');

// Define used MeshIm
mim=gf_mesh_im(m);
gf_mesh_im_set(mim, 'integ', gf_integ('IM_TRIANGLE(6)')); // Gauss methods on triangles

// Define used MeshFem
if (with_hardening == 1)
  mf_u=gf_mesh_fem(m,2); gf_mesh_fem_set(mf_u, 'fem',gf_fem('FEM_PK(2,2)'));
else
  mf_u=gf_mesh_fem(m,2); gf_mesh_fem_set(mf_u, 'fem',gf_fem('FEM_PK(2,1)'));
end
mf_data=gf_mesh_fem(m); gf_mesh_fem_set(mf_data, 'fem', gf_fem('FEM_PK_DISCONTINUOUS(2,0)'));
mf_data2=gf_mesh_fem(m,2); gf_mesh_fem_set(mf_data2, 'fem', gf_fem('FEM_PK_DISCONTINUOUS(2,0)'));
// mf_sigma=gf_mesh_fem(m,4); gf_mesh_fem_set(mf_sigma, 'fem',gf_fem('FEM_PK_DISCONTINUOUS(2,1)'));
mf_sigma=gf_mesh_fem(m,4); gf_mesh_fem_set(mf_sigma, 'fem',gf_fem('FEM_PK_DISCONTINUOUS(2,0)'));
mf_vm = gf_mesh_fem(m); set(mf_vm, 'fem', gf_fem('FEM_PK_DISCONTINUOUS(2,1)'));

// Find the border of the domain
P=gf_mesh_get(m, 'pts');
pidleft=find(abs(P(1,:))<1e-6); // Retrieve index of points which x near to 0
pidright=find(abs(P(1,:) - LX)<1e-6); // Retrieve index of points which x near to L
fleft =gf_mesh_get(m,'faces from pid',pidleft); 
fright=gf_mesh_get(m,'faces from pid',pidright);
gf_mesh_set(m,'boundary',1,fleft); // for Dirichlet condition
gf_mesh_set(m,'boundary',2,fright); // for Neumann condition

// Decomposed the mesh into 2 regions with different values of Lamé coeff
if (bi_material) separation = LY/2; else separation = 0; end
pidtop    = find(P(2,:)>=separation-1E-6); // Retrieve index of points of the top part
pidbottom = find(P(2,:)<=separation+1E-6); // Retrieve index of points of the bottom part
cvidtop   = gf_mesh_get(m, 'cvid from pid', pidtop);
cvidbottom= gf_mesh_get(m, 'cvid from pid', pidbottom);
CVtop     = gsort(gf_mesh_fem_get(mf_data, 'basic dof from cvid', cvidtop));
CVbottom  = gsort(gf_mesh_fem_get(mf_data, 'basic dof from cvid', cvidbottom));

// Definition of Lame coeff
lambda(CVbottom,1) = 121150; // Steel
lambda(CVtop,1) = 84605; // Iron
mu(CVbottom,1) = 80769; //Steel
mu(CVtop,1) = 77839; // Iron
// Definition of plastic threshold
von_mises_threshold(CVbottom) = 7000;
von_mises_threshold(CVtop) = 8000;
// Definition of hardening parameter
if (with_hardening)
  H = mu(1)/5;
else
  H = 0;
end

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


  
  
if (with_hardening)
  N = gf_mesh_get(m, 'dim');
  gf_model_set(md, 'add fem data', 'Previous_u', mf_u);
  mim_data = gf_mesh_im_data(mim, -1, [N, N]);
  gf_model_set(md, 'add im data', 'sigma', mim_data);
  
  // Declare that alpha is a data of the system 
 
  gf_model_set(md, 'add initialized data', 'alpha', [alpha]);
  gf_model_set(md, 'add initialized data', 'H', [H]);

  Is = 'Reshape(Id(meshdim*meshdim),meshdim,meshdim,meshdim,meshdim)';
  IxI = 'Id(meshdim)@Id(meshdim)';
  coeff_long = '((lambda)*(H))/((2*(mu)+(H))*(meshdim*(lambda)+2*(mu)+(H)))';
  B_inv = sprintf('((2*(mu)/(2*(mu)+(H)))*(%s) + (%s)*(%s))', Is, coeff_long, IxI);
  B = sprintf('((1+(H)/(2*(mu)))*(%s) - (((lambda)*(H))/(2*(mu)*(meshdim*(lambda)+2*(mu))))*(%s))', Is, IxI);
  ApH = sprintf('((2*(mu)+(H))*(%s) + (lambda)*(%s))', Is, IxI);
  Enp1 = '((Grad_u+Grad_u'')/2)';
  En = '((Grad_Previous_u+Grad_Previous_u'')/2)';
  
  //expression of sigma for Implicit Euler method
  //expr_sigma = strcat(['(', B_inv, '*(Von_Mises_projection((-(H)*', Enp1, ')+(', ApH, '*(',Enp1,'-',En,')) + (', B, '*sigma), von_mises_threshold) + H*', Enp1, '))']);
  
  //expression of sigma for generalized alpha algorithms
  expr_sigma = strcat(['(', B_inv, '*(Von_Mises_projection((',B,'*((1-alpha)*sigma))+(-(H)*(((1-alpha)*',En,')+(alpha*', Enp1, ')))+(alpha*', ApH, '*(',Enp1,'-',En,')) + (alpha*', ...
    B, '*sigma), von_mises_threshold) + (H)*(((1-alpha)*',En,')+(alpha*', Enp1, '))))']);
  
  gf_model_set(md, 'add nonlinear generic assembly brick', mim, expr_sigma + ':Grad_Test_u');
  // gf_model_set(md, 'add finite strain elasticity brick', mim, 'u', 'SaintVenant Kirchhoff', '[lambda; mu]');
else
    
  // Declare that sigma is a data of the system on mf_sigma
  gf_model_set(md, 'add fem data', 'sigma', mf_sigma);
  // Add plasticity brick on u
  gf_model_set(md, 'add elastoplasticity brick', mim, 'VM', 'u', 'lambda', 'mu', 'von_mises_threshold', 'sigma');
end

// Add homogeneous Dirichlet condition to u on the left hand side of the domain
gf_model_set(md, 'add Dirichlet condition with multipliers', mim, 'u', mf_u, 1);

// Add a source term to the system
gf_model_set(md,'add initialized fem data', 'VolumicData', mf_data2, gf_mesh_fem_get_eval(mf_data2, list(['f(1,1)*t(1)','f(2,1)*t(1)'])));
gf_model_set(md, 'add source term brick', mim, 'u', 'VolumicData', 2);

VM=zeros(1,gf_mesh_fem_get(mf_vm, 'nbdof'));

if (do_plot)
      h = scf();
      h.color_map = jetcolormap(255);
end

for step=1:size(t,2),
    disp(sprintf('step %d / %d, coeff = %g', step, size(t,2), t(step)));
    gf_model_set(md, 'variable', 'VolumicData', gf_mesh_fem_get_eval(mf_data2, list(['f(1,1)*t(step)','f(2,1)*t(step)'])));
    
    if (test_tangent_matrix)
      gf_model_get(md, 'test tangent matrix', 1E-8, 10, 0.000001);
    end;
   
    // Solve the system
    gf_model_get(md, 'solve', 'noisy', 'lsearch', 'simplest',  'alpha min', 0.8, 'max_iter', 100, 'max_res', 1e-6);
    // gf_model_get(md, 'solve', 'noisy', 'max_iter', 80);

    // Retrieve the solution U
    U = gf_model_get(md, 'variable', 'u', 0);
    
    // Compute new plasticity constraints used to compute 
    // the Von Mises or Tresca stress
    if (with_hardening)
      sigma_0 = gf_model_get(md, 'variable', 'sigma');
      sigma = gf_model_get(md, 'interpolation', expr_sigma, mim_data);
      U_0 = gf_model_get(md, 'variable', 'Previous_u');
      U_nalpha = alpha*U + (1-alpha)*U_0;
      
      M = gf_asm('mass matrix', mim, mf_vm);
      L = gf_asm('generic', mim, 1, 'sqrt(3/2)*Norm(Deviator(sigma))*Test_vm', -1, 'sigma', 0, mim_data, sigma, 'vm', 1, mf_vm, zeros(gf_mesh_fem_get(mf_vm, 'nbdof'),1));
      VM = (M \ L)';
      coeff1='-lambda/(2*mu*(meshdim*lambda+2*mu))';
      coeff2='1/(2*mu)';
      Ainv=sprintf('(%s)*(%s) + (%s)*(%s)', coeff1, IxI, coeff2, Is);
      Ep = sprintf('(Grad_u+Grad_u'')/2 - (%s)*sigma', Ainv);
      L = gf_asm('generic', mim, 1, sprintf('Norm(%s)*Test_vm', Ep), -1, 'sigma', 0, mim_data, sigma, 'u', 0, mf_u, U, 'vm', 1, mf_vm, zeros(gf_mesh_fem_get(mf_vm, 'nbdof'),1), 'mu', 0, mf_data, mu, 'lambda', 0, mf_data, lambda);
      plast = (M \ L)';
      
      gf_model_set(md, 'variable', 'u', U_nalpha);
      Epsilon_u = gf_model_get(md, 'interpolation', '((Grad_u+Grad_u'')/2)', mim_data);
      gf_model_set(md, 'variable', 'u', U);
      ind_gauss_pt = 22500;
      if (size(sigma, 2) <= N*(ind_gauss_pt + 1))
        ind_gauss_pt = floor(3*size(sigma, 2) / (4*N*N));
      end
      sigma_fig(1,step)=sigma(N*N*ind_gauss_pt + 1);
      Epsilon_u_fig(1,step)=Epsilon_u(N*N*ind_gauss_pt + 1);
      
      sigma = (sigma - (1-alpha)*sigma_0)/alpha;
      gf_model_set(md, 'variable', 'sigma', sigma);
      gf_model_set(md, 'variable', 'Previous_u', U);
    else
      gf_model_get(md, 'elastoplasticity next iter', mim, 'u', 'VM', 'lambda', 'mu', 'von_mises_threshold', 'sigma');
      plast = gf_model_get(md, 'compute plastic part', mim, mf_vm, 'u', 'VM', 'lambda', 'mu', 'von_mises_threshold', 'sigma');
      // Compute Von Mises or Tresca stress
      VM = gf_model_get(md, 'compute elastoplasticity Von Mises or Tresca', 'sigma', mf_vm, 'Von Mises');
    end
       
       
    if (do_plot)
      drawlater;
      clf();
      subplot(3,1,1);
      gf_plot(mf_vm,VM, 'deformation',U,'deformation_mf',mf_u,'refine', 4, 'deformation_scale',1, 'disp_options', 0); // 'deformed_mesh', 'on')
      colorbar(min(VM),max(VM));
      a = get("current_axes"); a.data_bounds = [-20 120 -20 40];
      // caxis([0 10000]);
      n = t(step);
      title(sprintf('Von Mises criterion for t = %d', step));
      
      subplot(3,1,2);
      gf_plot(mf_vm,plast, 'deformation',U,'deformation_mf',mf_u,'refine', 4, 'deformation_scale',1, 'disp_options', 0);  // 'deformed_mesh', 'on')
      colorbar(min(plast),max(plast));
      a = get("current_axes"); a.data_bounds = [-20 120 -20 40];
      // caxis([0 10000]);
      n = t(step);
      title(sprintf('Plastification for t = %d', step));
    
      if (with_hardening)
        subplot(3,1,3);
        plot(Epsilon_u_fig, sigma_fig,'r','LineWidth',2)
        xlabel('Strain');
        ylabel('Stress')
        a = get("current_axes"); a.data_bounds = [-0.1 0.35 -16000 16000];
      end;
      drawnow;
      sleep(1000);
    end
 
end;










