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

path = get_absolute_file_path('demo_wave2D.sce');

printf('demo wave2D started\n');

if getos()=='Windows' then
  // Under Windows, all the trace messages are available in the dos console
  // Under Linuxs, all the trace messages are redirected to the Scilab console
  consolebox('on');
end
gf_util('trace level',3);
gf_util('warning level',3);

gf_workspace('clear all');

disp('2D scalar wave equation (helmholtz) demonstration');
printf(' we present three approaches for the solution of the helmholtz problem\n')
printf(' - the first one is to use the new getfem ''model bricks''\n')
printf(' - the second one is to use the old getfem ''model bricks''\n')
printf(' - the third one is to use the ''low level'' approach, i.e. to assemble\n')
printf('   and solve the linear systems.\n')

printf('The result is the wave scattered by a disc, the incoming wave beeing a plane wave coming from the top\n');
printf(' \delta u + k^2 = 0\n');
printf(' u = -uinc              on the interior boundary\n');
printf(' \partial_n u + iku = 0 on the exterior boundary\n');

//PK = 10; gt_order = 6; k = 7; use_hierarchical = 0; load_the_mesh=0;
PK       = 3; 
gt_order = 3; 
k        = 1; 
use_hierarchical = 1; 
load_the_mesh    = 1;

if (use_hierarchical) then
  s = 'hierarchical'; 
else 
  s = 'classical'; 
end

disp(sprintf('using %s P%d FEM with geometric transformations of degree %d',s,PK,gt_order));

if (load_the_mesh) then
  disp('the mesh is loaded from a file, gt_order ignored');
end

if load_the_mesh == 0 then
  // a quadrangular mesh is generated, with a high degree geometric transformation
  // number of cells for the regular mesh
  Nt = 10; 
  Nr = 8;
  m  = gf_mesh('empty',2);
  dtheta  = 2*%pi*1/Nt; R=1+9*(0:Nr-1)/(Nr-1);
  gt      = gf_geo_trans(sprintf('GT_PRODUCT(GT_PK(1,%d),GT_PK(1,1))',gt_order));
  ddtheta = dtheta/gt_order;
  for i=1:Nt
    for j=1:Nr-1
      ti=(i-1)*dtheta:ddtheta:i*dtheta;
      X = [R(j)*cos(ti) R(j+1)*cos(ti)];
      Y = [R(j)*sin(ti) R(j+1)*sin(ti)];
      gf_mesh_set(m,'add convex',gt,[X;Y]);
    end
  end
  fem_u = gf_fem(sprintf('FEM_QK(2,%d)',PK));
  fem_d = gf_fem(sprintf('FEM_QK(2,%d)',PK));
  mfu = gf_mesh_fem(m,1);
  mfd = gf_mesh_fem(m,1);  
  gf_mesh_fem_set(mfu'fem',fem_u);
  gf_mesh_fem_set(mfd'fem',fem_d);
  sIM = sprintf('IM_GAUSS_PARALLELEPIPED(2,%d)',gt_order+2*PK);
  mim = gf_mesh_im(m, g_integ(sIM));
else
  // the mesh is loaded
  m = gf_mesh('import','gid',path + 'data/holed_disc_with_quadratic_2D_triangles.msh');
  if (use_hierarchical) then
    // hierarchical basis improve the condition number
    // of the final linear system
    fem_u = gf_fem(sprintf('FEM_PK_HIERARCHICAL(2,%d)',PK));
    //fem_u=gf_fem('FEM_HCT_TRIANGLE');
    //fem_u=gf_fem('FEM_HERMITE(2)');
  else
    fem_u = gf_fem(sprintf('FEM_PK(2,%d)',PK));
  end
  fem_d = gf_fem(sprintf('FEM_PK(2,%d)',PK));
  mfu   = gf_mesh_fem(m,1);
  mfd   = gf_mesh_fem(m,1);  
  gf_mesh_fem_set(mfu,'fem',fem_u);
  gf_mesh_fem_set(mfd,'fem',fem_d);
  mim = gf_mesh_im(m,gf_integ('IM_TRIANGLE(13)'));
end
nbdu = gf_mesh_fem_get(mfu,'nbdof');
nbdd = gf_mesh_fem_get(mfd,'nbdof');

// identify the inner and outer boundaries
P = gf_mesh_get(m,'pts'); // get list of mesh points coordinates
pidobj = find(sum(P.^2,'r') < 1*1+1e-6);
pidout = find(sum(P.^2,'r') > 10*10-1e-2);

// build the list of faces from the list of points
fobj = gf_mesh_get(m,'faces from pid',pidobj); 
fout = gf_mesh_get(m,'faces from pid',pidout);
gf_mesh_set(m,'boundary',1,fobj);
gf_mesh_set(m,'boundary',2,fout);

// expression of the incoming wave
disp(k)
wave_expr = sprintf('cos(%f*y+.2)+1*%%i*sin(%f*y+.2)',k,k);
Uinc      = gf_mesh_fem_get_eval(mfd,list(list(wave_expr)));

//
// we present two approaches for the solution of the Helmholtz problem
// - the first one is to use the new getfem "model bricks"
// - the second one is to use the "low level" approach, i.e. to assemble
//   and solve the linear systems.
if 1 then
  timer();
  // solution using new model bricks
  md = gf_model('complex');
  gf_model_set(md, 'add fem variable', 'u', mfu);
  gf_model_set(md, 'add initialized data', 'k', [k]);
  gf_model_set(md, 'add Helmholtz brick', mim, 'u', 'k');
  gf_model_set(md, 'add initialized data', 'Q', [1*%i*k]);
  gf_model_set(md, 'add Fourier Robin brick', mim, 'u', 'Q', 2);
  gf_model_set(md, 'add initialized fem data', 'DirichletData', mfd, Uinc);
  gf_model_set(md, 'add Dirichlet condition with multipliers', mim, 'u', mfd, 1, 'DirichletData');
  // gf_model_set(md, 'add Dirichlet condition with penalization', mim, 'u', 1e12, 1, 'DirichletData');

  gf_model_get(md, 'solve');
  U = gf_model_get(md, 'variable', 'u');
  disp(sprintf('solve done in %.2f sec', timer()));
else
  // solution using the "low level" approach
  [H,R] = gf_asm('dirichlet', 1, mim, mfu, mfd, gf_mesh_fem_get(mfd,'eval',1),Uinc);
  [_null,ud] = gf_spmat_get(H,'dirichlet nullspace', R);
  
  Qb2 = gf_asm('boundary qu term', 2, mim, mfu, mfd, ones(1,nbdd));
  M   = gf_asm('mass matrix',mim, mfu);
  L   = -gf_asm('laplacian',mim, mfu,mfd,ones(1,nbdd));

  // builds the matrix associated to
  // (\Delta u + k^2 u) inside the domain, and 
  // (\partial_n u + ik u) on the exterior boundary
  A = L + (k*k) * M + (1*%i*k)*Qb2;


  // eliminate dirichlet conditions and solve the system
  RF = _null'*(-A*ud(:));
  RK = _null'*A*_null;
  U  = _null*(RK\RF)+ud(:);
  U  = U(:).';
end

Ud = gf_compute(mfu,U,'interpolate on',mfd);

h = scf(); 
h.color_map = jetcolormap(255);
drawlater;
gf_plot(mfu,imag(U(:)'),'mesh','on','refine',32,'contour',0); 
colorbar(min(imag(U)),max(imag(U)));
h.color_map = jetcolormap(255);
drawnow;

h = scf(); 
h.color_map = jetcolormap(255);
drawlater;
gf_plot(mfd,abs(Ud(:)'),'mesh','on','refine',24,'contour',0.5); 
colorbar(min(abs(Ud)),max(abs(Ud)));
h.color_map = jetcolormap(255);
drawnow;

// compute the "exact" solution from its developpement 
// of bessel functions:
// by \Sum_n c_n H^(1)_n(kr)exp(i n \theta)
N     = 1000;
theta = 2*%pi*(0:N-1)/N;
y     = sin(theta); 
w     = eval(wave_expr);
fw    = fft(w); 
C     = fw/N;
S     = zeros(w);
S(:)  = C(1);
Nc    = 20;
for i=2:Nc
  n=i-1;  
  S = S + C(i)*exp(1*%i*n*theta) + C(N-(n-1))*exp(-1*%i*n*theta);
end
P = gf_mesh_fem_get(mfd,'basic dof nodes');
[T,R] = cart2pol(P(1,:),P(2,:));
Uex   = zeros(size(R));
nbes  = 1;
Uex   = besselh(0,nbes,k*R) * C(1)/besselh(0,nbes,k);
old_ieee = ieee();
ieee(2);
for i=2:Nc
  n   = i-1;  
  Uex = Uex + besselh(n,nbes,k*R) * C(i)/besselh(n,nbes,k) .* exp(1*%i*n*T);
  Uex = Uex + besselh(-n,nbes,k*R) * C(N-(n-1))/besselh(-n,nbes,k) .* exp(-1*%i*n*T);
end
ieee(old_ieee);

disp('the error won''t be less than ~1e-2 as long as a first order absorbing boundary condition will be used');
disp(sprintf('rel error ||Uex-U||_inf=%g',max(abs(Ud-Uex))/max(abs(Uex))));
disp(sprintf('rel error ||Uex-U||_L2=%g', gf_compute(mfd,Uex-Ud,'L2 norm',mim)/gf_compute(mfd,Uex,'L2 norm',mim)));
disp(sprintf('rel error ||Uex-U||_H1=%g', gf_compute(mfd,Uex-Ud,'H1 norm',mim)/gf_compute(mfd,Uex,'H1 norm',mim)));

h = scf();
h.color_map = jetcolormap(255);
// adjust the 'refine' parameter to enhance the quality of the picture
drawlater;
gf_plot(mfu,real(U(:)'),'mesh','on','refine',8);
h.color_map = jetcolormap(255);
drawnow;

printf('demo wave2D terminated\n');
