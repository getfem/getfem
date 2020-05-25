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

printf('\nThis is the ''legacy'' getfem-scilab demonstration.\n')
printf('This demo does not use the model bricks introduced with getfem 2.0.\n')
printf('Instead it shows how the linear system is built with direct calls\n')
printf('to the assembly routines.\n')

lines(0);
stacksize('max');

path = get_absolute_file_path('demo_tripod_alt.sce');

if getos()=='Windows' then
  // Under Windows, all the trace messages are available in the dos console
  // Under Linuxs, all the trace messages are redirected to the Scilab console
  consolebox('on');
end
gf_util('trace level',3);
gf_util('warning level',3);

gf_workspace('clear all');

printf('demo tripod_alt started\n');

// import the mesh
m   = gf_mesh('import','gid', path + '/data/tripod.GiD.msh');
mfu = gf_mesh_fem(m,3);     // mesh-fem supporting a 3D-vector field
mfd = gf_mesh_fem(m,1);     // scalar mesh_fem

// the mesh_im stores the integration methods for each tetrahedron
mim = gf_mesh_im(m,gf_integ('IM_TETRAHEDRON(5)'));

// we choose a P2 fem for the main unknown
gf_mesh_fem_set(mfu,'fem',gf_fem('FEM_PK(3,2)'));

// the material is homogeneous, hence we use a P0 fem for the data
gf_mesh_fem_set(mfd,'fem',gf_fem('FEM_PK(3,0)'));

// display some informations about the mesh
disp(sprintf('nbcvs=%d, nbpts=%d, nbdof=%d',gf_mesh_get(m,'nbcvs'), gf_mesh_get(m,'nbpts'),gf_mesh_fem_get(mfu,'nbdof')));
P = gf_mesh_get(m,'pts'); // get list of mesh points coordinates
pidtop = find(abs(P(2,:)-13)<1e-6); // find those on top of the object
pidbot = find(abs(P(2,:)+10)<1e-6); // find those on the bottom

// build the list of faces from the list of points
ftop = gf_mesh_get(m,'faces from pid',pidtop);
fbot = gf_mesh_get(m,'faces from pid',pidbot);

// assign boundary numbers
gf_mesh_set(m,'boundary',1,ftop);
gf_mesh_set(m,'boundary',2,fbot);

E  = 1e3;
nu = 0.3;
lambda = E*nu/((1+nu)*(1-2*nu));
mu     = E/(2*(1+nu));
nbd    = gf_mesh_fem_get(mfd, 'nbdof');
F = gf_asm('boundary_source', 1, mim, mfu, mfd, repmat([0;-10;0],1,nbd));
K = gf_asm('linear_elasticity', mim, mfu, mfd, lambda*ones(1,nbd),mu*ones(1,nbd));

// handle Dirichlet condition
[H,R]  = gf_asm('dirichlet', 2, mim, mfu, mfd, ones(1,1,nbd) .*. eye(3,3), zeros(3, nbd));
[N,U0] = gf_spmat_get(H, 'dirichlet_nullspace', R);

// N:        nnz(N) = 16341   size(N) = 16764 x 16341
// K:        nnz(K) = 1147742 size(K) = 16764 x 16764
// A = K*N:  nnz(A) = 1123597 size(A) = 16764 x 16341
// B = N'*A: nnz(B) = 1110396 size(B) = 16341 x 16341

// KK = N'*K*N; // This computation doesn't fit in the scilab stack. I must split it into parts
//K = K*N;
//K = N'*K;
K = N'*K*N;

F = N'*F;

// solve ...
//sleep(100); // bug with timer
t_start = timer();

disp('solving...'); 
lsolver = 3; // change this to compare the different solvers

if (lsolver == 1) then   // conjugate gradient
  P  = gf_precond('ildlt',K);
  UU = gf_linsolve('cg',K,F,P,'noisy','res',1e-9);
elseif (lsolver == 2) then // superlu
  UU = gf_linsolve('superlu',K,F);
elseif (lsolver == 3) then
  UU = umfpack(K, "\", F);
elseif (lsolver == 4) then
  UU = linsolve(K, F);
else                   // the scilab "slash" operator 
  UU = K\F;
end
t_end = timer();

disp(sprintf('linear system solved in %f sec', t_end-t_start));

U = (N*UU).'+U0;

// now that we have the solution, we want to compute the von mises stress
// first, we need to get the derivate of the solution:
mfdu = gf_mesh_fem(m,1);
// the P2 fem is not derivable across elements, hence we use a discontinuous
// fem for the derivative of U.
gf_mesh_fem_set(mfdu,'fem',gf_fem('FEM_PK_DISCONTINUOUS(3,1)'));

// on output size(DU)=[3,3,nbdof(mfdu)]
DU = gf_compute(mfu,U,'gradient',mfdu);

// from the derivative, we compute the von mises stress
VM = zeros(1,gf_mesh_fem_get(mfdu,'nbdof'));
N  = gf_mesh_get(m,'dim');
for i=1:size(DU,3)
  t     = DU(:,:,i);
  E     = (t+t')/2;
  VM(i) = sum(E(:).^2) - ((1) ./ N)*sum(diag(E),'r')^2;
end
VM = 4*mu^2*VM;

disp('plotting ... can also take some minutes!');

h = scf();

//r = [0.7 .7 .7]; l = r($,:); s=63; s1=20; s2=25; s3=48;s4=55; 
//for i=1:s
//  c1 = max(min((i-s1)/(s2-s1),1),0);
//  c2 = max(min((i-s3)/(s4-s3),1),0); 
//  r($+1,:)=(1-c2)*((1-c1)*l + c1*[1 0 0]) + c2*[1 .8 .2]; 
//end
//h.color_map = r;
h.color_map = jetcolormap(255);

// we plot the von mises on the deformed object, in superposition with the initial mesh.
drawlater;
gf_plot(mfdu,VM,'mesh','on', 'cvlst', gf_mesh_get(m, 'outer faces'), 'deformation',U,'deformation_mf',mfu);

h.children.rotation_angles = [135 75];
a = gca();
a.box = 'off';
a.axes_visible = 'off';
a.x_label.visible = 'off';
a.y_label.visible = 'off';
a.z_label.visible = 'off';
colorbar(min(VM),max(VM));
drawnow;

printf('demo tripod_alt terminated\n');
