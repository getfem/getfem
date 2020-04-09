% Copyright (C) 2005-2020 Julien Pommier.
%
% This file is a part of GetFEM
%
% GetFEM  is  free software;  you  can  redistribute  it  and/or modify it
% under  the  terms  of the  GNU  Lesser General Public License as published
% by  the  Free Software Foundation;  either version 3 of the License,  or
% (at your option) any later version along with the GCC Runtime Library
% Exception either version 3.1 or (at your option) any later version.
% This program  is  distributed  in  the  hope  that it will be useful,  but
% WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
% or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
% License and GCC Runtime Library Exception for more details.
% You  should  have received a copy of the GNU Lesser General Public License
% along  with  this program;  if not, write to the Free Software Foundation,
% Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.


disp('This is the "legacy" getfem-matlab demonstration')
disp('This demo does not use the model bricks')
disp('instead it show how the linear system is built with direct calls')
disp('to the assembly routines.')

gf_workspace('clear all');
% import the mesh
m=gfMesh('import','gid','../meshes/tripod.GiD.msh');
mfu=gfMeshFem(m,3);     % mesh-fem supporting a 3D-vector field
mfd=gfMeshFem(m,1);     % scalar mesh_fem
% the mesh_im stores the integration methods for each tetrahedron
mim=gfMeshIm(m,gf_integ('IM_TETRAHEDRON(5)'));
% we choose a P2 fem for the main unknown
set(mfu,'fem',gf_fem('FEM_PK(3,2)'));
% the material is homogeneous, hence we use a P0 fem for the data
set(mfd,'fem',gf_fem('FEM_PK(3,0)'));
% display some informations about the mesh
disp(sprintf('nbcvs=%d, nbpts=%d, nbdof=%d',get(m,'nbcvs'),...
             get(m,'nbpts'),get(mfu,'nbdof')));
P=get(m,'pts'); % get list of mesh points coordinates
pidtop=find(abs(P(2,:)-13)<1e-6); % find those on top of the object
pidbot=find(abs(P(2,:)+10)<1e-6); % find those on the bottom
% build the list of faces from the list of points
ftop=get(m,'faces from pid',pidtop); 
fbot=get(m,'faces from pid',pidbot);
% assign boundary numbers
set(m,'boundary',1,ftop);
set(m,'boundary',2,fbot);

E=1e3;
nu=0.3;
lambda = E*nu/((1+nu)*(1-2*nu));
mu =E/(2*(1+nu));
nbd=get(mfd, 'nbdof');
F = gf_asm('boundary_source', 1, mim, mfu, mfd, repmat([0;-10;0],1,nbd));
K = gf_asm('linear_elasticity', mim, mfu, mfd, ...
	   lambda*ones(1,nbd),mu*ones(1,nbd));

% handle Dirichlet condition
[H,R]=gf_asm('dirichlet', 2, mim, mfu, mfd, repmat(eye(3),[1,1,nbd]), zeros(3, nbd));
[N,U0]=gf_spmat_get(H, 'dirichlet_nullspace', R);
KK=N'*K*N;
FF=N'*F;
% solve ...
disp('solving...'); t0 = cputime;
lsolver = 1 % change this to compare the different solvers
if (lsolver == 1),     % conjugate gradient
  P=gfPrecond('ildlt',KK);
  UU=gf_linsolve('cg',KK,FF,P,'noisy','res',1e-9);
elseif (lsolver == 2), % superlu
  UU=gf_linsolve('superlu',KK,FF);
else                   % the matlab "slash" operator 
  UU=KK\FF;
end;
disp(sprintf('linear system solved in %.2f sec', cputime-t0));
U=(N*UU).'+U0;

% now that we have the solution, we want to compute the von mises stress
% first, we need to get the derivate of the solution:
mfdu=gfMeshFem(m,1);
% the P2 fem is not derivable across elements, hence we use a discontinuous
% fem for the derivative of U.
set(mfdu,'fem',gf_fem('FEM_PK_DISCONTINUOUS(3,1)'));

% on output size(DU)=[3,3,nbdof(mfdu)]
DU=gf_compute(mfu,U,'gradient',mfdu);

% from the derivative, we compute the von mises stress
VM=zeros(1,get(mfdu,'nbdof'));
N=gf_mesh_get(m,'dim');
for i=1:size(DU,3),
  t=DU(:,:,i);
  E=(t+t')/2;
  VM(i) = sum(E(:).^2) - (1./N)*sum(diag(E))^2;
end;
VM = 4*mu^2*VM;

disp('plotting ... can also take some minutes!');

% we plot the von mises on the deformed object, in superposition with the initial mesh.
gf_plot(mfdu,VM,'mesh','on', 'cvlst', get(m, 'outer faces'),...
	'deformation',U,'deformation_mf',mfu);
caxis([0 100]); colorbar; view(180,-50); camlight;

r=[0.7 .7 .7]; l = r(end,:); s=63; s1=20; s2=25; s3=48;s4=55; for i=1:s, c1 = max(min((i-s1)/(s2-s1),1),0);c2 = max(min((i-s3)/(s4-s3),1),0); r(end+1,:)=(1-c2)*((1-c1)*l + c1*[1 0 0]) + c2*[1 .8 .2]; end; colormap(r); colorbar;

