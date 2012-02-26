% Copyright (C) 2006-2012 Yves Renard, Julien Pommier.
%
% This file is a part of GETFEM++
%
% Getfem++  is  free software;  you  can  redistribute  it  and/or modify it
% under  the  terms  of the  GNU  Lesser General Public License as published
% by  the  Free Software Foundation;  either version 2.1 of the License,  or
% (at your option) any later version.
% This program  is  distributed  in  the  hope  that it will be useful,  but
% WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
% or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
% License for more details.
% You  should  have received a copy of the GNU Lesser General Public License
% along  with  this program;  if not, write to the Free Software Foundation,
% Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.


clear all
gf_workspace('clear all');
clc


L=100; H=20;
m=gfMesh('triangles grid',0:4:L, 0:2:H);


mim=gfMeshIm(m);  set(mim, 'integ',gfInteg('IM_TRIANGLE(6)'));
mfu=gfMeshFem(m,2); set(mfu, 'fem',gfFem('FEM_PK(2,1)'));
mfd=gfMeshFem(m); set(mfd, 'fem',gfFem('FEM_PK_DISCONTINUOUS(2,1)'));
mf0=gfMeshFem(m); set(mf0, 'fem',gfFem('FEM_PK(2,0)'));
mfdu=gfMeshFem(m); set(mfdu, 'fem',gfFem('FEM_PK_DISCONTINUOUS(2,1)'));
lambda=121150;
mu=80769;
von_mises_threshold=7000;

P=get(m, 'pts');
pidleft=find(abs(P(1,:))<1e-6);
pidright=find(abs(P(1,:) - L)<1e-6);

fleft =get(m,'faces from pid',pidleft); 
fright=get(m,'faces from pid',pidright);
% assign boundary numbers
set(m,'boundary',1,fleft);
set(m,'boundary',2,fright);

b0=gfMdBrick('small deformations plasticity',mim,mfu, von_mises_threshold);
set(b0, 'param','lambda',lambda);
set(b0, 'param','mu',mu);
b1=gfMdBrick('generalized dirichlet',b0,1);
b2=gfMdBrick('source term',b1,2);

mds=gfMdState(b2)


VM=zeros(1,get(mfdu, 'nbdof'));

%F=[0 -200; 0 -300; 0 0]';
F = [0 -330]';
t = [0 0.9032 1 1.1 1.3 1.5 1.7 1.74 1.7 1.5 1.3 1.1 1 0.9032 0.7 0.5 0.3 0.1 0];

%nbstep = size(F,2);
nbstep = size(t,2);

dd=get(mf0, 'basic dof from cvid');

for step=1:nbstep,

  %step = 1;
  %b2.set('param','source_term', mfd, get(mfd, 'eval',{0;-400*sin(step*pi/2)}));
  set(b2, 'param','source_term', mfd, get(mfd, 'eval',{F(1,1);F(2,1)*t(step)}));
  get(b2, 'solve', mds, 'very noisy', 'max_iter', 1000, 'max_res', 1e-6);
  
  
  T = get(mds, 'tangent_matrix');
  
  
  U=get(mds, 'state');
  U=U(1:get(mfu, 'nbdof'));
  VM = get(b0, 'von mises', mds, mfdu);
  max(abs(VM))
  
  subplot(2,1,1);
  gf_plot(mfdu,VM,'deformed_mesh','on', 'deformation',U,'deformation_mf',mfu,'refine', 4, 'deformation_scale',1); 
  colorbar;
  caxis([0 10000]);
  n = t(step);
    title(['Von Mises criterion for t = ', num2str(n)]);
  
  ERR=gf_compute(mfu,U,'error estimate', mim);
  E=ERR; E(dd)=ERR;
  subplot(2,1,2);
  gf_plot(mf0, E, 'mesh','on', 'refine', 1); colorbar;
  title('Error estimate');

    
  pause;
end;
