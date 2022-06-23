% Copyright (C) 2005-2020 Yves Renard, Julien Pommier.
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
% Example of automatic refinement of the mesh
% In this example, the refinement will focus on the
% transition between the Dirichlet and the Neumann boundary.

gf_workspace('clear all');
% clear all; clf;
L=100; H=22;
N=2;
draw = true;

asize =  size(who('automatic_var654'));
if (asize(1)) draw = false; end;


if (N == 2), % 2D beam
  m=gfMesh('regular simplices',0:10:L, 0:11:H);
  mim=gfMeshIm(m);    set(mim, 'integ',gfInteg('IM_TRIANGLE(6)'));
  mfu=gfMeshFem(m,N); set(mfu, 'fem',gfFem('FEM_PK(2,2)'));
  mfd=gfMeshFem(m);   set(mfd, 'fem',gfFem('FEM_PK(2,1)'));
  mf0=gfMeshFem(m);   set(mf0, 'fem',gfFem('FEM_PK(2,0)'));
  mfdu=gfMeshFem(m);  set(mfdu,'fem',gfFem('FEM_PK_DISCONTINUOUS(2,2)'));
else         % 3D beam
  m=gfMesh('regular simplices',0:10:L, 0:11:H, 0:11:H);
  mim=gfMeshIm(m);    set(mim, 'integ',gfInteg('IM_TETRAHEDRON(5)'));
  mfu=gfMeshFem(m,N); set(mfu, 'fem',gfFem('FEM_PK(3,2)'));
  mfd=gfMeshFem(m);   set(mfd, 'fem',gfFem('FEM_PK(3,1)'));
  mf0=gfMeshFem(m);   set(mf0, 'fem',gfFem('FEM_PK(3,0)'));
  mfdu=gfMeshFem(m);  set(mfdu,'fem',gfFem('FEM_PK_DISCONTINUOUS(3,1)'));
end;

lambda=121150; mu=80769;

P=get(m,'pts');
fleft =gf_mesh_get(m,'faces from pid',find(abs(P(1,:))<1e-6));
fright=gf_mesh_get(m,'faces from pid',find(abs(P(1,:) - L)<1e-6));
% assign boundary numbers
gf_mesh_set(m,'region',1,fleft);
gf_mesh_set(m,'region',2,fright);


F=zeros(N,1); F(2) = -20; % the external force


md=gf_model('real');
gf_model_set(md, 'add fem variable', 'u', mfu);
gf_model_set(md, 'add initialized data', 'lambda', [lambda]);
gf_model_set(md, 'add initialized data', 'mu', [mu]);
gf_model_set(md, 'add isotropic linearized elasticity brick', mim, 'u', 'lambda', 'mu');
gf_model_set(md, 'add initialized data', 'VolumicData', F);
gf_model_set(md, 'add source term brick', mim, 'u', 'VolumicData');
gf_model_set(md, 'add Dirichlet condition with multipliers', mim, 'u', mfu, 1);



for step=1:7,
  gf_model_get(md, 'solve');
  U = gf_model_get(md, 'variable', 'u');
  
  VM = gf_model_get(md, 'compute isotropic linearized Von Mises or Tresca', 'u', 'lambda', 'mu', mfdu);
  
  if (N==3) opt = {'cvlst', get(m,'outer_faces')}; 
  else opt = {}; end;
  
  if (draw)
    subplot(2,1,1);
    gf_plot(mfdu,VM,'deformed_mesh','on', 'deformation',U,...
	    'deformation_mf',mfu,'refine', 4, 'deformation_scale',1, opt{:}); 
    gf_colormap('chouette');
    caxis([0 1e7]); colorbar; 
    title('Von Mises stress');
  end
  
  dd=get(mf0, 'basic dof from cvid');
  ERR=gf_compute(mfu, U, 'error estimate', mim);
  E=ERR; E(dd)=ERR;
  
  if (draw)
    subplot(2,1,2);
    gf_plot(mf0, E, 'mesh','on', 'refine', 1, opt{:}); colorbar;
    title('Error estimate')
    pause(1.5);
  end
  set(m, 'refine', find(ERR > 2e-6));
  set(m, 'optimize structure', false);
  norm(E)
end;

if (norm(E) > 0.005)
   error('Refine test: final error to big');
end
