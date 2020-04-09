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


disp('This demo is an adaption of the original tripod demo')
disp('which uses the "brick" framework of getfem')
disp('You can easily switch between linear/non linear')
disp('compressible/incompressible elasticity!')

linear = true;
incompressible = false;


gf_workspace('clear all');
% import the mesh
m=gfMesh('import','gid','../meshes/tripod.GiD.msh');
mfu=gfMeshFem(m,3);     % mesh-fem supporting a 3D-vector field
mfd=gfMeshFem(m,1);     % scalar mesh_fem, for data fields.
% the mesh_im stores the integration methods for each tetrahedron
mim=gfMeshIm(m,gf_integ('IM_TETRAHEDRON(5)'));
% we choose a P2 fem for the main unknown
gf_mesh_fem_set(mfu,'fem',gf_fem('FEM_PK(3,2)'));
% the material is homogeneous, hence we use a P0 fem for the data
gf_mesh_fem_set(mfd,'fem',gf_fem('FEM_PK(3,0)'));
% display some informations about the mesh
disp(sprintf('nbcvs=%d, nbpts=%d, nbdof=%d',gf_mesh_get(m,'nbcvs'),...
             gf_mesh_get(m,'nbpts'),gf_mesh_fem_get(mfu,'nbdof')));
P=gf_mesh_get(m,'pts'); % get list of mesh points coordinates
pidtop=find(abs(P(2,:)-13)<1e-6); % find those on top of the object
pidbot=find(abs(P(2,:)+10)<1e-6); % find those on the bottom
% build the list of faces from the list of points
ftop=gf_mesh_get(m,'faces from pid',pidtop); 
fbot=gf_mesh_get(m,'faces from pid',pidbot);
% assign boundary numbers
gf_mesh_set(m,'boundary',1,ftop);
gf_mesh_set(m,'boundary',2,fbot);

E = 1e3; Nu = 0.3;
% set the Lame coefficients
lambda = E*Nu/((1+Nu)*(1-2*Nu));
mu = E/(2*(1+Nu));

% create a meshfem for the pressure field (used if incompressible ~= 0)
mfp=gfMeshFem(m); set(mfp, 'fem',gfFem('FEM_PK_DISCONTINUOUS(3,0)'));

md=gf_model('real');
gf_model_set(md, 'add fem variable', 'u', mfu);
if (linear)
  % the linearized elasticity , for small displacements
  gf_model_set(md, 'add initialized data', 'cmu', [mu]);
  gf_model_set(md, 'add initialized data', 'clambda', [lambda]);
  gf_model_set(md, 'add isotropic linearized elasticity brick', mim, 'u', 'clambda', 'cmu');

  if (incompressible)
    gf_model_set(md, 'add fem variable', 'p', mfp);
    gf_model_set(md, 'add linear incompressibility brick', mim, 'u', 'p');
  end;
else
  params = [lambda;mu];
  gf_model_set(md,'add initialized data','params', params);
  if (incompressible)
    lawname = 'Incompressible Mooney Rivlin';
    gf_model_set(md, 'add finite strain elasticity brick', mim, lawname, 'u', 'params');
    gf_model_set(md, 'add fem variable', 'p', mfp);
    gf_model_set(md, 'add finite strain incompressibility brick',  mim, 'u', 'p');
  else
    lawname = 'SaintVenant Kirchhoff';
    gf_model_set(md, 'add finite strain elasticity brick', mim, lawname, 'u', 'params');
  end;
end

% set a vertical force on the top of the tripod

gf_model_set(md, 'add initialized data', 'VolumicData', [0;-10;0]);
gf_model_set(md, 'add source term brick', mim, 'u', 'VolumicData');

% attach the tripod to the ground
gf_model_set(md, 'add Dirichlet condition with multipliers', mim, 'u', mfu, 2);

disp('running solve...')

t0=cputime; 

gf_model_get(md, 'solve', 'noisy', 'max iter', 1);
U = gf_model_get(md, 'variable', 'u');

disp(sprintf('solve done in %.2f sec', cputime-t0));

mfdu=gf_mesh_fem(m,1);
% the P2 fem is not derivable across elements, hence we use a discontinuous
% fem for the derivative of U.
gf_mesh_fem_set(mfdu,'fem',gf_fem('FEM_PK_DISCONTINUOUS(3,1)'));
if (linear)
  VM = gf_model_get(md, 'compute isotropic linearized Von Mises or Tresca', 'u', 'clambda', 'cmu', mfdu);
else
  VM = gf_model_get(md, 'compute finite strain elasticity Von Mises', lawname, 'u', 'params', mfdu);
end

disp('plotting ... ');

% we plot the von mises on the deformed object, in superposition
% with the initial mesh.
if (linear),
  gf_plot(mfdu,VM,'mesh','on', 'cvlst', get(m, 'outer faces'),...
	  'deformation',U,'deformation_mf',mfu);
else
  gf_plot(mfdu,VM,'mesh','on', 'cvlst', get(m, 'outer faces'),...
	  'deformation',U,'deformation_mf',mfu,'deformation_scale',1);
end;

caxis([0 100]);
colorbar; view(180,-50); camlight;
gf_colormap('tripod');

% the von mises stress is exported into a VTK file
% (which can be viewed with 'mayavi -d tripod.vtk -m BandedSurfaceMap')
% see http://mayavi.sourceforge.net/
gf_mesh_fem_get(mfdu,'export to vtk','tripod.vtk','ascii',VM,'vm')

