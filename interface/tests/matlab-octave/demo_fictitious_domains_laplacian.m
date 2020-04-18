% Copyright (C) 2006-2020 Yves Renard, Julien Pommier.
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


disp('This demo use levelset to impose (weakly) a Dirichlet condition on a part of an ');
disp('implicit boundary defined by the zero of the levelset and a Neumann condition on ');
disp('the remaining part of that boundary. A Poisson problem');

clear all;
gf_workspace('clear all');

  

NX= 10
N = 3
ls_degree = 1
R = 0.4;

if (N == 3)
  m=gf_mesh('cartesian', -.5:(1/NX):.5, -.5:(1/NX):.5, -.5:(1/NX):.5);
  %m=gf_Mesh('regular simplices', -.5:(1/NX):.5, -.5:(1/NX):.5, -.5:(1/NX):.5);
  mfu0=gfMeshFem(m,1);
  mf_mult=gfMeshFem(m,1);
  set(mfu0, 'fem', gf_fem('FEM_QK(3,1)'))
  %set(mfu0, 'fem', gf_fem('FEM_PK(3,1)'));
  set(mf_mult, 'fem', gf_fem('FEM_QK(3,1)'))
  %set(mf_mult, 'fem', gf_fem('FEM_PK(3,0)'));
  adapt_im = 'IM_TETRAHEDRON(6)';
elseif (N == 2)
  %m=gf_mesh('cartesian', -.5:(1/NX):.5, -.5:(1/NX):.5);
  m=gfMesh('triangles grid', -.5:(1/NX):.5, -.5:(1/NX):.5);
  mfu0=gfMeshFem(m,1);
  mf_mult=gfMeshFem(m,1);
  set(mfu0, 'fem', gf_fem('FEM_PK(2,3)'));
  set(mf_mult, 'fem', gf_fem('FEM_PK(2,1)'));
  adapt_im = 'IM_TRIANGLE(6)'
else 
  error('Wrong dimension');
end

ls=gf_levelset(m, ls_degree);
ls2s=gf_levelset(m, ls_degree, 'with_secondary');
ls2=gf_LevelSet(m, ls_degree, 'with_secondary');

mf_ls=gfObject(gf_levelset_get(ls, 'mf'));
mf_ls2=gfObject(gf_levelset_get(ls2, 'mf'));
P=get(mf_ls, 'basic dof nodes');
x = P(1,:); y = P(2,:);
if (N == 3)
  z = P(3,:);
else
  z = 0 * x;
end
ULS=1000*ones(1,numel(x));
ULS2=1000*ones(1,numel(x));
ULS = min(ULS,((x.^2 + y.^2 + z.^2).^1.5 - R^3));
%ULS2 = min(ULS,x);
if (N == 3)
  ULS2 = min(ULS2,z);%ULS2=-ULS2;
else
 ULS2 = min(ULS2,y);%ULS2=-ULS2;
end


gf_levelset_set(ls, 'values', ULS);
gf_levelset_set(ls2, 'values', ULS, ULS2);
%gf_levelset_set(lss, 'values', ULS);
gf_levelset_set(ls2s, 'values', ULS, -ULS2);

mls=gfMeshLevelSet(m);
set(mls, 'add', ls);
set(mls, 'add',ls2);
set(mls, 'adapt');

mls2=gfMeshLevelSet(m);
set(mls2, 'add', ls);
set(mls2, 'add',ls2s);
set(mls2, 'adapt');


mim_bounddown = gfMeshIm('levelset',mls,'boundary(b)', gf_integ(adapt_im));
%mim_bound2 = gfMeshIm('levelset',mls,'boundary(a)', gf_integ(adapt_im));
mim_boundup = gfMeshIm('levelset',mls2,'boundary(b)', gf_integ(adapt_im));
mim_int = gfMeshIm('levelset', mls, 'inside(a)', gf_integ(adapt_im));
set(mim_int, 'integ', 4);

% Some verifications
A1=gf_asm('volumic','V()+=comp()',mim_bounddown);
%A2=gf_asm('volumic','V()+=comp()',mim_bound2);
A2=gf_asm('volumic','V()+=comp()',mim_boundup);
V =gf_asm('volumic','V()+=comp()',mim_int);
if (N == 2)
  disp(sprintf('length : %g should be %g', A1, pi*R));
  disp(sprintf('length : %g should be %g', A2, 2*pi*R));
  disp(sprintf('area : %g should be %g', V, pi*R^2));
else
  disp(sprintf('area : %g should be %g', A1, 4*pi*R^2/2));
  disp(sprintf('area : %g should be %g', A2, 4*pi*R^2));
  disp(sprintf('volume : %g should be %g', V, 4*pi*R^3/3.));   
end

% partial mesh fem
dof_out = get(mfu0, 'dof from im', mim_int);
cv_out = get(mim_int, 'convex_index');
cv_in = setdiff(gf_mesh_get(m, 'cvid'), cv_out);
mfu = gfMeshFem('partial', mfu0, dof_out, cv_in);

% data
if (N == 2)
  Volumic_data = gf_mesh_fem_get(mfu0, 'eval', { '45*sqrt(x.^2+y.^2)' });
  surface_data = gf_mesh_fem_get(mfu0, 'eval', { '-15*(x.^2+y.^2)' });
  Sol_U = gf_mesh_fem_get(mfu0, 'eval', { sprintf('5*((%g)^3-(x.^2+y.^2).^1.5)', R) });
else
  Volumic_data = gf_mesh_fem_get(mfu0, 'eval', { '60*sqrt(x.^2+y.^2+z.^2)' });
  surface_data = gf_mesh_fem_get(mfu0, 'eval', { '-15*(x.^2+y.^2+z.^2)' });
  Sol_U = gf_mesh_fem_get(mfu0, 'eval', { sprintf('5*((%g)^3-(x.^2+y.^2+z.^2).^1.5)', R) });
end

% getfem model


md=gf_model('real');
gf_model_set(md, 'add fem variable', 'u', mfu);
gf_model_set(md, 'add Laplacian brick', mim_int, 'u');
gf_model_set(md, 'add initialized fem data', 'VolumicData', mfu0, Volumic_data);
gf_model_set(md, 'add source term brick', mim_int, 'u', 'VolumicData');

if 0,
B2=gf_asm('mass matrix', mim_boundup, mfu0, mfu);
Rh=gf_spmat_get(B2, 'tmult', surface_data);
gf_model_set(md, 'add explicit rhs', 'u', Rh);
else

gf_model_set(md, 'add initialized fem data', 'SurfaceData', mfu0, surface_data);
gf_model_set(md, 'add source term brick', mim_boundup, 'u', 'SurfaceData');   
    
% gf_model_set(md, 'add initialized fem data', 'SurfaceData', mfu0, surface_data);
% gf_model_set(md, 'add source term brick', mim_bound2, 'u', 'SurfaceData');
% gf_model_set(md, 'add initialized fem data', 'SurfaceData2', mfu0, -surface_data);
% gf_model_set(md, 'add source term brick', mim_bounddown, 'u', 'SurfaceData2');
end;

%range bases
BRBB=gf_asm('mass matrix', mim_bounddown, mf_mult, mf_mult);
gf_mesh_fem_set(mf_mult,'reduce meshfem', BRBB);

gf_model_set(md, 'add fem variable', 'Lambda', mf_mult);
B=gf_asm('mass matrix', mim_bounddown, mfu, mf_mult);
gf_model_set(md, 'add explicit matrix', 'u', 'Lambda', B, 1);

% gf_model_set(md, 'add multiplier', 'Lambda', mf_mult, 'u');
% gf_model_set(md, 'add Dirichlet condition with multipliers', ...
% 	     mim_bounddown, 'u', 'Lambda', -1);

gf_model_get(md, 'solve');
U = gf_model_get(md, 'variable', 'u');
Lambda = gf_model_get(md, 'variable', 'Lambda');

% Comparison with the exaxt solution
U0=gf_compute(mfu,U,'interpolate on',mfu0);
ERRL2 = gf_compute(mfu0, U0, 'L2 dist', mim_int, mfu0, Sol_U);
L2_Error=100*gf_compute(mfu, U, 'L2 dist', mim_int, mfu0, Sol_U)/gf_compute(mfu0, Sol_U,'L2 norm',mim_int);
ERRH1 = gf_compute(mfu0, U0, 'H1 semi dist', mim_int, mfu0, Sol_U);
H1_Error=100*gf_compute(mfu, U, 'H1 semi dist', mim_int, mfu0, Sol_U)/gf_compute(mfu0,Sol_U,'H1 norm',mim_int);

%Compute error on the Multiplier


map_Error=(abs(U0-Sol_U));
%LambdaE = gf_mesh_fem_get(mf_mult, 'eval', { '-15*(x.^2+y.^2+z.^2)' });
BA=gf_asm('mass matrix', mim_bounddown, mf_mult, mf_mult);
BLS=gf_asm('lsneuman matrix', mim_bounddown, mfu0, mf_mult, ls);
KA=gf_asm('nlsgrad matrix', mim_bounddown, mfu0, mfu0, ls);
%ERRL_mult= (Lambda*gf_spmat_get(BA, 'mult',Lambda)+ U0*gf_spmat_get(KA, 'mult',U0) + 2*Lambda*gf_spmat_get(BLS, 'mult',U0))/(U0*gf_spmat_get(KA, 'mult',U0))
ERRL2_mult= (Lambda*gf_spmat_get(BA, 'mult',Lambda)+ Sol_U*gf_spmat_get(KA, 'mult',Sol_U) + 2*Lambda*gf_spmat_get(BLS, 'mult',Sol_U))/(Sol_U*gf_spmat_get(KA, 'mult',Sol_U));
L2_error_mult=100* sqrt(abs(ERRL2_mult));

disp(sprintf('L2 norm %g\n H1 error = %g\n L2 norm mult=%g', L2_Error, H1_Error, L2_error_mult));

if (N == 2)
  figure(1);
  gf_plot(mfu, U, 'mesh','on', 'refine', 2);
  xlabel('x'); ylabel('y');
  title('Displacement solution');
  figure(2);
  gf_plot(mfu0, map_Error, 'mesh','on', 'refine', 2);
  xlabel('x'); ylabel('y');
  title('Map Error in displacement');
else
 % gf_plot(mfu, U, 'mesh','on', 'cvlst', gf_mesh_get(m, 'outer faces'), 'refine', 2);
 %Plot displacement
  figure(1);
  sl=gf_slice({'boundary',{'intersection',{'ball', -1,[0;0;0], R},{'planar',1,[0;0;0],[0;0;1]}}},mfu,3);
  %sl=gf_slice({'boundary',{'intersection',{'ball',-1,[0;0;0],R}}},mfu,1);
  Usl=gf_compute(mfu,U,'interpolate on', sl);
  gf_plot_slice(sl,'mesh_faces','on','mesh','off','data',Usl,'mesh_slice_edges','on');
  xlabel('x'); ylabel('y'); zlabel('z');
  title('Displacement solution');
  
  
  %plot exact solution
  figure(2);
  sl=gf_slice({'boundary',{'intersection',{'ball', -1,[0;0;0], R},{'planar',1,[0;0;0],[0;0;1]}}},mfu0,3);
  %sl=gf_slice({'boundary',{'intersection',{'ball',-1,[0;0;0],R}}},mfu0,10);
  Usl=gf_compute(mfu0,Sol_U,'interpolate on', sl);
  gf_plot_slice(sl,'mesh_faces','on','mesh','off','data',Usl,'mesh_slice_edges','on');
  xlabel('x'); ylabel('y'); zlabel('z');
  title('Exact solution');
  
  %Plot map eroor on u
  figure(3);
  sl=gf_slice({'boundary',{'intersection',{'ball', -1,[0;0;0], R},{'planar',-1,[0;0;0],[0;0;1]}}},mfu0,10);
  %sl=gf_slice({'boundary',{'intersection',{'ball',-1,[0;0;0],R}}},mfu0,10);
  Usl=gf_compute(mfu0,map_Error,'interpolate on', sl);
  gf_plot_slice(sl,'mesh_faces','on','mesh','off','data',Usl,'mesh_slice_edges','on');
  xlabel('x'); ylabel('y'); zlabel('z');
  title('Map Error in displacement');
  
  
  %Plot multiplier
  figure(4);
  %sl=gf_slice({'boundary',{'intersection',{'ball', 0,[0;0;0], R},{'planar',-1,[0;0;0],[0;0;1]}}},mf_mult,10);
  sl=gf_slice({'boundary',{'intersection',{'ball',0,[0;0;0],R}}},mf_mult,3);
  Usl=gf_compute(mf_mult, Lambda,'interpolate on', sl);
  gf_plot_slice(sl,'mesh_faces','on','mesh','off','data',Usl,'mesh_slice_edges','on');
  xlabel('x'); ylabel('y'); zlabel('z');
  title('Multiplier solution');
  end
gf_colormap('chouette');

