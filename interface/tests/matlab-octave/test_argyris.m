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


if (1) % draw the shape functions

gf_workspace('clear all');
m = gf_mesh('triangles grid',[0:0.5:1],[0:0.5:1]);
mf  = gfMeshFem(m,1);
gf_mesh_fem_set(mf,'fem',gf_fem('FEM_REDUCED_HCT_TRIANGLE'));
nbdof=gf_mesh_fem_get(mf, 'nbdof');
U = zeros(1, nbdof);

for i=13:nbdof
  U(i) = 0.5;
gf_plot(mf,U,'mesh','on','refine',10, 'zplot', 'on'); colorbar;
  title('shape function');
  pause;
gf_plot(mf,U,'mesh','on','refine',10); colorbar;
  title('shape function');
  
  U(i) = 0.0;
end;

return;

end;



% trace on;
gf_workspace('clear all');
NX=10
m = gf_mesh('triangles grid',[0:1/NX:1],[0:1/NX:1]);
%gf_mesh_set(m,'transform', [.3 .8; .8 -.2]);
%m=gfMesh('pt2d', [0 0; 1 0; 0 1]', [1 2 3]');
% create a mesh_fem of for a field of dimension 1 (i.e. a scalar field)
mf  = gfMeshFem(m,1);
mfl = gfMeshFem(m,1);
mflg = gfMeshFem(m,1);
mflh = gfMeshFem(m,1);
% assign the Q2 fem to all convexes of the mesh_fem,
%gf_mesh_fem_set(mf,'fem',gf_fem('FEM_PK(2,3)'));
%gf_mesh_fem_set(mf,'fem',gf_fem('FEM_ARGYRIS'));
gf_mesh_fem_set(mf,'fem',gf_fem('FEM_HCT_TRIANGLE'));
%gf_mesh_fem_set(mf,'fem',gf_fem('FEM_HERMITE(2)'));
gf_mesh_fem_set(mfl,'fem',gf_fem('FEM_PK(2,5)'));
gf_mesh_fem_set(mflg,'fem',gf_fem('FEM_PK_DISCONTINUOUS(2,4)'));
gf_mesh_fem_set(mflh,'fem',gf_fem('FEM_PK_DISCONTINUOUS(2,3)'));
% an exact integration will be used
%mim = gf_mesh_im(m, gf_integ('IM_TRIANGLE(13)'));
mim = gf_mesh_im(m, gf_integ('IM_HCT_COMPOSITE(IM_TRIANGLE(13))'));
% detect the border of the mesh
border = gf_mesh_get(m,'outer faces');
% mark it as boundary #1
gf_mesh_set(m, 'boundary', 1, border);
gf_plot_mesh(m, 'regions', [1]); % the boundary edges appears in red
pause(1);

                                                  % exact solution

expr_u = 'y.^5';
Uexact = gf_mesh_fem_get(mfl,'eval', { expr_u }); % interpolate the
M=gf_asm('mass matrix', mim, mf, mf);
F=gf_asm('volumic source', mim, mf, mfl, Uexact);
U=(M\F)';
  

						
Ul = gf_compute(mf,U,'interpolate on', mfl);

DUl = gf_compute(mfl, Ul, 'gradient',mflg);
D2Ul = gf_compute(mflg, DUl, 'gradient',mflh);
D2Ul2 = gf_compute(mfl,Ul, 'hessian',mflh);
nref=4

subplot(2,2,1); gf_plot(mfl,Ul,'mesh','on','refine',nref,'contour',.01:.01:.1); colorbar;title('computed solution');
subplot(2,2,2); gf_plot(mfl,Ul-Uexact,'mesh','on','refine',nref); colorbar;title('difference with exact solution');

subplot(2,2,3); gf_plot(mflg,DUl(1,:),'mesh','on', 'refine', nref); colorbar;title('gradx');
subplot(2,2,4); gf_plot(mflg,DUl(2,:),'mesh','on', 'refine', nref); colorbar;title('grady');

disp(sprintf('H1 norm of error: %g', gf_compute(mfl,Ul-Uexact,'H1 norm',mim)));
