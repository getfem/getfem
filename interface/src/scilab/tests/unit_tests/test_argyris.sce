// Copyright (C) 2009-2020 Yann Colette
// 
//  This file is a part of GetFEM++
// 
//  GetFEM++  is  free software;  you  can  redistribute  it  and/or modify it
//  under  the  terms  of the  GNU  Lesser General Public License as published
//  by  the  Free Software Foundation;  either version 3 of the License,  or
//  (at your option) any later version along with the GCC Runtime Library
//  Exception either version 3.1 or (at your option) any later version.
//  This program  is  distributed  in  the  hope  that it will be useful,  but
//  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
//  or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
//  License and GCC Runtime Library Exception for more details.
//  You  should  have received a copy of the GNU Lesser General Public License
//  along  with  this program;  if not, write to the Free Software Foundation,
//  Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.

clear pde;
gf_workspace('clear all');

NX = 10
m  = gf_mesh('triangles grid',[0:1/NX:1],[0:1/NX:1]);
//gf_mesh_set(m,'transform', [.3 .8; .8 -.2]);
//m = gf_mesh('pt2d', [0 0; 1 0; 0 1]', [1 2 3]');

// create a mesh_fem of for a field of dimension 1 (i.e. a scalar field)
mf   = gf_mesh_fem(m,1);
mfl  = gf_mesh_fem(m,1);
mflg = gf_mesh_fem(m,1);
mflh = gf_mesh_fem(m,1);

// assign the Q2 fem to all convexes of the mesh_fem,

//gf_mesh_fem_set(mf,'fem',gf_fem('FEM_PK(2,3)'));
//gf_mesh_fem_set(mf,'fem',gf_fem('FEM_ARGYRIS'));
gf_mesh_fem_set(mf,'fem',gf_fem('FEM_HCT_TRIANGLE'));
//gf_mesh_fem_set(mf,'fem',gf_fem('FEM_HERMITE(2)'));

gf_mesh_fem_set(mfl,'fem',gf_fem('FEM_PK(2,5)'));
gf_mesh_fem_set(mflg,'fem',gf_fem('FEM_PK_DISCONTINUOUS(2,4)'));
gf_mesh_fem_set(mflh,'fem',gf_fem('FEM_PK_DISCONTINUOUS(2,3)'));

// an exact integration will be used
//mim = gf_mesh_im(m, gf_integ('IM_TRIANGLE(13)'));
mim = gf_mesh_im(m, gf_integ('IM_HCT_COMPOSITE(IM_TRIANGLE(13))'));

// detect the border of the mesh
border = gf_mesh_get(m,'outer faces');

// mark it as boundary #1
gf_mesh_set(m, 'boundary', 1, border);

h = scf();
h.color_map = jetcolormap(255);
drawlater;
gf_plot_mesh(m, 'regions', [1]); // the boundary edges appears in red
h.color_map = jetcolormap(255);
drawnow;

pause;

// exact solution

if 0 then
  // setup a pde structure for gf_solve
  pde = init_pde();

  pde('type')   = 'laplacian'; 
  pde('lambda') = list(1);
  pde('mim')    = mim;
  pde('mf_u')   = mf;       // this does not copy whole objects, just their handles
  pde('mf_d')   = mfl;
  expr_u        = 'y.*(y-1).*x.*(x-1)+x.^5/10';
  expr_f        = '-(2*(x.^2+y.^2)-2*x-2*y+20*x.^3/10)';
  pde('F')      = list(expr_f);
  
  pde = add_empty_bound(pde);
  pde('bound')($)('type') = 'Dirichlet';
  pde('bound')($)('R')    = list(expr_u);

  U      = gf_solve(pde);
  Uexact = gf_mesh_fem_get_eval(mfl, list(list(expr_u)));
else
  expr_u = 'y.^5';
  Uexact = gf_mesh_fem_get_eval(mfl, list(list(expr_u)));
  M = gf_asm('mass matrix', mim, mf, mf);
  F = gf_asm('volumic source', mim, mf, mfl, Uexact);
  U = (M\F)';
end
						
Ul    = gf_compute(mf,U,'interpolate on', mfl);
DUl   = gf_compute(mfl, Ul, 'gradient',mflg);
D2Ul  = gf_compute(mflg, DUl, 'gradient',mflh);
D2Ul2 = gf_compute(mfl,Ul, 'hessian',mflh);
nref  = 4

h = scf();
h.color_map = jetcolormap(255);
drawlater;
subplot(2,2,1); 
gf_plot(mfl,Ul,'mesh','on','refine',nref,'contour',.01:.01:.1); 
colorbar(min(Ul),max(Ul));
title('computed solution');

subplot(2,2,2); 
gf_plot(mfl,Ul-Uexact,'mesh','on','refine',nref); 
colorbar(min(Ul-Uexact),max(Ul-Uexact));
title('difference with exact solution');

subplot(2,2,3); 
gf_plot(mflg,DUl(1,:),'mesh','on', 'refine', nref); 
colorbar(min(DUl(1,:)),max(DUl(1,:)));
title('gradx');

subplot(2,2,4); 
gf_plot(mflg,DUl(2,:),'mesh','on', 'refine', nref); 
colorbar(min(DUl(2,:)),max(DUl(2,:)));
title('grady');
drawnow;

disp(sprintf('H1 norm of error: %g', gf_compute(mfl,Ul-Uexact,'H1 norm',mim)));

