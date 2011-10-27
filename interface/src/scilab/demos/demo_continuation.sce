// Scilab GetFEM++ interface
// Copyright (C) 2011-2011 Tomas Ligursky, Yves Renard.
//
// This file is a part of GetFEM++
//
// GetFEM++  is  free software;  you  can  redistribute  it  and/or modify it
// under  the  terms  of the  GNU  Lesser General Public License as published
// by  the  Free Software Foundation;  either version 2.1 of the License,  or
// (at your option) any later version.
// This program  is  distributed  in  the  hope  that it will be useful,  but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
// License for more details.
// You  should  have received a copy of the GNU Lesser General Public License
// along  with  this program;  if not, write to the Free Software Foundation,
// Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.
//
// Simple example of a bifurcation problem: Delta(u) + u = lambda exp(u)
//
// This program is used to check that scilab-getfem is working. This is also
// a good example of use of GetFEM++.
//

lines(0);
stacksize('max');
// path = get_absolute_file_path('demo_continuation.sce');

gf_workspace('clear all');
lambda0 = 0;
direction = 1;
nbstep = 10;

noisy = 2;

// create a simple cartesian mesh
m = gf_mesh('cartesian', [0:.1:1]);

// create a mesh_fem for a field of dimension 1 (i.e. a scalar field)
mf = gf_mesh_fem(m, 1);
// assign the Q1 fem to all convexes of the mesh_fem,
gf_mesh_fem_set(mf, 'classical fem', 1);

// Integration which will be used
mim = gf_mesh_im(m, 4);


md = gf_model('real');
gf_model_set(md, 'add fem variable', 'u', mf);
gf_model_set(md, 'add Laplacian brick', mim, 'u');
lambda = lambda0;
gf_model_set(md, 'add initialized data', 'lambda', [lambda]);
gf_model_set(md, 'add basic nonlinear brick', mim, 'u', 'u-lambda*exp(u)', '1-lambda*exp(u)', 'lambda');


if noisy > 0 then printf('computing initial point\n'); end
gf_model_get(md, 'solve', 'noisy', 'max_iter', 100);
U = gf_model_get(md, 'variable', 'u');
//lambda = gf_model_get(md, 'variable', 'lambda');
[T_U, T_lambda, h] = gf_model_get(md, 'init Moore-Penrose continuation', 'lambda', U, lambda, direction);
printf('U = '); disp(U); printf('lambda = %e\n', lambda);
printf('lambda - U(1) * exp(-U(1)) = %e\n', lambda - U(1) * exp(-U(1)));

clf();
drawlater;
x  = gf_mesh_fem_get(mf, 'basic dof nodes'); plot(x, U, 'k+-');
//gf_plot(mf, U, 'mesh', 'on', 'contour', .01:.01:.1);  colorbar;
xtitle('lambda = '+ string(lambda), '$x$', '$u$');
a = get("current_axes");
a.title.font_size = 4; a.x_label.font_size= 4; a.y_label.font_size= 4; a.y_label.font_angle = 0;
drawnow;


for step = 1:nbstep
    sleep(1000);
    printf('\nbeginning of step %d\n', step);
    [T_U, T_lambda, h] = gf_model_get(md, 'Moore-Penrose continuation', 'lambda', U, lambda, T_U, T_lambda, h);
    U = gf_model_get(md, 'variable', 'u');
    lambda = gf_model_get(md, 'variable', 'lambda');
    printf('U = '); disp(U); printf('lambda = %e\n', lambda);
    printf('lambda - U(1) * exp(-U(1)) = %e\n', lambda - U(1) * exp(-U(1)));

    clf();
    drawlater;
    x  = gf_mesh_fem_get(mf, 'basic dof nodes'); plot(x, U, 'k+-');
    //gf_plot(mf, U, 'mesh', 'on', 'contour', .01:.01:.1);  colorbar;
    xtitle('lambda = '+ string(lambda), '$x$', '$u$');
    a = get("current_axes");
    a.title.font_size = 4;
    a.x_label.font_size= 4; a.y_label.font_size= 4; a.y_label.font_angle = 0;
    drawnow;
    printf('end of step n° %d', step); printf(' / %d\n', nbstep);
end