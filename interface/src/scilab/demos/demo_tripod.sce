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

path = get_absolute_file_path('demo_tripod.sce');

printf('\nThis demo is an adaption of the original tripod demo\n')
printf('which uses the new ''brick'' framework of getfem.\n')
printf('The code is shorter, faster and much more powerful.\n')
printf('You can easily switch between linear/non linear\n')
printf('compressible/incompressible elasticity!\n\n')

linear = 1;
incompressible = 0;

if getos()=='Windows' then
  // Under Windows, all the trace messages are available in the dos console
  // Under Linuxs, all the trace messages are redirected to the Scilab console
  consolebox('on');
end
gf_util('trace level',3);
gf_util('warning level',3);

gf_workspace('clear all');

// import the mesh
m   = gf_mesh('import','gid', path + '/data/tripod.GiD.msh');
mfu = gf_mesh_fem(m,3);     // mesh-fem supporting a 3D-vector field
mfd = gf_mesh_fem(m,1);     // scalar mesh_fem, for data fields.

// the mesh_im stores the integration methods for each tetrahedron
mim = gf_mesh_im(m,gf_integ('IM_TETRAHEDRON(5)'));

// we choose a P2 fem for the main unknown
gf_mesh_fem_set(mfu,'fem',gf_fem('FEM_PK(3,2)'));

// the material is homogeneous, hence we use a P0 fem for the data
gf_mesh_fem_set(mfd,'fem',gf_fem('FEM_PK(3,0)'));

// display some informations about the mesh
printf('nbcvs=%d, nbpts=%d, nbdof=%d\n',gf_mesh_get(m,'nbcvs'), gf_mesh_get(m,'nbpts'),gf_mesh_fem_get(mfu,'nbdof'));

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
Nu = 0.3;

// set the Lame coefficients
lambda = E*Nu/((1+Nu)*(1-2*Nu));
mu     = E/(2*(1+Nu));

// create a meshfem for the pressure field (used if incompressible ~= 0)
mfp = gf_mesh_fem(m); 
gf_mesh_fem_set(mfp,'fem',gf_fem('FEM_PK_DISCONTINUOUS(3,0)'));




md=gf_model('real');
gf_model_set(md, 'add fem variable', 'u', mfu);
if (linear)
  // the linearized elasticity , for small displacements
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

// set a vertical force on the top of the tripod

gf_model_set(md, 'add initialized data', 'VolumicData', [0;-10;0]);
gf_model_set(md, 'add source term brick', mim, 'u', 'VolumicData');

// attach the tripod to the ground
gf_model_set(md, 'add Dirichlet condition with multipliers', mim, 'u', mfu, 2);

disp('running solve...')

gf_model_get(md, 'solve', 'noisy', 'max iter', 1);
U = gf_model_get(md, 'variable', 'u');

mfdu=gf_mesh_fem(m,1);
// the P2 fem is not derivable across elements, hence we use a discontinuous
// fem for the derivative of U.
gf_mesh_fem_set(mfdu,'fem',gf_fem('FEM_PK_DISCONTINUOUS(3,1)'));
if (linear)
  VM = gf_model_get(md, 'compute isotropic linearized Von Mises or Tresca', 'u', 'clambda', 'cmu', mfdu);
else
  VM = gf_model_get(md, 'compute finite strain elasticity Von Mises', lawname, 'u', 'params', mfdu);
end



disp('plotting ...');

h = scf();
h.color_map = jetcolormap(255); //gf_colormap('tripod');

// we plot the von mises on the deformed object, in superposition
// with the initial mesh.
drawlater;
if (linear) then
  gf_plot(mfdu,VM,'mesh','on', 'mesh_edges_color',name2rgb('white'), 'cvlst', gf_mesh_get(m, 'outer faces'), 'deformation',U,'deformation_mf',mfu);
else
  gf_plot(mfdu,VM,'mesh','on', 'mesh_edges_color',name2rgb('white'), 'cvlst', gf_mesh_get(m, 'outer faces'), 'deformation',U,'deformation_mf',mfu);
end

h.children.rotation_angles = [135 75];
colorbar(min(VM),max(VM)); 
xlabel('');
ylabel('');
zlabel('');
a = gca();
a.box = 'off';
a.axes_visible = 'off';
drawnow;

printf('the von mises stress is exported into a VTK file\n');
printf('(which can be viewed with ''mayavi -d tripod.vtk -m BandedSurfaceMap'')\n');
printf('see http://mayavi.sourceforge.net/\n');

gf_mesh_fem_get(mfdu,'export to vtk', path + '/tripod.vtk','ascii',VM,'vm')

printf('demo tripod terminated\n');
