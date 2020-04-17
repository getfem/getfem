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

printf('demo nonlinear_elasticity started\n');

path = get_absolute_file_path('demo_nonlinear_elasticity.sce');

// Load the axrot_matrix macro
exec(path + 'axrot_matrix.sci');

if getos()=='Windows' then
  // Under Windows, all the trace messages are available in the dos console
  // Under Linuxs, all the trace messages are redirected to the Scilab console
  consolebox('on');
end
gf_util('trace level',3);
gf_util('warning level',3);

gf_workspace('clear all');

// set a custom colormap
r=[0.7 .7 .7]; l = r($,:); s=63; s1=20; s2=25; s3=48;s4=55; 
for i=1:s
  c1 = max(min((i-s1)/(s2-s1),1),0);
  c2 = max(min((i-s3)/(s4-s3),1),0); 
  r($+1,:)=(1-c2)*((1-c1)*l + c1*[1 0 0]) + c2*[1 .8 .2]; 
end

incompressible = 1;

lawname = 'Ciarlet Geymonat';
params  = [1;1;0.5];
params = [0;1];
if (incompressible) then
  lawname = 'Mooney Rivlin';
  params  = [1;1];
end

if 0 then
  h = 20;
  // import the mesh
  //m = gf_mesh('load', path + '/data/ladder.mesh');
  //m = gf_mesh('load', path + '/data/ladder_1500.mesh');
  m = gf_mesh('load', path + '/data/holed_bar.mesh');
  gf_mesh_set(m, 'transform', [1 0 0; 0 0 1; 0 1 0]);
  mfu = gf_mesh_fem(m,3);     // mesh-fem supporting a 3D-vector field
  mfd = gf_mesh_fem(m,1);     // scalar mesh_fem
  // the mesh_im stores the integration methods for each tetrahedron
  mim = gf_mesh_im(m,gf_integ('IM_TETRAHEDRON(5)'));
  // we choose a P2 fem for the main unknown
  gf_mesh_fem_set(mfu, 'fem',gf_fem('FEM_HERMITE(3)'));
  //gf_mesh_fem_set(mfu, 'fem',gf_fem('FEM_PK(3,2)'));
  mfdu = gf_mesh_fem(m,1);
  // the material is homogeneous, hence we use a P0 fem for the data
  gf_mesh_fem_set(mfd,'fem',gf_fem('FEM_PK(3,1)'));
  // the P2 fem is not derivable across elements, hence we use a discontinuous
  // fem for the derivative of U.
  gf_mesh_fem_set(mfdu,'fem',gf_fem('FEM_PK_DISCONTINUOUS(3,2)'));
else
  N1 = 1; 
  N2 = 4; 
  h  = 20;
  m   = gf_mesh('cartesian',(0:N1)/N1 - .5, (0:N2)/N2*h, ((0:N1)/N1 - .5)*3);
  mfu = gf_mesh_fem(m,3);     // mesh-fem supporting a 3D-vector field
  mfd = gf_mesh_fem(m,1);     // scalar mesh_fem
  // the mesh_im stores the integration methods for each tetrahedron
  mim = gf_mesh_im(m,gf_integ('IM_GAUSS_PARALLELEPIPED(3,6)'));
  // we choose a P2 fem for the main unknown
  gf_mesh_fem_set(mfu, 'fem',gf_fem('FEM_QK(3,2)'));
  mfdu = gf_mesh_fem(m,1);
  // the material is homogeneous, hence we use a P0 fem for the data
  gf_mesh_fem_set(mfd,'fem',gf_fem('FEM_QK(3,1)'));
  // the P2 fem is not derivable across elements, hence we use a discontinuous
  // fem for the derivative of U.
  gf_mesh_fem_set(mfdu,'fem',gf_fem('FEM_QK_DISCONTINUOUS(3,2)'));
end

m_char    = gf_mesh_get(m, 'char');
mfu_char  = gf_mesh_fem_get(mfu, 'char');
mfdu_char = gf_mesh_fem_get(mfdu, 'char');

// display some informations about the mesh
disp(sprintf('nbcvs=%d, nbpts=%d, nbdof=%d',gf_mesh_get(m,'nbcvs'), gf_mesh_get(m,'nbpts'),gf_mesh_fem_get(mfu,'nbdof')));
P = gf_mesh_get(m,'pts'); // get list of mesh points coordinates
//pidtop = find(abs(P(2,:)-13)<1e-6); // find those on top of the object
//pidbot = find(abs(P(2,:)+10)<1e-6); // find those on the bottom

pidtop = find(abs(P(2,:)-h)<1e-6); // find those on top of the object
pidbot = find(abs(P(2,:)-0)<1e-6); // find those on the bottom

// build the list of faces from the list of points
ftop = gf_mesh_get(m,'faces from pid',pidtop); 
fbot = gf_mesh_get(m,'faces from pid',pidbot);

// assign boundary numbers
gf_mesh_set(m,'boundary',1,ftop);
gf_mesh_set(m,'boundary',2,fbot);
gf_mesh_set(m,'boundary',3,[ftop fbot]);


md = gf_model('real');
gf_model_set(md, 'add fem variable', 'u', mfu);
gf_model_set(md, 'add initialized data', 'params', params);
gf_model_set(md, 'add nonlinear elasticity brick', mim, 'u', lawname, 'params');
if (incompressible) then
  mfp = gf_mesh_fem(m,1); 
  gf_mesh_fem_set(mfp, 'classical discontinuous fem', 1);
  gf_model_set(md, 'add fem variable', 'p', mfp);
  gf_model_set(md, 'add nonlinear incompressibility brick',  mim, 'u', 'p');
end
 
gf_model_set(md, 'add fem data', 'DirichletData', mfd, 3);
gf_model_set(md, 'add Dirichlet condition with penalization', mim, 'u', 1e10, 3, 'DirichletData');

VM = zeros(1,gf_mesh_fem_get(mfdu,'nbdof'));

reload = 0;

if (reload == 0) then
  UU     = [];
  VVM    = [];
  nbstep = 40;
else
  load(path + '/demo_nonlinear_elasticity_U.mat');
  nb_step = size(UU,1);
end
P = gf_mesh_fem_get(mfd, 'basic dof_nodes');
r = sqrt(P(1 ,:).^2 + P(3, :).^2);
theta = atan(P(3,:),P(1,:));

scf();

for step=1:nbstep
  w = 3*step/nbstep;
  //set(b2, 'param', 'R', [0;0;0]);

  if (~reload) then
    R = zeros(3, gf_mesh_fem_get(mfd, 'nbdof'));
    dtheta  = %pi;
    dtheta2 = %pi/2;
    
    i_top = gf_mesh_fem_get(mfd, 'basic dof on region', 1);
    i_bot = gf_mesh_fem_get(mfd, 'basic dof on region', 2);
    dd = max(P(1,i_top)*sin(w*dtheta));
    if (w < 1) then
      RT1 = axrot_matrix([0 h*.75 0], [0 h*.75 1], w*dtheta);
      RT2 = axrot_matrix([0 0 0], [0 1 0], sqrt(w)*dtheta2);
      RB1 = axrot_matrix([0 h*.25 0], [0 h*.25 1], -w*dtheta);
      RB2 = RT2';
    elseif (w < 2) then
      RT1 = axrot_matrix([0 h*.75 0], [0 h*.75 1], (2-w)*dtheta);
      RT2 = axrot_matrix([0 0 0], [0 1 0], w*dtheta2);
      RB1 = axrot_matrix([0 h*.25 0], [0 h*.25 1], -(2-w)*dtheta);
      RB2 = RT2';
    else
      RT1 = axrot_matrix([0 h*.75 0], [0 h*.75 1], 0);
      RT2 = axrot_matrix([0 0 0], [0 1 0], (3-w)*2*dtheta2);
      RB1 = axrot_matrix([0 h*.25 0], [0 h*.25 1], 0);
      RB2 = RT2';    
    end

    for i=i_top
      ro = RT1*RT2*[P(:,i);1];
      R(:, i) = ro(1:3) - P(:,i);
    end

    for i=i_bot
      ro = RB1*RB2*[P(:,i);1];
      R(:, i) = ro(1:3) - P(:,i);
    end

    gf_model_set(md, 'variable', 'DirichletData', R);
    gf_model_get(md, 'solve', 'very noisy', 'max_iter', 100, 'max_res', 1e-5, 'lsearch', 'simplest');
    // full(gf_model_get(md, 'tangent matrix'))
    U  = gf_model_get(md, 'variable', 'u');
    VM = gf_model_get(md, 'compute Von Mises or Tresca', 'u', lawname, 'params', mfdu);
    UU  = [UU;U]; 
    VVM = [VVM;VM];
    save(path + '/demo_nonlinear_elasticity_U.mat', 'UU', 'VVM', 'm_char', 'mfu_char', 'mfdu_char');
  else
    U  = UU(step,:);
    VM = VVM(step,:);
  end
  disp(sprintf('step %d/%d : |U| = %g',step,nbstep,norm(U)));

  drawlater;
  clf();
  h_graph = gcf();
  h_graph.color_map = jetcolormap(255);
  gf_plot(mfdu,VM,'mesh','off', 'cvlst',gf_mesh_get(mfdu,'outer faces'), 'deformation',U,'deformation_mf',mfu,'deformation_scale', 1, 'refine', 8);
  colorbar(min(U),max(U));
  h_graph.color_map = jetcolormap(255);
  drawnow;
  sleep(1000); 
  // save a picture..
  xs2png(h_graph.figure_id, path + sprintf('/torsion%03d.png',step));
end
  
printf('end of computations, you can now replay the animation with\n');
printf('exec demo_nonlinear_elasticity_anim.sce;\n');

printf('demo nonlinear_elasticity terminated\n');
