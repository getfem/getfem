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

gf_workspace('clear all');
gf_util('trace level', 1);

% set a custom colormap
r=[0.7 .7 .7]; l = r(end,:); s=63; s1=20; s2=25; s3=48;s4=55;
for i=1:s, c1 = max(min((i-s1)/(s2-s1),1),0);c2 = max(min((i-s3)/(s4-s3),1),0); r(end+1,:)=(1-c2)*((1-c1)*l + c1*[1 0 0]) + c2*[1 .8 .2]; end;
colormap(r);

dirichlet_version = 2; % 1 = simplification, 2 = penalisation
drawing = true;
test_tangent_matrix = false;
incompressible = false;

lawname = 'Ciarlet Geymonat';
params = [1;1;0.25];
% lawname = 'SaintVenant Kirchhoff';
% params = [1;1];
if (incompressible)
    lawname = 'Incompressible Mooney Rivlin';
    params = [1;1];
end

N1=2; N2=4; h=20;
m=gf_mesh('cartesian',(0:N1)/N1 - .5, (0:N2)/N2*h, ((0:N1)/N1 - .5)*3);
mfu=gf_mesh_fem(m,3);     % mesh-fem supporting a 3D-vector field
% the mesh_im stores the integration methods for each tetrahedron
mim=gf_mesh_im(m,gf_integ('IM_GAUSS_PARALLELEPIPED(3,4)'));
% we choose a P2 fem for the main unknown
gf_mesh_fem_set(mfu, 'fem',gf_fem('FEM_QK(3,2)'));
mfdu=gf_mesh_fem(m,1);
% the material is homogeneous, hence we use a P0 fem for the data
if (dirichlet_version == 1)
  mfd=mfu;
else
  mfd=gf_mesh_fem(m,1);     % scalar mesh_fem
  gf_mesh_fem_set(mfd,'fem',gf_fem('FEM_QK(3,1)'));
end
% the P2 fem is not derivable across elements, hence we use a discontinuous
% fem for the derivative of U.
gf_mesh_fem_set(mfdu,'fem',gf_fem('FEM_QK_DISCONTINUOUS(3,2)'));


m_char=gf_mesh_get(m, 'char');
mfu_char=gf_mesh_fem_get(mfu, 'char');
mfdu_char=gf_mesh_fem_get(mfdu, 'char');

% display some informations about the mesh
disp(sprintf('nbcvs=%d, nbpts=%d, nbdof=%d',gf_mesh_get(m,'nbcvs'),...
            gf_mesh_get(m,'nbpts'),gf_mesh_fem_get(mfu,'nbdof')));
P=gf_mesh_get(m,'pts'); % get list of mesh points coordinates
%pidtop=find(abs(P(2,:)-13)<1e-6); % find those on top of the object
%pidbot=find(abs(P(2,:)+10)<1e-6); % find those on the bottom

pidtop=find(abs(P(2,:)-h)<1e-6); % find those on top of the object
pidbot=find(abs(P(2,:)-0)<1e-6); % find those on the bottom


% build the list of faces from the list of points
ftop=gf_mesh_get(m, 'faces from pid', pidtop); 
fbot=gf_mesh_get(m, 'faces from pid', pidbot);
% assign boundary numbers
gf_mesh_set(m,'boundary', 1, ftop);
gf_mesh_set(m,'boundary', 2, fbot);
gf_mesh_set(m,'boundary', 3, [ftop fbot]);


md=gf_model('real');
gf_model_set(md, 'add fem variable', 'u', mfu);
gf_model_set(md,'add initialized data','params', params);
% gf_model_set(md, 'add nonlinear elasticity brick', mim, 'u', lawname, 'params');
gf_model_set(md, 'add finite strain elasticity brick', mim, lawname, 'u',  'params');
% gf_model_set(md, 'add nonlinear term', mim, ...
%             'sqr(Trace(Green_Lagrangian(Id(meshdim)+Grad_u)))/8 + Norm_sqr(Green_Lagrangian(Id(meshdim)+Grad_u))/4');
% gf_model_set(md, 'add nonlinear term', mim, ...
%            '((Id(meshdim)+Grad_u)*(params(1)*Trace(Green_Lagrangian(Id(meshdim)+Grad_u))*Id(meshdim)+2*params(2)*Green_Lagrangian(Id(meshdim)+Grad_u))):Grad_Test_u');
% gf_model_set(md, 'add nonlinear term', mim, 'Saint_Venant_Kirchhoff_potential(Grad_u,params)'); 
% gf_model_set(md, 'add nonlinear term', mim, ...
%           '((Id(meshdim)+Grad_u)*(Ciarlet_Geymonat_PK2(Grad_u,params))):Grad_Test_u');
% gf_model_set(md, 'add nonlinear term', mim, ...
%                'Ciarlet_Geymonat_potential(Grad_u,params)');
% gf_model_set(md, 'add nonlinear term', mim, ...
%         '((Id(meshdim)+Grad_u)*(Incompressible_Mooney_Rivlin_PK2(Grad_u,params))):Grad_Test_u');
% gf_model_set(md, 'add nonlinear term', mim, ...
%                  'Incompressible_Mooney_Rivlin_potential(Grad_u,params)');
% gf_model_set(md, 'add nonlinear term', mim, ...
%      '((Id(meshdim)+Grad_u)*(Saint_Venant_Kirchhoff_PK2(Grad_u,params))):Grad_Test_u');

    
if (incompressible)
  mfp = gf_mesh_fem(m,1); 
  % gf_mesh_fem_set(mfp, 'classical discontinuous fem', 1);
  gf_mesh_fem_set(mfp, 'classical fem', 1);
  gf_model_set(md, 'add fem variable', 'p', mfp);
  % gf_model_set(md, 'add nonlinear incompressibility brick',  mim, 'u', 'p');
  gf_model_set(md, 'add finite strain incompressibility brick',  mim, 'u', 'p');
  % gf_model_set(md, 'add nonlinear term', mim, ...
  %                'p*(1-Det(Id(meshdim)+Grad_u))');
  % gf_model_set(md, 'add nonlinear term', mim, ...
  %                 '-p*Det(Id(meshdim)+Grad_u)*(Inv(Id(meshdim)+Grad_u))'':Grad_Test_u + Test_p*(1-Det(Id(meshdim)+Grad_u))');
end

if (dirichlet_version == 1)
  gf_model_set(md, 'add fem data', 'DirichletData', mfu);
  gf_model_set(md, 'add Dirichlet condition with simplification', 'u', 3, 'DirichletData');
else
  gf_model_set(md, 'add fem data', 'DirichletData', mfd, 3);
  gf_model_set(md, 'add Dirichlet condition with penalization', mim, 'u', 1e4, 3, 'DirichletData');
end

VM=zeros(1,gf_mesh_fem_get(mfdu, 'nbdof'));

reload = 0;

if (reload == 0),
  UU=[];
  VVM=[];
  nbstep=40;
else
  load 'demo_nonlinear_elasticity_U.mat';
  nb_step = size(UU,1);
end;


P=gf_mesh_fem_get(mfd, 'basic dof_nodes');
r = sqrt(P(1 ,:).^2 + P(3, :).^2);
theta = atan2(P(3,:),P(1,:));

for step=1:nbstep,
  w = 3*step/nbstep;
  %set(b2, 'param', 'R', [0;0;0]);

  if (~reload)
    dtheta =  pi;
    dtheta2 = pi/2;
      
    if (dirichlet_version == 1)   
      R=zeros(gf_mesh_fem_get(mfd, 'nbdof'), 1);
    else
      R=zeros(3, gf_mesh_fem_get(mfd, 'nbdof'));
    end
    
    i_top = gf_mesh_fem_get(mfd, 'basic dof on region', 1);
    i_bot = gf_mesh_fem_get(mfd, 'basic dof on region', 2);
    
    
    dd = max(P(1,i_top)*sin(w*dtheta));
    if (w < 1), 
      RT1 = axrot_matrix([0 h*.75 0], [0 h*.75 1], w*dtheta);
      RT2 = axrot_matrix([0 0 0], [0 1 0], sqrt(w)*dtheta2);
      RB1 = axrot_matrix([0 h*.25 0], [0 h*.25 1], -w*dtheta);
      RB2 = RT2';
    elseif (w < 2),
      RT1 = axrot_matrix([0 h*.75 0], [0 h*.75 1], (2-w)*dtheta);
      RT2 = axrot_matrix([0 0 0], [0 1 0], w*dtheta2);
      RB1 = axrot_matrix([0 h*.25 0], [0 h*.25 1], -(2-w)*dtheta);
      RB2 = RT2';
    else
      RT1 = axrot_matrix([0 h*.75 0], [0 h*.75 1], 0);
      RT2 = axrot_matrix([0 0 0], [0 1 0], (3-w)*2*dtheta2);
      RB1 = axrot_matrix([0 h*.25 0], [0 h*.25 1], 0);
      RB2 = RT2';    
    end;
    if (dirichlet_version == 1)
      for i=i_top,
        ro = RT1*RT2*[P(:,i);1];
        R(i) = ro(1+mod(i-1,3)) - P(1+mod(i-1,3),i);
      end
      for i=i_bot,
        ro = RB1*RB2*[P(:,i);1];
        R(i) = ro(1+mod(i-1,3)) - P(1+mod(i-1,3),i);
      end 
    else
      for i=i_top,
        ro = RT1*RT2*[P(:,i);1];
        R(:, i) = ro(1:3) - P(:,i);
      end
      for i=i_bot,
        ro = RB1*RB2*[P(:,i);1];
        R(:, i) = ro(1:3) - P(:,i);
      end
    end
    
    gf_model_set(md, 'variable', 'DirichletData', R);
    gf_model_get(md, 'solve', 'noisy', 'max_iter', 100, 'max_res', 1e-5); %  , 'lsearch', 'simplest');
      
    if (test_tangent_matrix)
      gf_model_get(md, 'test tangent matrix', 1E-8, 10, 0.0001);
    end;
       
    U = gf_model_get(md, 'variable', 'u');
    % VM0 = gf_model_get(md, 'compute Von Mises or Tresca', 'u', lawname, 'params', mfdu);
    % sigma = gf_model_get(md, 'compute second Piola Kirchhoff tensor', 'u', lawname, 'params', mfdu);
    
    % Direct interpolation of the Von Mises stress
    % VM = gf_model_get(md, 'interpolation', '(sqrt(3/2)/Det(Id(meshdim)+Grad_u))*Norm((Id(meshdim)+Grad_u)*Saint_Venant_Kirchhoff_PK2(Grad_u,params)*(Id(meshdim)+Grad_u'') - Id(meshdim)*Trace((Id(meshdim)+Grad_u)*Saint_Venant_Kirchhoff_PK2(Grad_u,params)*(Id(meshdim)+Grad_u''))/meshdim)', mfdu);
    % VM = gf_model_get(md, 'interpolation', '(sqrt(3/2)/Det(Id(meshdim)+Grad_u))*Norm(Deviator((Id(meshdim)+Grad_u)*Saint_Venant_Kirchhoff_PK2(Grad_u,params)*(Id(meshdim)+Grad_u'')))', mfdu);
    % VM = gf_model_get(md, 'interpolation', 'sqrt(3/2)*Norm(Deviator(Cauchy_stress_from_PK2(Saint_Venant_Kirchhoff_PK2(Grad_u,params),Grad_u)))', mfdu);
    VM = gf_model_get(md, 'compute finite strain elasticity Von Mises', lawname, 'u', 'params', mfdu);
    % norm(VM-VM0)
    
    UU = [UU;U]; 
    VVM = [VVM;VM];
    save demo_nonlinear_elasticity_U.mat UU VVM m_char mfu_char mfdu_char;
  else
    U=UU(step,:);
    VM=VVM(step,:);
  end;
  disp(sprintf('step %d/%d : |U| = %g',step,nbstep,norm(U)));

  if (drawing)
    gf_plot(mfdu,VM,'mesh','off', 'cvlst',gf_mesh_get(mfdu,'outer faces'), 'deformation',U,'deformation_mf',mfu,'deformation_scale', 1, 'refine', 8); colorbar;
    axis([-3     6     0    20    -2     2]); caxis([0 .3]);
    view(30+20*w, 23+30*w);  
    campos([50 -30 80]);
    camva(8);
    camup
    camlight; 
    axis off;
    pause(1);
    % save a picture..
    %print(gcf, '-dpng', '-r150', sprintf('torsion%03d',step));
  end
end;
  
disp('end of computations, you can now replay the animation with')
disp('demo_nonlinear_elasticity_anim')

