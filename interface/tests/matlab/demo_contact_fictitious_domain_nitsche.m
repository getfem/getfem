% Copyright (C) 2006-2020 Mathieu Fabre.
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

 

disp('Resolution of a contact problem in 2D or 3D with two elastics bodies');
disp('with a fictitious domain method and Nitsche s method');
clear all;

draw_mesh = false;

% gf_workspace('clear all');

ref_sol = 0;% 0 Reference solution
            % 1 Error in L2 and H1 in Omega1 Omega2
            % 2 link between gamma0 and theta
            
N = 2; % 2 ou 3 dimensions
R = 0.25;
dirichlet_val = 0;
f_coeff = 0.1;

if (ref_sol == 0)
    method = [-1]; % theta
    gamma = [1/200]; % 1/200
    nxy = [60]; % 2D ->400 and 3D -> 30
    ls_degree = 2;
    penalty_parameter = 10E-8;
    vertical_force = -0.1;
end
if (ref_sol == 1)
    method = [-1];
    gamma = [1/200];
    if (N==2)
        nxy=[5 15 10 15 27 37 51 61 75 90 101 121];
    else
        nxy = [5 10 15 20 25 30];
    end
    ls_degree = 2;
    penalty_parameter = 10E-8;
    vertical_force = -0.1;  
end
if (ref_sol == 2)
    method = [0 1 -1];
    gamma = [400 200 100 50 25 10 1 1/10 1/25 1/50 1/100 1/200 1/400];
    nxy = [21];
    ls_degree = 2;
    penalty_parameter = 10E-8;
    vertical_force = -0.1;
end


for xx = 1:1:size(method,2)
for yy = 1:1:size(gamma,2)
for zz = 1:1:size(nxy,2)
theta = method(xx);
gamma0 = gamma(yy);
NX = nxy(zz);

% Definition of fictitious domain mesh with quadrangles and order 1 of level-set

if (N==2)
    m=gf_mesh('regular simplices', -.5:(1/NX):.5, -.5:(1/NX):.5);
    %m=gf_mesh('cartesian', -.5:(1/NX):.5, -.5:(1/NX):.5);
else
    m=gf_mesh('regular simplices', -.5:(1/NX):.5, -.5:(1/NX):.5, -.5:(1/NX):.5);
end


% gf_plot_mesh(m, 'convexes', 'on');
% pause;

ls1=gf_levelset(m, ls_degree);
ls2=gf_levelset(m, ls_degree);
mf_ls1=gfObject(gf_levelset_get(ls1, 'mf'));
mf_ls2=gfObject(gf_levelset_get(ls2, 'mf'));
mfu=gfMeshFem(m,N);
if (N==2)
    set(mfu, 'fem', gf_fem('FEM_PK(2,2)')); % set(mfu, 'fem', gf_fem('FEM_PK(2,2)'));
else
    set(mfu, 'fem', gf_fem('FEM_PK(3,2)')); % set(mfu, 'fem', gf_fem('FEM_PK(2,2)'));
end
mfvm=gfMeshFem(m,1);
if (N==2)
    set(mfvm, 'fem', gf_fem('FEM_PK_DISCONTINUOUS(2,1)')); % set(mfvm, 'fem', gf_fem('FEM_PK_DISCONTINUOUS(2,1)'));
else
    set(mfvm, 'fem', gf_fem('FEM_PK_DISCONTINUOUS(3,1)')); % set(mfvm, 'fem', gf_fem('FEM_PK_DISCONTINUOUS(2,1)'));
end

mls1=gfMeshLevelSet(m);
mls2=gfMeshLevelSet(m);

% definition of Omega 1 (circle)
 
P=get(mf_ls1, 'basic dof nodes');
if (N==2)
    x = P(1,:); y = P(2,:);
    ULS1=1000*ones(1,numel(x));
    ULS1 = min(ULS1, sqrt(x.^2 + y.^2) - R);
else
    x = P(1,:); y = P(2,:); z = P(3,:);
    ULS1=1000*ones(1,numel(x));
    ULS1 = min(ULS1, sqrt(x.^2 + y.^2 + z.^2) - R);
end
gf_levelset_set(ls1, 'values', ULS1);

% definition of Omega 2 (rectangle)

P=get(mf_ls2, 'basic dof nodes');
if (N==2)
    x = P(1,:); y = P(2,:);
    ULS2=1000*ones(1,numel(x));
    yc = -0.25; xc=0; 
    ULS2=min(ULS2,y-yc);
else
    x = P(1,:); y = P(2,:); z = P(3,:);
    ULS2=1000*ones(1,numel(x));
    zc = -0.25; xc=0; 
    ULS2=min(ULS2,z-zc);
end
gf_levelset_set(ls2, 'values', ULS2);

set(mls1, 'add', ls1);
set(mls1, 'adapt');

set(mls2, 'add', ls2);
set(mls2, 'adapt');

% Dirichlet's boundary
 
GAMMAC = 1; GAMMAD = 2;


border = gf_mesh_get(m,'outer faces');
normals = gf_mesh_get(m, 'normal of faces', border);
if (N==2)
    contact_boundary=border(:, find(normals(2, :) < -0.01)); % Normal vector is -e2
else
    contact_boundary=border(:, find(normals(3, :) < -0.01)); % Normal vector is -e3
end
gf_mesh_set(m, 'region', GAMMAD, contact_boundary);



%figure 1 : plot figure

if (draw_mesh)
  figure(1);  
  if (N==2)
    gf_plot_mesh(get(mls1,'cut mesh')); % ,'curved', 'on'
    hold on; gf_plot_mesh(get(mls2,'cut mesh')); hold off;
    hold on; gf_plot_mesh(m, 'regions', GAMMAD, 'convexes', 'on'); %plot de bord avec condition de type Dirichlet
    title('boundary with Dirichlet condition in red');hold off;
  else
   
   gf_plot_mesh(get(mls1,'cut mesh')); % ,'curved', 'on'
    hold on; gf_plot_mesh(get(mls2,'cut mesh')); hold off;
   hold on; gf_plot_mesh(m, 'regions', GAMMAD, 'convexes', 'on'); %plot de bord avec condition de type Dirichlet
   title('boundary with Dirichlet condition in red');hold off;
    xlabel('x'); ylabel('y'); zlabel('z');
   title('Displacement solution'); 
  end
end


%Finites elements' method on mls1 and mls2

if (N==2)
    mim = gfMeshIm('levelset', mls1,'all', gf_integ('IM_TRIANGLE(5)'));
    mim_bound = gfMeshIm('levelset', mls1, 'boundary', gf_integ('IM_TRIANGLE(5)'));
    mim1 = gfMeshIm('levelset', mls1, 'inside', gf_integ('IM_TRIANGLE(5)')); 
    mim2 = gfMeshIm('levelset', mls2, 'inside', gf_integ('IM_TRIANGLE(5)'));
else
    mim_bound = gfMeshIm('levelset', mls1, 'boundary', gf_integ('IM_TETRAHEDRON(5)'));
    mim = gfMeshIm('levelset', mls1, 'all', gf_integ('IM_TETRAHEDRON(5)')); 
    mim1 = gfMeshIm('levelset', mls1, 'inside', gf_integ('IM_TETRAHEDRON(5)')); 
    mim2 = gfMeshIm('levelset', mls2, 'inside', gf_integ('IM_TETRAHEDRON(5)'));
end

set(mim, 'integ', 4);
set(mim1, 'integ', 4);
set(mim2, 'integ', 4);

dof_out = get(mfu, 'dof from im', mim1);
cv_out = get(mim1, 'convex_index');
cv_in = setdiff(gf_mesh_get(m, 'cvid'), cv_out);
mfu1 = gfMeshFem('partial', mfu, dof_out, cv_in);

dof_out = get(mfu, 'dof from im', mim2);
cv_out = get(mim2, 'convex_index');
cv_in = setdiff(gf_mesh_get(m, 'cvid'), cv_out);
mfu2 = gfMeshFem('partial', mfu, dof_out, cv_in);

%Elastic model 

md=gf_model('real');
gf_model_set(md,'add fem variable', 'u1', mfu1);
gf_model_set(md,'add fem variable', 'u2', mfu2);
gf_model_set(md,'add initialized fem data', 'd1', mf_ls1, ULS1);
gf_model_set(md,'add initialized fem data', 'd2', mf_ls2, ULS2);
gf_model_set(md,'add initialized data', 'gamma0', gamma0);


gf_model_set(md, 'add initialized data', 'friction_coeff',[f_coeff]);
clambda = 1;           % Lame coefficient
cmu = 1;               % Lame coefficient
gf_model_set(md, 'add initialized data', 'cmu', [cmu]);
gf_model_set(md, 'add initialized data', 'clambda', [clambda]);
gf_model_set(md, 'add isotropic linearized elasticity brick', mim1, 'u1','clambda', 'cmu');
gf_model_set(md, 'add isotropic linearized elasticity brick', mim2, 'u2','clambda', 'cmu');

if (N==2)
    gf_model_set(md, 'add initialized data', 'Fdata', [0 vertical_force]);
else
    gf_model_set(md, 'add initialized data', 'Fdata', [0 0 vertical_force]);
end

gf_model_set(md, 'add source term brick', mim1, 'u1', 'Fdata');

if (N==2)
    Ddata = zeros(1, 2);
else
    Ddata = zeros(1, 3);
end

gf_model_set(md, 'add initialized data', 'Ddata', Ddata);
gf_model_set(md, 'add Dirichlet condition with simplification', 'u2', GAMMAD, 'Ddata'); 

if (N==2)
    cpoints = [0, 0,   0, 0.1]; % constrained points for 2d
    cunitv  = [1, 0,  1, 0];    % corresponding constrained directions for 2d, better with [0, 0.1]
else
    cpoints = [0, 0, 0,    0, 0, 0,   0, 0, 0.1]; % constrained points for 3d
    cunitv  = [1, 0, 0,   0, 1, 0,   0, 1, 0];    % corresponding constrained directions for 3d, better with [0, 0.1]
end

gf_model_set(md, 'add initialized data', 'cpoints', cpoints);
gf_model_set(md, 'add initialized data', 'cunitv', cunitv);
gf_model_set(md, 'add pointwise constraints with multipliers', 'u1', 'cpoints', 'cunitv');


gf_model_set(md, 'add initialized data', 'penalty_param1', [penalty_parameter]);
indmass = gf_model_set(md, 'add mass brick', mim1, 'u1', 'penalty_param1');

% gf_model_set(md, 'add initialized data', 'penalty_param2', [penalty_parameter]);
% indmass = gf_model_set(md, 'add mass brick', mim2, 'u2', 'penalty_param2');

gf_model_set(md,'add Nitsche fictitious domain contact brick', mim_bound, 'u1', 'u2', 'd1', 'd2', 'gamma0', theta, 'friction_coeff'); 

disp('solve');

% niter= 20;
% gf_model_get(md, 'test tangent matrix term', 'u1', 'u2', 1e-6, niter, 10.0);
% gf_model_get(md, 'test tangent matrix', 1e-6, niter, 10);
% gf_model_get(md, 'test tangent matrix', 1e-6, 20, 10);

niter= 100;
gf_model_get(md, 'solve', 'max_res', 1E-9, 'max_iter', niter, 'noisy');

U1 = gf_model_get(md, 'variable', 'u1');  UU1 = gf_model_get(md, 'variable', 'u1');
VM1 = gf_model_get(md, 'compute_isotropic_linearized_Von_Mises_or_Tresca', ...
                          'u1', 'clambda', 'cmu', mfvm);
U2 = gf_model_get(md, 'variable', 'u2');  UU2 = gf_model_get(md, 'variable', 'u2');
VM2 = gf_model_get(md, 'compute_isotropic_linearized_Von_Mises_or_Tresca', ...
                          'u2', 'clambda', 'cmu', mfvm);

% plot figure

if (ref_sol == 0)
    if (N==2)
        sl1=gf_slice({'isovalues', -1, mf_ls1, ULS1, 0}, m, 5);
        P1=gf_slice_get(sl1,'pts'); dP1=gf_compute(mfu1,U1,'interpolate on',sl1);
        gf_slice_set(sl1, 'pts', P1 + dP1);
        VMsl1=gf_compute(mfvm,VM1,'interpolate on',sl1);
        sl2=gf_slice({'isovalues', -1, mf_ls2, ULS2, 0}, m, 5);
        P2=gf_slice_get(sl2,'pts'); dP2=gf_compute(mfu2,U2,'interpolate on',sl2);
        gf_slice_set(sl2, 'pts', P2+dP2);
        VMsl2=gf_compute(mfvm,VM2,'interpolate on',sl2);
        gf_plot_slice(sl1,'mesh','off','mesh_slice_edges','off','data',VMsl1);
        hold on;
        gf_plot_slice(sl2,'mesh','off','mesh_slice_edges','off','data',VMsl2);
        hold off;   
    else
        sl1=gf_slice({'boundary',{'intersection',{'ball',-1,[0;0;0],R},{'planar',1,[0;0;0],[1;0;0]}}}, m, 5);
        sl2=gf_slice({'boundary',{'intersection',{'planar',-1,[0;0;zc],[0;0;1]},{'planar',1,[0;0;0],[1;0;0]}}}, m, 5); 
        P1=gf_slice_get(sl1,'pts');P2=gf_slice_get(sl2,'pts');
        dP1=gf_compute(mfu1,U1,'interpolate on',sl1);
        dP2=gf_compute(mfu2,U2,'interpolate on',sl2);
        gf_slice_set(sl1, 'pts', P1+dP1);
        gf_slice_set(sl2, 'pts', P2+dP2);           
        VMsl1=gf_compute(mfvm,VM1,'interpolate on',sl1);
        VMsl2=gf_compute(mfvm,VM2,'interpolate on',sl2);
        set(gcf,'renderer','zbuffer');
        h=gf_plot_slice(sl1,'mesh','off','mesh_faces','on','mesh_slice_edges','on','data',VMsl1);
        hold on;
        h=gf_plot_slice(sl2,'mesh','off','mesh_slice_edges','off','data',VMsl2);
        hold off;
        view(-55,10); axis on; camlight('headlight'); %gf_colormap('tank');
        xlabel('x'); ylabel('y'); zlabel('z');
        %title('3D');       
        
    end
end;  


%plot false fictious domain

%         m_fict=gf_mesh('regular simplices', -.5:(1/10):.5, -.5:(1/10):.5);
%         hold on; gf_plot_mesh(m_fict); %plot de bord avec condition de type Dirichlet
%         hold off;


        map(1,:)=[0 200 255]; % bleu foncé
        map(2,:)=[0 255 255];
        map(3,:)=[128 255 128];
        map(4,:)=[255 255 0];
        map(5,:)=[255 128 0];
        map(6,:)=[255 0 0];
        colormap(map./255)

% save the reference solution if res_sol= 0 and else errors in L2 and H1  
% Discrétisation de réf N > n/4 

if (ref_sol == 0)
    gf_mesh_fem_get(mfu, 'save', 'sol_ref_mesh_fem','with_mesh');
    save sol_de_reference1 UU1;
    save sol_de_reference2 UU2;
    dls = gf_levelset_get(ls1, 'degree');
    save degree_levelset dls;
else
    meshref = gf_mesh('load', 'sol_ref_mesh_fem');
  
    % Reconstruction of mfuref, mfu1ref, mfu2ref, min1ref, min2ref
  
    dlsref = load('degree_levelset', 'dls');
    ls1ref=gf_levelset(meshref, dlsref.dls);
    ls2ref=gf_levelset(meshref, dlsref.dls);
    mf_ls1ref=gfObject(gf_levelset_get(ls1ref, 'mf'));
    mf_ls2ref=gfObject(gf_levelset_get(ls2ref, 'mf'));
    mfuref=gfMeshFem(meshref,N);
    if(N==2)
        set(mfuref, 'fem', gf_fem('FEM_PK(2,2)'));
    else
        set(mfuref, 'fem', gf_fem('FEM_PK(3,2)'));
    end
    mls1ref=gfMeshLevelSet(meshref);
    mls2ref=gfMeshLevelSet(meshref);

    P=get(mf_ls1ref, 'basic dof nodes');
    if(N==2)
        x = P(1,:); y = P(2,:);
        ULS1ref=1000*ones(1,numel(x));
        ULS1ref = min(ULS1ref, sqrt(x.^2 + y.^2) - R);
    else
        x = P(1,:); y = P(2,:); z = P(3,:);
        ULS1ref=1000*ones(1,numel(x));
        ULS1ref = min(ULS1ref, sqrt(x.^2 + y.^2 + z.^2) - R);
    end
    gf_levelset_set(ls1ref, 'values', ULS1ref);

    P=get(mf_ls2ref, 'basic dof nodes');
    if(N==2)
        x = P(1,:); y = P(2,:);
        ULS2ref=1000*ones(1,numel(x));
        yc = -0.25; xc=0; 
        ULS2ref=min(ULS2ref,y-yc);
    else
        x = P(1,:); y = P(2,:); z = P(3,:);
        ULS2ref=1000*ones(1,numel(x));
        zc = -0.25; xc=0; 
        ULS2ref=min(ULS2ref,z-zc);
    end
    gf_levelset_set(ls2ref, 'values', ULS2ref); 

    set(mls1ref, 'add', ls1ref);
    set(mls1ref, 'adapt');

    set(mls2ref, 'add', ls2ref);
    set(mls2ref, 'adapt');
    if(N==2)
        mim1ref = gfMeshIm('levelset', mls1ref, 'inside', gf_integ('IM_TRIANGLE(5)')); 
        mim2ref = gfMeshIm('levelset', mls2ref, 'inside', gf_integ('IM_TRIANGLE(5)')); 
    else
        mim1ref = gfMeshIm('levelset', mls1ref, 'inside', gf_integ('IM_TETRAHEDRON(5)')); 
        mim2ref = gfMeshIm('levelset', mls2ref, 'inside', gf_integ('IM_TETRAHEDRON(5)')); 
    end
    set(mim1ref, 'integ', 4);
    set(mim2ref, 'integ', 4);


    dof_out = get(mfuref, 'dof from im', mim1ref);
    cv_out = get(mim1ref, 'convex_index');
    cv_in = setdiff(gf_mesh_get(meshref, 'cvid'), cv_out);
    mfu1ref = gfMeshFem('partial', mfuref, dof_out, cv_in);

    dof_out = get(mfuref, 'dof from im', mim2ref);
    cv_out = get(mim2ref, 'convex_index');
    cv_in = setdiff(gf_mesh_get(meshref, 'cvid'), cv_out);
    mfu2ref = gfMeshFem('partial', mfuref, dof_out, cv_in);

    U1ref = load('sol_de_reference1', 'UU1');
    U2ref = load('sol_de_reference2', 'UU2');

    bB1 = gf_mesh_fem_get(mfu1, 'extension matrix');
    U1e = gf_compute(mfu, U1*bB1', 'interpolate on', mfu1ref); 
    bB2 = gf_mesh_fem_get(mfu2, 'extension matrix');
    U2e = gf_compute(mfu, U2*bB2', 'interpolate on', mfu2ref);
    n_tot1 = gf_compute(mfu1ref, U1e-U1ref.UU1, 'L2 norm', mim1ref);
    n_tot2 = gf_compute(mfu2ref, U2e-U2ref.UU2, 'L2 norm', mim2ref);
    n_ref1 = gf_compute(mfu1ref, U1ref.UU1, 'L2 norm',  mim1ref);
    n_ref2 = gf_compute(mfu2ref, U2ref.UU2, 'L2 norm',  mim2ref);
    m_tot1 = gf_compute(mfu1ref, U1e-U1ref.UU1, 'H1 norm', mim1ref);
    m_tot2 = gf_compute(mfu2ref, U2e-U2ref.UU2, 'H1 norm', mim2ref);
    m_ref1 = gf_compute(mfu1ref, U1ref.UU1, 'H1 norm',  mim1ref);
    m_ref2 = gf_compute(mfu2ref, U2ref.UU2, 'H1 norm',  mim2ref);
    n1 = 100*n_tot1/n_ref1; 
    n2 = 100*n_tot2/n_ref2;
    m1 = 100*m_tot1/m_ref1;
    m2 = 100*m_tot2/m_ref2;
    %nddl(zz)= gf_model_get(md,'nbdof');
    Y11(yy,zz,xx)=n1;
    Y12(yy,zz,xx)=n2;
    Y21(yy,zz,xx)=m1;
    Y22(yy,zz,xx)=m2;
  end
end
end
end

if (ref_sol == 1 ) % Curve of error depending of h
    Y11(1,:,1)
    Y12(1,:,1)
    Y21(1,:,1)
    Y22(1,:,1)
    X=1./nxy % [1/10 1/16 1/22 1/28 1/34 1/40];

   
    figure(1);
    msize= size(X,2);
    loglog(X,Y11(1,:,1),'o-k', 'linewidth', 2, 'MarkerSize', 15 )
    hold on;
    loglog(X,Y12(1,:,1),'+--k', 'linewidth', 2, 'MarkerSize', 15);
    hold off;
    P1 = polyfit(log(X),log(Y11(1,:,1)),1); % the first and second are too bad;
    P2 = polyfit(log(X),log(Y12(1,:,1)),1);
    legend(strcat('norm on Omega 1  (slope=',num2str(P1(1)), ')'), ...
           strcat('norm on Omega 2 (slope=',num2str(P2(1)), ')'), ...
           'Location', 'NorthWest');
    grid on;
    axesobj = findobj('type', 'axes');
    set(axesobj, 'fontname', 'times'); set(axesobj, 'fontunits', 'points');
    set(axesobj, 'fontsize', 18); set(axesobj, 'fontweight', 'bold');
    set(axesobj, 'linewidth', 4);
    xlabel('h');
    ylabel('L^2 relative error (in %)');
    set(gca,'XTickLabel',{'0.01';'0.1';'1';'...'}) 
    

    figure(2);
 
    loglog(X,Y21(1,:,1),'o-k', 'linewidth', 2, 'MarkerSize', 15 )
    hold on;
    loglog(X,Y22(1,:,1),'+--k', 'linewidth', 2, 'MarkerSize', 15);
    hold off;
    P3 = polyfit(log(X),log(Y21(1,:,1)),1);
    P4 = polyfit(log(X),log(Y22(1,:,1)),1);
    legend(strcat('norm on Omega 1  (slope=',num2str(P3(1)), ')'), ...
           strcat('norm on Omega 2 (slope=',num2str(P4(1)), ')'), ...
           'Location', 'NorthWest');
    grid on;
    axesobj = findobj('type', 'axes');
    set(axesobj, 'fontname', 'times'); set(axesobj, 'fontunits', 'points');
    set(axesobj, 'fontsize', 18); set(axesobj, 'fontweight', 'bold');
    set(axesobj, 'linewidth', 4);
    xlabel('h');
    ylabel('H^1 relative error (in %)');
    set(gca,'XTickLabel',{'0.01';'0.1';'1';'...'}) 
theta
gamma0
N
       
end

if (ref_sol == 2 )  % Curve of error depending of gamma0
    
    % method = [0 1 -1];
    
    Y11(:,1,:)
    Y12(:,1,:)
    Y21(:,1,:)
    Y22(:,1,:)
    X = gamma 

   
    
    figure(1);
    

    loglog(X,Y11(:,1,1)','o-k', 'linewidth', 2, 'MarkerSize', 15 )
    hold on;
    loglog(X,Y11(:,1,2)','+-k', 'linewidth', 2, 'MarkerSize', 15 );
    loglog(X,Y11(:,1,3)','x-k', 'linewidth', 2, 'MarkerSize', 15 );
    hold off;
        
    P1 = polyfit(log(X),log(Y11(:,1,1)'),1); % the first and second are too bad;
    P2 = polyfit(log(X),log(Y11(:,1,2)'),1);
    P3 = polyfit(log(X),log(Y11(:,1,3)'),1);
    legend(strcat('norm for theta = 0 '), ...
           strcat('norm for theta = 1 '), ...
           strcat('norm for theta = -1'), ...
           'Location', 'NorthWest');
    grid on;
    axesobj = findobj('type', 'axes');
    set(axesobj, 'fontname', 'times'); set(axesobj, 'fontunits', 'points');
    set(axesobj, 'fontsize', 18); set(axesobj, 'fontweight', 'bold');
    set(axesobj, 'linewidth', 4);
    xlabel('gamma0');
    ylabel('L^2(Omega 1) relative error (in %)');
    % set(gca,'XTickLabel',{'0.01';'0.1';'1';'...'})  
    
    
    figure(2);
    

    loglog(X,Y12(:,1,1)','o-k', 'linewidth', 2, 'MarkerSize', 15 )
    hold on;
    loglog(X,Y12(:,1,2)','+-k', 'linewidth', 2, 'MarkerSize', 15 );
    loglog(X,Y12(:,1,3)','x-k', 'linewidth', 2, 'MarkerSize', 15 );
    hold off;
        
    P1 = polyfit(log(X),log(Y12(:,1,1)'),1); % the first and second are too bad;
    P2 = polyfit(log(X),log(Y12(:,1,2)'),1);
    P3 = polyfit(log(X),log(Y12(:,1,3)'),1);
    legend(strcat('norm for theta = 0 '), ...
           strcat('norm for theta = 1 '), ...
           strcat('norm for theta = -1'), ...
           'Location', 'NorthWest');
    grid on;
    axesobj = findobj('type', 'axes');
    set(axesobj, 'fontname', 'times'); set(axesobj, 'fontunits', 'points');
    set(axesobj, 'fontsize', 18); set(axesobj, 'fontweight', 'bold');
    set(axesobj, 'linewidth', 4);
    xlabel('gamma0');
    ylabel('L^2(Omega 2) relative error (in %)');
    % set(gca,'XTickLabel',{'0.01';'0.1';'1';'...'}) 
   
    
    
     figure(3);
    

    loglog(X,Y21(:,1,1)','o-k', 'linewidth', 2, 'MarkerSize', 15 )
    hold on;
    loglog(X,Y21(:,1,2)','+-k', 'linewidth', 2, 'MarkerSize', 15 );
    loglog(X,Y21(:,1,3)','x-k', 'linewidth', 2, 'MarkerSize', 15 );
    hold off;
        
    P1 = polyfit(log(X),log(Y21(:,1,1)'),1); % the first and second are too bad;
    P2 = polyfit(log(X),log(Y21(:,1,2)'),1);
    P3 = polyfit(log(X),log(Y21(:,1,3)'),1);
    legend(strcat('norm for theta = 0 '), ...
           strcat('norm for theta = 1 '), ...
           strcat('norm for theta = -1'), ...
           'Location', 'NorthWest');
    grid on;
    axesobj = findobj('type', 'axes');
    set(axesobj, 'fontname', 'times'); set(axesobj, 'fontunits', 'points');
    set(axesobj, 'fontsize', 18); set(axesobj, 'fontweight', 'bold');
    set(axesobj, 'linewidth', 4);
    xlabel('gamma0');
    ylabel('H^1(Omega 1) relative error (in %)');
    % set(gca,'XTickLabel',{'0.01';'0.1';'1';'...'}) 

    
    figure(4);
    

    loglog(X,Y22(:,1,1),'o-k', 'linewidth', 2, 'MarkerSize', 15 )
    hold on;
    loglog(X,Y22(:,1,2),'+-k', 'linewidth', 2, 'MarkerSize', 15 );
    loglog(X,Y22(:,1,3),'x-k', 'linewidth', 2, 'MarkerSize', 15 );
    hold off;
        
    P1 = polyfit(log(X),log(Y22(:,1,1)'),1); % the first and second are too bad;
    P2 = polyfit(log(X),log(Y22(:,1,2)'),1);
    P3 = polyfit(log(X),log(Y22(:,1,3)'),1);
    legend(strcat('norm for theta = 0 '), ...
           strcat('norm for theta = 1 '), ...
           strcat('norm for theta = -1'), ...
           'Location', 'NorthWest');
    grid on;
    axesobj = findobj('type', 'axes');
    set(axesobj, 'fontname', 'times'); set(axesobj, 'fontunits', 'points');
    set(axesobj, 'fontsize', 18); set(axesobj, 'fontweight', 'bold');
    set(axesobj, 'linewidth', 4);
    xlabel('gamma0');
    ylabel('H^1(Omega 2) relative error (in %)');
    % set(gca,'XTickLabel',{'0.01';'0.1';'1';'...'}) 
   
end










