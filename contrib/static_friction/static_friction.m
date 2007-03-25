% addpath ~/source++/getfem_toolbox
% addpath ~/source++/getfem++/contrib/static_friction/

gf_workspace('clear all');

option = 4; % 0 : reference solution
            % 1 : error
            % 2 : solution
            % 3 : contact stress
            % 4 : some convergence curves 

if (option == 0)
  mesh = gf_mesh('load', 'reference_sol2d.meshfem');
  mf = gf_mesh_fem('load', 'reference_sol2d.meshfem', mesh);
  mf_vm = gf_mesh_fem('load', 'reference_sol2d.meshfem_vm', mesh);
  U = load('reference_sol2d.U')';
  VM = load('reference_sol2d.VM')';
  gf_plot(mf_vm, VM, 'norm', 'on', 'refine', 1, 'deformation', U, 'deformation_mf', mf, 'deformed_mesh','off', 'deformation_scale', 1.0);
  gf_colormap('chouette');
elseif (option == 1)
  mesh = gf_mesh('load', 'reference_sol2d.meshfem');
  mf = gf_mesh_fem('load', 'reference_sol2d.meshfem', mesh);
  U = load('reference_sol2d.U')';
  UERR = load('reference_sol2d_error.U')';
  % UERR = log(abs(UERR));
  gf_plot(mf, UERR, 'norm', 'on', 'refine', 1, 'deformation', U, 'deformation_mf', mf, 'deformed_mesh','off', 'deformation_scale', 1.0);
  colorbar;
  gf_colormap('chouette');
elseif (option == 2)
  mesh = gf_mesh('load', 'static_friction.meshfem');
  N = gf_mesh_get(mesh, 'dim');
  mf = gf_mesh_fem('load', 'static_friction.meshfem', mesh);
  mf_vm = gf_mesh_fem('load', 'static_friction.meshfem_vm', mesh);
  U = load('static_friction.U')';
  VM = load('static_friction.VM')';
  if (N==2) 
    gf_plot(mf_vm, VM, 'norm', 'on', 'refine', 1, 'deformation', U, ...
      'deformation_mf', mf, 'deformed_mesh','on', 'deformation_scale', 1.0);
    gf_colormap('chouette');
  else

    sl1=gf_slice({'boundary',{'none'}},mesh,5);
    c=[0.1;0;20];x=[1;0;0];y=[0;1;0];z=[0;0;1];

    % trois plans de coupe:
    % sl2=gf_slice({'boundary',{'union',{'planar',+1,c,x},{'planar',+1,c,y},{'planar',+1,c,z}}},mesh,5);

    % deux plans de coupe:
    % sl2=gf_slice({'boundary',{'union',{'planar',+1,c,x},{'planar',+1,c,y}}},mesh,5);

    % un seul plan de coupe:
    sl2=gf_slice({'boundary',{'planar',+1,c,x}},mesh,5);


    P=gf_slice_get(sl2,'pts'); dP=gf_compute(mf,U,'interpolate on',sl2); gf_slice_set(sl2, 'pts', P+dP);

    VMsl=gf_compute(mf_vm,VM,'interpolate on',sl2);
    set(gcf,'renderer','zbuffer');
    figure(1); h=gf_plot_slice(sl2,'mesh','on','mesh_slice_edges','off','data',VMsl);
    view(-80,-15); axis on;
    camlight;
    gf_colormap('chouette');
    map=[1:-1/50:0]'*[1 1 1]; colormap(map); % for NB
    figure(2); h=gf_plot_slice(sl1,'mesh_faces','off','mesh','on'); view(-85,-15);
    axis on; camlight; set(h,'facecolor',[.8 0 0]);

    figure(3);
    gf_plot(mf_vm, VM, 'norm', 'on', 'refine', 4, 'deformation', U, ...
      'deformation_mf', mf, 'deformed_mesh','on', 'deformation_scale', ...
      1.0, 'cvlst', gf_mesh_get(mesh, 'outer faces'));
    view(-5,-10); camlight; colormap(map); xlabel('x'); ylabel('y'); zlabel('z');
    figure(1); 
    zlabel('z'); 
  end;
  xlabel('x'); ylabel('y');
  colorbar;
elseif (option == 3)
  mesh = gf_mesh('load', 'static_friction.meshfem');
  N = gf_mesh_get(mesh, 'dim');
  sll=gfSlice('load','static_friction.sl');
  LN=load('static_friction.LN')';
  P0=gf_slice_get(sll, 'pts');
  % [h1,h2,h3,h4]=gf_plot_slice(sll, 'tube','off', ...
  %                            'mesh_slice_edges_color',[.3 .3 .3]);
  hold on;
  gf_slice_set(sll,'pts',[P0 ; LN(N:N:size(LN,2))]);
  [hh1,hh2,hh3,hh4]=gf_plot_slice(sll, 'tube','off', ...
        'mesh_slice_edges_color','black','mesh_slice_edges_width',1.5);
  % sl=gfSlice('load','xfem_dirichlet.sl');
  gf_plot_mesh(mesh, 'edges_width', 1, 'curved', 'on', 'refine', 4, 'edges_color',[0.6 1 0.6] );
  
  npt = size(P0, 2);
  P0 = [P0;zeros(1,npt)];
  P1 = gf_slice_get(sll,'pts'); 
  lseg = gf_slice_get(sll,'splxs', 1);
  F=[lseg(1,:) lseg(2,:); lseg(2,:) npt+lseg(2,:); npt+lseg(1,:) npt+lseg(1,:)];
  % %F=[lseg; npt+lseg(2,:)];
  h=patch('Vertices',[P0 P1]', 'Faces', F');
  hold off;
  set(h,'FaceAlpha',0.3);
  set(h,'LineStyle','none');
  set(gcf,'renderer','opengl');
  set(gcf,'color','white');
  % axis off;
  view(3);
  campos([-145  -645  31])
  camtarget([-1.8 20 -1.98]);
  camva(4.45);
  camup([0.48 2.06 0.94]);
  % camzoom(1.2);
  % axis([-0.5000    0.5000   -0.5000    0.5000 -.5 .5]);
  % % print(gcf,'-dpng','-r300', 'lagrange_multipliers.png');
else
  H    = [0.075  0.15    0.25   0.5    1      2      3      4     5.5];
  % P1/P1 Non Augmenté
  L2_1 = [0.0020 0.0082  0.039  0.086  0.30   1.68   3.75   9.09];
  H1_1 = [0.042  0.072   0.136  0.284  0.64   2.04   4.20   9.50];
  % P1/P1 augmenté
  L2_2 = [0.0021 0.0079  0.039  0.087  0.31   1.68   3.75   9.09];
  H1_2 = [0.042  0.072   0.136  0.285  0.64   2.04   4.20   9.50];
  % P1/P0 augmenté
  L2_3 = [0.0027 0.010   0.021  0.092  0.45   2.54   7.28];
  H1_3 = [0.0418 0.072   0.132  0.29   0.71   2.79   7.55];
  % P1/P2 augmenté
  L2_4 = [0.003  0.0090  0.022  0.087  0.45   1.40   3.75   8.30   36.8];
  H1_4 = [0.0419 0.072   0.13   0.28   0.72   1.78   4.15   8.73   36.9];
  % P2/P1 augmenté
  L2_5 = [100    0.0045  0.0082 0.050  0.12   0.88   2.04   11.9];
  H1_5 = [100    0.0074  0.014  0.092  0.15   0.92   2.11   12.03];
  % P2/P0 augmenté
  L2_6 = [100    0.0039  0.013  0.048  0.18   0.59   1.74   6.92];
  H1_6 = [100    0.0097  0.021  0.059  0.22   0.71   1.91   7.01];
  % P2/P2 augmenté
  L2_7 = [100    0.00042 0.0022 0.010  0.048  0.13   0.58   4.06   4.89];
  H1_7 = [100    0.0052  0.010  0.047  0.10   0.30   0.79   4.13   5.15];

  loglog(H(1:7), L2_1(1:7), 'o-k');
  hold on;
  loglog(H(1:7), L2_2(1:7), 'x:k');
  loglog(H(1:7), L2_3(1:7), '+-.k');
  loglog(H(1:7), L2_4(1:7), '*--k');
  loglog(H(2:7), L2_5(2:7), 's-k');
  loglog(H(2:7), L2_6(2:7), 'd:k');
  loglog(H(2:7), L2_7(2:7), '<-.k');
  hold off;
  xlabel('h');
  ylabel('L^2(\Omega) error');
  P1 = polyfit(log(H(1:5)), log(L2_1(1:5)), 1);
  P2 = polyfit(log(H(1:5)), log(L2_2(1:5)), 1);
  P3 = polyfit(log(H(1:5)), log(L2_3(1:5)), 1);
  P4 = polyfit(log(H(1:5)), log(L2_4(1:5)), 1);
  P5 = polyfit(log(H(2:5)), log(L2_5(2:5)), 1);
  P6 = polyfit(log(H(2:5)), log(L2_6(2:5)), 1);
  P7 = polyfit(log(H(2:5)), log(L2_7(2:5)), 1);
  legend(strcat('P1/P1 org (slope=',num2str(P1(1)), ')'), ...
         strcat('P1/P1     (slope=',num2str(P2(1)), ')'), ...
         strcat('P1/P0     (slope=',num2str(P3(1)), ')'), ...
         strcat('P1/P2     (slope=',num2str(P4(1)), ')'), ...
         strcat('P2/P1     (slope=',num2str(P5(1)), ')'), ...
         strcat('P2/P0     (slope=',num2str(P6(1)), ')'), ...
         strcat('P2/P2     (slope=',num2str(P7(1)), ')'), ...
         'Location', 'NorthWest');
  grid on;
  axis([0.05 7 1e-4 10]);
  pause;
  loglog(H(1:7), H1_1(1:7), 'o-k');
  hold on;
  loglog(H(1:7), H1_2(1:7), 'x:k');
  loglog(H(1:7), H1_3(1:7), '+-.k');
  loglog(H(1:7), H1_4(1:7), '*--k');
  loglog(H(2:7), H1_5(2:7), 's-k');
  loglog(H(2:7), H1_6(2:7), 'd:k');
  loglog(H(2:7), H1_7(2:7), '<-.k');
  hold off;
  xlabel('h');
  ylabel('H^1(\Omega) error');
  P1 = polyfit(log(H(1:5)), log(H1_1(1:5)), 1);
  P2 = polyfit(log(H(1:5)), log(H1_2(1:5)), 1);
  P3 = polyfit(log(H(1:5)), log(H1_3(1:5)), 1);
  P4 = polyfit(log(H(1:5)), log(H1_4(1:5)), 1);
  P5 = polyfit(log(H(2:5)), log(H1_5(2:5)), 1);
  P6 = polyfit(log(H(2:5)), log(H1_6(2:5)), 1);
  P7 = polyfit(log(H(2:5)), log(H1_7(2:5)), 1);
  legend(strcat('P1/P1 org (slope=',num2str(P1(1)), ')'), ...
         strcat('P1/P1     (slope=',num2str(P2(1)), ')'), ...
         strcat('P1/P0     (slope=',num2str(P3(1)), ')'), ...
         strcat('P1/P2     (slope=',num2str(P4(1)), ')'), ...
         strcat('P2/P1     (slope=',num2str(P5(1)), ')'), ...
         strcat('P2/P0     (slope=',num2str(P6(1)), ')'), ...
         strcat('P2/P2     (slope=',num2str(P7(1)), ')'), ...
         'Location', 'NorthWest');
  grid on;
  axis([0.05 7 1e-3 10]);
  pause;
  GAMMA=[1    0.1  0.01   1e-3   1e-4   1e-5   1e-6   1e-7   1e-8   1e-9   1e-10  1e-11   1e-12   1e-13 ];
  H1_G =[7.5  1.45 0.2736 0.2717 0.2723 0.2724 0.2724 0.2724 0.2724 0.2724 0.2724 0.2724  0.2724  0.2724];
  condG=[1e6  2e5  6.45e4 6.58e4 6.59e4 6.59e4 2.3e5  2.3e6  2.3e7  2.3e8  2.22e9 2.25e10 2.25e11 2.22e12];

  loglog(GAMMA(1:14), H1_G(1:14), 'o-k');
  axis([1e-14 2 0.1 10]);
  xlabel('\gamma_0');
  ylabel('H^1(\Omega) error');
  grid on;
  pause;
  loglog(GAMMA(1:14), condG(1:14), 'o-k');
  axis([1e-14 2 1 3e22]);
  xlabel('\gamma_0');
  ylabel('Condition number');
  grid on;
  pause;
  H_3D  = [ 1    2    3   5   11  ];
  L2_3D_0=[ 3.53 12.8 35  167 1461 ];
  H1_3D_0=[ 5.60 14.8 37  168 1462 ];
  L2_3D_1=[ 4.33 13   40  173 1022 ];
  H1_3D_1=[ 6.20 15   42  174 1023 ];
  
  loglog(H_3D(1:4), L2_3D_0(1:4), 'o-k');
  hold on;
  loglog(H_3D(1:4), L2_3D_1(1:4), 'x:k');
  hold off;
  xlabel('h');
  ylabel('L^2(\Omega) error');
  P1 = polyfit(log(H_3D(1:4)), log(L2_3D_0(1:4)), 1);
  P2 = polyfit(log(H_3D(1:4)), log(L2_3D_1(1:4)), 1);
  legend(strcat('P1/P0 (slope=',num2str(P1(1)), ')'), ...
	 strcat('P1/P1 (slope=',num2str(P2(1)), ')'));
  axis([0.5 20 1 300]);
  grid on;
  pause;
  loglog(H_3D(1:4), H1_3D_0(1:4), 'o-k');
  hold on;
  loglog(H_3D(1:4), H1_3D_1(1:4), 'x:k');
  hold off;
  xlabel('h');
  ylabel('H^1(\Omega) error');
  P1 = polyfit(log(H_3D(1:3)), log(H1_3D_0(1:3)), 1);
  P2 = polyfit(log(H_3D(1:3)), log(H1_3D_1(1:3)), 1);
  legend(strcat('P1/P0 (slope=',num2str(P1(1)), ')'), ...
	 strcat('P1/P1 (slope=',num2str(P2(1)), ')'));
  axis([0.5 20 1 300]);
  grid on;
end;

% gf_plot(mf_vm, VM, 'norm', 'on', 'refine', 1, 'deformation', U, 'deformation_mf', mf, 'deformed_mesh','on', 'deformation_scale', 1.0);
%caxis([-0.04 0.22]);
% axis([-25 25 0 35]);
%colorbar;



% gf_plot_mesh(mesh, 'edges_width', 1, 'edges_color', [0 0 0], 'curved', 'on', 'refine', 4); % pour le dessin du maillage
