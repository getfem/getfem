% Copyright (C) 2008-2012 Yves Renard, Julien Pommier.
%
% This file is a part of GETFEM++
%
% Getfem++  is  free software;  you  can  redistribute  it  and/or modify it
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

% addpath ~/source++/getfem++/contrib/static_friction/

gf_workspace('clear all');

option = 3; % 0 : reference solution
            % 1 : error
            % 2 : solution
            % 3 : contact stress
            % 4 : some convergence curves 

if (option == 1)
  mesh = gf_mesh('load', 'reference_sol2d.meshfem');
  mf = gf_mesh_fem('load', 'reference_sol2d.meshfem', mesh);
  U = load('reference_sol2d.U')';
  UERR = load('reference_sol2d_error.U')';
  % UERR = log(abs(UERR));
  gf_plot(mf, UERR, 'norm', 'on', 'refine', 1, 'deformation', U, 'deformation_mf', mf, 'deformed_mesh','off', 'deformation_scale', 1.0);
  colorbar;
  gf_colormap('chouette');
  A = colormap; colormap(A(7:size(A,1),:));
elseif (option == 2 || option == 0)
  if (option == 0)
    mesh = gf_mesh('load', 'reference_sol2d.meshfem');
    mf = gf_mesh_fem('load', 'reference_sol2d.meshfem', mesh);
    mf_vm = gf_mesh_fem('load', 'reference_sol2d.meshfem_vm', mesh);
    U = load('reference_sol2d.U')';
    VM = load('reference_sol2d.VM')';
  else 
    mesh = gf_mesh('load', 'static_friction.meshfem');
    mf = gf_mesh_fem('load', 'static_friction.meshfem', mesh);
    mf_vm = gf_mesh_fem('load', 'static_friction.meshfem_vm', mesh);
    U = load('static_friction.U')';
    VM = load('static_friction.VM')';
  end;
  N = gf_mesh_get(mesh, 'dim');
  if (N==2) 
    gf_plot(mf_vm, VM, 'norm', 'on', 'refine', 1, 'deformation', U, ...
      'deformation_mf', mf, 'deformed_mesh','off', 'deformation_scale', 1.0);
    gf_colormap('chouette');
    A = colormap; colormap(A(7:size(A,1),:));
    xlabel('x'); ylabel('y');
    colorbar;
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
    h=gf_plot_slice(sl2,'mesh','on','mesh_slice_edges','off','data',VMsl);
    view(-80,-15); axis on;
    camlight;
    gf_colormap('chouette');
    % map=[1:-1/10:0]'*[1 1 1]; colormap(map); % for NB
    xlabel('x'); ylabel('y'); zlabel('z');  colorbar;
    pause;
    h=gf_plot_slice(sl1,'mesh_faces','off','mesh','on'); view(-85,-15);
    axis on; camlight; set(h,'facecolor',[.8 0 0]);
    xlabel('x'); ylabel('y'); zlabel('z');  colorbar;
    pause;
    gf_plot(mf_vm, VM, 'norm', 'on', 'refine', 4, 'deformation', U, ...
      'deformation_mf', mf, 'deformed_mesh','on', 'deformation_scale', ...
      1.0, 'cvlst', gf_mesh_get(mesh, 'outer faces'));
    view(-5,-10); camlight; colormap(map); xlabel('x'); ylabel('y'); zlabel('z');
    xlabel('x'); ylabel('y'); zlabel('z');
  end;
elseif (option == 3)
  mesh = gf_mesh('load', 'reference_sol2d.meshfem');
  sll=gfSlice('load','reference_sol2d.sl');
  LN=load('reference_sol2d.LN')';
  % mesh = gf_mesh('load', 'static_friction.meshfem');
  % sll=gfSlice('load','static_friction.sl');
  % LN=load('static_friction.LN')';
  N = gf_mesh_get(mesh, 'dim');
  P0=gf_slice_get(sll, 'pts');
  % [h1,h2,h3,h4]=gf_plot_slice(sll, 'tube','off', ...
  %                            'mesh_slice_edges_color',[.3 .3 .3]);
  hold on;
  gf_slice_set(sll,'pts',[P0 ; LN(N:N:size(LN,2))]);
  [hh1,hh2,hh3,hh4]=gf_plot_slice(sll, 'tube','off', ...
        'mesh_slice_edges_color','black','mesh_slice_edges_width',1.5);
  % sl=gfSlice('load','xfem_dirichlet.sl');
  % gf_plot_mesh(mesh, 'edges_width', 1, 'curved', 'off', 'refine', 1, 'edges_color',[0.6 1 0.6] );
  
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
  axis([-21 21 0.1 41 -5 0]);
  % disp('saving figure ...');
  % print(gcf,'-dpng','-r300', 'titi.png');
else
  H    = [0.075  0.15    0.25   0.5    1      2      3      4     5.5];
  % P1/P1 Non Augmenté
  L2_1 = [0.0020 0.0082  0.039  0.086  0.30   1.68   3.75   9.09] / 1.41;
  H1_1 = [0.042  0.072   0.136  0.284  0.64   2.04   4.20   9.50] / 1.41;
  L2C_1= [0.167  0.253   0.408  0.686  1.80   2.17   4.23   5.55] / 0.168;
  % P1/P1 augmenté
  L2_2 = [0.0021 0.0079  0.039  0.087  0.31   1.68   3.75   9.09] / 1.41;
  H1_2 = [0.042  0.072   0.136  0.285  0.64   2.04   4.20   9.50] / 1.41;
  L2C_2= [0.165  0.250   0.403  0.686  1.78   2.15   4.22   5.55] / 0.168;
  % P1/P0 augmenté
  L2_3 = [0.0027 0.010   0.021  0.092  0.45   2.54   7.28] / 1.41;
  H1_3 = [0.0418 0.072   0.132  0.29   0.71   2.79   7.55] / 1.41;
  L2C_3= [0.8    0.703   1.6   1.27   3.53   7.05   9.91] / 0.168;
  % P1/P2 augmenté
  L2_4 = [0.003  0.0090  0.022  0.087  0.45   1.40   3.75   8.30   36.8]/1.41;
  H1_4 = [0.0419 0.072   0.13   0.28   0.72   1.78   4.15   8.73   36.9]/1.41;
  L2C_4= [0.1    0.300   0.7   0.447   2.24   3.89   6.19   4.66   17.9]/0.168;
  % P2/P1 augmenté
  L2_5 = [100    0.0045  0.0082 0.050  0.12   0.88   2.04   11.9] / 1.41;
  H1_5 = [100    0.0074  0.014  0.092  0.15   0.92   2.11   12.03] / 1.41;
  L2C_5= [100    0.112   0.17   0.24   0.38   0.99   1.54   5.55] / 0.168;
  % P2/P0 augmenté
  L2_6 = [100    0.0039  0.013  0.048  0.18   0.59   1.74   6.92] / 1.41;
  H1_6 = [100    0.0097  0.021  0.059  0.22   0.71   1.91   7.01] / 1.41;
  L2C_6= [100    0.314   0.499  0.698  1.72   3.47   5.27   4.42] / 0.168;
  % P2/P2 augmenté
  L2_7 = [100    0.00042 0.0022 0.010  0.048  0.13   0.58   4.06   4.89]/1.41;
  H1_7 = [100    0.0052  0.010  0.047  0.10   0.30   0.79   4.13   5.15]/1.41;
  L2C_7= [100    0.086   0.119 0.224  0.490  1.314  2.23   1.87   5.57]/0.168;

  loglog(H(1:7), L2_1(1:7), 'o-k', 'linewidth', 2, 'MarkerSize', 15);
  hold on;
  loglog(H(1:7), L2_2(1:7), 'x-.k', 'linewidth', 2, 'MarkerSize', 15);
  loglog(H(1:7), L2_3(1:7), '+--k', 'linewidth', 2, 'MarkerSize', 15);
  loglog(H(1:7), L2_4(1:7), '*-k', 'linewidth', 2, 'MarkerSize', 15);
  loglog(H(2:7), L2_5(2:7), 's-.k', 'linewidth', 2, 'MarkerSize', 15);
  loglog(H(2:7), L2_6(2:7), 'd--k', 'linewidth', 2, 'MarkerSize', 15);
  loglog(H(2:7), L2_7(2:7), '<-k', 'linewidth', 2, 'MarkerSize', 15);
  hold off;
  xlabel('h');
  ylabel('L^2(\Omega) relative error in %');
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
  loglog(H(1:7), L2C_1(1:7), 'o-k', 'linewidth', 2, 'MarkerSize', 15);
  hold on;
  loglog(H(1:7), L2C_2(1:7), 'x-.k', 'linewidth', 2, 'MarkerSize', 15);
  loglog(H(1:7), L2C_3(1:7), '+--k', 'linewidth', 2, 'MarkerSize', 15);
  loglog(H(1:7), L2C_4(1:7), '*-k', 'linewidth', 2, 'MarkerSize', 15);
  loglog(H(2:7), L2C_5(2:7), 's-.k', 'linewidth', 2, 'MarkerSize', 15);
  loglog(H(2:7), L2C_6(2:7), 'd--k', 'linewidth', 2, 'MarkerSize', 15);
  loglog(H(2:7), L2C_7(2:7), '<-k', 'linewidth', 2, 'MarkerSize', 15);
  hold off;
  xlabel('h');
  ylabel('L^2(\Gamma_C) relative error in %');
  P1 = polyfit(log(H(1:7)), log(L2C_1(1:7)), 1);
  P2 = polyfit(log(H(1:7)), log(L2C_2(1:7)), 1);
  P3 = polyfit(log(H(1:7)), log(L2C_3(1:7)), 1);
  P4 = polyfit(log(H(1:7)), log(L2C_4(1:7)), 1);
  P5 = polyfit(log(H(2:7)), log(L2C_5(2:7)), 1);
  P6 = polyfit(log(H(2:7)), log(L2C_6(2:7)), 1);
  P7 = polyfit(log(H(2:7)), log(L2C_7(2:7)), 1);
  legend(strcat('P1/P1 org (slope=',num2str(P1(1)), ')'), ...
         strcat('P1/P1     (slope=',num2str(P2(1)), ')'), ...
         strcat('P1/P0     (slope=',num2str(P3(1)), ')'), ...
         strcat('P1/P2     (slope=',num2str(P4(1)), ')'), ...
         strcat('P2/P1     (slope=',num2str(P5(1)), ')'), ...
         strcat('P2/P0     (slope=',num2str(P6(1)), ')'), ...
         strcat('P2/P2     (slope=',num2str(P7(1)), ')'), ...
         'Location', 'SouthEast');
  grid on;
  axis([0.05 7 1e-1 100]);



  pause;
  loglog(H(1:7), H1_1(1:7), 'o-k', 'linewidth', 2, 'MarkerSize', 15);
  hold on;
  loglog(H(1:7), H1_2(1:7), 'x-.k', 'linewidth', 2, 'MarkerSize', 15);
  loglog(H(1:7), H1_3(1:7), '+--k', 'linewidth', 2, 'MarkerSize', 15);
  loglog(H(1:7), H1_4(1:7), '*-k', 'linewidth', 2, 'MarkerSize', 15);
  loglog(H(2:7), H1_5(2:7), 's-.k', 'linewidth', 2, 'MarkerSize', 15);
  loglog(H(2:7), H1_6(2:7), 'd--k', 'linewidth', 2, 'MarkerSize', 15);
  loglog(H(2:7), H1_7(2:7), '<-k', 'linewidth', 2, 'MarkerSize', 15);
  hold off;
  xlabel('h');
  ylabel('H^1(\Omega) relative error in %');
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
  H1_G =[7.5  1.45 0.2736 0.2717 0.2723 0.2724 0.2724 0.2724 0.2724 0.2724 0.2724 0.2724  0.2724  0.2724] / 1.41;
  condG=[1e6  2e5  6.45e4 6.58e4 6.59e4 6.59e4 2.3e5  2.3e6  2.3e7  2.3e8  2.22e9 2.25e10 2.25e11 2.22e12];

  loglog(GAMMA(1:14), H1_G(1:14), 'o-k', 'linewidth', 2, 'MarkerSize', 15);
  axis([1e-14 2 0.1 10]);
  xlabel('\gamma_0');
  ylabel('H^1(\Omega) relative error in %');
  grid on;
  pause;
  loglog(GAMMA(1:14), condG(1:14), 'o-k', 'linewidth', 2, 'MarkerSize', 15);
  axis([1e-14 2 1 3e22]);
  xlabel('\gamma_0');
  ylabel('Condition number');
  grid on;
  pause;
  H_3D  = [ 0.5  1    2    3   5   ];
  L2_3D_0=[ 0.91 3.52 12.8 35  142 ] / 9.97;
  H1_3D_0=[ 2.52 5.18 14.8 30  72  ] / 9.98;
  L2_3D_1=[ 0.98 4.36 13   40  162 ] / 9.97;
  H1_3D_1=[ 2.65 5.86 15   32  76  ] / 9.98;
  
  loglog(H_3D(1:5), L2_3D_0(1:5), 'o-k', 'linewidth', 2, 'MarkerSize', 15);
  hold on;
  loglog(H_3D(1:5), L2_3D_1(1:5), 'x:k', 'linewidth', 2, 'MarkerSize', 15);
  hold off;
  xlabel('h');
  ylabel('L^2(\Omega) relative error in %');
  P1 = polyfit(log(H_3D(1:5)), log(L2_3D_0(1:5)), 1);
  P2 = polyfit(log(H_3D(1:5)), log(L2_3D_1(1:5)), 1);
  legend(strcat('P1/P0 (slope=',num2str(P1(1)), ')'), ...
	 strcat('P1/P1 (slope=',num2str(P2(1)), ')'));
  axis([0.25 10 0.05 50]);
  grid on;
  pause;
  loglog(H_3D(1:5), H1_3D_0(1:5), 'o-k', 'linewidth', 2, 'MarkerSize', 15);
  hold on;
  loglog(H_3D(1:5), H1_3D_1(1:5), 'x:k', 'linewidth', 2, 'MarkerSize', 15);
  hold off;
  xlabel('h');
  ylabel('H^1(\Omega) relative error in %');
  P1 = polyfit(log(H_3D(1:3)), log(H1_3D_0(1:3)), 1);
  P2 = polyfit(log(H_3D(1:3)), log(H1_3D_1(1:3)), 1);
  legend(strcat('P1/P0 (slope=',num2str(P1(1)), ')'), ...
	 strcat('P1/P1 (slope=',num2str(P2(1)), ')'));
  axis([0.25 10 0.1 50]);
  grid on;
end;



% Pour mettre des fontes plus grosses.
% une commande
% get(findobj, 'type')
% renseigne sur les type d'objets à chercher.
% ensuite on recupère les handles par
% axesobj = findobj('type', 'axes')
% par exemple, puis on peut faire
% set(axesobj, 'fontunits', 'points');
% set(axesobj, 'fontsize', 15);
% set(axesobj, 'fontweight', 'bold');
% Il vaut mieux a la fin decouper les images avec gimp par exemple.

axesobj = findobj('type', 'axes');
set(axesobj, 'fontname', 'times');
set(axesobj, 'fontunits', 'points');
set(axesobj, 'fontsize', 15);
set(axesobj, 'fontweight', 'bold');


% Pour certains graphiques, il vaut mieux renommer les "ticks" par
%  set(gca,'XTickLabel',{'0.1';'1';'10';'...'})
%  set(gca,'YTickLabel',{'0.0001%';'0.001%';'0.01%';'0.1%';'1%';'10%'})     


% Pour sortir le graphique en png, faire par exemple :
% print(gcf,'-dpng','-r450', 'toto.png');
