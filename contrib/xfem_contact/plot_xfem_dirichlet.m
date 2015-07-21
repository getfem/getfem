% Copyright (C) 2008-2015 Yves Renard, Julien Pommier.
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


% addpath ~/source++/getfem++/contrib/xfem_contact/

gf_workspace('clear all');
mf = gf_mesh_fem('load', 'xfem_dirichlet_ls.mf');
lsU = load('xfem_dirichlet_ls.U')';
lsU1 = load('xfem_dirichlet_exact.U')';
mf1 = gf_mesh_fem('load', 'xfem_dirichlet.mfE');
nn=4;  % 0 : plot the exported mf
       % 1 : 
       % 2 : 
       % 3 : 
       % 4 : plotting the lagrange multipliers on the dirichlet boundary
       % 5 : the solution. 
       % 6 : plot some convergence curves. 


clf
if nn==0,
  disp('plot the exported mf');
  [hsur, hcont] = gf_plot(mf,lsU,'refine',2,'contour',0.,'mesh','on', 'pcolor','off');
  set(hcont{1}, 'LineWidth', 3);
  set(hcont{1}, 'Color', 'black');
elseif nn==1
  disp('plot the cut mesh');
  mc = gfMesh('load','cut.mesh');
  mfc = gfMeshFem(mc); 
  set(mfc,'classical_fem',2);
  
  lsUc = gf_compute(mf, lsU, 'interpolate on', mfc);
  
  [hsur, hcont] = gf_plot(mf, lsU, 'refine', 2, 'zplot', 'on');
  hold on;
  [hsur, hcont] = gf_plot(mfc, lsUc.*(lsUc<0), 'refine', 1, 'mesh','on', 'pcolor','on','zplot', 'on'); hold on;
  %colormap([.8 1 .8]);
  [hsur, hcont] = gf_plot(mf,lsU,'refine',1,'mesh','on','zplot', 'on','contour',0.,'pcolor','off','zplot', 'on');
  
  %set(hcont{1}, 'LineWidth', 2);
  %set(hcont{1}, 'Color', 'read');
  %axis('tight'); axis off;
elseif nn==2 || nn==3,
  disp('plot the solution, with the 0 isovalue');
  sl=gfSlice('load','xfem_dirichlet.sl');
  slU=-load('xfem_dirichlet.slU')';
  P=gf_slice_get(sl,'pts'); P=[P(1:2,:);slU];
  gf_slice_set(sl,'pts',P);
  gf_plot_slice(sl, 'data', slU, 'mesh','on','mesh_edges','off');
  
  

  m=gf_mesh_fem_get(mf, 'linked_mesh');
  slc=gf_Slice({'isovalues', 0, mf, lsU, 0}, m, 16);
  hold on;
  P2=gf_slice_get(slc, 'pts');
  gf_slice_set(slc, 'pts', [P2;0.1 * ones(1,size(P2,2))]);
  [h1,h2,h3,h4]=gf_plot_slice(slc, 'tube','off','mesh_slice_edges_color','black');
  set(h4, 'LineWidth', 2);
  
  if (nn == 2),
%    set(hcont{1}, 'Color', 'black');
    view(3); camlight; axis off;camzoom(1.8);
  else
    slc2=gfSlice('load', 'xfem_dirichlet.sl0');
    hold on;
    set(gcf,'renderer','zbuffer');
    [h1,h2,h3,h4]=gf_plot_slice(slc2, 'tube','off','mesh_slice_edges_color','white');
    set(h4, 'LineWidth', 4);
    view(3); 
    %caxis([-.2 .3]); gf_colormap('froid');
  end;
elseif nn==4,
  disp('plotting the lagrange multipliers on the dirichlet boundary');
  sll=gf_Slice('load','xfem_dirichlet.sll');
  slL=load('xfem_dirichlet.slL')';
  P0=gf_slice_get(sll, 'pts');
  [h1,h2,h3,h4]=gf_plot_slice(sll, 'tube','off','mesh_slice_edges_color','black');
  hold on;
  gf_slice_set(sll,'pts',[P0 ; max(slL,-100)*0.05]);
  [hh1,hh2,hh3,hh4]=gf_plot_slice(sll, 'tube','off','mesh_slice_edges_color','black','mesh_slice_edges_width',2,'showoptions','on');
  sl=gf_Slice('load','xfem_dirichlet.sl');
  gf_plot_slice(sl,'mesh','on');
  
  npt = size(P0, 2);
  P0 = [P0;zeros(1,npt)];
  P1 = gf_slice_get(sll,'pts'); 
  lseg = gf_slice_get(sll,'splxs', 1);
  F=[lseg(1,:) lseg(2,:); lseg(2,:) npt+lseg(2,:); npt+lseg(1,:) npt+lseg(1,:)];
  %F=[lseg; npt+lseg(2,:)];
  h=patch('Vertices',[P0 P1]', 'Faces', F');
  hold on;
  set(h,'FaceAlpha',0.3);
  set(h,'LineStyle','none');
  set(gcf,'renderer','zbuffer');
  set(gcf,'color','white');
  set(h,'facecolor',[.5 .5 .5]);
  axis off;
  view(3);
  camzoom(2.5);
  %axis([-0.5000    0.5000   -0.5000    0.5000 -.5 .5]);
  % print(gcf,'-dpng','-r300', 'lagrange_multipliers.png');
elseif nn==5,
  disp('plot the solution');
  lsU = max(-lsU, 0) *10;
  [hsur, hcont] = gf_plot(mf, lsU, 'refine',2,'mesh','on','zplot', 'on');
  %colormap([0.7 0.8 1.0; 0.8 1 0.8; 0.8 1 0.8; 0.8 1 0.8; 0.8 1 0.8; 0.8 1 0.8; 0.8 1 0.8; 0.8 1 0.8; 0.8 1 0.8; 0.8 1 0.8; 0.8 1 0.8; 0.8 1 0.8; 0.8 1 0.8; 0.8 1 0.8; 0.8 1 0.8; 0.8 1 0.8; 0.8 1 0.8; 0.8 1 0.8; 0.8 1 0.8; 0.8 1 0.8; 0.8 1 0.8; 0.8 1 0.8]);
  axis off;
  %camlight;
  % set(hcont{1}, 'LineWidth', 2);
  % set(hcont{1}, 'Color', 'black');
elseif nn==6,%plot multiplier
  slL=load('xfem_dirichlet.slL')';
  sll=gf_Slice('load','xfem_dirichlet.sll');
  gf_plot_slice(sll,'mesh','on','mesh_slice_edges_color','black','data',slL,'showoptions','on');
  slL=load('xfem_dirichlet.slL')';
  %sll=gf_Slice('load','xfem_dirichlet.sl0');
  %slL=load('xfem_dirichlet.slU0')';
  %gf_plot_slice(sll,'mesh','on','mesh_slice_edges_color','black','data',slL,'showoptions','on');
  axis on;
elseif nn==7,%plot displacement at the bondary
  slU=load('xfem_dirichlet.slU')';  
  sll=gf_Slice('load','xfem_dirichlet.sl');
  gf_plot_slice(sll, 'mesh','on','mesh_slice_edges_color','black','data',slU,'showoptions','on');
  axis on;
  
elseif nn==8,%plot displacement at the half sphere 
 % sll=gf_Slice('load','xfem_dirichlet.sl0');
 
 mfs = gf_mesh_fem('load', 'xfem_dirichlet.mf');
lsUs = load('xfem_dirichlet_ls.U')';

 sl=gf_slice({'boundary',{'intersection',{'ball',-1,[0;0],0.4}}},mfs,9);
 Usl=gf_compute(mfs,lsUs,'interpolate on', sl);
 gf_plot_slice(sl,'mesh_faces','on','mesh','on','data',Usl,'mesh_slice_edges','on');

elseif nn==9,

  % Without stabilization, FEM_RHS = 'FEM_PK(2,3)'; LEVEL_SET_DEGREE = 2;
  H    = [1/320   1/160  1/80   1/40  1/20  1/10  1/5];
  % P1/P0 Non stabilis�
  L2_1 = [2.74    7.4    11.7   26    48.9  70    82 ];
  H1_1 = [15.28   24     31     44    59.5  77    89 ];
  L2C_1= [41700   31300  18300  9560  5000  6620  323];
  % P1+/P0 Non stabilis�
  L2_2 = [0.021   0.083  0.31   1.0   5.5   14    29];
  H1_2 = [1.73    3.5    6.89   13.4  28    48    64];
  L2C_2= [435     1320   659    384   868   230   92];
  % P2/P1 Non stabilis�
  L2_3 = [0.0011  0.0016 0.0049 0.038 0.43  4     12];
  H1_3 = [0.017   0.05   0.204  0.773 2.9   11    32];
  L2C_3= [0.23    0.89   0.798  1.78  4.5   10    38];
  % Q1/Q0 Non stabilis�
  L2_4 = [0.018   0.066  0.24   0.95  3.5   10    38];
  H1_4 = [1.55    3.1    6.12   12    23    43    70];
  L2C_4= [12.42   10.5   3.70   6.36  13.7  45    44];
  % Q2/Q1 Non stabilis�
  L2_5 = [0.0059  0.013  0.031  0.11  0.41  5     11];
  H1_5 = [0.037   0.11   0.19   0.68  2.10  7.3   24];
  L2C_5= [1.07    1.03   1.06   2.79  4.3   5.9   21];
  

  loglog(H(1:7), L2_1(1:7), 'o-k', 'linewidth', 2, 'MarkerSize', 15);
  hold on;
  loglog(H(1:7), L2_2(1:7), 'x-.k', 'linewidth', 2, 'MarkerSize', 15);
  loglog(H(1:7), L2_3(1:7), '+--k', 'linewidth', 2, 'MarkerSize', 15);
  loglog(H(1:7), L2_4(1:7), '*-k', 'linewidth', 2, 'MarkerSize', 15);
  loglog(H(1:7), L2_5(1:7), 's-.k', 'linewidth', 2, 'MarkerSize', 15);
  hold off;
  P1 = polyfit(log(H(1:5)), log(L2_1(1:5)), 1);
  P2 = polyfit(log(H(1:5)), log(L2_2(1:5)), 1);
  P3 = polyfit(log(H(1:5)), log(L2_3(1:5)), 1);
  P4 = polyfit(log(H(1:5)), log(L2_4(1:5)), 1);
  P5 = polyfit(log(H(2:5)), log(L2_5(2:5)), 1);
  legend(strcat('P1/P0  (slope=',num2str(P1(1)), ')'), ...
         strcat('P1+/P0 (slope=',num2str(P2(1)), ')'), ...
         strcat('P2/P1  (slope=',num2str(P3(1)), ')'), ...
         strcat('Q1/Q0  (slope=',num2str(P4(1)), ')'), ...
         strcat('Q2/Q1  (slope=',num2str(P5(1)), ')'), ...
         'Location', 'NorthWest');
  grid on;
  axesobj = findobj('type', 'axes');
  set(axesobj, 'fontname', 'times'); set(axesobj, 'fontunits', 'points');
  set(axesobj, 'fontsize', 18); set(axesobj, 'fontweight', 'bold');
  set(axesobj, 'linewidth', 4);
  xlabel('h');
  ylabel('L^2(\Omega) relative error (in %)');
  set(gca,'XTickLabel',{'0.001';'0.01';'0.1';'1';'...'}) 
  % axis([0.05 7 1e-4 10]);
  pause;

  loglog(H(1:7), H1_1(1:7), 'o-k', 'linewidth', 2, 'MarkerSize', 15);
  hold on;
  loglog(H(1:7), H1_2(1:7), 'x-.k', 'linewidth', 2, 'MarkerSize', 15);
  loglog(H(1:7), H1_3(1:7), '+--k', 'linewidth', 2, 'MarkerSize', 15);
  loglog(H(1:7), H1_4(1:7), '*-k', 'linewidth', 2, 'MarkerSize', 15);
  loglog(H(1:7), H1_5(1:7), 's-.k', 'linewidth', 2, 'MarkerSize', 15);
  hold off;
  P1 = polyfit(log(H(1:5)), log(H1_1(1:5)), 1);
  P2 = polyfit(log(H(1:5)), log(H1_2(1:5)), 1);
  P3 = polyfit(log(H(1:5)), log(H1_3(1:5)), 1);
  P4 = polyfit(log(H(1:5)), log(H1_4(1:5)), 1);
  P5 = polyfit(log(H(1:5)), log(H1_5(1:5)), 1);
  legend(strcat('P1/P0  (slope=',num2str(P1(1)), ')'), ...
         strcat('P1+/P0 (slope=',num2str(P2(1)), ')'), ...
         strcat('P2/P1  (slope=',num2str(P3(1)), ')'), ...
         strcat('Q1/Q0  (slope=',num2str(P4(1)), ')'), ...
         strcat('Q2/Q1  (slope=',num2str(P5(1)), ')'), ...
         'Location', 'NorthWest');
  grid on;
  axesobj = findobj('type', 'axes');
  set(axesobj, 'fontname', 'times'); set(axesobj, 'fontunits', 'points');
  set(axesobj, 'fontsize', 18); set(axesobj, 'fontweight', 'bold');
  set(axesobj, 'linewidth', 2);
  xlabel('h');
  ylabel('H^1(\Omega) relative error (in %)');
  set(gca,'XTickLabel',{'0.001';'0.01';'0.1';'1';'...'}) 
  % axis([0.05 7 1e-4 10]);
  pause;


  loglog(H(1:7), L2C_1(1:7), 'o-k', 'linewidth', 2, 'MarkerSize', 15);
  hold on;
  loglog(H(1:7), L2C_2(1:7), 'x-.k', 'linewidth', 2, 'MarkerSize', 15);
  loglog(H(1:7), L2C_3(1:7), '+--k', 'linewidth', 2, 'MarkerSize', 15);
  loglog(H(1:7), L2C_4(1:7), '*-k', 'linewidth', 2, 'MarkerSize', 15);
  loglog(H(1:7), L2C_5(1:7), 's-.k', 'linewidth', 2, 'MarkerSize', 15);
  hold off;
  P1 = polyfit(log(H(1:5)), log(L2C_1(1:5)), 1);
  P2 = polyfit(log(H(1:5)), log(L2C_2(1:5)), 1);
  P3 = polyfit(log(H(1:5)), log(L2C_3(1:5)), 1);
  P4 = polyfit(log(H(1:5)), log(L2C_4(1:5)), 1);
  P5 = polyfit(log(H(1:5)), log(L2C_5(1:5)), 1);
  legend(strcat('P1/P0  (slope=',num2str(P1(1)), ')'), ...
         strcat('P1+/P0 (slope=',num2str(P2(1)), ')'), ...
         strcat('P2/P1  (slope=',num2str(P3(1)), ')'), ...
         strcat('Q1/Q0  (slope=',num2str(P4(1)), ')'), ...
         strcat('Q2/Q1  (slope=',num2str(P5(1)), ')'), ...
         'Location', 'NorthWest');
  grid on;
  axesobj = findobj('type', 'axes');
  set(axesobj, 'fontname', 'times'); set(axesobj, 'fontunits', 'points');
  set(axesobj, 'fontsize', 18); set(axesobj, 'fontweight', 'bold');
  set(axesobj, 'linewidth', 2);
  xlabel('h');
  ylabel('L^2(\Gamma_D) relative error (in %)');
  set(gca,'XTickLabel',{'0.001';'0.01';'0.1';'1';'...'}) 
  % axis([0.05 7 1e-4 10]);
  pause;


  % With BB stabilization, gamma0 = 0.1, FEM_RHS = 'FEM_PK(2,3)';
  % LEVEL_SET_DEGREE = 2;
  H    = [1/320     1/160   1/80   1/40  1/20  1/10  1/5];
  % P1/P0 stabilis�
  L2_1 = [0.022     0.086   0.34   1.3   5.2   16    41];
  H1_1 = [2.0       3.7     7.3    14    28    50    74];
  L2C_1= [5.8       4.2     9      16    35    59    57];
  % P1+/P0 stabilis�
  L2_2 = [0.02      0.078   0.30   1.14  4.5   19    51];
  H1_2 = [1.73      3.4     6.7    13    26    80    62];
  L2C_2= [1.81      4.14    5.4    9     28    60    38];
  % P2/P1 stabilis�
  L2_3 = [0.000062  0.0005  0.0037 0.033 0.34  2.44  9.9];
  H1_3 = [0.012     0.054   0.58   1.12  3.26  10.5  33];
  L2C_3= [0.017     0.061   1      1.03  1.93  6.3   19];
  % Q1/Q0 stabilis�
  L2_4 = [0.014     0.05    0.22   0.87  3.23  9.9   23];
  H1_4 = [1.65      3.10    6.14   12.08 23.43 44    69];
  L2C_4= [0.9       1.63    2.82   6.55  15.19 34    33];
  % Q2/Q1 stabilis�
  L2_5 = [0.000035  0.00033 0.0029 0.018 0.13  1.1   17];
  H1_5 = [0.0093    0.051   0.15   0.58  1.9   6.9   65];
  L2C_5= [0.019     0.079   0.17   0.38  0.8   4.2   17];


  loglog(H(1:7), L2_1(1:7), 'o-k', 'linewidth', 2, 'MarkerSize', 15);
  hold on;
  loglog(H(1:7), L2_2(1:7), 'x-.k', 'linewidth', 2, 'MarkerSize', 15);
  loglog(H(1:7), L2_3(1:7), '+--k', 'linewidth', 2, 'MarkerSize', 15);
  loglog(H(1:7), L2_4(1:7), '*-k', 'linewidth', 2, 'MarkerSize', 15);
  loglog(H(1:7), L2_5(1:7), 's-.k', 'linewidth', 2, 'MarkerSize', 15);
  hold off;
  P1 = polyfit(log(H(1:5)), log(L2_1(1:5)), 1);
  P2 = polyfit(log(H(1:5)), log(L2_2(1:5)), 1);
  P3 = polyfit(log(H(1:5)), log(L2_3(1:5)), 1);
  P4 = polyfit(log(H(1:5)), log(L2_4(1:5)), 1);
  P5 = polyfit(log(H(2:5)), log(L2_5(2:5)), 1);
  legend(strcat('P1/P0  (slope=',num2str(P1(1)), ')'), ...
         strcat('P1+/P0 (slope=',num2str(P2(1)), ')'), ...
         strcat('P2/P1  (slope=',num2str(P3(1)), ')'), ...
         strcat('Q1/Q0  (slope=',num2str(P4(1)), ')'), ...
         strcat('Q2/Q1  (slope=',num2str(P5(1)), ')'), ...
         'Location', 'NorthWest');
  grid on;
  axesobj = findobj('type', 'axes');
  set(axesobj, 'fontname', 'times'); set(axesobj, 'fontunits', 'points');
  set(axesobj, 'fontsize', 18); set(axesobj, 'fontweight', 'bold');
  set(axesobj, 'linewidth', 2);
  xlabel('h');
  ylabel('L^2(\Omega) relative error (in %)');
  set(gca,'XTickLabel',{'0.001';'0.01';'0.1';'1';'...'}) 
  % axis([0.05 7 1e-4 10]);
  pause;

  loglog(H(1:7), H1_1(1:7), 'o-k', 'linewidth', 2, 'MarkerSize', 15);
  hold on;
  loglog(H(1:7), H1_2(1:7), 'x-.k', 'linewidth', 2, 'MarkerSize', 15);
  loglog(H(1:7), H1_3(1:7), '+--k', 'linewidth', 2, 'MarkerSize', 15);
  loglog(H(1:7), H1_4(1:7), '*-k', 'linewidth', 2, 'MarkerSize', 15);
  loglog(H(1:7), H1_5(1:7), 's-.k', 'linewidth', 2, 'MarkerSize', 15);
  hold off;
  P1 = polyfit(log(H(1:5)), log(H1_1(1:5)), 1);
  P2 = polyfit(log(H(1:5)), log(H1_2(1:5)), 1);
  P3 = polyfit(log(H(1:5)), log(H1_3(1:5)), 1);
  P4 = polyfit(log(H(1:5)), log(H1_4(1:5)), 1);
  P5 = polyfit(log(H(1:5)), log(H1_5(1:5)), 1);
  legend(strcat('P1/P0  (slope=',num2str(P1(1)), ')'), ...
         strcat('P1+/P0 (slope=',num2str(P2(1)), ')'), ...
         strcat('P2/P1  (slope=',num2str(P3(1)), ')'), ...
         strcat('Q1/Q0  (slope=',num2str(P4(1)), ')'), ...
         strcat('Q2/Q1  (slope=',num2str(P5(1)), ')'), ...
         'Location', 'NorthWest');
  grid on;
  axesobj = findobj('type', 'axes');
  set(axesobj, 'fontname', 'times'); set(axesobj, 'fontunits', 'points');
  set(axesobj, 'fontsize', 18); set(axesobj, 'fontweight', 'bold');
  set(axesobj, 'linewidth', 2);
  xlabel('h');
  ylabel('H^1(\Omega) relative error (in %)');
  set(gca,'XTickLabel',{'0.001';'0.01';'0.1';'1';'...'}) 
  % axis([0.05 7 1e-4 10]);
  pause;


  loglog(H(1:7), L2C_1(1:7), 'o-k', 'linewidth', 2, 'MarkerSize', 15);
  hold on;
  loglog(H(1:7), L2C_2(1:7), 'x-.k', 'linewidth', 2, 'MarkerSize', 15);
  loglog(H(1:7), L2C_3(1:7), '+--k', 'linewidth', 2, 'MarkerSize', 15);
  loglog(H(1:7), L2C_4(1:7), '*-k', 'linewidth', 2, 'MarkerSize', 15);
  loglog(H(1:7), L2C_5(1:7), 's-.k', 'linewidth', 2, 'MarkerSize', 15);
  hold off;
  P1 = polyfit(log(H(1:5)), log(L2C_1(1:5)), 1);
  P2 = polyfit(log(H(1:5)), log(L2C_2(1:5)), 1);
  P3 = polyfit(log(H(1:5)), log(L2C_3(1:5)), 1);
  P4 = polyfit(log(H(1:5)), log(L2C_4(1:5)), 1);
  P5 = polyfit(log(H(1:5)), log(L2C_5(1:5)), 1);
  legend(strcat('P1/P0  (slope=',num2str(P1(1)), ')'), ...
         strcat('P1+/P0 (slope=',num2str(P2(1)), ')'), ...
         strcat('P2/P1  (slope=',num2str(P3(1)), ')'), ...
         strcat('Q1/Q0  (slope=',num2str(P4(1)), ')'), ...
         strcat('Q2/Q1  (slope=',num2str(P5(1)), ')'), ...
         'Location', 'NorthWest');
  grid on;
  axesobj = findobj('type', 'axes');
  set(axesobj, 'fontname', 'times'); set(axesobj, 'fontunits', 'points');
  set(axesobj, 'fontsize', 18); set(axesobj, 'fontweight', 'bold');
  set(axesobj, 'linewidth', 2);
  xlabel('h');
  ylabel('L^2(\Gamma_D) relative error (in %)');
  set(gca,'XTickLabel',{'0.001';'0.01';'0.1';'1';'...'}) 
  % axis([0.05 7 1e-4 10]);
  pause;



  % With BB stabilization and stabilized normal derivative, gamma0 = 0.1,
  % FEM_RHS = 'FEM_PK(2,3)'; LEVEL_SET_DEGREE = 2; MINIMAL_ELT_RATIO = 0.01;
  H    = [1/320     1/160    1/80   1/40  1/20  1/10  1/5];
  % P1/P0 fully stabilised
  L2_1 = [0.022     0.086    0.34   1.31  5.2   16    41];
  H1_1 = [2.27      3.69     7.3    14    29    51    74];
  L2C_1= [5.5       4.15     10     16    44    65    57];   
  % P1+/P0 fully stabilised
  L2_2 = [0.02      0.078    0.30   1.14  4.55  13    24];
  H1_2 = [1.88      3.4      6.7    13.2  25.6  46    61];
  L2C_2= [1.81      2.6      5.3    8.8   27.9  45    35];
  % P2/P1 fully stabilised
  L2_3 = [0.000063  0.0005   0.004  0.032 0.34  2.44  9.9];
  H1_3 = [0.012     0.054    0.6    1.11  3.26  10.5  32];
  L2C_3= [0.046     0.056    1.1    0.92  1.9   5.8   18];
  % Q1/Q0 fully stabilised
  L2_4 = [0.014     0.057    0.22   0.87  3.22  9.86  23];
  H1_4 = [1.6       3.1      6.14   12.1  24    44.4  70];
  L2C_4= [0.87      1.37     2.45   5.8   11.7  21.1  33];
  % Q2/Q1 fully stabilised
  L2_5 = [0.000035  0.000033 0.0033 0.019 0.13  1.09  16];
  H1_5 = [0.0096    0.051    0.16   0.58  1.92  6.63  64];
  L2C_5= [0.019     0.079    0.18   0.38  0.79  4.04  17];


  loglog(H(1:7), L2_1(1:7), 'o-k', 'linewidth', 2, 'MarkerSize', 15);
  hold on;
  loglog(H(1:7), L2_2(1:7), 'x-.k', 'linewidth', 2, 'MarkerSize', 15);
  loglog(H(1:7), L2_3(1:7), '+--k', 'linewidth', 2, 'MarkerSize', 15);
  loglog(H(1:7), L2_4(1:7), '*-k', 'linewidth', 2, 'MarkerSize', 15);
  loglog(H(1:7), L2_5(1:7), 's-.k', 'linewidth', 2, 'MarkerSize', 15);
  hold off;
  P1 = polyfit(log(H(1:5)), log(L2_1(1:5)), 1);
  P2 = polyfit(log(H(1:5)), log(L2_2(1:5)), 1);
  P3 = polyfit(log(H(1:5)), log(L2_3(1:5)), 1);
  P4 = polyfit(log(H(1:5)), log(L2_4(1:5)), 1);
  P5 = polyfit(log(H(2:5)), log(L2_5(2:5)), 1);
  legend(strcat('P1/P0  (slope=',num2str(P1(1)), ')'), ...
         strcat('P1+/P0 (slope=',num2str(P2(1)), ')'), ...
         strcat('P2/P1  (slope=',num2str(P3(1)), ')'), ...
         strcat('Q1/Q0  (slope=',num2str(P4(1)), ')'), ...
         strcat('Q2/Q1  (slope=',num2str(P5(1)), ')'), ...
         'Location', 'NorthWest');
  grid on;
  axesobj = findobj('type', 'axes');
  set(axesobj, 'fontname', 'times'); set(axesobj, 'fontunits', 'points');
  set(axesobj, 'fontsize', 18); set(axesobj, 'fontweight', 'bold');
  set(axesobj, 'linewidth', 2);
  xlabel('h');
  ylabel('L^2(\Omega) relative error (in %)');
  set(gca,'XTickLabel',{'0.001';'0.01';'0.1';'1';'...'}) 
  % axis([0.05 7 1e-4 10]);
  pause;

  loglog(H(1:7), H1_1(1:7), 'o-k', 'linewidth', 2, 'MarkerSize', 15);
  hold on;
  loglog(H(1:7), H1_2(1:7), 'x-.k', 'linewidth', 2, 'MarkerSize', 15);
  loglog(H(1:7), H1_3(1:7), '+--k', 'linewidth', 2, 'MarkerSize', 15);
  loglog(H(1:7), H1_4(1:7), '*-k', 'linewidth', 2, 'MarkerSize', 15);
  loglog(H(1:7), H1_5(1:7), 's-.k', 'linewidth', 2, 'MarkerSize', 15);
  hold off;
  P1 = polyfit(log(H(1:5)), log(H1_1(1:5)), 1);
  P2 = polyfit(log(H(1:5)), log(H1_2(1:5)), 1);
  P3 = polyfit(log(H(1:5)), log(H1_3(1:5)), 1);
  P4 = polyfit(log(H(1:5)), log(H1_4(1:5)), 1);
  P5 = polyfit(log(H(1:5)), log(H1_5(1:5)), 1);
  legend(strcat('P1/P0  (slope=',num2str(P1(1)), ')'), ...
         strcat('P1+/P0 (slope=',num2str(P2(1)), ')'), ...
         strcat('P2/P1  (slope=',num2str(P3(1)), ')'), ...
         strcat('Q1/Q0  (slope=',num2str(P4(1)), ')'), ...
         strcat('Q2/Q1  (slope=',num2str(P5(1)), ')'), ...
         'Location', 'NorthWest');
  grid on;
  axesobj = findobj('type', 'axes');
  set(axesobj, 'fontname', 'times'); set(axesobj, 'fontunits', 'points');
  set(axesobj, 'fontsize', 18); set(axesobj, 'fontweight', 'bold');
  set(axesobj, 'linewidth', 2);
  xlabel('h');
  ylabel('H^1(\Omega) relative error (in %)');
  set(gca,'XTickLabel',{'0.001';'0.01';'0.1';'1';'...'}) 
  % axis([0.05 7 1e-4 10]);
  pause;


  loglog(H(1:7), L2C_1(1:7), 'o-k', 'linewidth', 2, 'MarkerSize', 15);
  hold on;
  loglog(H(1:7), L2C_2(1:7), 'x-.k', 'linewidth', 2, 'MarkerSize', 15);
  loglog(H(1:7), L2C_3(1:7), '+--k', 'linewidth', 2, 'MarkerSize', 15);
  loglog(H(1:7), L2C_4(1:7), '*-k', 'linewidth', 2, 'MarkerSize', 15);
  loglog(H(1:7), L2C_5(1:7), 's-.k', 'linewidth', 2, 'MarkerSize', 15);
  hold off;
  P1 = polyfit(log(H(1:5)), log(L2C_1(1:5)), 1);
  P2 = polyfit(log(H(1:5)), log(L2C_2(1:5)), 1);
  P3 = polyfit(log(H(1:5)), log(L2C_3(1:5)), 1);
  P4 = polyfit(log(H(1:5)), log(L2C_4(1:5)), 1);
  P5 = polyfit(log(H(1:5)), log(L2C_5(1:5)), 1);
  legend(strcat('P1/P0  (slope=',num2str(P1(1)), ')'), ...
         strcat('P1+/P0 (slope=',num2str(P2(1)), ')'), ...
         strcat('P2/P1  (slope=',num2str(P3(1)), ')'), ...
         strcat('Q1/Q0  (slope=',num2str(P4(1)), ')'), ...
         strcat('Q2/Q1  (slope=',num2str(P5(1)), ')'), ...
         'Location', 'NorthWest');
  grid on;
  axesobj = findobj('type', 'axes');
  set(axesobj, 'fontname', 'times'); set(axesobj, 'fontunits', 'points');
  set(axesobj, 'fontsize', 18); set(axesobj, 'fontweight', 'bold');
  set(axesobj, 'linewidth', 2);
  xlabel('h');
  ylabel('L^2(\Gamma_D) relative error (in %)');
  set(gca,'XTickLabel',{'0.001';'0.01';'0.1';'1';'...'}) 
  % axis([0.05 7 1e-4 10]);
  pause;

end;

% Pour mettre des fontes plus grosses.
% une commande
% get(findobj, 'type')
% renseigne sur les type d'objets � chercher.
% ensuite on recup�re les handles par
% axesobj = findobj('type', 'axes')
% par exemple, puis on peut faire
% set(axesobj, 'fontunits', 'points');
% set(axesobj, 'fontsize', 15);
% set(axesobj, 'fontweight', 'bold');
% Il vaut mieux a la fin decouper les images avec gimp par exemple.

axesobj = findobj('type', 'axes'); set(axesobj, 'fontname', 'times'); set(axesobj, 'fontunits', 'points'); set(axesobj, 'fontsize', 18); set(axesobj, 'fontweight', 'bold'); set(axesobj, 'linewidth', 2);


% Pour certains graphiques, il vaut mieux renommer les "ticks" par
%  set(gca,'XTickLabel',{'0.1';'1';'10';'...'})
%  set(gca,'YTickLabel',{'0.0001%';'0.001%';'0.01%';'0.1%';'1%';'10%'})     


% Pour sortir le graphique en png, faire par exemple :
% print(gcf,'-dpng','-r450', 'toto.png');

