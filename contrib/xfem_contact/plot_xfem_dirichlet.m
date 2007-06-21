% addpath ~/source++/getfem++/contrib/xfem_contact/

gf_workspace('clear all');
mf = gfMeshFem('load', 'xfem_dirichlet_ls.mf');
lsU = -load('xfem_dirichlet_ls.U')';

nn=5 % 0 : plot the exported mf
     % 1 : 
     % 2 : 
     % 3 : 
     % 4 : 
     % 5 : plot some convergence curves. 


clf
if nn==0,
  disp('plot the exported mf');
  [hsur, hcont] = gf_plot(mf,lsU,'refine',16,'contour',0,'mesh','on', 'pcolor','off');
  set(hcont{1}, 'LineWidth', 2);
  set(hcont{1}, 'Color', 'black');
elseif nn==1
  disp('plot the cut mesh');
  mc = gfMesh('load','cut.mesh');
  mfc = gfMeshFem(mc); 
  set(mfc,'classical_fem',2);
  
  lsUc = gf_compute(mf, lsU, 'interpolate on', mfc);
  
  %[hsur, hcont] = gf_plot(mf, lsU, 'refine', 24, 'zplot', 'on');
  %hold on;
  [hsur, hcont] = gf_plot(mfc, lsUc.*(lsUc>0), 'refine', 8, 'zplot', 'on', 'mesh','on', 'pcolor','on'); hold on;
  colormap([.8 1 .8]);
  %[hsur, hcont] = gf_plot(mf,lsU,'refine',4,'contour',0,'pcolor','off');
  
  %set(hcont{1}, 'LineWidth', 2);
  %set(hcont{1}, 'Color', 'black');
  %axis('tight'); axis off;
elseif nn==2 || nn==3,
  disp('plot the solution, with the 0 isovalue');
  sl=gfSlice('load','xfem_dirichlet.sl');
  slU=load('xfem_dirichlet.slU')';
  P=gf_slice_get(sl,'pts'); P=[P(1:2,:);slU];
  gf_slice_set(sl,'pts',P);
  gf_plot_slice(sl, 'data', slU, 'mesh','on','mesh_edges','on');
  
  

  m=get(mf, 'linked_mesh');
  slc=gfSlice({'isovalues', 0, mf, lsU, 0}, m, 16);
  hold on;
  P2=gf_slice_get(slc, 'pts');
  gf_slice_set(slc, 'pts', [P2;0.1 * ones(1,size(P2,2))]);
  [h1,h2,h3,h4]=gf_plot_slice(slc, 'tube','off','mesh_slice_edges_color','black');
  set(h4, 'LineWidth', 2);
  
  if (nn == 2),
    %set(hcont{1}, 'Color', 'black');
    view(3); camlight; view(2);
  else
    slc2=gfSlice('load', 'xfem_dirichlet.sl0');
    hold on;
    [h1,h2,h3,h4]=gf_plot_slice(slc2, 'tube','off','mesh_slice_edges_color','white');
    set(h4, 'LineWidth', 4);
    caxis([-.2 .3]); gf_colormap('froid');
  end;
elseif nn==4,
  disp('plotting the lagrange multipliers on the dirichlet boundary');
  sll=gfSlice('load','xfem_dirichlet.sll');
  slL=load('xfem_dirichlet.slL')';
  P0=gf_slice_get(sll, 'pts');
  [h1,h2,h3,h4]=gf_plot_slice(sll, 'tube','off','mesh_slice_edges_color',[.3 .3 .3]);
  hold on;
  gf_slice_set(sll,'pts',[P0 ; max(slL,-100)*0.05]);
  [hh1,hh2,hh3,hh4]=gf_plot_slice(sll, 'tube','off','mesh_slice_edges_color','black','mesh_slice_edges_width',1.5);
  sl=gfSlice('load','xfem_dirichlet.sl');
  gf_plot_slice(sl,'mesh','on');
  
  npt = size(P0, 2);
  P0 = [P0;zeros(1,npt)];
  P1 = gf_slice_get(sll,'pts'); 
  lseg = gf_slice_get(sll,'splxs', 1);
  F=[lseg(1,:) lseg(2,:); lseg(2,:) npt+lseg(2,:); npt+lseg(1,:) npt+lseg(1,:)];
  %F=[lseg; npt+lseg(2,:)];
  h=patch('Vertices',[P0 P1]', 'Faces', F');
  hold off;
  set(h,'FaceAlpha',0.3);
  set(h,'LineStyle','none');
  set(gcf,'renderer','opengl');
  set(gcf,'color','white');
  axis off;
  view(3);
  camzoom(1.7);
  axis([-0.5000    0.5000   -0.5000    0.5000 -.5 .5]);
  % print(gcf,'-dpng','-r300', 'lagrange_multipliers.png');
elseif nn==5,

  % Without stabilization, FEM_RHS = 'FEM_PK(2,3)'; LEVEL_SET_DEGREE = 2;
  H    = [1/320   1/160  1/80   1/40    1/20    1/10    1/5];
  % P1/P0 Non stabilisé
  L2_1 = [3.93    5.69   14.68  47.8    89.673  393.07  425.8 ]; % a refaire
  H1_1 = [15.44   19.82  27.21  49.23   56.57   127.73  197.5 ]; % a refaire
  L2C_1= [3.83e6  2.29e6 673394 1.078e6 56154   8.92e6 174884]; % a refaire
  % P1+/P0 Non stabilisé
  L2_2 = [0.89    0.19   0.98   2.95    12.72   42.75   221]; % a refaire
  H1_2 = [2.80    5.54   11.03  21.69   41.19   73.51   103]; % a refaire
  L2C_2= [1244    1436   3517   2636    1828    759     192]; % a refaire
  % P2/P1 Non stabilisé
  L2_3 = [0.00994 0.0417 0.203  1.00    7.42    42.1    109];
  H1_3 = [0.061   0.223  0.797  2.74    9.76    27.41   57];
  L2C_3= [0.222   1.22   3.03   3.59    16.3    10.3    74.5];
  % Q1/Q0 Non stabilisé
  L2_4 = [0.0434  0.61   0.79   1.963   8.63    52.16   131]; % a refaire
  H1_4 = [2.28    4.55   8.97   17.44   33.25   61.84   95.64]; % a refaire
  L2C_4= [203     205    208    214     237     267     348]; % a refaire
  % Q2/Q1 Non stabilisé
  L2_5 = [0.0752  0.67   0.405  1.636   7.916   61.91   106.2]; % a refaire
  H1_5 = [0.0148  0.3174 0.562  1.98    7.43    30.06   66.56]; % a refaire
  L2C_5= [201     202    225    202     258     201     136]; % a refaire
  

  loglog(H(1:7), L2_1(1:7), 'o-k', 'linewidth', 2, 'MarkerSize', 15);
  hold on;
  loglog(H(1:7), L2_2(1:7), 'x-.k', 'linewidth', 2, 'MarkerSize', 15);
  loglog(H(1:7), L2_3(1:7), '+--k', 'linewidth', 2, 'MarkerSize', 15);
  loglog(H(1:7), L2_4(1:7), '*-k', 'linewidth', 2, 'MarkerSize', 15);
  loglog(H(1:7), L2_5(1:7), 's-.k', 'linewidth', 2, 'MarkerSize', 15);
  hold off;
  xlabel('h');
  ylabel('L^2(\Omega) relative error in %');
  P1 = polyfit(log(H(1:5)), log(L2_1(1:5)), 1);
  P2 = polyfit(log(H(1:5)), log(L2_2(1:5)), 1);
  P3 = polyfit(log(H(1:5)), log(L2_3(1:5)), 1);
  P4 = polyfit(log(H(1:5)), log(L2_4(1:5)), 1);
  P5 = polyfit(log(H(2:5)), log(L2_5(2:5)), 1);
  legend(strcat('P1/P0     (slope=',num2str(P1(1)), ')'), ...
         strcat('P1+/P0 cubic bubble (slope=',num2str(P2(1)), ')'), ...
         strcat('P2/P1     (slope=',num2str(P3(1)), ')'), ...
         strcat('Q1/Q0     (slope=',num2str(P4(1)), ')'), ...
         strcat('Q2/Q1     (slope=',num2str(P5(1)), ')'), ...
         'Location', 'NorthWest');
  grid on;
  % axis([0.05 7 1e-4 10]);
  pause;

  loglog(H(1:7), H1_1(1:7), 'o-k', 'linewidth', 2, 'MarkerSize', 15);
  hold on;
  loglog(H(1:7), H1_2(1:7), 'x-.k', 'linewidth', 2, 'MarkerSize', 15);
  loglog(H(1:7), H1_3(1:7), '+--k', 'linewidth', 2, 'MarkerSize', 15);
  loglog(H(1:7), H1_4(1:7), '*-k', 'linewidth', 2, 'MarkerSize', 15);
  loglog(H(1:7), H1_5(1:7), 's-.k', 'linewidth', 2, 'MarkerSize', 15);
  hold off;
  xlabel('h');
  ylabel('H^1(\Omega) relative error in %');
  P1 = polyfit(log(H(1:5)), log(H1_1(1:5)), 1);
  P2 = polyfit(log(H(1:5)), log(H1_2(1:5)), 1);
  P3 = polyfit(log(H(1:5)), log(H1_3(1:5)), 1);
  P4 = polyfit(log(H(1:5)), log(H1_4(1:5)), 1);
  P5 = polyfit(log(H(1:5)), log(H1_5(1:5)), 1);
  legend(strcat('P1/P0     (slope=',num2str(P1(1)), ')'), ...
         strcat('P1+/P0 cubic bubble (slope=',num2str(P2(1)), ')'), ...
         strcat('P2/P1     (slope=',num2str(P3(1)), ')'), ...
         strcat('Q1/Q0     (slope=',num2str(P4(1)), ')'), ...
         strcat('Q2/Q1     (slope=',num2str(P5(1)), ')'), ...
         'Location', 'NorthWest');
  grid on;
  % axis([0.05 7 1e-4 10]);
  pause;


  loglog(H(1:7), L2C_1(1:7), 'o-k', 'linewidth', 2, 'MarkerSize', 15);
  hold on;
  loglog(H(1:7), L2C_2(1:7), 'x-.k', 'linewidth', 2, 'MarkerSize', 15);
  loglog(H(1:7), L2C_3(1:7), '+--k', 'linewidth', 2, 'MarkerSize', 15);
  loglog(H(1:7), L2C_4(1:7), '*-k', 'linewidth', 2, 'MarkerSize', 15);
  loglog(H(1:7), L2C_5(1:7), 's-.k', 'linewidth', 2, 'MarkerSize', 15);
  hold off;
  xlabel('h');
  ylabel('L^2(\Gamma_D) relative error in %');
  P1 = polyfit(log(H(1:5)), log(L2C_1(1:5)), 1);
  P2 = polyfit(log(H(1:5)), log(L2C_2(1:5)), 1);
  P3 = polyfit(log(H(1:5)), log(L2C_3(1:5)), 1);
  P4 = polyfit(log(H(1:5)), log(L2C_4(1:5)), 1);
  P5 = polyfit(log(H(1:5)), log(L2C_5(1:5)), 1);
  legend(strcat('P1/P0     (slope=',num2str(P1(1)), ')'), ...
         strcat('P1+/P0 cubic bubble (slope=',num2str(P2(1)), ')'), ...
         strcat('P2/P1     (slope=',num2str(P3(1)), ')'), ...
         strcat('Q1/Q0     (slope=',num2str(P4(1)), ')'), ...
         strcat('Q2/Q1     (slope=',num2str(P5(1)), ')'), ...
         'Location', 'NorthWest');
  grid on;
  % axis([0.05 7 1e-4 10]);
  pause;


  % With BB stabilization, gamma0 = 0.01, FEM_RHS = 'FEM_PK(2,3)';
  % LEVEL_SET_DEGREE = 2;
  H    = [1/320   1/160  1/80   1/40    1/20    1/10    1/5];
  % P1/P0 stabilisé
  L2_1 = [0.058   0.231  0.971  3.64    18.81   73.78   67]; % a refaire
  H1_1 = [2.78    5.54   11.04  21.85   41.55   71.92   98]; % a refaire
  L2C_1= [18.57   34.86  65.91  120.1   152.7   133.1   175]; % a refaire
  % P1+/P0 stabilisé
  L2_2 = [0.11    0.18   0.98   2.73    11.8    39.7    218]; % a refaire
  H1_2 = [2.79    5.49   10.86  21.2    40.5    72.7    103]; % a refaire
  L2C_2= [206     207    207    208     212     214     191]; % a refaire
  % P2/P1 stabilisé
  L2_3 = [0.0094  0.041  0.20   1.00    7.42    42.2    113];
  H1_3 = [0.0613  0.223  0.79   2.74    9.76    27.4    57];
  L2C_3= [0.087   0.114  0.31   1.18    14      9.80    52];
  % Q1/Q0 stabilisé
  L2_4 = [0.72    0.72   0.72   0.72    0.72    0.72    0.72]; % a refaire
  H1_4 = [0.72    0.72   0.72   0.72    0.72    0.72    0.72]; % a refaire
  L2C_4= [300     300    300    300     300     300     300];  % a refaire
  % Q2/Q1 stabilisé
  L2_5 = [0.72    0.72   0.72   0.72    0.72    0.72    0.72]; % a refaire
  H1_5 = [0.72    0.72   0.72   0.72    0.72    0.72    0.72]; % a refaire
  L2C_5= [300     300    300    300     300     300     300];  % a refaire


  loglog(H(1:7), L2_1(1:7), 'o-k', 'linewidth', 2, 'MarkerSize', 15);
  hold on;
  loglog(H(1:7), L2_2(1:7), 'x-.k', 'linewidth', 2, 'MarkerSize', 15);
  loglog(H(1:7), L2_3(1:7), '+--k', 'linewidth', 2, 'MarkerSize', 15);
  loglog(H(1:7), L2_4(1:7), '*-k', 'linewidth', 2, 'MarkerSize', 15);
  loglog(H(1:7), L2_5(1:7), 's-.k', 'linewidth', 2, 'MarkerSize', 15);
  hold off;
  xlabel('h');
  ylabel('L^2(\Omega) relative error in %');
  P1 = polyfit(log(H(1:5)), log(L2_1(1:5)), 1);
  P2 = polyfit(log(H(1:5)), log(L2_2(1:5)), 1);
  P3 = polyfit(log(H(1:5)), log(L2_3(1:5)), 1);
  P4 = polyfit(log(H(1:5)), log(L2_4(1:5)), 1);
  P5 = polyfit(log(H(2:5)), log(L2_5(2:5)), 1);
  legend(strcat('P1/P0     (slope=',num2str(P1(1)), ')'), ...
         strcat('P1+/P0 cubic bubble (slope=',num2str(P2(1)), ')'), ...
         strcat('P2/P1     (slope=',num2str(P3(1)), ')'), ...
         strcat('Q1/Q0     (slope=',num2str(P4(1)), ')'), ...
         strcat('Q2/Q1     (slope=',num2str(P5(1)), ')'), ...
         'Location', 'NorthWest');
  grid on;
  % axis([0.05 7 1e-4 10]);
  pause;

  loglog(H(1:7), H1_1(1:7), 'o-k', 'linewidth', 2, 'MarkerSize', 15);
  hold on;
  loglog(H(1:7), H1_2(1:7), 'x-.k', 'linewidth', 2, 'MarkerSize', 15);
  loglog(H(1:7), H1_3(1:7), '+--k', 'linewidth', 2, 'MarkerSize', 15);
  loglog(H(1:7), H1_4(1:7), '*-k', 'linewidth', 2, 'MarkerSize', 15);
  loglog(H(1:7), H1_5(1:7), 's-.k', 'linewidth', 2, 'MarkerSize', 15);
  hold off;
  xlabel('h');
  ylabel('H^1(\Omega) relative error in %');
  P1 = polyfit(log(H(1:5)), log(H1_1(1:5)), 1);
  P2 = polyfit(log(H(1:5)), log(H1_2(1:5)), 1);
  P3 = polyfit(log(H(1:5)), log(H1_3(1:5)), 1);
  P4 = polyfit(log(H(1:5)), log(H1_4(1:5)), 1);
  P5 = polyfit(log(H(1:5)), log(H1_5(1:5)), 1);
  legend(strcat('P1/P0     (slope=',num2str(P1(1)), ')'), ...
         strcat('P1+/P0 cubic bubble (slope=',num2str(P2(1)), ')'), ...
         strcat('P2/P1     (slope=',num2str(P3(1)), ')'), ...
         strcat('Q1/Q0     (slope=',num2str(P4(1)), ')'), ...
         strcat('Q2/Q1     (slope=',num2str(P5(1)), ')'), ...
         'Location', 'NorthWest');
  grid on;
  % axis([0.05 7 1e-4 10]);
  pause;


  loglog(H(1:7), L2C_1(1:7), 'o-k', 'linewidth', 2, 'MarkerSize', 15);
  hold on;
  loglog(H(1:7), L2C_2(1:7), 'x-.k', 'linewidth', 2, 'MarkerSize', 15);
  loglog(H(1:7), L2C_3(1:7), '+--k', 'linewidth', 2, 'MarkerSize', 15);
  loglog(H(1:7), L2C_4(1:7), '*-k', 'linewidth', 2, 'MarkerSize', 15);
  loglog(H(1:7), L2C_5(1:7), 's-.k', 'linewidth', 2, 'MarkerSize', 15);
  hold off;
  xlabel('h');
  ylabel('L^2(\Gamma_D) relative error in %');
  P1 = polyfit(log(H(1:5)), log(L2C_1(1:5)), 1);
  P2 = polyfit(log(H(1:5)), log(L2C_2(1:5)), 1);
  P3 = polyfit(log(H(1:5)), log(L2C_3(1:5)), 1);
  P4 = polyfit(log(H(1:5)), log(L2C_4(1:5)), 1);
  P5 = polyfit(log(H(1:5)), log(L2C_5(1:5)), 1);
  legend(strcat('P1/P0     (slope=',num2str(P1(1)), ')'), ...
         strcat('P1+/P0 cubic bubble (slope=',num2str(P2(1)), ')'), ...
         strcat('P2/P1     (slope=',num2str(P3(1)), ')'), ...
         strcat('Q1/Q0     (slope=',num2str(P4(1)), ')'), ...
         strcat('Q2/Q1     (slope=',num2str(P5(1)), ')'), ...
         'Location', 'NorthWest');
  grid on;
  % axis([0.05 7 1e-4 10]);
  pause;

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

