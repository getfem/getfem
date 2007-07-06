% addpath ~/source++/getfem++/contrib/xfem_contact/

gf_workspace('clear all');
mf = gfMeshFem('load', 'xfem_dirichlet_ls.mf');
lsU = -load('xfem_dirichlet_ls.U')';

nn=5 % 0 : plot the exported mf
     % 1 : 
     % 2 : 
     % 3 : 
     % 4 : plotting the lagrange multipliers on the dirichlet boundary
     % 5 : plot some convergence curvesthe solution. 
     % 6 : plot some convergence curves. 


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
  disp('plot the solution');
  lsU = max(lsU, 0) * 10;
  [hsur, hcont] = gf_plot(mf, lsU, 'refine',16,'mesh','on', 'zplot', 'on');
  colormap([0.7 0.8 1.0; 0.8 1 0.8; 0.8 1 0.8; 0.8 1 0.8; 0.8 1 0.8; 0.8 1 0.8; 0.8 1 0.8; 0.8 1 0.8; 0.8 1 0.8; 0.8 1 0.8; 0.8 1 0.8; 0.8 1 0.8; 0.8 1 0.8; 0.8 1 0.8; 0.8 1 0.8; 0.8 1 0.8; 0.8 1 0.8; 0.8 1 0.8; 0.8 1 0.8; 0.8 1 0.8; 0.8 1 0.8; 0.8 1 0.8]);
  axis off;
  camlight;
  % set(hcont{1}, 'LineWidth', 2);
  % set(hcont{1}, 'Color', 'black');
elseif nn==6,

  % Without stabilization, FEM_RHS = 'FEM_PK(2,3)'; LEVEL_SET_DEGREE = 2;
  H    = [1/320   1/160  1/80   1/40    1/20    1/10    1/5];
  % P1/P0 Non stabilisé
  L2_1 = [3.93    8.93   16.65  25.11   77.88   118     95 ];
  H1_1 = [15.44   23.89  34.83  45.18   62.96   85      102];
  L2C_1= [41700   31300  18300  9560    5000    6620    323];
  % P1+/P0 Non stabilisé
  L2_2 = [0.14    0.25   1.26   3.67    23      73.17   133];
  H1_2 = [2.99    5.28   10.7   20.73   41.2    66.67   82];
  L2C_2= [3870    816    271    2351    5466    176     66];
  % P2/P1 Non stabilisé
  L2_3 = [0.0095  0.0417 0.203  1.00    7.42    42.1    109];
  H1_3 = [0.061   0.223  0.797  2.74    9.76    27.41   57];
  L2C_3= [0.222   1.22   3.03   3.59    16.3    10.3    74.5];
  % Q1/Q0 Non stabilisé
  L2_4 = [0.0383  0.154  0.63   2.86    15.04   39.57   251];
  H1_4 = [2.29    4.57   9.07   18.05   34.81   63.54   149];
  L2C_4= [18.64   10.9   25.29  591     259     136     95];
  % Q2/Q1 Non stabilisé
  L2_5 = [0.00442 0.0199 0.086  0.43    2.47    17.92   68.52];
  H1_5 = [0.0451  0.3224 0.596  2.07    6.90    20.99   57.85];
  L2C_5= [0.132   15.51  1.22   1.69    2.07    5.89    40.43];
  

  loglog(H(1:7), L2_1(1:7), 'o-k', 'linewidth', 2, 'MarkerSize', 15);
  hold on;
  loglog(H(1:7), L2_2(1:7), 'x-.k', 'linewidth', 2, 'MarkerSize', 15);
  loglog(H(1:7), L2_3(1:7), '+--k', 'linewidth', 2, 'MarkerSize', 15);
  loglog(H(1:7), L2_4(1:7), '*-k', 'linewidth', 2, 'MarkerSize', 15);
  loglog(H(1:7), L2_5(1:7), 's-.k', 'linewidth', 2, 'MarkerSize', 15);
  hold off;
  P1 = polyfit(log(H(1:7)), log(L2_1(1:7)), 1);
  P2 = polyfit(log(H(1:7)), log(L2_2(1:7)), 1);
  P3 = polyfit(log(H(1:7)), log(L2_3(1:7)), 1);
  P4 = polyfit(log(H(1:7)), log(L2_4(1:7)), 1);
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
  ylabel('L^2(\Omega) relative error');
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
  P1 = polyfit(log(H(1:7)), log(H1_1(1:7)), 1);
  P2 = polyfit(log(H(1:7)), log(H1_2(1:7)), 1);
  P3 = polyfit(log(H(1:7)), log(H1_3(1:7)), 1);
  P4 = polyfit(log(H(1:7)), log(H1_4(1:7)), 1);
  P5 = polyfit(log(H(1:7)), log(H1_5(1:7)), 1);
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
  ylabel('H^1(\Omega) relative error');
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
  P1 = polyfit(log(H(1:7)), log(L2C_1(1:7)), 1);
  P2 = polyfit(log(H(1:7)), log(L2C_2(1:7)), 1);
  P3 = polyfit(log(H(1:7)), log(L2C_3(1:7)), 1);
  P4 = polyfit(log(H(1:7)), log(L2C_4(1:7)), 1);
  P5 = polyfit(log(H(1:7)), log(L2C_5(1:7)), 1);
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
  ylabel('L^2(\Gamma_D) relative error');
  set(gca,'XTickLabel',{'0.001';'0.01';'0.1';'1';'...'}) 
  % axis([0.05 7 1e-4 10]);
  pause;


  % With BB stabilization, gamma0 = 0.1, FEM_RHS = 'FEM_PK(2,3)';
  % LEVEL_SET_DEGREE = 2;
  H    = [1/320   1/160  1/80   1/40    1/20    1/10    1/5];
  % P1/P0 stabilisé
  L2_1 = [0.058   0.231  0.971  3.74    17.23   69.08   71];
  H1_1 = [2.77    5.54   10.98  21.74   41.40   71.06   96.8];
  L2C_1= [3.40    5.86   10.98  32.42   42.34   48.59   70.91];
  % P1+/P0 stabilisé
  L2_2 = [0.054   0.214  0.88   3.43    16.27   64.18   83.37];
  H1_2 = [2.58    5.14   10.21  20.0    38.67   64.91   82.69];
  L2C_2= [3.28    3.86   7.31   21.93   31.59   42.60   42.07];
  % P2/P1 stabilisé
  L2_3 = [0.0095  0.041  0.20   1.01    7.40    42.6    120];
  H1_3 = [0.0613  0.358  0.80   2.75    10.55   27.4    72];
  L2C_3= [0.030   0.801  1.12   1.35    17.85   7.78    53];
  % Q1/Q0 stabilisé
  L2_4 = [0.036   0.1451 0.577  2.33    9.56    37.43   52.72];
  H1_4 = [2.29    4.58   9.09   17.95   35.02   61.98   94.34];
  L2C_4= [0.87    3.36   3.35   6.68    14.86   31.82   38.15];
  % Q2/Q1 stabilisé
  L2_5 = [0.0044  0.0262 0.087  0.439   2.51    18.44   137];
  H1_5 = [0.127   0.433  0.61   2.40    8.15    36.89   74];
  L2C_5= [0.172   0.22   0.25   0.578   2.39    11.32   65];


  loglog(H(1:7), L2_1(1:7), 'o-k', 'linewidth', 2, 'MarkerSize', 15);
  hold on;
  loglog(H(1:7), L2_2(1:7), 'x-.k', 'linewidth', 2, 'MarkerSize', 15);
  loglog(H(1:7), L2_3(1:7), '+--k', 'linewidth', 2, 'MarkerSize', 15);
  loglog(H(1:7), L2_4(1:7), '*-k', 'linewidth', 2, 'MarkerSize', 15);
  loglog(H(1:7), L2_5(1:7), 's-.k', 'linewidth', 2, 'MarkerSize', 15);
  hold off;
  P1 = polyfit(log(H(1:7)), log(L2_1(1:7)), 1);
  P2 = polyfit(log(H(1:7)), log(L2_2(1:7)), 1);
  P3 = polyfit(log(H(1:7)), log(L2_3(1:7)), 1);
  P4 = polyfit(log(H(1:7)), log(L2_4(1:7)), 1);
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
  ylabel('L^2(\Omega) relative error');
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
  P1 = polyfit(log(H(1:7)), log(H1_1(1:7)), 1);
  P2 = polyfit(log(H(1:7)), log(H1_2(1:7)), 1);
  P3 = polyfit(log(H(1:7)), log(H1_3(1:7)), 1);
  P4 = polyfit(log(H(1:7)), log(H1_4(1:7)), 1);
  P5 = polyfit(log(H(1:7)), log(H1_5(1:7)), 1);
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
  ylabel('H^1(\Omega) relative error');
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
  P1 = polyfit(log(H(1:7)), log(L2C_1(1:7)), 1);
  P2 = polyfit(log(H(1:7)), log(L2C_2(1:7)), 1);
  P3 = polyfit(log(H(1:7)), log(L2C_3(1:7)), 1);
  P4 = polyfit(log(H(1:7)), log(L2C_4(1:7)), 1);
  P5 = polyfit(log(H(1:7)), log(L2C_5(1:7)), 1);
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
  ylabel('L^2(\Gamma_D) relative error');
  set(gca,'XTickLabel',{'0.001';'0.01';'0.1';'1';'...'}) 
  % axis([0.05 7 1e-4 10]);
  pause;



  % With BB stabilization and stabilized normal derivative, gamma0 = 0.1,
  % FEM_RHS = 'FEM_PK(2,3)'; LEVEL_SET_DEGREE = 2; MINIMAL_ELT_RATIO = 0.01;
  H    = [1/320   1/160  1/80   1/40    1/20    1/10    1/5];
  % P1/P0 fully stabilised
  L2_1 = [0.058   0.231  0.955  3.73    17.24   69.08   71];
  H1_1 = [2.77    5.53   10.98  21.60   40.95   71      96.8];
  L2C_1= [3.27    5.74   10.70  22.36   30.67   48      71];   
  % P1+/P0 fully stabilised
  L2_2 = [0.054   0.21   0.88   3.43    16.36   64.17   83.43];
  H1_2 = [2.62    5.15   10.20  19.9    38.12   64.81   82.64];
  L2C_2= [2.39    3.83   7.24   15.73   30.32   41.64   41.85];
  % P2/P1 fully stabilised
  L2_3 = [0.0094  0.041  0.20   1.01    7.38    42.6    119];
  H1_3 = [0.06    0.22   0.79   2.75    10.38   27.8    71.8];
  L2C_3= [0.030   0.064  0.182  0.66    9.01    10.7    41];
  % Q1/Q0 fully stabilised
  L2_4 = [0.036   0.1451 0.577  2.33    9.56    37.80   52.20];
  H1_4 = [2.29    4.58   9.09   17.95   35.02   61.96   94.39];
  L2C_4= [0.78    3.31   3.16   6.67    14.34   30.16   39.01];
  % Q2/Q1 fully stabilised
  L2_5 = [0.0044  0.0262 0.087  0.439   2.51    18.22   137];
  H1_5 = [0.045   0.533  0.61   2.40    8.32    23.22   74];
  L2C_5= [0.021   0.32   0.27   0.56    2.63    6.06    65];


  loglog(H(1:7), L2_1(1:7), 'o-k', 'linewidth', 2, 'MarkerSize', 15);
  hold on;
  loglog(H(1:7), L2_2(1:7), 'x-.k', 'linewidth', 2, 'MarkerSize', 15);
  loglog(H(1:7), L2_3(1:7), '+--k', 'linewidth', 2, 'MarkerSize', 15);
  loglog(H(1:7), L2_4(1:7), '*-k', 'linewidth', 2, 'MarkerSize', 15);
  loglog(H(1:7), L2_5(1:7), 's-.k', 'linewidth', 2, 'MarkerSize', 15);
  hold off;
  P1 = polyfit(log(H(1:7)), log(L2_1(1:7)), 1);
  P2 = polyfit(log(H(1:7)), log(L2_2(1:7)), 1);
  P3 = polyfit(log(H(1:7)), log(L2_3(1:7)), 1);
  P4 = polyfit(log(H(1:7)), log(L2_4(1:7)), 1);
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
  ylabel('L^2(\Omega) relative error');
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
  P1 = polyfit(log(H(1:7)), log(H1_1(1:7)), 1);
  P2 = polyfit(log(H(1:7)), log(H1_2(1:7)), 1);
  P3 = polyfit(log(H(1:7)), log(H1_3(1:7)), 1);
  P4 = polyfit(log(H(1:7)), log(H1_4(1:7)), 1);
  P5 = polyfit(log(H(1:7)), log(H1_5(1:7)), 1);
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
  ylabel('H^1(\Omega) relative error');
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
  P1 = polyfit(log(H(1:7)), log(L2C_1(1:7)), 1);
  P2 = polyfit(log(H(1:7)), log(L2C_2(1:7)), 1);
  P3 = polyfit(log(H(1:7)), log(L2C_3(1:7)), 1);
  P4 = polyfit(log(H(1:7)), log(L2C_4(1:7)), 1);
  P5 = polyfit(log(H(1:7)), log(L2C_5(1:7)), 1);
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
  ylabel('L^2(\Gamma_D) relative error');
  set(gca,'XTickLabel',{'0.001';'0.01';'0.1';'1';'...'}) 
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

axesobj = findobj('type', 'axes'); set(axesobj, 'fontname', 'times'); set(axesobj, 'fontunits', 'points'); set(axesobj, 'fontsize', 18); set(axesobj, 'fontweight', 'bold'); set(axesobj, 'linewidth', 2);


% Pour certains graphiques, il vaut mieux renommer les "ticks" par
%  set(gca,'XTickLabel',{'0.1';'1';'10';'...'})
%  set(gca,'YTickLabel',{'0.0001%';'0.001%';'0.01%';'0.1%';'1%';'10%'})     


% Pour sortir le graphique en png, faire par exemple :
% print(gcf,'-dpng','-r450', 'toto.png');

