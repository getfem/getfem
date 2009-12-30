% addpath ~/source++/getfem++/contrib/mixed_dynamic_friction

if (0)
  
  %
  % last data produced
  %
      
  A = load('mixed_scalar_hyperbolic.data');

  % energy curves

  plot(A(:, 1), A(:, 2), '-k', 'linewidth', 2, 'MarkerSize', 15);
  % axis([0 0.7 0 0.1]);
  title('energy');
  pause;

  % displacement curves

  plot(A(:, 1), A(:, 3), '-k', 'linewidth', 2, 'MarkerSize', 15);
  % axis([0 0.7 -0.01 0.06]);
  title('displacement');
  pause;
 
  % contact stress curves

  plot(A(:, 1), A(:, 4), '-k', 'linewidth', 2, 'MarkerSize', 15);
  % axis([0 0.7 -0.4 0.02]);
  title('stress');
  pause;

end;

%
% Convergence for h and dt decreasing and P1+/P0 method
%

if (1)

  A1 = load('mixed_scalar_hyperbolic_P1P0_NX4_DT01.data');
  A2 = load('mixed_scalar_hyperbolic_P1P0_NX10_DT005.data');
  A3 = load('mixed_scalar_hyperbolic_P1P0_NX20_DT0025.data');
  A4 = load('mixed_scalar_hyperbolic_P1P0_NX40_DT00125.data');
  A5 = load('mixed_scalar_hyperbolic_P1P0_NX80_DT000625.data');

  % energy curves

  plot(A1(:, 1), A1(:, 2), '-k', 'linewidth', 2, 'MarkerSize', 15);
  hold on;
  plot(A2(:, 1), A2(:, 2), '--k', 'linewidth', 2, 'MarkerSize', 15);
  plot(A3(:, 1), A3(:, 2), '-.k', 'linewidth', 2, 'MarkerSize', 15);
  plot(A4(:, 1), A4(:, 2), '-.k', 'linewidth', 2, 'MarkerSize', 15);
  plot(A5(:, 1), A5(:, 2), '-.k', 'linewidth', 2, 'MarkerSize', 15);
  hold off;
  % axis([0 0.7 0 0.1]);
  axis([0 0.7 0 0.02]);
  xlabel('t');
  ylabel('total energy');
  legend('dt = 0.08', 'dt = 0.04', 'dt = 0.02', 'dt = 0.01', 'dt = 0.005', 'Location', 'SouthWest');
  axesobj = findobj('type', 'axes'); set(axesobj, 'fontname', 'times'); set(axesobj, 'fontunits', 'points'); set(axesobj, 'fontsize', 24); set(axesobj, 'fontweight', 'bold'); set(axesobj, 'linewidth', 2);
  pause;
  print(gcf,'-deps','-r450', 'energy.eps');

  % displacement curves

  plot(A1(:, 1), A1(:, 3), '-k', 'linewidth', 2, 'MarkerSize', 15);
  hold on;
  plot(A2(:, 1), A2(:, 3), '--k', 'linewidth', 2, 'MarkerSize', 15);
  plot(A3(:, 1), A3(:, 3), '-.k', 'linewidth', 2, 'MarkerSize', 15);
  plot(A4(:, 1), A4(:, 3), '-.k', 'linewidth', 2, 'MarkerSize', 15);
  plot(A5(:, 1), A5(:, 3), '-.k', 'linewidth', 2, 'MarkerSize', 15);
  hold off;
  % axis([0 0.7 -0.01 0.06]);
  axis([0 0.7 -0.01 0.06]);
  xlabel('t');
  ylabel('center point displacement');
  legend('dt = 0.01', 'dt = 0.001', 'dt = 0.0001', 'Location', 'SouthWest');
  axesobj = findobj('type', 'axes'); set(axesobj, 'fontname', 'times'); set(axesobj, 'fontunits', 'points'); set(axesobj, 'fontsize', 24); set(axesobj, 'fontweight', 'bold'); set(axesobj, 'linewidth', 2);
  pause;
  print(gcf,'-deps','-r450', 'displacement.eps');

  % contact stress curves

  plot(A1(:, 1), A1(:, 4), '-k', 'linewidth', 2, 'MarkerSize', 15);
  hold on;
  plot(A2(:, 1), A2(:, 4), '--k', 'linewidth', 2, 'MarkerSize', 15);
  plot(A3(:, 1), A3(:, 4), '-.k', 'linewidth', 2, 'MarkerSize', 15);
  plot(A4(:, 1), A4(:, 4), '-.k', 'linewidth', 2, 'MarkerSize', 15);
  plot(A5(:, 1), A5(:, 4), '-.k', 'linewidth', 2, 'MarkerSize', 15);
  hold off;
  axis([0 0.7 -0.6 0.02]);
  xlabel('t');
  ylabel('center point contact stress');
  legend('dt = 0.01', 'dt = 0.001', 'dt = 0.0001', 'Location', 'SouthWest');
  axesobj = findobj('type', 'axes'); set(axesobj, 'fontname', 'times'); set(axesobj, 'fontunits', 'points'); set(axesobj, 'fontsize', 24); set(axesobj, 'fontweight', 'bold'); set(axesobj, 'linewidth', 2);
  pause;
  print(gcf,'-deps','-r450', 'stress');


end;



%
% Convergence in energy curves
%

if (1)
  % P1+/P0 method
  A = [0.1    0.04   0.01     0.004     0.001     0.0004    0.0001;  % dt
       165.9  1.129  0.3874  8.083e-2  9.166e-3  8.333e-4  9.166e-5; % NX=4
       217.8  138.0  7.498   3.064     0.4933    2.666e-2  1.666e-3; % NX=10
       208.6  185.1  27.69   3.251     0.3483    6.0e-2    4.999e-3; % NX=20
       207    208    96.94   22.84499  0.8558    0.2424    0.019999; % NX=40
       207    211    75.33   168       4.30      0.5608    0.06833   % NX=80
       ];
  loglog(A(1,:), A(2, :), '-k', 'linewidth', 2, 'MarkerSize', 15);
  hold on;
  loglog(A(1,:), A(3, :), '--k', 'linewidth', 2, 'MarkerSize', 15);
  loglog(A(1,:), A(4, :), '-.k', 'linewidth', 2, 'MarkerSize', 15);
  loglog(A(1,:), A(5, :), '-.k', 'linewidth', 2, 'MarkerSize', 15);
  loglog(A(1,:), A(6, :), '-.k', 'linewidth', 2, 'MarkerSize', 15);
  hold off;
  xlabel('dt');
  ylabel('relative error on total energy (%)');
  legend('h = 0.25', 'h = 0.1', 'h = 0.05', 'h = 0.025', 'h = 0.0125', 'Location', 'SouthEast');
  axesobj = findobj('type', 'axes'); set(axesobj, 'fontname', 'times'); set(axesobj, 'fontunits', 'points'); set(axesobj, 'fontsize', 24); set(axesobj, 'fontweight', 'bold'); set(axesobj, 'linewidth', 2);
  pause;
  print(gcf,'-deps','-r450', 'toto');



  % P2/P1 method
  A = [0.1    0.04   0.01     0.004     0.001     0.0004    0.0001;     % dt
       42     3.869  1.7191   0.1124    0.01199   0.0008333 0.00006666; % NX=4
       88     20.43  0.6866   0.02666   0.03999   0.0008333 0.0008333;  % NX=10
       167    69     1.29     0.41      0.012499  0.0041666 0.0003333;  % NX=20
       195    166    6.0674   2.44666   0.0641666 0.032499  0.00166666; % NX=40
       204    198    23.15    5.53160   0.45083   0.1816    0.0033333   % NX=80
       ];
  loglog(A(1,:), A(2, :), '-k', 'linewidth', 2, 'MarkrSize', 15);
  hold on;
  loglog(A(1,:), A(3, :), '--k', 'linewidth', 2, 'MarkerSize', 15);
  loglog(A(1,:), A(4, :), '-.k', 'linewidth', 2, 'MarkerSize', 15);
  loglog(A(1,:), A(5, :), '-.k', 'linewidth', 2, 'MarkerSize', 15);
  loglog(A(1,:), A(6, :), '-.k', 'linewidth', 2, 'MarkerSize', 15);
  hold off;
  xlabel('dt');
  ylabel('relative error on total energy (%)');
  legend('h = 0.25', 'h = 0.1', 'h = 0.05', 'h = 0.025', 'h = 0.0125', 'Location', 'SouthEast');
  axesobj = findobj('type', 'axes'); set(axesobj, 'fontname', 'times'); set(axesobj, 'fontunits', 'points'); set(axesobj, 'fontsize', 24); set(axesobj, 'fontweight', 'bold'); set(axesobj, 'linewidth', 2);
  pause;
  print(gcf,'-deps','-r450', 'toto');

end;





%
% P1plusP0
%

if (0)

  A1 = load('msh_dt0.01.data');
  A2 = load('msh_dt0.001.data');
  A3 = load('msh_dt0.0001.data');

  % energy curves

  plot(A1(:, 1), A1(:, 2), '-k', 'linewidth', 2, 'MarkerSize', 15);
  hold on;
  plot(A2(:, 1), A2(:, 2), '--k', 'linewidth', 2, 'MarkerSize', 15);
  plot(A3(:, 1), A3(:, 2), '-.k', 'linewidth', 2, 'MarkerSize', 15);
  hold off;
  % axis([0 0.7 0 0.1]);
  axis([0 0.7 0 0.02]);
  xlabel('t');
  ylabel('total energy');
  legend('dt = 0.01', 'dt = 0.001', 'dt = 0.0001', 'Location', 'SouthWest');
  axesobj = findobj('type', 'axes'); set(axesobj, 'fontname', 'times'); set(axesobj, 'fontunits', 'points'); set(axesobj, 'fontsize', 24); set(axesobj, 'fontweight', 'bold'); set(axesobj, 'linewidth', 2);
  pause;
  print(gcf,'-deps','-r450', 'energy.eps');

  % displacement curves

  plot(A1(:, 1), A1(:, 3), '-k', 'linewidth', 2, 'MarkerSize', 15);
  hold on;
  plot(A2(:, 1), A2(:, 3), '--k', 'linewidth', 2, 'MarkerSize', 15);
  plot(A3(:, 1), A3(:, 3), '-.k', 'linewidth', 2, 'MarkerSize', 15);
  hold off;
  % axis([0 0.7 -0.01 0.06]);
  axis([0 0.7 -0.01 0.06]);
  xlabel('t');
  ylabel('center point displacement');
  legend('dt = 0.01', 'dt = 0.001', 'dt = 0.0001', 'Location', 'SouthWest');
  axesobj = findobj('type', 'axes'); set(axesobj, 'fontname', 'times'); set(axesobj, 'fontunits', 'points'); set(axesobj, 'fontsize', 24); set(axesobj, 'fontweight', 'bold'); set(axesobj, 'linewidth', 2);
  pause;
  print(gcf,'-deps','-r450', 'displacement.eps');

  % contact stress curves

  plot(A1(:, 1), A1(:, 4), '-k', 'linewidth', 2, 'MarkerSize', 15);
  hold on;
  plot(A2(:, 1), A2(:, 4), '--k', 'linewidth', 2, 'MarkerSize', 15);
  plot(A3(:, 1), A3(:, 4), '-.k', 'linewidth', 2, 'MarkerSize', 15);
  hold off;
  axis([0 0.7 -0.6 0.02]);
  xlabel('t');
  ylabel('center point contact stress');
  legend('dt = 0.01', 'dt = 0.001', 'dt = 0.0001', 'Location', 'SouthWest');
  axesobj = findobj('type', 'axes'); set(axesobj, 'fontname', 'times'); set(axesobj, 'fontunits', 'points'); set(axesobj, 'fontsize', 24); set(axesobj, 'fontweight', 'bold'); set(axesobj, 'linewidth', 2);
  pause;
  print(gcf,'-deps','-r450', 'stress');


end;


  

%  loglog(H(1:7), L2_1(1:7), 'o-k', 'linewidth', 2, 'MarkerSize', 15);
%  hold on;
%  loglog(H(1:7), L2_2(1:7), 'x-.k', 'linewidth', 2, 'MarkerSize', 15);
%  loglog(H(1:7), L2_3(1:7), '+--k', 'linewidth', 2, 'MarkerSize', 15);
%  loglog(H(1:7), L2_4(1:7), '*-k', 'linewidth', 2, 'MarkerSize', 15);
%  loglog(H(1:7), L2_5(1:7), 's-.k', 'linewidth', 2, 'MarkerSize', 15);
%  hold off;
%  P1 = polyfit(log(H(1:7)), log(L2_1(1:7)), 1);
%  P2 = polyfit(log(H(1:7)), log(L2_2(1:7)), 1);
%  P3 = polyfit(log(H(1:7)), log(L2_3(1:7)), 1);
%  P4 = polyfit(log(H(1:7)), log(L2_4(1:7)), 1);
%  P5 = polyfit(log(H(2:5)), log(L2_5(2:5)), 1);
%  legend(strcat('P1/P0  (slope=',num2str(P1(1)), ')'), ...
%         strcat('P1+/P0 (slope=',num2str(P2(1)), ')'), ...
%         strcat('P2/P1  (slope=',num2str(P3(1)), ')'), ...
%         strcat('Q1/Q0  (slope=',num2str(P4(1)), ')'), ...
%         strcat('Q2/Q1  (slope=',num2str(P5(1)), ')'), ...
%         'Location', 'NorthWest');
%  grid on;
%  axesobj = findobj('type', 'axes');
%  set(axesobj, 'fontname', 'times'); set(axesobj, 'fontunits', 'points');
%  set(axesobj, 'fontsize', 18); set(axesobj, 'fontweight', 'bold');
%  set(axesobj, 'linewidth', 2);
%  xlabel('h');
%  ylabel('L^2(\Omega) relative error');
%  set(gca,'XTickLabel',{'0.001';'0.01';'0.1';'1';'...'}) 
%  % axis([0.05 7 1e-4 10]);
%  pause;

 

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

% axesobj = findobj('type', 'axes'); set(axesobj, 'fontname', 'times'); set(axesobj, 'fontunits', 'points'); set(axesobj, 'fontsize', 18); set(axesobj, 'fontweight', 'bold'); set(axesobj, 'linewidth', 2);


% Pour certains graphiques, il vaut mieux renommer les "ticks" par
%  set(gca,'XTickLabel',{'0.1';'1';'10';'...'})
%  set(gca,'YTickLabel',{'0.0001%';'0.001%';'0.01%';'0.1%';'1%';'10%'})     


% Pour sortir le graphique en png, faire par exemple :
% print(gcf,'-dpng','-r450', 'toto.png');

