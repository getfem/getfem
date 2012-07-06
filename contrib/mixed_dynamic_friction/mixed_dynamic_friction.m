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


% addpath ~/source++/getfem++/contrib/mixed_dynamic_friction

%
% current
%

if (0)

A = load('mixed_dynamic_friction.data');

% energy curves

plot(A(:, 1), A(:, 2), '-k', 'linewidth', 2, 'MarkerSize', 15);
% axis([0 0.02 0 0.1]);
pause;

% displacement curves

plot(A(:, 1), A(:, 3), '-k', 'linewidth', 2, 'MarkerSize', 15);
% axis([0 0.02 -0.01 0.06]);
pause;

% contact stress curves

plot(A(:, 1), A(:, 4), '-k', 'linewidth', 2, 'MarkerSize', 15);
% axis([0 0.02 -0.4 0.02]);
pause;

else


%
% P1plusP0
%

A1 = load('mdfconv_P2P1_0.5.data');
A2 = load('mdfconv_P2P1_0.25.data');
A3 = load('mdfconv_P2P1_0.075.data');

% energy curves

plot(A1(:, 1), A1(:, 2), '-k', 'linewidth', 2, 'MarkerSize', 15);
hold on;
plot(A2(:, 1), A2(:, 2), '--k', 'linewidth', 2, 'MarkerSize', 15);
plot(A3(:, 1), A3(:, 2), '-.k', 'linewidth', 2, 'MarkerSize', 15);
hold off;
axis([0 0.02 730 770]);
% axis([0 0.7 0 0.02]);
xlabel('t');
ylabel('total energy');
legend('dt = 2\times10^{-4}', 'dt = 10^{-4}', 'dt = 2.5\times10^{-5}', 'Location', 'SouthWest');
axesobj = findobj('type', 'axes'); set(axesobj, 'fontname', 'times'); set(axesobj, 'fontunits', 'points'); set(axesobj, 'fontsize', 24); set(axesobj, 'fontweight', 'bold'); set(axesobj, 'linewidth', 2);
pause;
print(gcf,'-deps','-r450', 'energy.eps');
% print(gcf,'-dpng','-r450', 'energy.png');

% contact stress curves

plot(A1(:, 1), A1(:, 3)*2./0.5, '-k', 'linewidth', 2, 'MarkerSize', 15);
hold on;
plot(A2(:, 1), A2(:, 3)*2./0.25 -50, '--k', 'linewidth', 2, 'MarkerSize', 15);
plot(A3(:, 1), A3(:, 3)*2./0.075 -100, '-.k', 'linewidth', 2, 'MarkerSize', 15);
hold off;
axis([0 0.02 -400 0.01]);
% axis([0 0.7 -0.01 0.06]);
xlabel('t');
ylabel('point A contact stress');
legend('dt = 2\times10^{-4}', 'dt = 10^{-4}', 'dt = 2.5\times10^{-5}', 'Location', 'SouthWest');
axesobj = findobj('type', 'axes'); set(axesobj, 'fontname', 'times'); set(axesobj, 'fontunits', 'points'); set(axesobj, 'fontsize', 24); set(axesobj, 'fontweight', 'bold'); set(axesobj, 'linewidth', 2);
pause;
print(gcf,'-deps','-r450', 'stress.eps');
% print(gcf,'-dpng','-r450', 'stress.png');

% displacement curves

plot(A1(:, 1), A1(:, 4), '-k', 'linewidth', 2, 'MarkerSize', 15);
hold on;
plot(A2(:, 1), A2(:, 4), '--k', 'linewidth', 2, 'MarkerSize', 15);
plot(A3(:, 1), A3(:, 4), '-.k', 'linewidth', 2, 'MarkerSize', 15);
hold off;
axis([0 0.02 -0.1 10]);
xlabel('t');
ylabel('point A normal displacement');
legend('dt = 10^{-3}', 'dt = 10^{-4}', 'dt = 10^{-5}', 'Location', 'SouthWest');
axesobj = findobj('type', 'axes'); set(axesobj, 'fontname', 'times'); set(axesobj, 'fontunits', 'points'); set(axesobj, 'fontsize', 24); set(axesobj, 'fontweight', 'bold'); set(axesobj, 'linewidth', 2);
pause;
print(gcf,'-deps','-r450', 'displacement.eps');
% print(gcf,'-dpng','-r450', 'displacement.png');



end;


if (1)

H = [5e-5, 1e-4, 2e-4, 4e-4, 8e-4];
H1 = [0.43, 1.42, 2.57, 4.41, 5.86];
loglog(H(2:5), H1(2:5), 'o-k', 'linewidth', 2, 'MarkerSize', 15);
P1 = polyfit(log(H(2:5)), log(H1(2:5)), 1);
legend(strcat('P2/P1  (slope=',num2str(P1(1)), ')'))
grid on;
axesobj = findobj('type', 'axes');
set(axesobj, 'fontname', 'times'); set(axesobj, 'fontunits', 'points');
set(axesobj, 'fontsize', 24); set(axesobj, 'fontweight', 'bold');
set(axesobj, 'linewidth', 2);
xlabel('dt');
ylabel('H^1(\Omega) relative error');
%set(gca,'XTickLabel',{'0.001';'0.01';'0.1';'1';'...'}) 
set(gca,'YTickLabel',{'1%', '10%'}) 
axis([1e-4 1e-3 0.5 10]);
pause;

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

