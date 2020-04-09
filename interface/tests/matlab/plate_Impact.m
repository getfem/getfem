% Copyright (C) 2011-2020 Cedric POZZOLINI
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
%
%
% Newmark-Dumont-Paoli for a Kirchhoff-Love plate in dynamics with
% obstacles.
%

clear all
gf_workspace('clear all');
NX=2; NY=2;
longX= 0.4;
longY= 1.2;
deltaX= longX/NX;
deltaY= longY/NY;
tic
%%%%%%% Create a simple cartesian mesh
m=gf_mesh('regular simplices',0:deltaX:longX,0:deltaY:longY);
nddl = 3*(NX+1)*(NY+1);
nbnoeud=(NX+1)*(NY+1);
nelt=(NX)*(NY);

%%%%%%% Physical parameters
alpha =  0;  %1e-5;%%%%%%coeff ammortissement
beta  = 1/2; %%%%%% beta Newmark parameter - pas d'ammortissement en masse
Thickness = 0.01; % Plate thickness
%E     = 6.9*10^(10); % module young
E     = 21*10^(10); % module young acier
rho   = 7770; % densit??? acier 
%rho   = 5700; % densit??? alu
NU=0.3;% coeff poisson alu
D = (E*(Thickness)^3)/(12*(1-NU^2)*rho*Thickness);%% rigidit??? flexion/masse*epaisseur
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
amp = 0; %  amp:    amplitude of excitation force
dt=2*1e-5;% pas de temps
tmax= 3; % dur???e max
Nmax=tmax/dt; % nombre de pas de temps
omega = 10;% pulsation
e=0;


%%%%%%% FEMs and integration methods
 
mf = gf_mesh_fem(m,1); %%%%%% create a mesh_fem of for a field of dimension 1 (i.e. a scalar field)

mim=gfMeshIm(m);  %%%%%%%% hold integration methods over a mesh
mfu=gfMeshFem(m); %%%%%%% mesh for the main unknow 
mfd=gfMeshFem(m); %%%%%%% mesh for the data
mfred=gfMeshFem(m); %%%%%%% mesh for the velocity
 
set(mim, 'integ',gfInteg('IM_TRIANGLE(10)'));
set(mfu, 'fem',gfFem('FEM_ARGYRIS')); %%%% for the main unknow 
set(mfd, 'fem',gfFem('FEM_PK(2,3)')); %%%%%%%%%%% for the data
set(mfred, 'fem',gfFem('FEM_PK(2,0)')); %%%%%%%%%%% for the velocity

nddl = gf_mesh_fem_get(mfu,'nbdof');
nddllag = gf_mesh_fem_get(mfd,'nbdof');

CoordMesh = gf_mesh_get(m, 'pts'); %%%Return the list of point coordinates of the mesh
% Coordonn???es des noeuds du maillage
X=CoordMesh(1,:);
Y=CoordMesh(2,:);
nelt = gf_mesh_get(m,'nbcvs');
    




%%%%%%%%%%% boundaries and normal vector definition 
% flst = get(m, 'outer_faces');
% normale = get(m, 'normal of faces', flst);
border = gf_mesh_get(m,'outer faces');
normale = gf_mesh_get(m, 'normal of faces', border);
ftop     = border(:,find(abs(normale(1,:)-1) < 1e-5));
fbottom  = border(:,find(abs(normale(1,:)+1) < 1e-5));
fleft    = border(:,find(abs(normale(2,:)+1) < 1e-5));
fright   = border(:,find(abs(normale(2,:)-1) < 1e-5));
CLAMPED_BOUNDARY_NUM = 2;
gf_mesh_set(m, 'region', CLAMPED_BOUNDARY_NUM, [fleft]);

figure
hold on
gf_mesh_set(m, 'region', 42, [fleft]); %%%%%%,ftop, fbottom% create the region #42
gf_plot_mesh(m, 'regions', [42]); % the boundary edges appears in red
gf_mesh_set(m, 'region', 43, [fright]); %%%%%%,ftop, fbottom% create the region #42
gf_plot_mesh(m, 'regions', [43]); % the boundary edges appears in red
gf_plot_mesh(mfu, 'vertices', 'on', 'convexes', 'on', 'dof','on');
% gf_mesh_set(m, 'region', 42, border); 
% gf_plot_mesh(m, 'regions', [42]); % the boundary edges appears in red 
%nbre = gf_mesh_fem_get( mfu, 'nbdof'); %%% Return the number of degrees of freedom (dof) of the mesh_fem
[Points,INDx]=gf_mesh_get(m, 'pid from cvid'); %%% Return the number of nodes attached to each convex
[DOFs, IDx] = gf_mesh_fem_get(mfu, 'basic dof from cvid'); %%% Return the degrees of freedom attached to each convex of the mesh
%%%IDx is a row vector, length(IDx) = length(CVids)+1. DOFs
%%%is a row vector containing the concatenated list of dof of each convex in CVids. 
%%%Each entry of IDx is the position of the corresponding convex point list in DOFs. Hence, for example, 
%%%the list of points of the second convex is DOFs(IDx(2):IDx(3)-1).
 
dofdisp=zeros(1,nddl);
for j= 1:nelt
    dofdisp(DOFs(21*(j-1)+1))=1;
    dofdisp(DOFs(21*(j-1)+7))=1;
    dofdisp(DOFs(21*(j-1)+13))=1;
end

DOFsbdlibre = gf_mesh_fem_get(mfu, 'basic dof on region',43);
DOFsbd = gf_mesh_fem_get(mfu, 'basic dof on region',42);

  ddlcoingche=DOFsbdlibre(1);
  ddlcoindt= DOFsbdlibre(length(DOFsbdlibre)-6);

%%%%%%%%% Assembly of the mass matrix 
%%%%%%%%% 
%Mass=gf_asm('mass matrix', mim, mfu , mfu); % build the mass matrix

MassB=gf_asm('mass matrix', mim, mfu , mfred); % build the mass matrixB
MassC=gf_asm('mass matrix', mim, mfred , mfred); % build the mass matrixB
Mass= MassB*inv(MassC)*MassB';

MassProjLH=gf_asm('mass matrix', mim, mfu , mfd); % build the mass matrix proj lagrange/ hermite

Massreg=gf_asm('mass matrix', mim, mfu , mfu); % build the mass matrix
Kass=  gf_asm('bilaplacian KL', mim, mfu, mfd, D*ones(1, nddllag), NU*ones(1, nddllag));
%%%%%%%%% Assembly of the matrix for bilaplacian problem

 

% for i=DOFs
%     if dofdisp(i)==1
%   
%       Kass(i,:)=zeros(1,nddl);
%       Kass(:,i)=zeros(nddl,1);
%       Kass(i,i)=1;
%       Mass(i,:)=zeros(1,nddl);
%       Mass(:,i)=zeros(nddl,1);
%       Mass(i,i)=1;
%     %  Force_dom (i)=0;
%     %  Force_dom_ini (i)=0;
%     %  Force_sin (i)=0;
%       
%     end
%     
% end
% 
% A = Mass +(((dt)^2)*beta +alpha*((dt)/2))*Kass;

%% Static problem solve


md=gf_model('real');
gf_model_set(md, 'add fem variable', 'u', mfu);

gf_model_set(md, 'add initialized data', 'D', [D*dt^2*beta]);
gf_model_set(md, 'add initialized data', 'nu', [NU]);
gf_model_set(md, 'add Kirchhoff-Love plate brick', mim, 'u', 'D', 'nu');


Force_dom_ini = dt^2*beta * get(mfd, 'eval', {'8600'});
Force_dom_ini=Force_dom_ini/(rho*Thickness);
gf_model_set(md, 'add initialized fem data', 'VolumicData', mfd, Force_dom_ini);

gf_model_set(md, 'add source term brick', mim, 'u', 'VolumicData');
 
gf_model_set(md, ...
 	     'add normal derivative Dirichlet condition with penalization', ...
 	       mim, 'u', 1e10, CLAMPED_BOUNDARY_NUM);
 
gf_model_set(md, 'add Dirichlet condition with penalization', ...
    	     mim, 'u', 1e10, CLAMPED_BOUNDARY_NUM);

gf_model_get(md, 'solve', 'noisy', 'max_res', 1e-14);
U1 = (gf_model_get(md, 'variable', 'u'))';

%gf_plot(mfu, U1', 'mesh','on');

colorbar;
gf_plot(mfu, U1', 'zplot', 'on', 'deformed_mesh','on');
%pause;


%% Dynamic problem solve



gf_model_set(md, 'variable', 'VolumicData', get(mfd, 'eval', {'0'}));
gf_model_set(md, 'add explicit matrix', 'u', 'u', Mass);
ind_rhs = gf_model_set(md, 'add explicit rhs', 'u', zeros(nddl, 1)); % -1 ?? enlever ... bug interface

%%% obstacle definition


  
B = sparse(0,nddl);
Val = zeros(nddl,1);
gap = zeros(0, 0);
nbconstraints = 0;

for i=DOFs
    if dofdisp(i)==1 && Val(i) == 0
        Val(i) = 1;
        nbconstraints = nbconstraints + 1;
        B(nbconstraints, i) = -1;
        gap(nbconstraints, 1) = 0.1;
    end
end
obstacle2=zeros(nddl,1);
    for j=1:nddl
        obstacle2(j)=inf;
    end
for i=DOFs%%%DOFsbd(1):DOFsbd(lenght(DOFsbd))%%%%%%%%%  "|" est le ou 
  
    if dofdisp(i)==1 %| (Y(inoeud(j))==3)
      
              obstacle2(i)= 0.05;%%%%%% 1xamp aux noeuds fleche ??? la base
    end
    if dofdisp(i)==0 %| (Y(inoeud(j))==3)
      
              obstacle2(i)= inf ;
    end
end



gf_model_set(md, 'add variable', 'lambda', nbconstraints);
gf_model_set(md, 'add initialized data', 'r', [1]);
gf_model_set(md, 'add initialized data', 'gap', gap);
gf_model_set(md, 'add initialized data', 'alpha', ones(nbconstraints, 1));
gf_model_set(md, 'add basic contact brick', 'u', 'lambda', 'r', B, 'gap', 'alpha', 1);

%Ud1= gf_mesh_fem_get(mfd, 'eval', { '0 ' })';
    %Ud1= gf_mesh_fem_get(mfd, 'eval', { ' 6*10^(-2)*(y.^2)' })'; %5*(y.^2)*10^(-2)%%%%%% deplacement impose initial
    %U1= gf_compute(mf, Ud1', 'extrapolate on', mfu)'; 
    %U1= inv(Mass)*MassProjLH*Ud1;
 %      U1 = Kass\Force_dom_ini;
Ud2= gf_mesh_fem_get(mfd, 'eval', { '0 ' })';
U2 = U1 - dt* Massreg \ (MassProjLH*Ud2);
% Force_sin = zeros(nddl,1);
%Force_sin = (omega^2)*Mass*U_nstat-Kass*U_nstat;
U3=U2;
%%%%%%%%%%%%%%%%%%%%% Declare source term (Force)
Force = gf_mesh_fem_get(mfd, 'eval', { '0' }); %5*(x.^2+y.^2)*
%%%%%%%%%%%%%%%%%%%%% Assembly of volumic source term (Force)
Force_dom = gf_asm('volumic source', mim, mfu, mfd, Force);


%%%%%%%%%%%%%%%%%%%%% Declare source term intial (Force)
%Force_ini = gf_mesh_fem_get(mfd, 'eval', { '9600' });%8600 %5*(x.^2+y.^2)*
%%%%%%%%%%%%%%%%%%%%% Assembly of volumic source term initial (Force)
%Force_dom_ini = gf_asm('volumic source', mim, mfu, mfd, Force_ini);
%Force_dom_ini=Force_dom_ini/(rho*Thickness);

for t = dt:dt:tmax
    
    n=round(t/dt);
    
    
   %G_n=Force_sin*sin(omega*t) +(Force_dom)/(rho*Thickness);
       
      
   F_n=(2*Mass -((dt)^2)*(1-2*beta)*Kass)*U2 - (Mass+(((dt)^2)*beta-alpha*((dt)/2))*Kass)*U1;%+((dt)^2)*(G_n);
   gf_model_set(md, 'set private rhs', ind_rhs, F_n);

   
  t
  gf_model_get(md, 'solve', 'noisy','max_res', 1e-14);%
  Q_n =(gf_model_get(md, 'variable', 'u') )';%
 % U3= A\F_n;%

  % gf_plot(mfu, U3','mesh','on');
   %colorbar;
   %pause;
   



%%%%%%%%%%%%%%%%%%%%%%% Test de la contrainte convexe 
%   Q_n= (U3(:)+e*U1(:))/(1+e);
%     
%   contact=0;
% for i=DOFs %&& e >=0
%     if dofdisp(i)==1 
%         if Q_n(i)<=-obstacle2(i)
%            
%            contact=1
%        end
%     end    
% end
%        
%       
%        
%        
%        
%    
%  if contact==1 %&& e>=0
% % 
%         F_ne=(F_n+ e*A*U1(:))/(1+e);
%         
%          gf_model_set(md, 'set private rhs', ind_rhs, F_ne);
%           
%    gf_model_get(md, 'solve', 'noisy', 'max_res', 1e-11);
%    Q_n = (gf_model_get(md, 'variable', 'u'))';
%     
%  end  
%   
   U3(:)= (1+e)*Q_n(:)-e*U1(:); 
%         
   


     U1=U2;
     U2=U3;
     % format long; U3(ddlcoindt)
     Ucoind(n)=U3(ddlcoindt);
     Ucoing(n)=U3(ddlcoingche);  
    
     %Ucentre(n)=U1(502);
  
  
       %%%% Evaluation de l'???nergie totale classique avec G_n
        if n > 1
    %ETOT(n)    ETOT(n)=(0.125/dt^2)*(U3(:)-U1(:))'*Mass*(U3(:)-U1(:)) + (0.5)*(U2(:))'*Kass*(U2(:));
    ETOT(n)=(0.5/dt^2)*(U2(:)-U1(:))'*Mass*(U2(:)-U1(:)) + (0.5)*(U1(:))'*Kass*(U1(:));% - (G_n)*U2(:);  
    
    %ETOT(n)=(0.125/dt^2)*(U(:,n+1)-U(:,n-1))'*Mg*(U(:,n+1)-U(:,n-1)) + (0.5)*(U(:,n))'*Kg*(U(:,n)) - (G_n)'*U(:,n);
      %%%% Evaluation de l'???nergie totale forme DP - 2008 splines 
     % ETOTDP(n)=0.5*((U(:,n)-U(:,n-1))'*((1/dt^2)*Mg + beta*Kg)*(U(:,n)-U(:,n-1))  + (U(:,n))'*Kg*(U(:,n-1))) - (G_n)'*U(:,n) ;
    % ETOT(n)=(0.125/dt^2)*(U3(:)-U1(:))'*Mass*(U3(:)-U1(:)) + (0.5)*(U2(:))'*Kass*(U2(:)) - (G_n)'*U2(:);
 
     ETOTDP(n)=0.5*((U2(:)-U1(:))'*((1/dt^2)*Mass + beta*Kass)*(U2(:)-U1(:))  + (U2(:))'*Kass*(U1(:))); % - (G_n)'*U2(:) ;
        end
%%%%%%%%% Static Equilibrium  
%U = inv(Kass)*Force_dom;
%%%%%%%%%%%%% plot  dynamic
%figure
  %gf_plot(mfu, U3', 'zplot', 'on', 'deformed_mesh','on');
 %axis([0 longX 0 longY -0.1 0.1]);
 %%%caxis auto   
  %caxis([-0.1 0.1]);
  %colorbar;
 %%% mov = avifile('NDP_Impact_plaque.avi');
   %Fr(n) = getframe(gcf);
   %%%caxis auto
   %%hold on
   %%%% gf_plot(mfu, Obstacle2, 'zplot', 'on', 'mesh','on');
   % drawnow; hold off %pause(.01)
         %%%  mov = addframe(mov,Fr);
         
         
%-------------------------------------------------------------------------%
%--- /!\ ici on sauve dans le GIF ----------------------------------------%
%-------------------------------------------------------------------------%
       % [RGB,badmap] = frame2im(Fr(n)); %on la convertie en image de type 'true-color'
       % [IND,map] = rgb2ind(RGB, 255); %on convertie en couleur ind???x???es. 255 est le nombre de couleur.
       % if isfirst
        %    imwrite(IND,map,'NDP_Impact_plaque.gif','gif','LoopCount',100); %---- premi???re image du fichier GIF
        %    isfirst=false;
        %else
        %    imwrite(IND,map,'NDP_Impact_plaque.gif','gif','WriteMode','append','DelayTime',0.09); %---- les images suivantes
        %end
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
end

toc

%%%%movie2avi(Fr,'NDP_Impact_plaque.avi','compression','None')
%close(gcf)          %---- fermeture du handle figure
%%%%%%mov = close(mov);   %---- fermeture du handle vid???o


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
%%%%%%%%%%% Trac??? des r???sultats
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
   t = linspace(0,tmax,Nmax); % variable TEMPS
   plot(t,Ucoing(:),'b');
  hold on
   t = linspace(0,tmax,Nmax); % variable TEMPS
   plot(t,Ucoind(:),'--black');
   hold on
   %   plot(t,Ucentre(:),'r');
   title('Displacement of the free corners of a plate impacting flat obstacles : Newmark-Dumont-Paoli Sing. Argyris method beta=1/2 ');
   xlabel('time dt= 10^-3'),ylabel('disp.');
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
  %t = linspace(0,tmax,Nmax); % variable TEMPS
   plot(t,ETOT(:),'--black');
  % xlabel('time dt= 10^-3'),ylabel('Total Energy');
    t = linspace(0,tmax,Nmax); % variable TEMPS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hold on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
title('Energy of a plate impacting pointwise obstacles : Newmark-Dumont-Paoli method, Argyris sing. 96 elts, beta=1/2 ');
   plot(t,ETOTDP(:),'r');
   xlabel('time dt= 10^-3'),ylabel('Total Energy');



return;








BilA= D*ones(1,nddllag);
%gf_mesh_fem_get(mf, 'eval', { '1 ' });
poisson= Poisson*ones(1,nddllag);
Kass=  gf_asm('bilaplacian KL', mim, mfu, mfd, BilA, poisson);

%%%%%%%%%%%%%%%%%%%%% Declare source term (Force)
Force = gf_mesh_fem_get(mfd, 'eval', { '0' }); %5*(x.^2+y.^2)*
%%%%%%%%%%%%%%%%%%%%% Assembly of volumic source term (Force)
Force_dom = gf_asm('volumic source', mim, mfu, mfd, Force);


%%%%%%%%%%%%%%%%%%%%% Declare source term intial (Force)
Force_ini = gf_mesh_fem_get(mfd, 'eval', { '8600' });%8600 %5*(x.^2+y.^2)*
%%%%%%%%%%%%%%%%%%%%% Assembly of volumic source term initial (Force)
Force_dom_ini = gf_asm('volumic source', mim, mfu, mfd, Force_ini);
Force_dom_ini=Force_dom_ini/(rho*Thickness);

%%%%%%%%%%%%%%%%%%%%%%%%
% Conditions aux limites sur les matrices de masse et de rigidit???
%%%%%%%%%%%%%%%%%%%%%%%%


 

%%%%%%%%%%% boundaries and normal vector definition 
% flst = get(m, 'outer_faces');
% normale = get(m, 'normal of faces', flst);
border = gf_mesh_get(m,'outer faces');
normale = gf_mesh_get(m, 'normal of faces', border);
ftop     = border(:,find(abs(normale(1,:)-1) < 1e-5));
fbottom  = border(:,find(abs(normale(1,:)+1) < 1e-5));
fleft    = border(:,find(abs(normale(2,:)+1) < 1e-5));
fright   = border(:,find(abs(normale(2,:)-1) < 1e-5));
figure
hold on
gf_mesh_set(m, 'region', 42, [fleft]); %%%%%%,ftop, fbottom% create the region #42
gf_plot_mesh(m, 'regions', [42]); % the boundary edges appears in red
gf_mesh_set(m, 'region', 43, [fright]); %%%%%%,ftop, fbottom% create the region #42
gf_plot_mesh(m, 'regions', [43]); % the boundary edges appears in red
gf_plot_mesh(mfu, 'vertices', 'on', 'convexes', 'on', 'dof','on');
% gf_mesh_set(m, 'region', 42, border); 
% gf_plot_mesh(m, 'regions', [42]); % the boundary edges appears in red 
%nbre = gf_mesh_fem_get( mfu, 'nbdof'); %%% Return the number of degrees of freedom (dof) of the mesh_fem
[Points,INDx]=gf_mesh_get(m, 'pid from cvid'); %%% Return the number of nodes attached to each convex
[DOFs, IDx] = gf_mesh_fem_get(mfu, 'basic dof from cvid'); %%% Return the degrees of freedom attached to each convex of the mesh
%%%IDx is a row vector, length(IDx) = length(CVids)+1. DOFs
%%%is a row vector containing the concatenated list of dof of each convex in CVids. 
%%%Each entry of IDx is the position of the corresponding convex point list in DOFs. Hence, for example, 
%%%the list of points of the second convex is DOFs(IDx(2):IDx(3)-1).
 
%%%%%%% selection des ddl bord libre / bord encastre

DOFsbdlibre = gf_mesh_fem_get(mfu, 'basic dof on region',43);
DOFsbd = gf_mesh_fem_get(mfu, 'basic dof on region',42);

%%%%%%%%%%% tableau des ddl deplacement
dofdisp=zeros(1,nddl);
for j= 1:nelt
    dofdisp(DOFs(21*(j-1)+1))=1;
    dofdisp(DOFs(21*(j-1)+7))=1;
    dofdisp(DOFs(21*(j-1)+13))=1;
end
  

% % % 
  ddlcoingche=DOFsbdlibre(1);
  ddlcoindt= DOFsbdlibre(length(DOFsbdlibre)-6);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%% si sinus ??? la base
 U_nstat = zeros(nddl,1);
for i=DOFsbd%%%DOFsbd(1):DOFsbd(lenght(DOFsbd))%%%%%%%%%  "|" est le ou 
  
    if dofdisp(i)==1 %| (Y(inoeud(j))==3)
      
              U_nstat(i)= amp; %%%%%%% 1xamp aux noeuds fleche ??? la base
    end
end

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%% si sinus ??? la base CLNH
 Force_sin = zeros(nddl,1);
 Force_sin = (omega^2)*Mass*U_nstat-Kass*U_nstat;


for i=DOFsbd %%%%%%%%%  "|" est le ou 
   
      Kass(i,:)=zeros(1,nddl);
      Kass(:,i)=zeros(nddl,1);
      Kass(i,i)=1;
      Mass(i,:)=zeros(1,nddl);
      Mass(:,i)=zeros(nddl,1);
      Mass(i,i)=1;
      Force_dom (i)=0;
      Force_dom_ini (i)=0;
      Force_sin (i)=0;
      
    
    
end




% Boucle sur tous les noeuds de d???placement impos??? nul
% for i=1:nelt %%%%%%%%%  "|" est le ou 
%     idof=DOFs(IDx(i):IDx(i+1)-1);
%     inoeud=Points(INDx(i):INDx(i+1)-1);
%     for j=1:3
%     if (Y(inoeud(j))==0) %| (Y(inoeud(j))==3)
%         ii=3*j-3;
%         for k=1:3
%             iik=idof(ii+k);
%       Kass(iik,:)=zeros(1,nddl);
%       Kass(:,iik)=zeros(nddl,1);
%       Kass(iik,iik)=1;
%       Mass(iik,:)=zeros(1,nddl);
%       Mass(:,iik)=zeros(nddl,1);
%       Mass(iik,iik)=1;
%       Force_dom (iik)=0;
%       Force_dom_ini (iik)=0;
%       Force_sin (iik)=0;
%         end
%     end
%     end
% end



%%%%%%%%%%%%%%%%%%%%% Declare obstacles 

obstacle2=zeros(nddl,1);
    for j=1:nddl
        obstacle2(j)=inf;
    end

for i=DOFs%%%DOFsbd(1):DOFsbd(lenght(DOFsbd))%%%%%%%%%  "|" est le ou 
  
    if dofdisp(i)==1 %| (Y(inoeud(j))==3)
      
              obstacle2(i)= inf ; %0.1%%%%%% 1xamp aux noeuds fleche ??? la base
    end
    if dofdisp(i)==0 %| (Y(inoeud(j))==3)
      
              obstacle2(i)= inf ; %%%%%%% 1xamp aux noeuds fleche ??? la base
    end
end

    
    gap=0.1;%
    
    obstacle2(ddlcoindt)=gap;
    obstacle2(ddlcoingche)=gap; 
 
 %%%% Declare that u is an unknown of the system on the finite element method 
 md=gf_model('real'); %%%%%%%%%% Declare a real unknown
 gf_model_set(md, 'add fem variable', 'U1', mfu);
 
 gf_model_set(md, 'add fem variable', 'U2', mfu);
 
 gf_model_set(md, 'add fem variable', 'U3', mfu);
 
 
%%%%% Find the boundary of the domain, in order to set a Dirichlet condition. 
%figure
%border = gf_mesh_get(m,'outer faces');
% gf_mesh_set(m, 'region', 42, border); % create the region #42
% gf_plot_mesh(m, 'regions', [42]); % the boundary edges appears in red
 
 
 %%%%%% Dirichlet condition on the domain boundary



 %%%% Def A 1er membre de (P^{n+1}_{jbeta}) definie positive sym
    A = Mass +(((dt)^2)*beta +alpha*((dt)/2))*Kass;
    B=inv(A);
    Ad=A;%zeros(nddl,nddl);
    Ad(:,ddlcoindt)=0;
    Ad(ddlcoindt,:)=0;
    Ad(ddlcoindt,ddlcoindt)= 1;
    Ad(:,ddlcoingche)=0;
    Ad(ddlcoingche,:)=0;
    Ad(ddlcoingche,ddlcoingche)= 1;

   Bd=inv(Ad);
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Initialization de U_1 et U_2 (CI) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%U=zeros(nddl,Nmax);

  Ud1= gf_mesh_fem_get(mfd, 'eval', { '0 ' })';
  %  Ud1= gf_mesh_fem_get(mfd, 'eval', { ' 6.9*10^(-2)*(y.^2)' })'; %%%%%%% deplacement impose initial
    %U1= gf_compute(mfd, Ud1', 'extrapolate on', mfu)'; 
   %U1= inv(Mass)*MassProjLH*Ud1;
       U1 = Kass\Force_dom_ini;
  Ud2= gf_mesh_fem_get(mfd, 'eval', { '0 ' })';
  %Ud2=gf_mesh_fem_get(mfd, 'eval', { ' y*5*10^(-1)  ' })'; 
    %Ud2= Ud2;%%%%%%% deplacement impose initial 2ieme instant
    %U2= U1-dt*gf_compute(mfd, Ud2', 'interpolate on', mfu)';
  U2= U1-dt*inv(Massreg)*MassProjLH*Ud2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure
% gf_plot(mfu, U1', 'zplot', 'on', 'mesh','on');
% colorbar;
% Fr(1) = getframe;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%  energie totale clasique
ETOT=zeros(1,Nmax);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ETOTDP=zeros(1,Nmax);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%-------------------------------------------------------------------------%
%---- initialisation vid???o et GIF ----------------------------------------%
%-------------------------------------------------------------------------%
%TAILLE  = [150 100 800 800];
%TAILLE2 = [0.0 0.0 1.0 1.0];
%f       = 15;
%a       = 2;   
%mov     = avifile('evolution.avi','compression','none','Quality',100);%,,'fps',f);
% isfirst = true;     %---- variable d'initialisation du GIF !!!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%      DEBUT BOUCLE EN TEMPS        %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
G_n = zeros(nddl,1);


    Ucoind=zeros(1,Nmax);
    Ucoind(1)=U1(ddlcoindt);
    
    Ucoing=zeros(1,Nmax);
    Ucoing(1)=U1(ddlcoingche);
tic
for n=2:Nmax; 

    
%%%%%%%%%%%%%%%%%%% si sinus ??? la base + force dom + dep initial impos???
G_n=Force_sin*sin(omega*n*dt) +(Force_dom)/(rho*Thickness);
       
      
F_n=(2*Mass -((dt)^2)*(1-2*beta)*Kass)*U2 - (Mass+(((dt)^2)*beta-alpha*((dt)/2))*Kass)*U1+((dt)^2)*(G_n);
      
     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% prediction sol bilaterale  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% D???placement sans obstacles %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
U3= B*F_n;
%U3(:)= quadprog(A,-F_n,[],[],[],[],-gap,gap);
%%%%%%%%%%%%%%%%%%%%%%% Test de la contrainte convexe 
  Q_n= (U3(:)+e*U1(:))/(1+e);
    
  contact=0;
% %      contact=0;
% %     
% %     for i=DOFs%%%DOFsbd(1):DOFsbd(lenght(DOFsbd))%%%%%%%%%  "|" est le ou 
% %   
% %     if dofdisp(i)==1 %| (Y(inoeud(j))==3)
% %        if Q_n(i)<=-obstacle2(i)
% %            
% %            contact=1;
% %        end
% %     end    
% %        
% %        
% %    end
       
       if Q_n(ddlcoindt)<=-gap
 
           
           %nddlcontact= ddlcoindt;
     % %%%%%%%%%%Correction
 % %%%%%%%%%%%%%%%%%%%%%%%%  QUADPROG LQP
       contact=1;
       
       end
       
        if Q_n(ddlcoingche)<=-gap
 
           
          % nddlcontact= ddlcoingche;
     % %%%%%%%%%%Correction
 % %%%%%%%%%%%%%%%%%%%%%%%%  QUADPROG LQP
       contact=1;
       
        end
       
        % if Q_n(ddlcoingche)>=gap
 
           
         %  nddlcontact= ddlcoingche;
     % %%%%%%%%%%Correction
 % %%%%%%%%%%%%%%%%%%%%%%%%  QUADPROG LQP
       %contact=2;
       
       %end
  
   
 if contact==1;
% 
        F_ne=(F_n+ e*A*U1(:))/(1+e);
        Ud=zeros(nddl,1);
        Ud(ddlcoindt,1)= -gap;
        Ud(ddlcoingche,1)= -gap;
        Q_n(:)= Bd*(F_ne-A*Ud);
        Q_n(ddlcoindt) = -gap;
        Q_n(ddlcoingche) = -gap;
        U3(:)= (1+e)*Q_n(:)-e*U1(:); 
        
            Q_n(:)=(U3(:)+e*U1(:))/(1+e);%U3(:);
 end  
  
  
     
%   if contact==1;
% 
%     
%     
%     
%     %%%%%%%%%%Correction
% %%%%%%%%%%%%%%%%%%%%%%%%  QUADPROG LQP
%    F_ne= (F_n+ e*A*U1(:))/(1+e);
%    
%            options=optimset('LargeScale','on','PrecondBandWidth',Inf,'Display','final','TolX',1e-14,'MaxIter',Inf,'TolFun',1e-14);%
% 
%    Q_n(:)= quadprog(A,-F_ne,[],[],[],[],-obstacle2,[],U1(:),options);
%    U3(:)= (1+e)*Q_n(:)-e*U1(:);  
%end

%Q_n(:)=(U3(:)+e*U1(:))/(1+e);


 for i=DOFsbd%%%DOFsbd(1):DOFsbd(lenght(DOFsbd))%%%%%%%%%  "|" est le ou 
  
    if dofdisp(i)==1 %| (Y(inoeud(j))==3)
      
              U3(i)= amp*sin(omega*n*dt) ; %%%%%%% 1xamp aux noeuds fleche ??? la base
    end
    if dofdisp(i)==0 %| (Y(inoeud(j))==3)
      
              U3(i)= 0 ; %%%%%%% 1xamp aux noeuds fleche ??? la base
    end
end

     U1=U2;
     U2=U3;
     Ucoind(n)=Q_n(ddlcoindt);
     Ucoing(n)=Q_n(ddlcoingche);  
    
     Ucentre(n)=Q_n(502);
  
  
       %%%% Evaluation de l'???nergie totale classique avec G_n
        if n > 1
    ETOT(n)=(0.125/dt^2)*(U3(:)-U1(:))'*Mass*(U3(:)-U1(:)) + (0.5)*(U2(:))'*Kass*(U2(:)) - (G_n)*U2(:);    
    %ETOT(n)=(0.125/dt^2)*(U(:,n+1)-U(:,n-1))'*Mg*(U(:,n+1)-U(:,n-1)) + (0.5)*(U(:,n))'*Kg*(U(:,n)) - (G_n)'*U(:,n);
      %%%% Evaluation de l'???nergie totale forme DP - 2008 splines 
     % ETOTDP(n)=0.5*((U(:,n)-U(:,n-1))'*((1/dt^2)*Mg + beta*Kg)*(U(:,n)-U(:,n-1))  + (U(:,n))'*Kg*(U(:,n-1))) - (G_n)'*U(:,n) ;
    % ETOT(n)=(0.125/dt^2)*(U3(:)-U1(:))'*Mass*(U3(:)-U1(:)) + (0.5)*(U2(:))'*Kass*(U2(:)) - (G_n)'*U2(:);
 
     ETOTDP(n)=0.5*((U2(:)-U1(:))'*((1/dt^2)*Mass + beta*Kass)*(U2(:)-U1(:))  + (U2(:))'*Kass*(U1(:))) - (G_n)'*U2(:) ;
        end
%%%%%%%%% Static Equilibrium  
%U = inv(Kass)*Force_dom;
%%%%%%%%%%%%% plot  dynamic
  figure
 gf_plot(mfu, U3', 'zplot', 'on', 'deformed_mesh','on');
% axis([0 longX 0 longY -0.1 0.1]);
 %%%caxis auto   
 caxis([-0.1 0.1]);
  %colorbar;
 %%% mov = avifile('NDP_Impact_plaque.avi');
   %Fr(n) = getframe(gcf);
   %%%caxis auto
   %%hold on
   %%%% gf_plot(mfu, Obstacle2, 'zplot', 'on', 'mesh','on');
   drawnow; hold off %pause(.01)
         %%%  mov = addframe(mov,Fr);
         
         
%-------------------------------------------------------------------------%
%--- /!\ ici on sauve dans le GIF ----------------------------------------%
%-------------------------------------------------------------------------%
       % [RGB,badmap] = frame2im(Fr(n)); %on la convertie en image de type 'true-color'
       % [IND,map] = rgb2ind(RGB, 255); %on convertie en couleur ind???x???es. 255 est le nombre de couleur.
       % if isfirst
        %    imwrite(IND,map,'NDP_Impact_plaque.gif','gif','LoopCount',100); %---- premi???re image du fichier GIF
        %    isfirst=false;
        %else
        %    imwrite(IND,map,'NDP_Impact_plaque.gif','gif','WriteMode','append','DelayTime',0.09); %---- les images suivantes
        %end
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
end

toc

%%%%movie2avi(Fr,'NDP_Impact_plaque.avi','compression','None')
%close(gcf)          %---- fermeture du handle figure
%%%%%%mov = close(mov);   %---- fermeture du handle vid???o


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
%%%%%%%%%%% Trac??? des r???sultats
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
   t = linspace(0,tmax,Nmax); % variable TEMPS
   plot(t,Ucoing(:),'b');
  hold on
   t = linspace(0,tmax,Nmax); % variable TEMPS
   plot(t,Ucoind(:),'--black');
   hold on
      plot(t,Ucentre(:),'r');
   title('Displacement of the free corners of a plate impacting flat obstacles : Newmark-Dumont-Paoli Sing. Argyris method beta=1/2 ');
   xlabel('time dt= 10^-3'),ylabel('disp.');
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
   t = linspace(0,tmax,Nmax); % variable TEMPS
   plot(t,ETOT(:),'--black');
   xlabel('time dt= 10^-3'),ylabel('Total Energy');
    t = linspace(0,tmax,Nmax); % variable TEMPS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hold on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
title('Energy of a plate impacting pointwise obstacles : Newmark-Dumont-Paoli method, Argyris sing. 96 elts, beta=1/2 ');
   plot(t,ETOTDP(:),'r');
   xlabel('time dt= 10^-3'),ylabel('Total Energy');
