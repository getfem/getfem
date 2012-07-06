% Copyright (C) 2005-2012 Yves Renard, Julien Pommier.
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

gf_workspace('clear all'); clear all;
N = 2;
NX=10; NY=14;
m=gfMesh('regular simplices',0:0.4/NX:0.4, 0:1.2/NY:1.2);
% m=gfMesh('cartesian',0:1/NX:1, 0:1/NY:1);
% m=gfMesh('cartesian',0:0.4/NX:0.4, 0:1.2/NY:1.2);


useKL=1; % use the Kirchhoff-Love plate model, or just a pure
         % bilaplacian problem

D=1.0;   % Flexion modulus
newbricks=1; % Use new bricks or old ones.

if useKL, NU=0.3; end; % poisson ratio (0 <= NU <= 1)

mim=gfMeshIm(m); 
mfu=gfMeshFem(m); 
mfd=gfMeshFem(m);
set(mim, 'integ',gfInteg('IM_TRIANGLE(13)'));
set(mfu, 'fem',gfFem('FEM_ARGYRIS'));
set(mfd, 'fem',gfFem('FEM_PK(2,5)'));

% set(mim, 'integ', gfInteg('IM_GAUSS_PARALLELEPIPED(2,10)'));
% set(mfu, 'fem', gfFem('FEM_REDUCED_QUADC1_COMPOSITE'));
% set(mfd, 'fem', gfFem('FEM_QK(2,3)'));

flst = get(m, 'outer_faces');
n = get(m, 'normal of faces', flst);
ftop     = flst(:,find(abs(n(1,:)-1) < 1e-5));
fbottom  = flst(:,find(abs(n(1,:)+1) < 1e-5));
fleft    = flst(:,find(abs(n(2,:)+1) < 1e-5));
fright   = flst(:,find(abs(n(2,:)-1) < 1e-5));

FORCE_BOUNDARY_NUM=1;
MOMENTUM_BOUNDARY_NUM=2;
SIMPLE_SUPPORT_BOUNDARY_NUM=3;
CLAMPED_BOUNDARY_NUM=4;

set(m, 'region', FORCE_BOUNDARY_NUM, fright);
set(m, 'region', SIMPLE_SUPPORT_BOUNDARY_NUM, [fleft ftop fbottom]);
set(m, 'region', CLAMPED_BOUNDARY_NUM, [fleft ftop fbottom]);
set(m, 'region', MOMENTUM_BOUNDARY_NUM, [ftop fbottom]);

FT=2.;
sol_u=get(mfd, 'eval',{sprintf('sin(%g*(x+y))',FT)});
sol_f=sol_u*FT*FT*FT*FT*N*N*D;
sol_lapl_u=-FT*FT*sol_u*N;

if (newbricks) % uses new bricks
  
  md=gf_model('real');
  gf_model_set(md, 'add fem variable', 'u', mfu);


  if useKL
    gf_model_set(md, 'add initialized data', 'D', [D]);
    gf_model_set(md, 'add initialized data', 'nu', [NU]);
    gf_model_set(md, 'add Kirchhoff-Love plate brick', mim, 'u', 'D', 'nu');
    M = zeros(N,N, get(mfd,'nbdof'));
  else
    gf_model_set(md, 'add initialized data', 'D', [D]);
    gf_model_set(md, 'add bilaplacian brick', mim, 'u', 'D');
    M = zeros(1, get(mfd, 'nbdof'));
  end;

  gf_model_set(md, 'add initialized fem data', 'VolumicData', mfd, ...
     	       get(mfd, 'eval', {'1-(x-y).^2'}));

  gf_model_set(md, 'add source term brick', mim, 'u', 'VolumicData');

  
  gf_model_set(md, 'add initialized fem data', 'M', mfd, M);
  gf_model_set(md, 'add normal derivative source term brick', mim, 'u', ...
    	       'M', MOMENTUM_BOUNDARY_NUM);

  if (useKL) 
    H = zeros(N, N, get(mfd, 'nbdof'));
    F = zeros(N, get(mfd, 'nbdof'));
    gf_model_set(md, 'add initialized fem data', 'H', mfd, H);
    gf_model_set(md, 'add initialized fem data', 'F', mfd, F);
    gf_model_set(md, 'add Kirchhoff-Love Neumann term brick', mim, 'u', ...
		         'H', 'F', FORCE_BOUNDARY_NUM);
  else
    F = zeros(1, N, get(mfd, 'nbdof'));
    gf_model_set(md, 'add initialized fem data', 'F', mfd, F);
    gf_model_set(md, 'add normal source term brick', mim, 'u', 'F', ...
        		 FORCE_BOUNDARY_NUM);
  end;
 
  gf_model_set(md, ...
 	     'add normal derivative Dirichlet condition with penalization', ...
 	       mim, 'u', 1e10, CLAMPED_BOUNDARY_NUM);
 
  gf_model_set(md, 'add Dirichlet condition with penalization', ...
    	       mim, 'u', 1e10, SIMPLE_SUPPORT_BOUNDARY_NUM);

  t0=cputime; 
  gf_model_get(md, 'solve', 'noisy');
  U = gf_model_get(md, 'variable', 'u');
  disp(sprintf('solve done in %.2f sec', cputime-t0));

else % uses old bricks

  if useKL
    b0 = gfMdBrick('bilaplacian', mim, mfu, 'Kirchhoff-Love');
    set(b0, 'param','D', D);
    set(b0, 'param','nu', NU);
    M = zeros(1, get(mfd,'nbdof'));
  else
    b0 = gfMdBrick('bilaplacian', mim, mfu);
    set(b0, 'param','D', D);
    M = zeros(1, get(mfd, 'nbdof'));
  end;

  b1 = gfMdBrick('source term', b0);
  set(b1, 'param', 'source_term', mfd, get(mfd, 'eval', {'1-(x-y).^2'}));

  b2 = gfMdBrick('normal derivative source term',b1,MOMENTUM_BOUNDARY_NUM);
  set(b2, 'param', 'source_term', mfd, M);

  if (useKL) 
    H = zeros(N*N, get(mfd, 'nbdof'));
    F = zeros(N, get(mfd, 'nbdof'));
    b3 = gfMdBrick('neumann Kirchhoff-Love source term',b2,FORCE_BOUNDARY_NUM);
    set(b3, 'param', 'M', mfd, H);
    set(b3, 'param', 'divM', mfd, F);
  else
    F = zeros(1, N, get(mfd, 'nbdof'));
    b3 = gfMdBrick('normal source term', b2, FORCE_BOUNDARY_NUM);
    set(b3, 'param', 'normal_source_term', mfd, F);
  end;

  b4 = gfMdBrick('dirichlet on normal derivative', b3, mfd, CLAMPED_BOUNDARY_NUM, 'penalized');


  b5 = gfMdBrick('dirichlet', b4, SIMPLE_SUPPORT_BOUNDARY_NUM, mfd, 'penalized');

  mds=gfMdState(b5);
  disp('running solve... ');
  t0=cputime; 

  get(b5, 'solve', mds, 'noisy');
  disp(sprintf('solve done in %.2f sec', cputime-t0));

  U=get(mds, 'state');
end


gf_plot(mfu,U,'mesh','on');
colorbar;

disp(sprintf('H2 norm of the solution: %g', gf_compute(mfu,U,'H2 norm', mim)));

%err=gf_compute(mfu,U,'interpolate on',mfd) - sol_u;

%disp(sprintf('H1 norm of the error: %g', gf_compute(mfd,err,'H1 norm', mim)));
%disp(sprintf('H2 norm of the error: %g', gf_compute(mfd,err,'H2 norm', mim)));
