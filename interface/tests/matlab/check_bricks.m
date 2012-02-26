% Copyright (C) 2005-2012 Julien Pommier.
%
% This file is a part of GETFEM++
%
% Getfem++  is  free software;  you  can  redistribute  it  and/or modify it
% under  the  terms  of the  GNU  Lesser General Public License as published
% by  the  Free Software Foundation;  either version 2.1 of the License,  or
% (at your option) any later version.
% This program  is  distributed  in  the  hope  that it will be useful,  but
% WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
% or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
% License for more details.
% You  should  have received a copy of the GNU Lesser General Public License
% along  with  this program;  if not, write to the Free Software Foundation,
% Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.

function check_bricks(iverbose,idebug)
  global gverbose;
  global gdebug;  
  gf_workspace('clear all');
  if (nargin >= 1),
    gverbose = iverbose;
    if (nargin == 2),
      gdebug = idebug;
    else 
      gdebug = 0; end;
  else 
    gverbose = 0;
  end;
  m = gf_mesh('cartesian', 0:.1:1, 0:.1:1);
  mf = gf_mesh_fem(m); gf_mesh_fem_set(mf, 'classical fem', 2);
  mfq = gf_mesh_fem(m,2); gf_mesh_fem_set(mfq, 'classical fem', 2);
  mim = gf_mesh_im(m); gf_mesh_im_set(mim, 'integ', 6);
  
  P=gf_mesh_get(m,'pts'); % get list of mesh points coordinates
  pidtop=find(abs(P(2,:)-1)<1e-6); % find those on top of the object
  pidbot=find(abs(P(2,:)-0)<1e-6); % find those on the bottom
  ftop=gf_mesh_get(m,'faces from pid',pidtop); 
  fbot=gf_mesh_get(m,'faces from pid',pidbot);
  % assign boundary numbers
  gf_mesh_set(m,'boundary',1,ftop);
  gf_mesh_set(m,'boundary',2,fbot);
  gf_mesh_set(m,'boundary',3,gf_mesh_get(m,'outer faces'));
  
  
  helm = gf_mdbrick('helmholtz', mim, mf, 'complex');
  lapl = gf_mdbrick('generic elliptic', mim, mf);
  mass = gf_mdbrick('mass matrix', mim, mf);
  elas = gf_mdbrick('isotropic linearized elasticity', mim, mfq);
  elasnl = gf_mdbrick('nonlinear elasticity',mim, mfq, 'Ciarlet Geymonat');
  
  rms  = gf_mdstate('real');
  cms  = gf_mdstate(helm);
  
  foo=gf_mdstate_get(cms,'is_complex'); gfassert('foo==1');

  gf_mdbrick_set(helm, 'param', 'wave_number', 2);
  wv=gf_mdbrick_get(helm, 'param', 'wave_number');

  gf_mdbrick_set(mass, 'param', 'rho', mf, 0.1*ones(gf_mesh_fem_get(mf, ...
						  'nbdof'),1));
  l=gf_mdbrick_get(elas, 'param', 'lambda');
  gf_mdbrick_set(elas, 'param', 'lambda', l*3);
  
  for b=[lapl, elas, helm, mass, elasnl],
    disp_info(b);
  end;
  
  
function disp_info(b)
  disp(sprintf('subclass %s', gf_mdbrick_get(b, 'subclass')));
  disp(sprintf(' dim            %d', gf_mdbrick_get(b, 'dim')));
  disp(sprintf(' nbdof          %d', gf_mdbrick_get(b, 'nbdof')));
  disp(sprintf(' nb_constraints %d', gf_mdbrick_get(b, 'nb_constraints')));
  disp(sprintf(' is_linear      %d', gf_mdbrick_get(b, 'is_linear')));
  disp(sprintf(' is_symmetric   %d', gf_mdbrick_get(b, 'is_symmetric')));
  disp(sprintf(' is_coercive    %d', gf_mdbrick_get(b, 'is_coercive')));
  disp(sprintf(' is_complex     %d', gf_mdbrick_get(b, 'is_complex')));
  l = gf_mdbrick_get(b, 'param_list');
  s = ''; for il=1:numel(l), s = [s  ' '  l{il}]; end;
  disp(sprintf(' paramlist     %s', s)); 

  dir = gf_mdbrick('generalized dirichlet', b, 2);
  R=gf_mdbrick_get(dir, 'param', 'R'); R = 0*R+1; 
  gf_mdbrick_set(dir, 'param', 'R', R);

  disp('gogogogo')
  state = gf_mdstate(dir);
  gf_mdbrick_get(dir, 'solve', state, 'noisy');
  res = gf_mdstate_get(state, 'residual');
  U=gf_mdstate_get(state, 'state');
  disp(sprintf('solve done: |residual|=%g, |U|=%g',norm(res),norm(U)));
  T=gf_mdstate_get(state, 'tangent matrix');
  nnz(T)
  
  gf_workspace('stats');
  
  gf_delete(dir);

  
