function [Q,G,H,R,F]=gf_asm_pdetoolbc(mf_u, mf_d, b, e, f_expr)
%  FUNCTION [Q,G,H,R,F]=gf_asm_pdetoolbc(mf_u, mf_d, b, e, f_expr)
%  'pdetool style' assembling of boundary conditions
%  See gf_asm
%  Copyright (C) 1999-2017 Yves Renard
%
%  This file is a part of GetFEM++
%
%  GetFEM++  is  free software;  you  can  redistribute  it  and/or modify it
%  under  the  terms  of the  GNU  Lesser General Public License as published
%  by  the  Free Software Foundation;  either version 3 of the License,  or
%  (at your option) any later version along with the GCC Runtime Library
%  Exception either version 3.1 or (at your option) any later version.
%  This program  is  distributed  in  the  hope  that it will be useful,  but
%  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
%  or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
%  License and GCC Runtime Library Exception for more details.
%  You  should  have received a copy of the GNU Lesser General Public License
%  along  with  this program;  if not, write to the Free Software Foundation,
%  Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.
% $Id$  
  if (nargin < 4), error('not enough input arguments'); end;
  nb_boundaries = size(b,2);
  xyeval = gf_mesh_fem_get(mf_d, 'dof nodes');
  nbdof = size(xyeval,2);
  N=b(1,1); % dimension of the system (1->scalar, 2->vector 2D)
  qdim = gf_mesh_fem_get(mf_u, 'qdim');

  
  if (qdim ~= N),
    error(sprintf(['the boundary condition b was generated for a %d-D problem, '...
		   'while the Qdim of the mesh_fem is %d'], N, qdim));
  end;
  if (nargin >= 5 && qdim ~= size(f_expr,1)),
    error('the qdim of the mesh fem and the size of f (the volumic source term) do not match');
  end;
  if (gf_mesh_fem_get(mf_d, 'qdim') ~= 1),
    error('the Qdim of the data mesh_fem should always be 1');
  end;
  
  for bnum=1:nb_boundaries,
    %ignores
    if (b(1,bnum)==0), continue; end;
    
    %select edges which belong to
    %boundary 'bnum'
    E=e(:, find(e(5,:)==bnum));
    
    EC = gf_mesh_get(mf_d, 'faces from pid', [E(1,:) E(2,:)]);
%    EC=[];
%    for i=1:size(E,2)
%      EC = [EC gf_mesh_get(mf_d, 'faces from pid', E(1:2,i))];
%    end;
    gf_mesh_fem_set(mf_d, 'boundary',bnum, EC);
    gf_mesh_fem_set(mf_u, 'boundary',bnum, EC);

    M = b(2,bnum); % number of dirichlet conditions (0,1,..N)
    if (M>N), disp('invalid geometry matrix'); return; end;
    
    clear qexpr gexpr;
    pos_len = 3;
    pos = 3+N*N+N+M*N+M;
    
    % reading Q expressions
    for j=1:N,
      for i=1:N,
	len=b(pos_len,bnum);
	qexpr{i,j}=char(b(pos:(pos+len-1), bnum)');
	pos_len = pos_len+1;
	pos = pos+len;
      end;
    end;
    
    
    % reading G expressions
    for i=1:N,
      len=b(pos_len,bnum);
      gexpr{i}=char(b(pos:(pos+len-1), bnum)');
      pos_len = pos_len+1;
      pos = pos+len;
    end;      
    
    
    % reading H expression
    for j=1:N,
      for i=1:M,
	len=b(pos_len,bnum);
	hexpr{i,j}=char(b(pos:(pos+len-1), bnum)');
	pos_len = pos_len+1;
	pos = pos+len;	
      end;
    end;
    
    
    % reading R expressions
    for i=1:M,
      len=b(pos_len,bnum);
      rexpr{i}=char(b(pos:(pos+len-1), bnum)');
      pos_len = pos_len+1;
      pos = pos+len;	
    end;
    
    
    % computations of expressions on the dof
    vQ = zeros(N,N,nbdof);
    for i=1:N,
      for j=1:N,
	vQ(i,j,:) = eval_expr(xyeval, qexpr{i,j});
      end;
    end;
    vG = zeros(N,nbdof);
    for i=1:N,
      vG(i,:) = eval_expr(xyeval, gexpr{i});
    end;

    vH = zeros(N,N,nbdof);
    for i=1:M,
      for j=1:N,
	vH(i,j,:) = eval_expr(xyeval, hexpr{i,j});
      end;
    end;
    
    vR = zeros(N,nbdof);
    for i=1:M,
      vR(i,:) = eval_expr(xyeval, rexpr{i});
    end;
    bQ = gf_asm('boundary qu term', bnum, mf_u, mf_d, reshape(vQ,N*N,nbdof));
%    bH = gf_asm('boundary qu term', bnum, mf_u, mf_d, reshape(vH,N*N,nbdof));
    bG = gf_asm('boundary source', bnum, mf_u, mf_d, reshape(vG,N,nbdof));
%    bR = gf_asm('neumann', bnum, mf_u, mf_d, reshape(vR,N,nbdof));
    [bH,bR] = gf_asm('dirichlet',bnum, mf_u, mf_d, reshape(vH,N*N,nbdof), reshape(vR,N,nbdof));
    if (bnum ~= 1),
      Q=Q+bQ; G=G+bG; H=H+bH; R=R+bR;
    else
      Q=bQ; G=bG; H=bH; R=bR;
    end
  end;

  % check for volumic source term
  if (nargin == 5 && nargout == 5),
    if (isstr(f_expr)),
      Fd = zeros(N,nbdof);
      for i=1:N
	Fd(i,:) = eval_expr(xyeval, f_expr(i,:));
      end;
    else
      Fd = f_expr;
    end;
    F=gf_asm('volumic source', mf_u, mf_d, Fd);
    F=F(:);
  end;
  
function V=eval_expr(xypos, expr)
  %disp(['expr=' expr]);
  V=zeros(1,size(xypos,2));
  x=xypos(1,:); 
  if (size(xypos,1) >= 2), y=xypos(2,:); else y=0; end;
  if (size(xypos,1) >= 3), z=xypos(3,:); else z=0; end;
  eval(['e=' expr ';']); V(:) = e(:);
  
