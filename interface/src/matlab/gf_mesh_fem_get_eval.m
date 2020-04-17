function X=gf_mesh_fem_get_eval(mf, what, dof)
% gf_mesh_fem_get_eval : see the help in gf_mesh_fem_get(mf,'eval')
%  Copyright (C) 1999-2020 Yves Renard
%
%  This file is a part of GetFEM
%
%  GetFEM  is  free software;  you  can  redistribute  it  and/or modify it
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
  
  if (nargin < 2) error('not enough input arguments'); end;
  qdim=gf_mesh_fem_get(mf, 'qdim');
  nbdof=gf_mesh_fem_get(mf,'nbdof');
  if (nargin==2) dof=1:qdim:nbdof; end;
  % --- TODO --- only test the dof, not whole mesh
  if (~gf_mesh_fem_get(mf, 'is lagrangian')),
    error('interpolating on a non-lagrangian mesh fem');
  end;
%  if (qdim ~= 1),
%    dof = dof(1:qdim:nbdof);
%  end;
  if (find(mod(dof-1,qdim)))
    error(['when qdim is different of 1, only dofs 1,qdim+1,',...
	   '2*qdim+1,... are authorized']);
  end;
  dxy = gf_mesh_fem_get(mf, 'basic dof nodes',dof);
  
  if (size(what, 2) == nbdof & isnumeric(what)),
    X = what;
    return;
  elseif (ischar(what))
    error(['string expressions must be enclosed in a cell array: try with { ',...
           'your_expression }']);
  elseif (size(what,2) ~= qdim)
    error(sprintf(['wrong dimensions for the expression: should have ',...
		   '%d (=Qdim) columns instead of %d'],qdim,size(what,2)));
  end;
  
  X=zeros(size(what,1),nbdof);
  if (isnumeric(what)),
    X(dof) = repmat(what, 1, nbdof/qdim);
    return;
  elseif iscell(what),
    m=size(what,1);
    xpos = dxy(1,:);
    if (size(dxy,1)>=2),
      ypos = dxy(2,:);
    else ypos = zeros(size(xpos)); end;
    if (size(dxy,1)>=3),
      zpos = dxy(3,:);
    else zpos = zeros(size(xpos)); end;
      
    for i=1:m,
      for j=1:qdim
	if (isnumeric(what{i,j})),
	  if (numel(what{i,j}) ~= 1) error('numeric values should be scalar'); end;
	  X(i,dof+j-1)=what{i,j};
	elseif (ischar(what{i,j})),
	  x=xpos; y=ypos; z=zpos;
	  X(i,dof+j-1)=eval(what{i,j});
	elseif (isa(what{i,j},'function_handle'))
	  X(i,dof+j-1)=feval(what{i,j}, xpos, ypos, zpos);
	else
	  error(['sorry, don''t know how to eval a ' class(what{i,j}),...
		' expression, only function handles, numeric constants and ',...
		'string expressions are handled']);
	end;
      end;
    end;
  else
    error(['can''t evaluate on mesh fem: argument is neither a numeric ',...
	   'constant nor a cell array of (strings|constants|function handles)']);
  end;
  
