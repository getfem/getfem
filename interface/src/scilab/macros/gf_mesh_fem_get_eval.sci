// Copyright (C) 2009-2020 Yann Colette
// 
//  This file is a part of GetFEM++
// 
//  GetFEM++  is  free software;  you  can  redistribute  it  and/or modify it
//  under  the  terms  of the  GNU  Lesser General Public License as published
//  by  the  Free Software Foundation;  either version 3 of the License,  or
//  (at your option) any later version along with the GCC Runtime Library
//  Exception either version 3.1 or (at your option) any later version.
//  This program  is  distributed  in  the  hope  that it will be useful,  but
//  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
//  or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
//  License and GCC Runtime Library Exception for more details.
//  You  should  have received a copy of the GNU Lesser General Public License
//  along  with  this program;  if not, write to the Free Software Foundation,
//  Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.

function X=gf_mesh_fem_get_eval(mf, _what, dof)
// gf_mesh_fem_get_eval : see the help in gf_mesh_fem_get(mf,'eval')

[nargout,nargin] = argn();

if (nargin < 2) then error('not enough input arguments'); end;

qdim  = gf_mesh_fem_get(mf, 'qdim');
nbdof = gf_mesh_fem_get(mf,'nbdof');

if (nargin==2) then dof=1:qdim:nbdof; end;

// --- TODO --- only test the dof, not whole mesh
if (~gf_mesh_fem_get(mf, 'is lagrangian')) then
  error('interpolating on a non-lagrangian mesh fem');
end

//  if (qdim ~= 1) then
//    dof = dof(1:qdim:nbdof);
//  end

if (find(modulo(dof-1,qdim))) then
  error(['when qdim is different of 1, only dofs 1,qdim+1,',...
         '2*qdim+1,... are authorized']);
end

dxy = gf_mesh_fem_get(mf, 'basic dof nodes',dof);

if typeof(_what)~='list' then
  if (size(_what,2) == nbdof & typeof(_what)=='constant') then
    X = _what;
    return;
  elseif (typeof(_what)=='string') then
    error(['string expressions must be enclosed in a list: try with list( ',...
           'your_expression )']);
  elseif (size(_what,2) ~= qdim) then
    error(sprintf(['wrong dimensions for the expression: should have ',...
                   '%d (=Qdim) columns instead of %d'],qdim,size(_what,2)));
  end

  X = zeros(size(_what,1),nbdof);
end


if (typeof(_what)=='constant') then
  X(dof) = repmat(_what, 1, nbdof/qdim);
  return;
end

if typeof(_what)=='list' then
  X = zeros(size(_what),nbdof);
  m = size(_what); 

  xpos = dxy(1,:);
  if (size(dxy,1)>=2) then
    ypos = dxy(2,:);
  else 
    ypos = zeros(size(xpos)); 
  end
  if (size(dxy,1)>=3) then
    zpos = dxy(3,:);
  else 
    zpos = zeros(size(xpos)); 
  end
      
  for i=1:m
    for j=1:qdim
      if (typeof(_what(i)(j))=='constant') then
        if (length(_what(i)(j)) ~= 1) then error('numeric values should be scalar'); end;
        X(i,dof+j-1) = _what(i)(j);
      elseif (typeof(_what(i)(j))=='string') then
        x = xpos;
        y = ypos;
        z = zpos;
        //X(i,dof+j-1) = eval(_what(i)(j));
        X(i,dof+j-1) = evstr(_what(i)(j));
      elseif ((typeof(_what(i)(j))=='function')|(typeof(_what(i)(j))=='fptr')) then
        X(i,dof+j-1) = feval(xpos,ypos,zpos,_what(i)(j));
      else
        error(['sorry, don''t know how to eval a ' typeof(_what(i)(j)),...
               ' expression, only function handles, numeric constants and ',...
               'string expressions are handled']);
      end
    end
  end
else
  error(['can''t evaluate on mesh fem: argument is neither a numeric ',...
         'constant nor a cell array of (strings|constants|function handles)']);
end
endfunction
