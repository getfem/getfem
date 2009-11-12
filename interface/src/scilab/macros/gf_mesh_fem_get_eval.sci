function X=gf_mesh_fem_get_eval(mf, _what, dof)
// gf_mesh_fem_get_eval : see the help in gf_mesh_fem_get(mf,'eval')

[nargout,nargin] = argn();

if (nargin < 2) then error('not enough input arguments'); end;
  
qdim  = gf_mesh_fem_get(mf,'qdim');
nbdof = gf_mesh_fem_get(mf,'nbdof'); 

if (nargin==2) then 
  dof=1:qdim:nbdof; 
end

// --- TODO --- only test the dof, not whole mesh
if (~gf_mesh_fem_get(mf, 'is lagrangian')) then
  error('interpolating on a non-lagrangian mesh fem');
end
//  if (qdim ~= 1),
//    dof = dof(1:qdim:nbdof);
//  end;
if (find(modulo(dof-1,qdim))) then
  error('when qdim is different of 1, only dofs 1,qdim+1, 2*qdim+1,... are authorized');
end
dxy = gf_mesh_fem_get(mf, 'basic dof nodes',dof);

if isnumeric(_what) then
  if size(_what, 2) == nbdof then
    X = _what;
    return;
  end
  if (size(_what,2) ~= qdim) then
    error(sprintf('wrong dimensions for the expression: should have %d (=Qdim) columns instead of %d',qdim,size(_what,2)));
  end
elseif (typeof(_what)=='string') then
  error('string expressions must be enclosed in a list: try with list(your_expression)');
end;

X = [];
if (isnumeric(_what)) then
  X = zeros(length(_what),nbdof);
  X(:,dof) = repmat(_what, 1, nbdof/qdim);
  return;
elseif typeof(_what)=='hypermat' then
  m = prod(size((_what)));
  X = zeros(m,nbdof);
  xpos = dxy(1,:);
  if (size(dxy,1)>=2),
    ypos = dxy(2,:);
  else 
    ypos = zeros(xpos); 
  end
  if (size(dxy,1)>=3),
    zpos = dxy(3,:);
  else 
    zpos = zeros(xpos); 
  end
    
  for i=1:m
    for j=1:qdim
      if (isnumeric(_what(i,j))) then
        if (length(_what(i,j)) ~= 1) then error('numeric values should be scalar'); end;
        X(i,dof+j-1) = _what(i,j);
      elseif (typeof(_what(i,j))=='string') then
        x = xpos;
        y = ypos;
        z = zpos;
        X(i,dof+j-1)=eval(_what(i,j));
      elseif (type(_what(i,j))==11 | type(_what(i,j))==13) then
        X(i,dof+j-1) = feval(_what(i,j), xpos, ypos, zpos);
      else
        error('sorry, don''t know how to eval a ' + typeof(_what(i,j)) + ...
              ' expression, only function handles, numeric constants and ' + ...
              'string expressions are handled');
      end
    end
  end
elseif typeof(_what)=='list' then
  xpos = dxy(1,:);
  if (size(dxy,1)>=2),
    ypos = dxy(2,:);
  else 
    ypos = zeros(xpos); 
  end
  if (size(dxy,1)>=3),
    zpos = dxy(3,:);
  else 
    zpos = zeros(xpos); 
  end
  
  m = length(_what);
  n = length(_what(1));
  X = zeros(m,nbdof);
  for i=1:m
    for j=1:qdim
      for k=1:n
        if (isnumeric(_what(i)(k))) then
          if (length(_what(i)(k)) ~= 1) then error('numeric values should be scalar'); end;
          X(i,dof+j-1) = _what(i)(k);
        elseif (typeof(_what(i)(k))=='string') then
          x = xpos;
          y = ypos;
          z = zpos;
          X(i,dof+j-1) = evstr(_what(i)(k));
        elseif (typeof(_what(i)(k))=='fptr' | typeof(_what(i)(k))=='function') then
          tmp = evstr(_what(i)(k));
          X(i,dof+j-1) = feval(tmp, xpos, ypos, zpos);
        end
      end
    end
  end
else
  error('can''t evaluate on mesh fem: argument is neither a numeric ' + ...
        'constant nor a cell array of (strings|constants|function handles)');
end
endfunction

