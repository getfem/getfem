function b=gfObject(a)
% converts a getfem struct to an object such as gfMesh, gfFem etc..    
  if (isempty(a))
    b=[];
  elseif (isobject(a))
    b=a;
  elseif (isstruct(a) & numel(a))
    sclass = gf_workspace('class name',a(1));
    b=feval(sclass,a);
  else error('neither an object nor a struct');
  end;