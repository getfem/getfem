function [m,b]=gfCvStruct(a,varargin)
% gfCvStruc class pseudo-constructor 
  this_class = 'gfCvStruct';
  if (nargin==0) error('can''t create an empty convex_structure reference'); end;
  if (isa(a,this_class)),
    m=a;
  else    
    if (isstruct(a) & isfield(a,'id') & isfield(a,'cid'))
      cname = gf_workspace('class name',a);
    else
      cname = class(a);
    end;
    if (strcmp(cname, this_class))
      m = a;
    else
      error(['can''t create a ' this_class ' object from a ' cname ...
	     ' object']);
    end;
    m = class(m, this_class);
  end;
 