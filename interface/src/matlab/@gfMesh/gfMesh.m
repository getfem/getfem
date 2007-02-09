function [m,b]=gfMesh(a,varargin)
% gfMesh class constructor (see the function gf_mesh for the list of arguments)

  this_class = 'gfMesh';
  if (nargin==0) error('can''t create an empty mesh reference'); end;
  if (isa(a,this_class)),
    m=a;
  else
    if (isstruct(a) & isfield(a,'id') & isfield(a,'cid'))
      cname = gf_workspace('class name',a);
      if (strcmp(cname, this_class))
	m = a;
      else 
	error(['can''t create a ' this_class ' object from a ' cname ...
	       ' object']);
      end;
    else
      m = gf_mesh(a,varargin{:});
    end;
    m.txt = '';
    m = class(m, this_class);
  end;
