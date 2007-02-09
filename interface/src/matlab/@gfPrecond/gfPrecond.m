function [m,b]=gfPrecond(a,varargin)
% gfPrecond class constructor (see the function gf_precond for the
% list of arguments)

  this_class = 'gfPrecond';
  if (nargin==0) error('can''t create an empty precond reference'); end;
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
      m = gf_precond(a,varargin{:});
    end;
    m.txt = '';
    m = class(m, this_class);
  end;
