function m=gfFem(a,varargin)
% gfFem class constructor for FEMs. (see the function gf_fem for the list of
% arguments)
  
  this_class = 'gfFem';
  if (nargin==0) error('can''t create an empty Finite Element Method reference'); end;
  if (isa(a,this_class)),
    if (nargin == 1)
      m=a;
    else error('too much input arguments'); end;
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
      m = gf_fem(a,varargin{:});
    end;
    m = class(m, this_class);
  end;
