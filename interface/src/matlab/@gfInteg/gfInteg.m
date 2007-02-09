function m=gfInteg(a,varargin)
% gfInteg class constructor for integration methods. (see the function
% gf_integ for the list of arguments)
  
  this_class = 'gfInteg';
  if (nargin==0) error('can''t create an empty integration method reference'); end;
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
      m = gf_integ(a,varargin{:});
    end;
    m = class(m, this_class);
  end;
