function [m,b]=gfModel(a,varargin)
% gfModel class constructor (see the function gf_model for argument
% list)
  
  this_class = 'gfModel';
  if (nargin==0) error('can''t create an empty model reference'); end;
  if (isa(a,this_class) & nargin==1),
    m=a;
  else
    if (isstruct(a) & isfield(a,'id') & isfield(a,'cid'))
      cname = gf_workspace('class name',a);
    else
      cname = class(a);
    end;
    if (strcmp(cname, this_class) & nargin == 1)
      m = a;
    else
      m = gf_model(a,varargin{:});
    end;
    m.txt = '';
    m = class(m, this_class);
  end;
  
