function [m,b]=gfMeshLevelSet(a,varargin)
% gfLevelSet class constructor (see the function gf_mesh_levelset for argument
% list)
  
  this_class = 'gfMeshLevelSet';
  if (nargin==0) error('can''t create an empty MeshLevelSet reference'); end;
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
      m = gf_mesh_levelset(a,varargin{:});
    end;
    m.txt = '';
    m = class(m, this_class);
  end;
  
