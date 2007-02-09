function [varargout] = subsref(obj,index)
% gfMeshLevelSet/SUBSREF
  nout = max(nargout,1); cnt=1;
  FGET = @gf_mesh_levelset_get;
  FSET = @gf_mesh_levelset_set;
  switch index(1).type    
   case '{}'
    error('Cell array indexing not supported by gfMeshLevelSet objects')
   case '()'
    error('array indexing not supported by gfMeshLevelSet objects')
   case '.'
    switch index(1).subs
     case 'id'
      [varargout{1:nout}] = obj.id;
     case 'set'
      if (nargout) 
	[varargout{1:nargout}] = FSET(obj,index(2).subs{:});
      else
	FSET(obj,index(2).subs{:});
      end;
      return;
     case 'get'
      if (nargout) 
        [varargout{1:nargout}] = FGET(obj,index(2).subs{:});
      else
	FGET(obj,index(2).subs{:});
	if (exist('ans','var') == 1), varargout{1}=ans; end;
      end;
      return;
     otherwise
      [varargout{1:nargout}] = FGET(obj,index(1).subs);
    end
  end
  % if there are others indexes, let matlab do its work
  if (numel(index) > cnt)
    for i=1:nout,
      varargout{i} = subsref(varargout{i},index((cnt+1):end));
    end;
  end;  
