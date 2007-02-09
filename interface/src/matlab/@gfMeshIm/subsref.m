function [varargout] = subsref(obj,index)
% gfMeshIm/SUBSREF Define field name indexing for gfMesh objects
  nout = max(nargout,1); cnt=1;
  FGET = @gf_mesh_im_get;
  FSET = @gf_mesh_im_set;
  switch index(1).type    
   case '.'
    switch index(1).subs
     case {'integ'}
      if (numel(index)>1 & index(2).type == '()' & numel(index(2).subs)==1)
        [varargout{1:nout}] = ...
            gfObject(FGET(obj,index(1).subs,index(2).subs{1}));
        cnt=cnt+1;
      else
        [varargout{1:nout}] = gfObject(FGET(obj,index(1).subs));
      end;
     case {'mesh','linked_mesh'}
      % extraction of the linked mesh
      [varargout{1:nout}] = gfMesh(FGET(obj, 'linked mesh'));
     case 'id'
      [varargout{1:nout}] = obj.id;
     case 'set'
      if (nargout) 
	[varargout{1:nargout}] = FSET(obj,index(2).subs{:});
      else
	FSET(obj,index(2).subs{:});
	if (exist('ans','var') == 1), varargout{1}=ans; end;
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
      [varargout{1:nout}] = FGET(obj,index(1).subs);
    end
  end
  % if there are others indexes, let matlab do its work
  if (numel(index) > cnt)
    for i=1:nout,
      varargout{i} = subsref(varargout{i},index((cnt+1):end));
    end;
  end;  
