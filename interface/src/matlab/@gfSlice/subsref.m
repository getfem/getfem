function [varargout] = subsref(obj,index)
% gfSlice/SUBSREF Define field name indexing for gfSlice objects
% you can use the following (pseudo) object members:
%   .nbpts
%   .dim
%   .id
  nout = max(nargout,1); cnt=1;
  FGET = @gf_slice_get;
  FSET = @gf_slice_set;
  switch index(1).type    
   case '{}'
    error('Cell array indexing not supported by gfSlice objects')
   case '.'
    switch index(1).subs
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
