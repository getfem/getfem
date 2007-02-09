function [varargout] = subsref(obj,index)
% gfMesh/SUBSREF Define field name indexing for gfMesh objects
%  accessible methods:
%    m.nbpts
%    m.nbcvs
%    m.dim
%    m.id
%    m.cvs ..  
%  m(i) returns the coordinate of point i
%  m{i} return the list of PIDs of convex i  
  FGET=@gf_mesh_get;
  FSET=@gf_mesh_set;
  nout = max(nargout,1);  
  cnt=1;
  switch index(1).type
   case '{}'
    error('Cell array indexing not supported by gfMesh objects')
   case {'()'}
    varargout{1}=obj(index(1).subs{1});
   case '.'
    switch index(1).subs
     % return cvstruct and geotrans arrays
     % by default, calls gf_mesh_get(m, 'cvstruct')
     % but if something like m.cvstruct(1:2) was asked,
     % the call is then gf_mesh_get(m, 'cvstruct',1:2)
     case {'cvstruct','geotrans'}
      if (numel(index)>1 & index(2).type == '()' & numel(index(2).subs)==1)
        [varargout{1:nout}] = gfObject(FGET(obj,index(1).subs,index(2).subs{1}));
        cnt=cnt+1;
      else
        [varargout{1:nout}] = gfObject(FGET(obj,index(1).subs));
      end;
     case {'pts','pid','cvid','pid_from_cvid','pid_from_coords',...
	   'cvid_from_pid','faces_from_pid','outer_faces'}
      if (numel(index)>1 & index(2).type == '()' & numel(index(2).subs)==1)
        [varargout{1:nout}] = FGET(obj,index(1).subs,index(2).subs{1});
        cnt=cnt+1;
      else
        [varargout{1:nout}] = FGET(obj,index(1).subs);
      end;
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
	varargout{1}=FGET(obj,index(2).subs{:});
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
