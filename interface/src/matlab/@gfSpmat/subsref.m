function [varargout] = subsref(obj,index)
% gfSpmat/SUBSREF Define field name indexing for gfSpmat objects
%  accessible methods:
%    m.nc
%    m.nr
%    m.id
%    m.cvs ..  
%  m(i) returns the coordinate of point i
%  m{i} return the list of PIDs of convex i  
  nout = max(nargout,1); cnt=1;
  FGET = @gf_spmat_get;
  FSET = @gf_spmat_set;
  switch index(1).type
   case '{}'
    error('Cell array indexing not supported by gfSpmat objects')
   case {'()'}
    ri=index(1).subs{1}; rj = index(1).subs{2}; sz = FGET(obj, 'size');
    if (ri == ':') ri = 1:double(sz(1)); end;
    if (rj == ':') rj = 1:double(sz(2)); end;
    varargout{1}=FGET(obj, 'full', ri, rj);
   case '.'
    switch index(1).subs
     case 'id'
      [varargout{1:nout}] = obj.id;
     case 'nr'
      sz = FGET(obj, 'size');
      [varargout{1:nout}] = sz(1);
     case 'nc'
      sz = FGET(obj, 'size');
      [varargout{1:nout}] = sz(2);
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
