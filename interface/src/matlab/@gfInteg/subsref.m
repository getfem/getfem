function [varargout] = subsref(obj,index)
% gfInteg/SUBSREF Define field name indexing for gfInteg objects
%  accessible methods:
%    im.id  
%    im.is_exact
%    im.nbpts
%    im.pts
%    im.coeffs
%    ...etc
  FGET = @gf_integ_get;
  nout = max(nargout,1);  
  switch index(1).type
   case '()'
    error('indexes not supported');
   case '{}'
    error('cell access not supported');
   case '.'
    switch index(1).subs
     case 'id'
      varargout{1} = obj.id;
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
  if (numel(index) > 1)
    for i=1:nargout,
      varargout{i} = subsref(varargout{i},index(2:end));
    end;
  end;
