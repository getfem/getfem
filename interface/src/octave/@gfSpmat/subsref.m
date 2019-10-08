function [varargout]=subsref(obj, index)
% gfSpmat/subsref
  nout = max(nargout, 1); cnt=1;
  FGET = @gf_spmat_get;
  FSET = @gf_spmat_set;
  switch index(1).type
    case '{}'
      error('Cell array indexing not supported by gfSpmat objects')
    case '()'
      error('array indexing not supported by gfSpmat objects')
    case '.'
      switch index(1).subs
        case 'id'
          [varargout{1:nout}] = obj.id;
        case 'clear'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'clear', index(2).subs{:});
          else
            FSET(obj, 'clear', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'scale'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'scale', index(2).subs{:});
          else
            FSET(obj, 'scale', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'transpose'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'transpose', index(2).subs{:});
          else
            FSET(obj, 'transpose', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'conjugate'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'conjugate', index(2).subs{:});
          else
            FSET(obj, 'conjugate', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'transconj'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'transconj', index(2).subs{:});
          else
            FSET(obj, 'transconj', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'to_csc'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'to_csc', index(2).subs{:});
          else
            FSET(obj, 'to_csc', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'to_wsc'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'to_wsc', index(2).subs{:});
          else
            FSET(obj, 'to_wsc', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'to_complex'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'to_complex', index(2).subs{:});
          else
            FSET(obj, 'to_complex', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'set_diag'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'diag', index(2).subs{:});
          else
            FSET(obj, 'diag', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'assign'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'assign', index(2).subs{:});
          else
            FSET(obj, 'assign', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'add'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'add', index(2).subs{:});
          else
            FSET(obj, 'add', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'set'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, index(2).subs{:});
          else
            FSET(obj,index(2).subs{:});
            if (exist('ans', 'var') == 1)
              h=ans;
              if (isstruct(h) & isfield(h,'id') & isfield(h,'cid')), h = gfObject(h); end;
              varargout{1}=h;
            end;
          end;
          return;
        case 'get'
          if (nargout) 
            h = FGET(obj, index(2).subs{:});
            if (isstruct(h) & isfield(h,'id') & isfield(h,'cid')), h = gfObject(h); end;
            [varargout{1:nargout}] = h;
          else
	     FGET(obj,index(2).subs{:});
            if (exist('ans', 'var') == 1)
              h=ans;
              if (isstruct(h) & isfield(h,'id') & isfield(h,'cid')), h = gfObject(h); end;
              varargout{1}=h;
            end;
          end;
          return;
        otherwise
          if ((numel(index) > 1) && (strcmp(index(2).type, '()')))
            h = FGET(obj,index(1).subs, index(2).subs{:});
            if (isstruct(h) & isfield(h,'id') & isfield(h,'cid')), h = gfObject(h); end;
            [varargout{1:nargout}] = h;
            cnt = cnt + 1;
          else
            h = FGET(obj, index(1).subs);
            if (isstruct(h) & isfield(h,'id') & isfield(h,'cid')), h = gfObject(h); end;
            [varargout{1:nargout}] = h;
          end
      end
  end
  % if there are others indexes, let matlab do its work
  if (numel(index) > cnt)
    for i=1:nout,
      varargout{i} = subsref(varargout{i}, index((cnt+1):end));
    end;
  end;
% autogenerated mfile;
