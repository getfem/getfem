function [varargout]=subsref(obj, index)
% gfMesh/subsref
  nout = max(nargout, 1); cnt=1;
  FGET = @gf_mesh_get;
  FSET = @gf_mesh_set;
  switch index(1).type
    case '{}'
      error('Cell array indexing not supported by gfMesh objects')
    case '()'
      error('array indexing not supported by gfMesh objects')
    case '.'
      switch index(1).subs
        case 'id'
          [varargout{1:nout}] = obj.id;
        case 'set_pts'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'pts', index(2).subs{:});
          else
            FSET(obj, 'pts', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'add_point'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'add_point', index(2).subs{:});
          else
            FSET(obj, 'add_point', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'del_point'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'del_point', index(2).subs{:});
          else
            FSET(obj, 'del_point', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'add_convex'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'add_convex', index(2).subs{:});
          else
            FSET(obj, 'add_convex', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'del_convex'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'del_convex', index(2).subs{:});
          else
            FSET(obj, 'del_convex', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'del_convex_of_dim'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'del_convex_of_dim', index(2).subs{:});
          else
            FSET(obj, 'del_convex_of_dim', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'translate'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'translate', index(2).subs{:});
          else
            FSET(obj, 'translate', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'transform'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'transform', index(2).subs{:});
          else
            FSET(obj, 'transform', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'set_boundary'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'boundary', index(2).subs{:});
          else
            FSET(obj, 'boundary', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'set_region'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'region', index(2).subs{:});
          else
            FSET(obj, 'region', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'extend_region'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'extend_region', index(2).subs{:});
          else
            FSET(obj, 'extend_region', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'region_intersect'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'region_intersect', index(2).subs{:});
          else
            FSET(obj, 'region_intersect', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'region_merge'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'region_merge', index(2).subs{:});
          else
            FSET(obj, 'region_merge', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'region_subtract'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'region_subtract', index(2).subs{:});
          else
            FSET(obj, 'region_subtract', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'delete_boundary'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'delete_boundary', index(2).subs{:});
          else
            FSET(obj, 'delete_boundary', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'delete_region'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'delete_region', index(2).subs{:});
          else
            FSET(obj, 'delete_region', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'merge'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'merge', index(2).subs{:});
          else
            FSET(obj, 'merge', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'optimize_structure'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'optimize_structure', index(2).subs{:});
          else
            FSET(obj, 'optimize_structure', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'refine'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'refine', index(2).subs{:});
          else
            FSET(obj, 'refine', index(2).subs{:});
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
