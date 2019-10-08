function [varargout]=subsref(obj, index)
% gfMeshFem/subsref
  nout = max(nargout, 1); cnt=1;
  FGET = @gf_mesh_fem_get;
  FSET = @gf_mesh_fem_set;
  switch index(1).type
    case '{}'
      error('Cell array indexing not supported by gfMeshFem objects')
    case '()'
      error('array indexing not supported by gfMeshFem objects')
    case '.'
      switch index(1).subs
        case 'id'
          [varargout{1:nout}] = obj.id;
        case 'set_fem'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'fem', index(2).subs{:});
          else
            FSET(obj, 'fem', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'set_classical_fem'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'classical_fem', index(2).subs{:});
          else
            FSET(obj, 'classical_fem', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'set_classical_discontinuous_fem'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'classical_discontinuous_fem', index(2).subs{:});
          else
            FSET(obj, 'classical_discontinuous_fem', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'set_qdim'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'qdim', index(2).subs{:});
          else
            FSET(obj, 'qdim', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'reduction_matrices'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'reduction_matrices', index(2).subs{:});
          else
            FSET(obj, 'reduction_matrices', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'reduction'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'reduction', index(2).subs{:});
          else
            FSET(obj, 'reduction', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'reduce_meshfem'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'reduce_meshfem', index(2).subs{:});
          else
            FSET(obj, 'reduce_meshfem', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'set_dof_partition'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'dof_partition', index(2).subs{:});
          else
            FSET(obj, 'dof_partition', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'set_partial'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'set_partial', index(2).subs{:});
          else
            FSET(obj, 'set_partial', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'adapt'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'adapt', index(2).subs{:});
          else
            FSET(obj, 'adapt', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'set_enriched_dofs'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'set_enriched_dofs', index(2).subs{:});
          else
            FSET(obj, 'set_enriched_dofs', index(2).subs{:});
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
