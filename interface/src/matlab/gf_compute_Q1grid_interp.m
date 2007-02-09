function [U2,Iq,MF2]=gf_compute_Q1grid_interp(MF1,U1,varargin)
% See help on gf_compute
%  $Id$
  if (nargin < 3), error('not enough input arguments'); end;
    
  gf_workspace('push', 'gf_compute_Q1grid_interp');
  
  
  
  meshpts = gf_mesh_get(MF1, 'pts');
  zmin = min(meshpts,[],2); zmax = max(meshpts,[],2);
  ndim = length(zmin);
  if (ndim > 3) error('this number of dim is not supported (patch me)'); end;
  X = cell(ndim,1);
  
  try
    switch varargin{1}
     case 'regular h'
      if (length(varargin) ~= 2), error('wrong number of arguments'); end;
      if (length(varargin{2}) ~= ndim), error('invalid dimension'); end;
      for i=1:ndim,
	if (varargin{2}(i) <= 0), error('invalid step value'); end;
	X{i} = zmin(i):(varargin{2}(i)):zmax(i);
      end;
     case 'regular N'
      if (length(varargin) ~= 2), error('wrong number of arguments'); end;
      if (length(varargin{2}) ~= ndim), error('invalid dimension'); end;
      for i=1:ndim,
	if (varargin{2}(i) <= 0), error('invalid number of cells'); end;
	h = (zmax(i) - zmin(i))/(varargin{2}(i));
	X{i} = zmin(i):h:zmax(i);
      end;
     otherwise
      X = varargin{1};
      if (~iscell(X)), error('grid points should be stored in a cell array of size nbdim'); end;
      if (length(X) ~= ndim) error('wrong number of dimension in the grid points argument'); end;
    end

    Q=gf_mesh_fem_get(MF1,'qdim');
    M = gf_mesh('cartesian', X{:});
    MF2 = gf_mesh_fem(M,Q);
    gf_mesh_fem_set(MF2, 'classical fem', 1); % Q1 fem
    mfU2 = gf_compute(MF1,U1, 'interpolate on', MF2);

    PTS = gf_mesh_fem_get(MF2, 'dof nodes');

    PTS = PTS(end:-1:1,1:Q:end);   % (x,y,z)->(z,y,x) and remove duplicate dof
    [PTS,I] = sortrows(PTS'); % sort points, by z then by y then by x etc..
    I = Q*(I-1) + 1;
    sz = Q;
    for i=1:numel(X) sz = [sz length(X{i})]; end;
    Iq=zeros(Q,numel(I));
    for q=1:Q,
      Iq(q,:) = I'+(q-1);
    end;
    Iq = Iq(:);
    U2 = reshape(mfU2(Iq),sz);
    if (nargout == 3),
      gf_workspace('keep', MF2);
    end;
  catch,
    gf_workspace('pop');
    error(lasterr);
  end
  gf_workspace('pop');
  