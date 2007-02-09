function varargout=get(gt,varargin)
  % gfGeoTrans/get
  n = max(nargout,1);
  [varargout{1:n}]=gf_geotrans_get(gt,varargin{:});
