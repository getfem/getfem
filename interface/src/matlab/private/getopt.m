function o=getopt(opt,v)
  o = opt;
  if (nargin==1) return; end;
  if (mod(length(v),2) ~= 0) error('odd number of property/value pairs'); end;
  for i=1:2:length(v),
    optname=v{i};
    optval =v{i+1};
    if (~ischar(optname)) error(['expecting a property name, found a ' class(optname)]); end;
    if (~isfield(opt,optname)) error(['unknown property :"' optname '"']); end;
    o = setfield(o, optname, optval);
  end;
