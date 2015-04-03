function opts = build_options_list(varargin)
opts = init_param();
for i=1:length(varargin)/2
  opts = add_param(opts,varargin(2*(i-1)+1),varargin(2*(i-1)+2));
end
endfunction

