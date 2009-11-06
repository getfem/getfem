///////////////
// has_field //
///////////////

function ok = has_field(pde,varargin)
ok = 0;
for i=1:length(varargin),
  if (~or(getfield(1,pde)==varargin(i))) then
    return;
  end
end
ok = 1;
endfunction
