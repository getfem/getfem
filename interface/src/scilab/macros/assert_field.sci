//////////////////
// assert_field //
//////////////////

function assert_field(pde,varargin)
for i=1:length(varargin),
  if (~or(getfield(1,pde)==varargin(i))) then
    error('no member ' + varargin(i) + ' in mlist pde!'); 
  end
end
endfunction
