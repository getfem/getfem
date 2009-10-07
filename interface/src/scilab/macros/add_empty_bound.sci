function pde = add_empty_bound(pde)
if isempty(pde('bound')) then
  pde('bound') = list();
end

pde('bound')($+1) = mlist(['bound','type','R','H']);
pde('bound')($)('type') = 'Dirichlet';
pde('bound')($)('R') = [];
pde('bound')($)('H') = [];
endfunction

