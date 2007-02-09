function display(o)
% displays a gfLevelSet object
  s = sprintf(['gfLevelSet object: ID=%u [%d bytes]'],...
	       double(o.id),...
	       gf_levelset_get(o,'memsize'));
  disp(s);

