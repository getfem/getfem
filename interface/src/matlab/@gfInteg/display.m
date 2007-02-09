function display(im)
% gfInteg/display displays a short summary for a gfInteg object
% see gfInteg/char for an exhaustive description of the object.
  
  if (gf_integ_get(im,'is_exact')),
    s = 'exact integration method';
  else
    n = gf_integ_get(im,'nbpts');
    s2 = '';
    if (numel(n)>1) 
      s2=['+[' sprintf('%d',n(2)) sprintf(',%d',n(3:end)) ']']; 
    end;    
    s = ['approximate integration method with ' ...
	 sprintf('%d',n(1)) s2 ' points'];
  end;
  disp(sprintf(['gfInteg object ID=%u dim=%d, %s\n -> %s'],...
	      double(im.id),...
	      gf_integ_get(im,'dim'),s,...
	      gf_integ_get(im,'char')));
  