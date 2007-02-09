function display(em)
  for i=1:numel(em),
    disp(sprintf(['gfEltm (elementary matrix type) object ID=%u '], ...
		 double(em(i).id)));
  end;
  
