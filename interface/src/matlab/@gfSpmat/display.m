function display(m)
% gfSpmat/display displays a short summary for a gfSpmat object
  for i=1:numel(m)
    disp(sprintf('gfSpmat object ID=%u %s', double(m(i).id), gf_spmat_get(m(i), 'info')));
  end;
