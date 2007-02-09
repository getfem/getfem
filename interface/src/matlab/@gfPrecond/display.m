function display(m)
% gfPrecond/display displays a short summary for a gfPrecond object
  for i=1:numel(m)
    disp(sprintf('gfPrecond object ID=%u %s', double(m(i).id), gf_precond_get(m(i), 'info')));
  end;
