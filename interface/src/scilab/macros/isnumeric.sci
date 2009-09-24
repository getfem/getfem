function res = isnumeric(param)
  if ~isdef('param','local') then param = ''; end
  res = or(type(param) == [1 5 8]);
  res = res | (typeof(param)=='hypermat');
endfunction

