function res = isnumeric(parameter)
  res = or(type(parameter) == [1 5 8]);
  res = res | (typeof(parameter)=='hypermat');
endfunction

