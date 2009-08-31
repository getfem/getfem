function assert(sx)
global gverbose;
global gdebug;
//x = evalin('caller',sx);
x = eval(sx);
if (~and(x(:))) then
  if (gverbose) then
    dbstack;      
  end
  if (gdebug) then
    disp('Assertion failed:' + sx);
    pause;
  else
    error('AssertionFailed' + sx);
  end
end
endfunction

