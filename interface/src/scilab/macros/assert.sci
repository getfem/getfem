function assert(sx)
global gverbose;
global gdebug;

execstr('x = ' + sx);
if (~and(x(:))) then
  if (gverbose) then
    dbstack;      
  end
  if (gdebug) then
    disp('Assertion failed: ' + sx);
    pause;
  else
    error('Assertion failed: ' + sx);
  end
end
endfunction

