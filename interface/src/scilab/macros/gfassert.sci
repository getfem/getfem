function gfassert(sx)
  global gverbose;
  global gdebug;
  x = eval(sx);
  if (~all(x(:))) then
    if (gverbose) then
      whereami();      
    end
    if (gdebug) then
      disp(['Assertion failed:' sx]);
      printf('enter ''continue'' to continue\n');
      pause;
    else
      printf('AssertionFailed:');
      disp(sx);
    end
  end
endfunction

