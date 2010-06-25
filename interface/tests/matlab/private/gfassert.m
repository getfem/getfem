function gfassert(sx)
  global gverbose;
  global gdebug;
  x = evalin('caller',sx);
  if (~all(x(:))),
    if (gverbose)
      dbstack;      
    end;
    if (gdebug)      
      disp(['Assertion failed:' sx]);
      keyboard;
    else
      error('AssertionFailed',['Assertion failed' sx]);
    end;
  end;
  
