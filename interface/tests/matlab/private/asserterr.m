function asserterr(sx)
  global gverbose;
  global gdebug;

  % ASSERTERR DOES NOT WORK FOR ASSIGNMENTS (i.e. sx='x=1')
  % ONLY FOR EXPRESSIONS
  x = evalin('caller',sx, '[''catched'']');
  [LASTMSG, LASTID] = lasterr;
  if (~strcmp(x,'catched') | ...
      strfind(lasterr,'internal error') | ...
      strfind(lasterr,'Internal')),
    if (gverbose)
      dbstack;      
    end;
    if (gdebug)      
      disp(['Error triggering test failed:' sx]);
      keyboard;
    else
      error('AssertionFailed',['Error triggering test failed' sx]);
    end;
  elseif (strcmp(LASTID, 'MATLAB:m_assignments_do_not_produce_results')),
    disp('WRONG USE OF asserterr');
    dbstack
    keyboard
  end;
