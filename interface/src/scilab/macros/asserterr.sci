function asserterr(sx)
global gdebug;

// ASSERTERR DOES NOT WORK FOR ASSIGNMENTS (i.e. sx='x=1')
// ONLY FOR EXPRESSIONS
//x = evalin('caller',sx, '[''catched'']');
ierr = execstr(sx, 'errcatch');
[LASTMSG, LASTID] = lasterror();
if (~ierr) then
  if (gdebug)      
    disp('Error triggering test failed:' + sx);
    pause;
  else
    error(sprintf('Error triggering test failed: %s', sx));
  end
end
endfunction

