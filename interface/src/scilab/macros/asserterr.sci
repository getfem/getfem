function asserterr(sx)
global gdebug;

// ASSERTERR DOES NOT WORK FOR ASSIGNMENTS (i.e. sx='x=1')
// ONLY FOR EXPRESSIONS
//x = evalin('caller',sx, '[''catched'']');
ierr = execstr(sx, 'errcatch');
[str,n,line,func]=lasterror();
[LASTMSG, LASTID] = lasterror();
if (~ierr) then
  if (gdebug)      
    disp('Error triggering test failed: ' + sx);
    disp('error: ' + str);
    disp(sprintf('error %d in %s at line %d\n', n, func, line));
    pause;
  else
    disp('Error triggering test failed: ' + sx);
    disp('error: ' + str);
    disp(sprintf('error %d in %s at line %d\n', n, func, line));
    error('');
  end
end
endfunction

