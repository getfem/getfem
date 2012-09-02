function result = rebuild_lib_windows(path, lib_name, machine, vc_version)
  // path: current working path
  // lib_name: scilab dll name (without extension) to be reconstructed as a lib 
  // vc_version: version of Visual studio (10.0 by default)
  
  if ~isdef('vc_version') then
    vc_version = '10.0';
  end

  if ~isdef('machine') then
    vc_version = 'X86';
  end

  if getos()=='Windows' then
    try
      value = winqueryreg('HKEY_LOCAL_MACHINE', 'SOFTWARE\Wow6432Node\Microsoft\VisualStudio\' + vc_version + '\Setup\VC\', 'ProductDir');
    catch
      printf('Error: can''t find Visual Studio %s\n', vc_version);
      result = %f;
    end
    msvc_dir = """" + value + 'bin' + filesep();
    
    filename = SCI + filesep() + 'bin' + filesep() + lib_name + '.dll';
    
    dllinfolist = dllinfo(filename,'exports');
    
    fd = mopen(path + lib_name + '.def','w');
    mputl('EXPORTS',fd);
    mputl(dllinfolist(2), fd);
    mclose(fd);
    
    [output,bOK,result] = dos(msvc_dir + 'vcvars32.bat""','-echo')
    [output,bOK,result] = dos(msvc_dir + 'lib.exe"" /machine:' + machine + ' /def:' + path + lib_name + '.def /out:' + path + lib_name + '.lib','-echo')
    result = (result==0);
  else
    result = %f;
  end
endfunction

