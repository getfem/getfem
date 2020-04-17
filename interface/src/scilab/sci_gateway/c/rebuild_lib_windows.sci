// Copyright (C) 2009-2020 Yann Colette
// 
//  This file is a part of GetFEM++
// 
//  GetFEM++  is  free software;  you  can  redistribute  it  and/or modify it
//  under  the  terms  of the  GNU  Lesser General Public License as published
//  by  the  Free Software Foundation;  either version 3 of the License,  or
//  (at your option) any later version along with the GCC Runtime Library
//  Exception either version 3.1 or (at your option) any later version.
//  This program  is  distributed  in  the  hope  that it will be useful,  but
//  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
//  or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
//  License and GCC Runtime Library Exception for more details.
//  You  should  have received a copy of the GNU Lesser General Public License
//  along  with  this program;  if not, write to the Free Software Foundation,
//  Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.

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
		  if win64() then
        value = winqueryreg('HKEY_LOCAL_MACHINE', 'SOFTWARE\Microsoft\VisualStudio\' + vc_version + '\Setup\VC\', 'ProductDir');
			else
        value = winqueryreg('HKEY_LOCAL_MACHINE', 'SOFTWARE\Wow6432Node\Microsoft\VisualStudio\' + vc_version + '\Setup\VC\', 'ProductDir');
		  end
    catch
      printf('Error: can''t find Visual Studio %s\n', vc_version);
      result = %f;
			return;
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
