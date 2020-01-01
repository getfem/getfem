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

errcnt=0;
t = 'check_integ [integration methods]            ';
try
exec('check_integ.sce');
disp(['== ' t ': SUCCESS']);
catch
errcnt=errcnt+1; disp(['== ' t ': FAILURE']);
end;
t = 'check_fem [finite element methods]           ';
try
exec('check_fem.sce');
disp(['== ' t ': SUCCESS']);
catch
errcnt=errcnt+1; disp(['== ' t ': FAILURE']);
end;
t = 'check_geotrans [geometric transformations]   ';
try
exec('check_geotrans.sce');
disp(['== ' t ': SUCCESS']);
catch
errcnt=errcnt+1; disp(['== ' t ': FAILURE']);
end;
t = 'check_workspace [objects management]         ';
try
exec('check_workspace.sce');
disp(['== ' t ': SUCCESS']);
catch
errcnt=errcnt+1; disp(['== ' t ': FAILURE']);
end;
t = 'check_mesh_fem [mesh_fem manipulations]      ';
try
exec('check_mesh_fem.sce');
disp(['== ' t ': SUCCESS']);
catch
errcnt=errcnt+1; disp(['== ' t ': FAILURE']);
end;
t = 'check_asm [assembly routines]                ';
try
exec('check_asm.sce');
disp(['== ' t ': SUCCESS']);
catch
errcnt=errcnt+1; disp(['== ' t ': FAILURE']);
end;
t = 'check_slices [mesh slicing functions]        ';
try
exec('check_slices.sce');
disp(['== ' t ': SUCCESS']);
catch
errcnt=errcnt+1; disp(['== ' t ': FAILURE']);
end;
t = 'check_spmat [sparse matrix functions]        ';
try
exec('check_spmat.sce');
disp(['== ' t ': SUCCESS']);
catch
errcnt=errcnt+1; disp(['== ' t ': FAILURE']);
end;
if (errcnt),
  printf('\n\n== %d/11 tests FAILED\n', errcnt);
else
  printf('\n\n== All tests succeeded\n');
end;
disp('end of check_all..');
