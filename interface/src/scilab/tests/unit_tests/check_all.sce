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
t = 'check_bricks [model bricks functions]        ';
try
exec('check_bricks.sce');
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
