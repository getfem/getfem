errcnt=0;
t = 'check_integ [integration methods]            ';
try
  check_integ;
  disp(['== ' t ': SUCCESS']);
catch
  errcnt=errcnt+1; disp(['== ' t ': FAILURE']);
end;

t = 'check_fem [finite element methods]           ';
try
  check_fem;
  disp(['== ' t ': SUCCESS']);
catch
  errcnt=errcnt+1; disp(['== ' t ': FAILURE']);
end;

t = 'check_geotrans [geometric transformations]   ';
try
  check_geotrans;
  disp(['== ' t ': SUCCESS']);
catch
  errcnt=errcnt+1; disp(['== ' t ': FAILURE']);
end;

t = 'check_workspace [objects management]         ';
try
  check_workspace;
  disp(['== ' t ': SUCCESS']);
catch
  errcnt=errcnt+1; disp(['== ' t ': FAILURE']);
end;

t = 'check_mesh_fem [mesh_fem manipulations]      ';
try
  check_mesh_fem;
  disp(['== ' t ': SUCCESS']);
catch
  errcnt=errcnt+1; disp(['== ' t ': FAILURE']);
end;

t = 'check_asm [assembly routines]                ';
try
  check_asm;
  disp(['== ' t ': SUCCESS']);
catch
  errcnt=errcnt+1; disp(['== ' t ': FAILURE']);
end;

t = 'check_solve                                  ';
try
  check_solve;
  disp(['== ' t ': SUCCESS']);
catch
  errcnt=errcnt+1; disp(['== ' t ': FAILURE']);
end;

t = 'check_oo [pseudo object oriented interface]  ';
try
  check_oo;
  disp(['== ' t ': SUCCESS']);
catch
  errcnt=errcnt+1; disp(['== ' t ': FAILURE']);
end;

t = 'check_slices [mesh slicing functions]        ';
try
  check_slices;
  disp(['== ' t ': SUCCESS']);
catch
  errcnt=errcnt+1; disp(['== ' t ': FAILURE']);
end;

t = 'check_spmat [sparse matrix functions]        ';
try
  check_spmat;
  disp(['== ' t ': SUCCESS']);
catch
  errcnt=errcnt+1; disp(['== ' t ': FAILURE']);
end;

t = 'check_bricks [model bricks functions]        ';
try
  check_bricks;
  disp(['== ' t ': SUCCESS']);
catch
  errcnt=errcnt+1; disp(['== ' t ': FAILURE']);
end;

if (errcnt),
  disp(sprintf('\n\n== %d/11 tests FAILED\n', errcnt));
else
  disp(sprintf('\n\n== All tests succeeded\n'));
end;
disp('end of check_all..');
