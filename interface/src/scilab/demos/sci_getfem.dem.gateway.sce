// ====================================================================
// Yann COLLETTE
// Copyright 2009
// This file is released into the public domain
// ====================================================================
demopath = get_absolute_file_path("sci_getfem.dem.gateway.sce");

subdemolist = ["demo scilab_sum", "scilab_sum.dem.sce"; ..
               "demo cpp_find",   "cpp_find.dem.sce" ; ];

subdemolist(:,2) = demopath + subdemolist(:,2);
// ====================================================================
