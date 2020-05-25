// ====================================================================
// Yann COLLETTE
// Copyright 2009-2020
// This file is released into the public domain
// ====================================================================

path = get_absolute_file_path('buildmacros.sce');

genlib('sci_getfemlib',path,%f,%t);
genlib('sci_getfemoverloadlib',path + filesep() + 'overload',%f,%t);

clear path;
