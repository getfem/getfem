// ====================================================================
// Yann COLLETTE
// Copyright 2009
// This file is released into the public domain
// ====================================================================
mode(-1);
lines(0);
try
 v = getversion('scilab');
catch
 error(gettext('Scilab 5.3.x or more is required.'));  
end;
if v(1) < 6 then
 error(gettext('Scilab 5.3.x or more is required.'));  
end
// ====================================================================
if ~with_module('development_tools') then
  error(msprintf(gettext('%s module not installed.'),'development_tools'));
end
// ====================================================================
TOOLBOX_NAME = 'sci_getfem';
TOOLBOX_TITLE = 'SciGetFem';
// ====================================================================
toolbox_dir = get_absolute_file_path('builder.sce');

// Under Windows, configure is not launched.
if getos()=='Windows' then
  copyfile(pwd() + '/sci_gateway/c/builder_gateway_c.sce.in',pwd() + '/sci_gateway/c/builder_gateway_c.sce');
end

tbx_builder_macros(toolbox_dir);
tbx_builder_src(toolbox_dir);
tbx_builder_gateway(toolbox_dir);
tbx_builder_help(toolbox_dir);
tbx_build_loader(TOOLBOX_NAME, toolbox_dir);
tbx_build_cleaner(TOOLBOX_NAME, toolbox_dir);

clear toolbox_dir TOOLBOX_NAME TOOLBOX_TITLE;
// ====================================================================
