// ====================================================================
// Yann COLLETTE
// Copyright 2009
// This file is released into the public domain
// ====================================================================


jdk_opts = "-Djdk.xml.totalEntitySizeLimit=1000000 -Djdk.xml.maxGeneralEntitySizeLimit=1000000 -Djdk.xml.entityExpansionLimit=0";
if getenv('_JAVA_OPTIONS') == [] then
  setenv('_JAVA_OPTIONS', jdk_opts);
else
  setenv('_JAVA_OPTIONS', getenv('_JAVA_OPTIONS') + " " + jdk_opts);
end

help_dir = get_absolute_file_path('builder_help.sce');
tbx_builder_help_lang("en_US", help_dir);
//tbx_builder_help_lang("fr_FR", help_dir);

clear help_dir;
