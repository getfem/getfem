#!/bin/bash
# run this file in cygwin to generate the M-files, and pack them with the gf_matlab.mexw32 dll 
# so that they can be distributed..
set -e # exit on error

if [ -d getfem_matlab_toolbox ]; then 
  rm -rf getfem_matlab_toolbox 
fi
mkdir getfem_matlab_toolbox 
cd getfem_matlab_toolbox 
../../bin/extract_doc ../../interface/src matlab-com
cp -p ../../interface/src/matlab/*.m .
mkdir private
cp -p ../../interface/src/matlab/private/*.m ./private
cp ../Release/gf_matlab.mexw32 .
cd ..
if [ -f getfem_matlab_toolbox.zip ]; then 
  rm getfem_matlab_toolbox.zip
fi
zip -r getfem_matlab_toolbox.zip getfem_matlab_toolbox
