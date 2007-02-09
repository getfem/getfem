#!/bin/sh
echo "running checks..."
echo "  srcdir='$srcdir'"

export MATLABPATH=$(pwd)/../../src/matlab:$(pwd)/../../tests/matlab:${srcdir}:${srcdir}/../../src/matlab

echo "  export MATLAB_ROOT=$MATLAB_ROOT"
echo "  export MATLABPATH=$MATLABPATH"
#echo "  setenv LD_LIBRARY_PATH $LD_LIBRARY_PATH"
#echo "  setenv PATH $PATH"
#export PATH;

# the MATLABPATH is set here again in order to override a 
# possible ~/matlab/startup.m pointing to an old getfem release


if test x$MATLAB_ROOT = x; then
  MLAB=matlab;
else
  MLAB=${MATLAB_ROOT}/bin/matlab
fi

#echo "s=getenv('MATLABPATH'); while (length(s)), [a,s]=strtok(s,':'); addpath(a); end; disp(pwd); check_all; pause(1)" | ${MLAB} -nodesktop -nojvm -nosplash;
s=$(echo "s=getenv('MATLABPATH'); while (length(s)), [a,s]=strtok(s,':'); addpath(a); end; disp(pwd); check_all; pause(1)" | ${MLAB} -nodesktop -nojvm -nosplash 2>&1);

k=`echo "$s" | grep "All tests succeeded"`;
if test x"$k" = x""; then
  echo "failure: "
  echo "$s" | grep FAIL
  exit 1;
else
  echo "$s" | grep -i "succe";
fi;

