all:
	if [ -d getfem_reference ]; then "echo cleaning.."; rm -rf ./getfem_reference; fi;
	( cd ../.. && doxygen doc/doxygen/Doxyfile > /dev/null ) && (mv html getfem_reference && tar czvf html_getfem_reference.tar.gz getfem_reference )