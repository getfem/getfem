all:
	rm -rf getfem_reference && mkdir getfem_reference
	time (( cd ../.. && doxygen doc/doxygen/Doxyfile > /dev/null ) && tar czvf html_getfem_reference.tar.gz getfem_reference && cp html_getfem_reference.tar.gz ../../../getfem_html/ )

remake:
	time (( cd ../.. && doxygen doc/doxygen/Doxyfile > /dev/null ) && tar czvf html_getfem_reference.tar.gz getfem_reference && cp html_getfem_reference.tar.gz ../../../getfem_html/ )