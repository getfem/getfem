all: clean remake

clean:
	rm -rf getfem_reference && mkdir getfem_reference

remake:
	( cd ../.. && doxygen doc/doxygen/Doxyfile > /dev/null ) && export srcdir=../.. && ../../bin/upload_documentation --delete getfem_reference


#( cd ../.. && doxygen doc/doxygen/Doxyfile > /dev/null ) && tar czvf html_getfem_reference.tar.gz getfem_reference && cp html_getfem_reference.tar.gz ../../../getfem_html/ 
