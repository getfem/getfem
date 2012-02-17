
all : getfem_im_list.h
check : getfem_im_list.h

EXTRA_DIST = $(IM_METHODS_LOC) make_getfem_im_list getfem_im_list.h

getfem_im_list.h :  $(IM_METHODS)
	$(top_srcdir)/cubature/make_getfem_im_list $(top_srcdir)/cubature
	cp ./getfem_im_list.h ../src/getfem/getfem_im_list.h



CLEANFILES = getfem_im_list.h
