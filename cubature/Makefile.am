
IM_METHODS=

all : im_list.mak getfem_im_list.h

im_list.mak :
	echo IM_METHODS=`ls $(top_srcdir)/cubature/*.IM` > im_list.mak

# DEP_FILES = im_list.mak
@INCLUDE_ANYWAY@ im_list.mak

EXTRA_DIST = $(IM_METHODS) getfem_im_list.h make_getfem_im_list im_list.mak
getfem_im_list.h :  $(IM_METHODS)
	$(top_srcdir)/cubature/make_getfem_im_list $(top_srcdir)/cubature
	cp ./getfem_im_list.h ../src/getfem_im_list.h
	echo IM_METHODS=`ls $(top_srcdir)/cubature/*.IM` > im_list.mak