#
# Copyright (c) 2003, The Regents of the University of California, through
# Lawrence Berkeley National Laboratory (subject to receipt of any required 
# approvals from U.S. Dept. of Energy) 
#
# All rights reserved. 
#
# The source code is distributed under BSD license, see the file License.txt
#


#  Copyright (C) 2004-2020 Yves Renard
#
#  This file is a part of GetFEM++
#
#  GetFEM++  is  free software;  you  can  redistribute  it  and/or modify it
#  under  the  terms  of the  GNU  Lesser General Public License as published
#  by  the  Free Software Foundation;  either version 3 of the License,  or
#  (at your option) any later version along with the GCC Runtime Library
#  Exception either version 3.1 or (at your option) any later version.
#  This program  is  distributed  in  the  hope  that it will be useful,  but
#  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
#  or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
#  License and GCC Runtime Library Exception for more details.
#  You  should  have received a copy of the GNU Lesser General Public License
#  along  with  this program;  if not, write to the Free Software Foundation,
#  Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.

noinst_HEADERS=\
	slu_Cnames.h slu_dcomplex.h slu_scomplex.h slu_util.h\
	supermatrix.h colamd.h slu_cdefs.h slu_ddefs.h\
	slu_sdefs.h slu_zdefs.h



if USEBLASLITE
BLASLITEFILES=BLAS.c BLAS_f2c.h f2c_lite.c
else
BLASLITEFILES=
endif

SRC = ccolumn_bmod.c \
	ccolumn_dfs.c \
	ccopy_to_ucol.c \
	cgscon.c \
	cgsequ.c \
	cgsrfs.c \
	cgssv.c \
	cgssvx.c \
	cgstrf.c \
	cgstrs.c \
	clacon.c \
	clangs.c \
	claqgs.c \
	cmemory.c \
	cmyblas2.c \
	colamd.c \
	cpanel_bmod.c \
	cpanel_dfs.c \
cpivotgrowth.c \
cpivotL.c \
cpruneL.c \
creadhb.c \
csnode_bmod.c \
csnode_dfs.c \
csp_blas2.c \
csp_blas3.c \
cutil.c \
dcolumn_bmod.c \
dcolumn_dfs.c \
dcomplex.c \
dcopy_to_ucol.c \
dgscon.c \
dgsequ.c \
dgsrfs.c \
dgssv.c \
dgssvx.c \
dgstrf.c \
dgstrs.c \
dlacon.c \
dlamch.c \
dlangs.c \
dlaqgs.c \
dmemory.c \
dmyblas2.c \
dpanel_bmod.c \
dpanel_dfs.c \
dpivotgrowth.c \
dpivotL.c \
dpruneL.c \
dreadhb.c \
dsnode_bmod.c \
dsnode_dfs.c \
dsp_blas2.c \
dsp_blas3.c \
dutil.c \
dzsum1.c \
get_perm_c.c \
heap_relax_snode.c \
icmax1.c \
izmax1.c \
lsame.c \
memory.c \
mmd.c \
relax_snode.c \
scolumn_bmod.c \
scolumn_dfs.c \
scomplex.c \
scopy_to_ucol.c \
scsum1.c \
sgscon.c \
sgsequ.c \
sgsrfs.c \
sgssv.c \
sgssvx.c \
sgstrf.c \
sgstrs.c \
slacon.c \
slamch.c \
slangs.c \
slaqgs.c \
smemory.c \
smyblas2.c \
spanel_bmod.c \
spanel_dfs.c \
sp_coletree.c \
sp_ienv.c \
spivotgrowth.c \
spivotL.c \
sp_preorder.c \
spruneL.c \
sreadhb.c \
ssnode_bmod.c \
ssnode_dfs.c \
ssp_blas2.c \
ssp_blas3.c \
superlu_timer.c \
sutil.c \
util.c \
zcolumn_bmod.c \
zcolumn_dfs.c \
zcopy_to_ucol.c \
zgscon.c \
zgsequ.c \
zgsrfs.c \
zgssv.c \
zgssvx.c \
zgstrf.c \
zgstrs.c \
zlacon.c \
zlangs.c \
zlaqgs.c \
zmemory.c \
zmyblas2.c \
zpanel_bmod.c \
zpanel_dfs.c \
zpivotgrowth.c \
zpivotL.c \
zpruneL.c \
zreadhb.c \
zsnode_bmod.c \
zsnode_dfs.c \
zsp_blas2.c \
zsp_blas3.c \
zutil.c $(BLASLITEFILES)

#vire: xerbla.c 


noinst_LTLIBRARIES = libsuperlu.la
libsuperlu_la_SOURCES = $(SRC)
#libsuperlu_la_LDFLAGS = ${LIBTOOL_VERSION_INFO}
libsuperlu_la_CPPFLAGS = @SUPERLU_CPPFLAGS@

CLEANFILES = ii_files/* *.o.d

EXTRA_DIST=License.txt \
BLAS/License.txt \
BLAS/caxpy.f \
BLAS/ccopy.f \
BLAS/cdotc.f \
BLAS/cdotu.f \
BLAS/cgbmv.f \
BLAS/cgemm.f \
BLAS/cgemv.f \
BLAS/cgerc.f \
BLAS/cgeru.f \
BLAS/chbmv.f \
BLAS/chemm.f \
BLAS/chemv.f \
BLAS/cher2.f \
BLAS/cher2k.f \
BLAS/cher.f \
BLAS/cherk.f \
BLAS/chpmv.f \
BLAS/chpr2.f \
BLAS/chpr.f \
BLAS/crotg.f \
BLAS/cscal.f \
BLAS/csrot.f \
BLAS/csscal.f \
BLAS/cswap.f \
BLAS/csymm.f \
BLAS/csyr2k.f \
BLAS/csyrk.f \
BLAS/ctbmv.f \
BLAS/ctbsv.f \
BLAS/ctpmv.f \
BLAS/ctpsv.f \
BLAS/ctrmm.f \
BLAS/ctrmv.f \
BLAS/ctrsm.f \
BLAS/ctrsv.f \
BLAS/dasum.f \
BLAS/daxpy.f \
BLAS/dcabs1.f \
BLAS/dcopy.f \
BLAS/ddot.f \
BLAS/dgbmv.f \
BLAS/dgemm.f \
BLAS/dgemv.f \
BLAS/dger.f \
BLAS/dnrm2.f \
BLAS/drot.f \
BLAS/drotg.f \
BLAS/drotm.f \
BLAS/drotmg.f \
BLAS/dsbmv.f \
BLAS/dscal.f \
BLAS/dsdot.f \
BLAS/dspmv.f \
BLAS/dspr2.f \
BLAS/dspr.f \
BLAS/dswap.f \
BLAS/dsymm.f \
BLAS/dsymv.f \
BLAS/dsyr2.f \
BLAS/dsyr2k.f \
BLAS/dsyr.f \
BLAS/dsyrk.f \
BLAS/dtbmv.f \
BLAS/dtbsv.f \
BLAS/dtpmv.f \
BLAS/dtpsv.f \
BLAS/dtrmm.f \
BLAS/dtrmv.f \
BLAS/dtrsm.f \
BLAS/dtrsv.f \
BLAS/dzasum.f \
BLAS/dznrm2.f \
BLAS/icamax.f \
BLAS/idamax.f \
BLAS/isamax.f \
BLAS/izamax.f \
BLAS/lsame.f \
BLAS/sasum.f \
BLAS/saxpy.f \
BLAS/scabs1.f \
BLAS/scasum.f \
BLAS/scnrm2.f \
BLAS/scopy.f \
BLAS/sdot.f \
BLAS/sdsdot.f \
BLAS/sgbmv.f \
BLAS/sgemm.f \
BLAS/sgemv.f \
BLAS/sger.f \
BLAS/snrm2.f \
BLAS/srot.f \
BLAS/srotg.f \
BLAS/srotm.f \
BLAS/srotmg.f \
BLAS/ssbmv.f \
BLAS/sscal.f \
BLAS/sspmv.f \
BLAS/sspr2.f \
BLAS/sspr.f \
BLAS/sswap.f \
BLAS/ssymm.f \
BLAS/ssymv.f \
BLAS/ssyr2.f \
BLAS/ssyr2k.f \
BLAS/ssyr.f \
BLAS/ssyrk.f \
BLAS/stbmv.f \
BLAS/stbsv.f \
BLAS/stpmv.f \
BLAS/stpsv.f \
BLAS/strmm.f \
BLAS/strmv.f \
BLAS/strsm.f \
BLAS/strsv.f \
BLAS/xerbla_array.f \
BLAS/xerbla.f \
BLAS/zaxpy.f \
BLAS/zcopy.f \
BLAS/zdotc.f \
BLAS/zdotu.f \
BLAS/zdrot.f \
BLAS/zdscal.f \
BLAS/zgbmv.f \
BLAS/zgemm.f \
BLAS/zgemv.f \
BLAS/zgerc.f \
BLAS/zgeru.f \
BLAS/zhbmv.f \
BLAS/zhemm.f \
BLAS/zhemv.f \
BLAS/zher2.f \
BLAS/zher2k.f \
BLAS/zher.f \
BLAS/zherk.f \
BLAS/zhpmv.f \
BLAS/zhpr2.f \
BLAS/zhpr.f \
BLAS/zrotg.f \
BLAS/zscal.f \
BLAS/zswap.f \
BLAS/zsymm.f \
BLAS/zsyr2k.f \
BLAS/zsyrk.f \
BLAS/ztbmv.f \
BLAS/ztbsv.f \
BLAS/ztpmv.f \
BLAS/ztpsv.f \
BLAS/ztrmm.f \
BLAS/ztrmv.f \
BLAS/ztrsm.f \
BLAS/ztrsv.f
