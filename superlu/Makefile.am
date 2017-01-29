# Copyright (c) 2003, The Regents of the University of California, through
# Lawrence Berkeley National Laboratory (subject to receipt of any required 
# approvals from U.S. Dept. of Energy) 
# 
# All rights reserved. 
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met: 
# 
# (1) Redistributions of source code must retain the above copyright notice,
# this list of conditions and the following disclaimer. 
# (2) Redistributions in binary form must reproduce the above copyright notice,
# this list of conditions and the following disclaimer in the documentation
# and/or other materials provided with the distribution. 
# (3) Neither the name of Lawrence Berkeley National Laboratory, U.S. Dept. of
# Energy nor the names of its contributors may be used to endorse or promote
# products derived from this software without specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
# IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
# PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. 




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

EXTRA_DIST=License.txt
