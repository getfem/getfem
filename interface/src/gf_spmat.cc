// -*- c++ -*- (enables emacs c++ mode)
//===========================================================================
//
// Copyright (C) 2006-2008 Yves Renard, Julien Pommier.
//
// This file is a part of GETFEM++
//
// Getfem++  is  free software;  you  can  redistribute  it  and/or modify it
// under  the  terms  of the  GNU  Lesser General Public License as published
// by  the  Free Software Foundation;  either version 2.1 of the License,  or
// (at your option) any later version.
// This program  is  distributed  in  the  hope  that it will be useful,  but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
// License for more details.
// You  should  have received a copy of the GNU Lesser General Public License
// along  with  this program;  if not, write to the Free Software Foundation,
// Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.
//
//===========================================================================

#include <getfemint_gsparse.h>
#include <getfemint_gsparse_misc.h>
#include <getfemint_workspace.h>
#include <gmm/gmm_inoutput.h>

using namespace getfemint;

template <typename TA, typename TB> static void
gf_spmat_add(gsparse &res, gsparse &A, gsparse &B, TA, TB) {
  switch (B.storage()) {
  case gsparse::CSCMAT: gmm::copy(B.csc(TB()), res.wsc(TB())); break;
  case gsparse::WSCMAT: gmm::copy(B.wsc(TB()), res.wsc(TB())); break;
  default: THROW_INTERNAL_ERROR;
  }
  switch (A.storage()) {
  case gsparse::CSCMAT: gmm::add(A.csc(TA()),res.wsc(TB())); break;
  case gsparse::WSCMAT: gmm::add(A.wsc(TA()),res.wsc(TB())); break;
  default: THROW_INTERNAL_ERROR;
  }
}

#if 0
template <typename T> static void
gf_spmat_Dirichlet_nullspace(gsparse &gsp,
			     getfemint::mexargs_in& in, getfemint::mexargs_out& out, T) {
  garray<T> R            = in.pop().to_garray(T());
  
  gf_real_sparse_by_col      NS(gmm::mat_ncols(H), gmm::mat_nrows(H));
    std::vector<double>   Ud(H.ncols());
    
    size_type nl = getfem::Dirichlet_nullspace(H, NS, 
					       R.to_vector<std::vector<double> >(), Ud);
    gmm::resize(NS,gmm::mat_nrows(NS),nl); /* remove unused columns */
    //NS.resize(nl);
    out.pop().from_sparse(NS);
    out.pop().from_dcvector(Ud);
}

void
gf_spmat_Dirichlet_nullspace(gsparse &gsp,
				 getfemint::mexargs_in& in, getfemint::mexargs_out& out) {
  if (gsp.is_complex())
    gf_spmat_Dirichlet_nullspace(gsp, in, out, complex_type());
  else gf_spmat_get_Dirichlet_nullspace(gsp, in, out, scalar_type());      
}
#endif

template <typename T> static void
copy_spmat(gsparse &src, gsparse &dest, mexargs_in &in, T) {
  if (!in.remaining()) {
    dest.allocate(src.nrows(), src.ncols(), src.storage(), T());
    switch (src.storage()) {
      case gsparse::WSCMAT: gmm::copy(src.wsc(T()), dest.wsc(T())); break;
      case gsparse::CSCMAT: gmm::copy(src.csc(T()), dest.csc_w(T())); break;
      default: THROW_INTERNAL_ERROR;
    }
  } else {
    sub_index ii = in.pop().to_sub_index().check_range(src.nrows());
    sub_index jj = in.remaining() ? 
      in.pop().to_sub_index().check_range(src.ncols()) : ii.check_range(src.ncols());
    size_type m = ii.size(), n = jj.size();
    dest.allocate(m, n, src.storage(), T());
    switch (src.storage()) {
      case gsparse::WSCMAT: gmm::copy(gmm::sub_matrix(src.wsc(T()), ii, jj), 
				      dest.wsc(T())); break;
      case gsparse::CSCMAT: gmm::copy(gmm::sub_matrix(src.csc(T()), ii, jj), 
				      dest.csc_w(T())); break;
      default: THROW_INTERNAL_ERROR;
    }
  }
}

void load_spmat(mexargs_in& in, gsparse &gsp) {
  std::string mt = in.pop().to_string();
  std::string fname = in.pop().to_string();
  if (cmd_strmatch(mt, "hb") || 
      cmd_strmatch(mt, "harwell-boeing")) {
    gmm::HarwellBoeing_IO h; h.open(fname.c_str());
    if (h.is_complex()) {
      gmm::csc_matrix<complex_type> cscH;
      h.read(cscH);
      gsp.destructive_assign(cscH);
    } else {
      gmm::csc_matrix<double> cscH;
      h.read(cscH);
      gsp.destructive_assign(cscH);
    }
  } else if (cmd_strmatch(mt, "mm") || 
	     cmd_strmatch(mt, "matrix-market")) {
    gmm::MatrixMarket_IO h; h.open(fname.c_str());
    if (h.is_complex()) {
      gf_cplx_sparse_by_col H; h.read(H); 
      gsp.destructive_assign(H);
    } else {
      gf_real_sparse_by_col H; h.read(H);
      gsp.destructive_assign(H);
    }
  } else THROW_BADARG("unknown sparse matrix file-format : " << mt);
}

/*MLABCOM
  FUNCTION F=gf_spmat(command, args)
  @TEXT SPMAT:INIT('SPMAT_init')
MLABCOM*/

/*@TEXT SPMAT:INIT('SPMAT_init')
  General constructor for getfem sparse matrices @MATLAB{(i.e. sparse matrices
  which are stored in the getfem workspace, not the matlab sparse
  matrices)}.<Par>

  These sparse matrix can be stored as CSC (compressed column sparse),
  which is the format used by Matlab, or they can be stored as WSC
  (internal format to getfem). The CSC matrices are not writable (it
  would be very inefficient), but they are optimized for
  multiplication with vectors, and memory usage. The WSC are writable,
  they are very fast with respect to random read/write
  operation. However their memory overhead is higher than CSC
  matrices, and they are a little bit slower for matrix-vector
  multiplications.<Par>

  By default, all newly created matrices are build as WSC
  matrices. This can be changed later with SPMAT:SET('to_csc'), or may
  be changed automatically by getfem (for example ::LINSOLVE()
  converts the matrices to CSC).<Par>

  The matrices may store REAL or COMPLEX values.<Par>

  * SPMAT:INIT('empty', @tint m [, @tint n])<par>
  Create a new empty (i.e. full of zeros) sparse matrix, of dimensions
  m x n.  If n is omitted, the matrix dimension is m x m.<Par>

  * SPMAT:INIT('identity', @tint n)<par>
  Create a n x n identity matrix.<Par>

  * SPMAT:SET('diag', @dmat D [, @ivec E [, @int n [,@int m]]])
  Create a diagonal matrix. If E is given, D might be a matrix and
  each column of E will contain the sub-diagonal number that will be
  filled with the corresponding column of D.<Par>

  * SPMAT:INIT('copy', @spmat K [,I [,J]])<par>

  Duplicate a matrix K (which might be a gfSpmat @MATLAB{or a native
  matlab sparse matrix}). If I and/or J are given, the matrix M will
  be a submatrix of K. For example<par> 
  M = gf_spmat('copy', sprand(50,50,.1), 1:40, [6 7 8 3 10])<par>
  will return a 40x5 matrix.<Par>

  * SPMAT:INIT('mult', @spmat A, @spmat B)<par>
  Create a sparse matrix as the product of the sparse matrices A and
  B.  It requires that A and B be both real or both complex, you may
  have to use SPMAT:SET('to_complex')<Par>

  * SPMAT:INIT('add', @spmat A, @spmat B)<par>
  Create a sparse matrix as the sum of the sparse matrices A and
  B. Adding a real matrix with a complex matrix is possible.<Par>

  * SPMAT:INIT('load','hb'|'harwell-boeing', filename)<par>
  Read a sparse matrix from an Harwell-Boeing file. 
  @MATLAB{See also ::UTIL('load matrix').}<Par>

  * SPMAT:INIT('mm', filename)<par>
  * SPMAT:INIT('matrix-market', filename)<par>
  Read a sparse matrix from a Matrix-Market file. 
  @MATLAB{See also ::UTIL('load matrix').}
  @*/


void gf_spmat(getfemint::mexargs_in& in, getfemint::mexargs_out& out)
{
  if (in.narg() < 1) {
    THROW_BADARG( "Wrong number of input arguments");
  }
  gsparse &gsp = out.pop().create_gsparse();
  std::string cmd = in.pop().to_string();
  if (check_cmd(cmd, "empty", in, out, 1, 2, 0, 1)) {
    size_type m = in.pop().to_integer(1,INT_MAX), n = m;
    if (in.remaining()) n = in.pop().to_integer(1, INT_MAX);
    gsp.allocate(m, n, gsparse::WSCMAT, gsparse::REAL);
  } else if (check_cmd(cmd, "copy", in, out, 1, 3, 0, 1)) {
    dal::shared_ptr<gsparse> A = in.pop().to_sparse();
    if (!A->is_complex()) {
      copy_spmat(*A, gsp, in, scalar_type());
    } else {
      copy_spmat(*A, gsp, in, complex_type());
    }
  } else if (check_cmd(cmd, "identity", in, out, 1, 1, 0, 1)) {
    size_type n = in.pop().to_integer(1, INT_MAX);
    gsp.real_wsc(new gsparse::t_wscmat_r(n,n));
    gmm::copy(gmm::identity_matrix(), gsp.real_wsc());
  } else if (check_cmd(cmd, "mult", in, out, 2, 2, 0, 1)) {
    dal::shared_ptr<gsparse> A = in.pop().to_sparse();
    dal::shared_ptr<gsparse> B = in.pop().to_sparse();
    size_type m = A->nrows(), n = B->ncols();
      
    if (A->is_complex() != B->is_complex()) 
      THROW_BADARG("cannot multiply a complex matrix with a real one, use to_complex()");
    if (!A->is_complex()) gsp.real_wsc(new gsparse::t_wscmat_r(m,n));
    else gsp.cplx_wsc(new gsparse::t_wscmat_c(m,n));
    if (A->storage() == gsparse::CSCMAT && B->storage() == gsparse::CSCMAT) {
      if (!A->is_complex()) gmm::mult(A->real_csc(),B->real_csc(), gsp.real_wsc());
      else                  gmm::mult(A->cplx_csc(),B->cplx_csc(), gsp.cplx_wsc());
    } else if (A->storage() == gsparse::CSCMAT && B->storage() == gsparse::WSCMAT) {
      if (!A->is_complex()) gmm::mult(A->real_csc(),B->real_wsc(), gsp.real_wsc());
      else                  gmm::mult(A->cplx_csc(),B->cplx_wsc(), gsp.cplx_wsc());
    } else if (A->storage() == gsparse::WSCMAT && B->storage() == gsparse::CSCMAT) {
      if (!A->is_complex()) gmm::mult(A->real_wsc(),B->real_csc(), gsp.real_wsc());
      else                  gmm::mult(A->cplx_wsc(),B->cplx_csc(), gsp.cplx_wsc());
    } else if (A->storage() == gsparse::WSCMAT && B->storage() == gsparse::WSCMAT) {
      if (!A->is_complex()) gmm::mult(A->real_wsc(),B->real_wsc(), gsp.real_wsc());
      else                  gmm::mult(A->cplx_wsc(),B->cplx_wsc(), gsp.cplx_wsc());
    } else THROW_INTERNAL_ERROR;
  } else if (check_cmd(cmd, "add", in, out, 2, 2, 0, 1)) {
    dal::shared_ptr<gsparse> A = in.pop().to_sparse();
    dal::shared_ptr<gsparse> B = in.pop().to_sparse();
    size_type m = A->nrows(), n = A->ncols();
    if (A->is_complex() != B->is_complex()) {
      gsp.cplx_wsc(new gsparse::t_wscmat_c(m,n));
      if (A->is_complex()) gf_spmat_add(gsp, *B, *A, scalar_type(), complex_type());
      else gf_spmat_add(gsp, *A, *B, scalar_type(), complex_type());
    } else if (A->is_complex()) {
      gsp.cplx_wsc(new gsparse::t_wscmat_c(m,n));
      gf_spmat_add(gsp, *A, *B, complex_type(), complex_type());
    } else {
      gsp.real_wsc(new gsparse::t_wscmat_r(m,n));
      gf_spmat_add(gsp, *A, *B, scalar_type(), scalar_type());
    }
  } else if (check_cmd(cmd, "diag", in, out, 1, 4, 0, 1)) {
    spmat_set_diag(gsp, in, true);
  } else if (check_cmd(cmd, "load", in, out, 2, 2, 1, 1)) {
    load_spmat(in, gsp);
  } else bad_cmd(cmd);
}

