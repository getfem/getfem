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
#include <getfemint_workspace.h>
#include <getfem/getfem_assembling.h>
#include <gmm/gmm_inoutput.h>

using namespace getfemint;

template <typename MAT>
void copydiags(const MAT &M, const std::vector<size_type> &v,
	       garray<typename MAT::value_type> &w) {
  size_type m = gmm::mat_nrows(M), n = gmm::mat_ncols(M);
  for (unsigned ii=0; ii < v.size(); ++ii) {
    int d = int(v[ii]), i,j;
    if (d < 0) { i = -d; j = 0; } else { i = 0; j = d; }
    cout << "m=" << m << "n=" << n << ", d=" << d << ", i=" << i << ", j=" << j << "\n";
    for (; i < int(m) && j < int(n); ++i,++j)
      w(i,ii) = M(i,j);
  }
}

template <typename T> static void
gf_spmat_get_full(gsparse &gsp, getfemint::mexargs_in& in, getfemint::mexargs_out& out, T) {
  gmm::dense_matrix<T> ww;
  size_type n,m;
  if (!in.remaining()) {
    m = gsp.nrows(); n = gsp.ncols();
    gmm::resize(ww, m, n);
    switch (gsp.storage()) {
    case gsparse::CSCMAT: gmm::copy(gsp.csc(T()), ww); break;
    case gsparse::WSCMAT: gmm::copy(gsp.wsc(T()), ww); break;
    default: THROW_INTERNAL_ERROR;
    }
  } else {
    sub_index ii = in.pop().to_sub_index().check_range(gsp.nrows());
    sub_index jj = in.remaining() ?
      in.pop().to_sub_index().check_range(gsp.ncols()) : ii.check_range(gsp.ncols());
    m = ii.size(); n = jj.size();
    gmm::resize(ww, m, n);
    switch (gsp.storage()) {
    case gsparse::CSCMAT: gmm::copy(gmm::sub_matrix(gsp.csc(T()),ii,jj), ww); break;
    case gsparse::WSCMAT: gmm::copy(gmm::sub_matrix(gsp.wsc(T()),ii,jj), ww); break;
    default: THROW_INTERNAL_ERROR;
    }
  }
  std::copy(ww.begin(), ww.end(), out.pop().create_array(int(m),int(n),T()).begin());
}

template <typename T> static void
gf_spmat_mult_or_tmult(gsparse &gsp,
		       getfemint::mexargs_in& in, getfemint::mexargs_out& out,
		       bool tmult, T) {
  size_type nj = gsp.ncols(), ni = gsp.nrows();
  if (tmult) std::swap(ni,nj);
  //cout << "NJ=" << nj << "NI=" << ni << ", tmaul=" << tmult << "\n";
  garray<T> v = in.pop().to_garray(int(nj),T());
  garray<T> w = out.pop().create_array_v(unsigned(ni), T());
  gsp.mult_or_transposed_mult(v,w,tmult);
  /*switch (gsp.storage()) {
  case gsparse::CSCMAT:
    if (!tmult) gmm::mult(gsp.csc(T()), v, w);
    else gmm::mult(gmm::conjugated(gsp.csc(T())), v, w);
    break;
  case gsparse::WSCMAT:
    if (!tmult) gmm::mult(gsp.wsc(T()), v, w);
    else gmm::mult(gmm::conjugated(gsp.wsc(T())), v, w);
    break;
  default: THROW_INTERNAL_ERROR;
  }*/

}

template <typename T> static void
gf_spmat_get_diag(gsparse &gsp,
		  getfemint::mexargs_in& in, getfemint::mexargs_out& out, T) {
  std::vector<size_type> v;
  if (in.remaining()) {
    iarray vv = in.pop().to_iarray(-1);
    for (unsigned i=0; i < vv.size(); ++i)
      v.push_back(vv[i]);
  } else v.push_back(0);
  garray<T> w = out.pop().create_array(unsigned(std::min(gsp.nrows(), gsp.ncols())), unsigned(v.size()), T());
  switch (gsp.storage()) {
  case gsparse::CSCMAT: copydiags(gsp.csc(T()), v, w); break;
  case gsparse::WSCMAT: copydiags(gsp.wsc(T()), v, w); break;
  default: THROW_INTERNAL_ERROR;
  }
}

template <typename T> static void
gf_spmat_get_data(gmm::csc_matrix_ref<const T*, const unsigned int *, const unsigned int *> M,
		  getfemint::mexargs_out& out, int which) {
  size_type nz = M.jc[M.nc];
  if (which == 0) {
    iarray w = out.pop().create_iarray_h(unsigned(M.nc+1));
    for (unsigned i=0; i < M.nc+1; ++i)
      { w[i] = M.jc[i] + config::base_index(); }
    if (out.remaining()) {
      w = out.pop().create_iarray_h(unsigned(nz));
      for (unsigned i=0; i < nz; ++i)
	{ w[i] = M.ir[i] + config::base_index(); }
    }
  } else {
    garray<T> w = out.pop().create_array_h(unsigned(nz), T());
    for (unsigned i=0; i < M.nc+1; ++i) { w[i] = M.pr[i]; }
  }
}

template <typename T> static void
gf_spmat_get_Dirichlet_nullspace(gsparse &H, getfemint::mexargs_in& in, getfemint::mexargs_out& out, T) {
  garray<T> R            = in.pop().to_garray(T());
  gmm::col_matrix<gmm::wsvector<T> > NS(H.ncols(), H.nrows());
  std::vector<T> Ud(H.ncols()), R2(R.begin(), R.end());
  size_type nl;
  switch (H.storage()) {
    case gsparse::WSCMAT: nl = getfem::Dirichlet_nullspace(H.wsc(T()), NS, R2, Ud); break;
    case gsparse::CSCMAT: nl = getfem::Dirichlet_nullspace(H.csc(T()), NS, R2, Ud); break;
    default: THROW_INTERNAL_ERROR;
  }
  gmm::resize(NS,gmm::mat_nrows(NS),nl); /* remove unused columns */
  out.pop().from_sparse(NS);
  out.pop().from_dcvector(Ud);
}

/*MLABCOM
  FUNCTION [...]=gf_spmat_get(M, args)

  General getfem sparse matrix inquiry function. M might also be a
  native matlab sparse matrix.

  @GET SPMAT:GET('size')
  @GET SPMAT:GET('nnz')
  @GET SPMAT:GET('is_complex')
  @GET SPMAT:GET('storage')
  @GET SPMAT:GET('full')
  @GET SPMAT:GET('mult')
  @GET SPMAT:GET('tmult')
  @GET SPMAT:GET('diag')
  @GET SPMAT:GET('csc_ind')
  @GET SPMAT:GET('csc_val')
  @GET SPMAT:GET('dirichlet nullspace')
  @GET SPMAT:GET('info')
  @GET SPMAT:GET('save')
MLABCOM*/

void gf_spmat_get(getfemint::mexargs_in& in, getfemint::mexargs_out& out)
{
  if (in.narg() < 2) {
    THROW_BADARG( "Wrong number of input arguments");
  }
  dal::shared_ptr<gsparse> pgsp = in.pop().to_sparse();
  gsparse &gsp = *pgsp;

  std::string cmd = in.pop().to_string();

  if (check_cmd(cmd, "nnz", in, out, 0, 0, 0, 1)) {
    /*@GET n = SPMAT:GET('nnz')
    Return the number of non-null values stored in the sparse matrix.@*/
    out.pop().from_integer(int(gsp.nnz()));
  } else if (check_cmd(cmd, "full", in, out, 0, 2, 0, 1)) {
    /*@GET Sm = SPMAT:GET('full'[, @list I[, @list J]])
    Return a full (sub-)matrix.

    The optional arguments `I` and `J`, are the sub-intervals for the
    rows and columns that are to be extracted.@*/
    if (gsp.is_complex()) gf_spmat_get_full(gsp, in, out, complex_type());
    else gf_spmat_get_full(gsp, in,out,scalar_type());
  } else if (check_cmd(cmd, "mult", in, out, 1, 1, 0, 1)) {
    /*@GET MV = SPMAT:GET('mult',@vec V)
    Product of the sparse matrix `M` with a vector `V`.

    For matrix-matrix multiplications, see SPMAT:INIT('mult').@*/
    if (!gsp.is_complex())
      gf_spmat_mult_or_tmult(gsp, in, out, false, scalar_type());
    else gf_spmat_mult_or_tmult(gsp, in, out, false, complex_type());
  } else if (check_cmd(cmd, "tmult", in, out, 1, 1, 0, 1)) {
    /*@GET MtV = SPMAT:GET('tmult',@vec V)
    Product of `M` transposed (conjugated if `M` is complex) with the vector `V`.@*/
    if (!gsp.is_complex())
      gf_spmat_mult_or_tmult(gsp, in, out, true, scalar_type());
    else gf_spmat_mult_or_tmult(gsp, in, out, true, complex_type());
  } else if (check_cmd(cmd, "diag", in, out, 0, 1, 0, 1)) {
    /*@GET D = SPMAT:GET('diag'[, @list E])
    Return the diagonal of `M` as a vector.

    If `E` is used, return the sub-diagonals whose ranks are given in E.@*/
    if (!gsp.is_complex())
      gf_spmat_get_diag(gsp, in, out, scalar_type());
    else gf_spmat_get_diag(gsp, in, out, complex_type());
  } else if (check_cmd(cmd, "storage", in, out, 0, 0, 0, 1)) {
    /*@GET s = SPMAT:GET('storage')
    Return the storage type currently used for the matrix.

    The storage is returned as a string, either 'CSC' or 'WSC'.@*/
    out.pop().from_string(gsp.name());
  } else if (check_cmd(cmd, "size", in, out, 0, 0, 0, 1)) {
    /*@GET @CELL{ni,nj} = SPMAT:GET('size')
    Return a vector where `ni` and `nj` are the dimensions of the matrix.@*/
    iarray sz = out.pop().create_iarray_h(2);
    sz[0] = int(gsp.nrows());
    sz[1] = int(gsp.ncols());
  } else if (check_cmd(cmd, "is_complex", in, out, 0, 0, 0, 1)) {
    /*@GET b = SPMAT:GET('is_complex')
    Return 1 if the matrix contains complex values.@*/
    out.pop().from_integer(gsp.is_complex());
  } else if (check_cmd(cmd, "csc_ind", in, out, 0, 0, 0, 2)) {
    /*@GET @CELL{JC, IR} = SPMAT:GET('csc_ind')
    Return the two usual index arrays of CSC storage.

    If `M` is not stored as a CSC matrix, it is converted into CSC.@*/
    gsp.to_csc();
    if (!gsp.is_complex()) gf_spmat_get_data(gsp.csc(scalar_type()),  out, 0);
    else                   gf_spmat_get_data(gsp.csc(complex_type()), out, 0);
  } else if (check_cmd(cmd, "csc_val", in, out, 0, 0, 0, 1)) {
    /*@GET V = SPMAT:GET('csc_val')
    Return the array of values of all non-zero entries of `M`.

    If `M` is not stored as a CSC matrix, it is converted into CSC.@*/
    gsp.to_csc();
    if (!gsp.is_complex()) gf_spmat_get_data(gsp.csc(scalar_type()),  out, 1);
    else                   gf_spmat_get_data(gsp.csc(complex_type()), out, 1);
  } else if (check_cmd(cmd, "dirichlet nullspace", in, out, 1, 1, 2, 2)) {
    /*@GET @CELL{N, U0} = SPMAT:GET('dirichlet nullspace',@vec R)
    Solve the dirichlet conditions `M.U=R`.

    A solution `U0` which has a minimum L2-norm is returned, with a
    sparse matrix `N` containing an orthogonal basis of the kernel of
    the (assembled) constraints matrix `M` (hence, the PDE linear system
    should be solved on this subspace): the initial problem<Par>

    `K.U = B` with constraints `M.U = R`<Par>

    is replaced by<Par>

    `(N'.K.N).UU = N'.B` with `U = N.UU + U0`@*/
    if (pgsp->is_complex())
      gf_spmat_get_Dirichlet_nullspace(*pgsp, in, out, complex_type());
    else
      gf_spmat_get_Dirichlet_nullspace(*pgsp, in, out, scalar_type());
  } else if (check_cmd(cmd, "info", in, out, 0, 1)) {
    /*@GET s = SPMAT:GET('info')
    Return a string contains a short summary on the sparse matrix (dimensions, filling, ...).@*/
    std::stringstream ss;
    size_type ncc = gsp.nrows()*gsp.ncols();
    ss << gsp.nrows() << "x" << gsp.ncols() << " "
       << (gsp.is_complex() ? "COMPLEX" : "REAL") << " " << gsp.name()
       << ", NNZ=" << gsp.nnz() << " (filling="
       << 100.*double(gsp.nnz())/(double(ncc == 0 ? 1 : ncc)) << "%)";
    out.pop().from_string(ss.str().c_str());
  } else if (check_cmd(cmd, "save", in, out, 2, 2, 0, 0)) {
    /*@GET SPMAT:GET('save',@str format, @str filename)
    Export the sparse matrix.

    the format of the file may be 'hb' for Harwell-Boeing, or 'mm'
    for Matrix-Market.@*/
    std::string fmt = in.pop().to_string();
    int ifmt;
    if (cmd_strmatch(fmt, "hb") || cmd_strmatch(fmt, "harwell-boeing")) ifmt = 0;
    else if (cmd_strmatch(fmt, "mm") || cmd_strmatch(fmt, "matrix-market")) ifmt = 1;
    else THROW_BADARG("unknown sparse matrix file-format : " << fmt);
    std::string fname = in.pop().to_string();
    gsp.to_csc();
    if (!gsp.is_complex()) {
      if (ifmt == 0) gmm::Harwell_Boeing_save(fname.c_str(), gsp.csc(scalar_type()));
      else           gmm::MatrixMarket_save(fname.c_str(), gsp.csc(scalar_type()));
    } else {
      if (ifmt == 0) gmm::Harwell_Boeing_save(fname.c_str(), gsp.csc(complex_type()));
      else           gmm::MatrixMarket_save(fname.c_str(), gsp.csc(complex_type()));
    }
  } else bad_cmd(cmd);
}
