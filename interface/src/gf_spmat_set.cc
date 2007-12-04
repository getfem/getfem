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

using namespace getfemint;

typedef enum { TRANSP, CONJ, TRANSCONJ } gf_spmat_set_transpose_enum;

template <typename T> static void
gf_spmat_set_transpose(gsparse &gsp, 
		       gf_spmat_set_transpose_enum op, T) {
  size_type ni = gsp.nrows(), nj = gsp.ncols();
  if (op == CONJ) std::swap(ni,nj);
  gmm::row_matrix<gmm::rsvector<T> > m(nj,ni);
  switch (gsp.storage()) {
  case gsparse::CSCMAT: 
    switch (op) {
    case TRANSP:   gmm::copy(gmm::transposed(gsp.csc(T())), m); break;
    case TRANSCONJ:     gmm::copy(gmm::conjugated(gsp.csc(T())), m); break;
    case CONJ:gmm::copy(gmm::transposed(gmm::conjugated(gsp.csc(T()))), m); break;
    }
    gmm::resize(gsp.csc_w(T()), nj, ni);
    gmm::copy(m, gsp.csc_w(T()));
    break;
  case gsparse::WSCMAT: 
    switch (op) {
    case TRANSP:   gmm::copy(gmm::transposed(gsp.wsc(T())), m); break;
    case TRANSCONJ:     gmm::copy(gmm::conjugated(gsp.wsc(T())), m); break;
    case CONJ:gmm::copy(gmm::transposed(gmm::conjugated(gsp.wsc(T()))), m); break;
    }
    gmm::resize(gsp.wsc(T()), nj, ni);
    gmm::copy(m, gsp.wsc(T()));
    break;	
  default: THROW_INTERNAL_ERROR;
  }
}

template <typename T> static void
spmat_set_or_add_sub_matrix(gsparse &gsp, getfemint::mexargs_in& in,
			    gmm::sub_index ii, gmm::sub_index jj, bool do_add, T) {
  if (gsp.storage() != gsparse::WSCMAT)
    THROW_BADARG("cannot write to a CSC matrix (would be too inefficient). "
		 "Use to_wsc first");
  size_type m = ii.size(), n = jj.size();
  if (in.front().is_sparse()) {
    dal::shared_ptr<gsparse> src = in.pop().to_sparse();
    switch (src->storage()) {
      case gsparse::WSCMAT: 
	if (do_add) gmm::add (src->wsc(T()), gmm::sub_matrix(gsp.wsc(T()), ii, jj));
	else        gmm::copy(src->wsc(T()), gmm::sub_matrix(gsp.wsc(T()), ii, jj)); 
	break;
      case gsparse::CSCMAT: 
	if (do_add) gmm::add (src->csc(T()), gmm::sub_matrix(gsp.wsc(T()), ii, jj));
	else        gmm::copy(src->csc(T()), gmm::sub_matrix(gsp.wsc(T()), ii, jj)); break;
      default: THROW_INTERNAL_ERROR;
    }
  } else {
    garray<T> v = in.pop().to_garray(m,n,T());
    gmm::dense_matrix<T> M(m,n); std::copy(v.begin(), v.end(), M.begin());
    if (do_add) gmm::add (M, gmm::sub_matrix(gsp.wsc(T()), ii, jj));
    else        gmm::copy(M, gmm::sub_matrix(gsp.wsc(T()), ii, jj));
  }
}

void
spmat_set_or_add_sub_matrix(gsparse &gsp, getfemint::mexargs_in& in, bool do_add) {
  sub_index ii = in.pop().to_sub_index().check_range(gsp.nrows());
  sub_index jj = in.remaining() ? 
    in.pop().to_sub_index().check_range(gsp.ncols()) : ii.check_range(gsp.ncols());
  if (!gsp.is_complex()) {
    if (in.front().is_complex()) gsp.to_complex();
    else spmat_set_or_add_sub_matrix(gsp, in, ii, jj, do_add, scalar_type());
  } 
  if (gsp.is_complex()) {
    spmat_set_or_add_sub_matrix(gsp, in, ii, jj, do_add, complex_type());
  }
}

template <typename T, typename SUBI> static void
spmat_do_clear(gsparse &gsp, SUBI &ii, SUBI &jj, T) {
  if (gsp.storage() == gsparse::CSCMAT) THROW_BADARG("cannot not clear a CSC matrix, convert to WSC first");
  gmm::clear(gmm::sub_matrix(gsp.wsc(T()), ii, jj));
}

/*MLABCOM
  FUNCTION gf_spmat_set(M, args)

  Modification of the content of a getfem sparse matrix.

  @SET SPMAT:SET('clear')
  @SET SPMAT:SET('scale')
  @SET SPMAT:SET('transpose')
  @SET SPMAT:SET('conjugate')
  @SET SPMAT:SET('transconj')
  @SET SPMAT:SET('to_csc')
  @SET SPMAT:SET('to_wsc')
  @SET SPMAT:SET('to_complex')
  @SET SPMAT:SET('diag')
  @SET SPMAT:SET('assign')
  @SET SPMAT:SET('add')
 MLABCOM*/

void gf_spmat_set(getfemint::mexargs_in& in, getfemint::mexargs_out& out)
{
  if (in.narg() < 1) {
    THROW_BADARG( "Wrong number of input arguments");
  }
  getfemint_gsparse *pgsp = in.pop().to_getfemint_gsparse();
  gsparse &gsp = pgsp->sparse();

  std::string cmd = in.pop().to_string();

  if (check_cmd(cmd, "clear", in, out, 0, 2, 0, 0)) {
    /*@SET SPMAT:SET('clear' [, I[, J]])
      Erase the non-zero entries of the matrix.

      The optional arguments I and J may be specified to clear a
      sub-matrix instead of the entire matrix. @*/
    if (in.remaining()) {
      sub_index ii = in.pop().to_sub_index().check_range(gsp.nrows());
      sub_index jj = in.remaining() ? 
	in.pop().to_sub_index().check_range(gsp.ncols()) : ii.check_range(gsp.ncols());
      if (!gsp.is_complex()) spmat_do_clear(gsp, ii, jj, scalar_type());
      else                   spmat_do_clear(gsp, ii, jj, complex_type());
    } else {
      gmm::sub_interval ii(0,gsp.nrows());
      gmm::sub_interval jj(0,gsp.ncols());
      if (!gsp.is_complex()) spmat_do_clear(gsp, ii, jj, scalar_type());
      else                   spmat_do_clear(gsp, ii, jj, complex_type());      
    }
  } else if (check_cmd(cmd, "scale", in, out, 1, 1, 0, 0)) {
    /*@SET SPMAT:SET('scale', V)
      Multiplies the matrix by a scalar value V.
      @*/
    gsp.to_wsc();
    if (!gsp.is_complex() && in.front().is_complex()) gsp.to_complex();
    if (!gsp.is_complex()) {
      gmm::scale(gsp.real_wsc(), in.pop().to_scalar());
    } else {
      gmm::scale(gsp.cplx_wsc(), in.pop().to_scalar(complex_type()));
    }
  } else if (check_cmd(cmd, "transpose", in, out, 0, 0, 0, 0)) {    
    /*@SET SPMAT:SET('transpose')
      Transpose the matrix.
      @*/
    if (!gsp.is_complex())
      gf_spmat_set_transpose(gsp, TRANSP, scalar_type());
    else gf_spmat_set_transpose(gsp, TRANSP, complex_type());
  } else if (check_cmd(cmd, "conjugate", in, out, 0, 0, 0, 0)) {    
    /*@SET SPMAT:SET('conjugate')
      Conjugate each element of the matrix.
      @*/
    if (!gsp.is_complex())
      gf_spmat_set_transpose(gsp, CONJ, scalar_type());
    else gf_spmat_set_transpose(gsp, CONJ, complex_type());
  } else if (check_cmd(cmd, "transconj", in, out, 0, 0, 0, 0)) {    
    /*@SET SPMAT:SET('transconj')
      Transpose and conjugate the matrix.
      @*/
    if (!gsp.is_complex())
      gf_spmat_set_transpose(gsp, TRANSCONJ, scalar_type());
    else gf_spmat_set_transpose(gsp, TRANSCONJ, complex_type());
  } else if (check_cmd(cmd, "to_csc", in, out, 0, 0, 0, 0)) {
    /*@SET SPMAT:SET('to_csc')
      Convert the matrix to CSC storage.

      CSC storage is recommended for matrix-vector multiplications.
      @*/
    gsp.to_csc();
  } else if (check_cmd(cmd, "to_wsc", in, out, 0, 0, 0, 0)) {
    /*@SET SPMAT:SET('to_wsc')
      Convert the matrix to WSC storage.

      Read and write operation are quite fast with WSC storage.
      @*/
    gsp.to_wsc();
  } else if (check_cmd(cmd, "to_complex", in, out, 0, 0, 0, 0)) {
    /*@SET SPMAT:SET('to_complex')
      Store complex numbers.
      @*/
    gsp.to_complex();
  } else if (check_cmd(cmd, "diag", in, out, 1, 2, 0, 0)) {
    /*@SET SPMAT:SET('diag', @dmat D [, @ivec E])
      Change the diagonal (or sub-diagonals) of the matrix.

      If E is given, D might be a matrix and each column of E will
      contain the sub-diagonal number that will be filled with the
      corresponding column of D.
      @*/
    spmat_set_diag(gsp, in, false);
  } else if (check_cmd(cmd, "assign", in, out, 3, 3, 0, 0)) { 
    /*@SET SPMAT:SET('assign', @ivec I, @ivec J, V)
      Copy V into the sub-matrix M(I,J).

      V might be a sparse matrix or a full matrix.
      @*/
    spmat_set_or_add_sub_matrix(gsp, in, false);
  } else if (check_cmd(cmd, "add", in, out, 3, 3, 0, 0)) {
    /*@SET SPMAT:SET('add', I, J, V)
      Add V to the sub-matrix M(I,J).

      V might be a sparse matrix or a full matrix.
      @*/
    spmat_set_or_add_sub_matrix(gsp, in, true);
  } else bad_cmd(cmd);
}
