/*===========================================================================

 Copyright (C) 2006-2020 Yves Renard, Julien Pommier.

 This file is a part of GetFEM

 GetFEM  is  free software;  you  can  redistribute  it  and/or modify it
 under  the  terms  of the  GNU  Lesser General Public License as published
 by  the  Free Software Foundation;  either version 3 of the License,  or
 (at your option) any later version along with the GCC Runtime Library
 Exception either version 3.1 or (at your option) any later version.
 This program  is  distributed  in  the  hope  that it will be useful,  but
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 License and GCC Runtime Library Exception for more details.
 You  should  have received a copy of the GNU Lesser General Public License
 along  with  this program;  if not, write to the Free Software Foundation,
 Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.

===========================================================================*/

#include <getfemint_gsparse.h>
#include <gmm/gmm_inoutput.h>

namespace getfemint {

  gsparse::gsparse(size_type m, size_type n, storage_type s_, value_type v_) 
    : pwscmat_r(0), pwscmat_c(0), pcscmat_r(0), pcscmat_c(0), gfimat(0) {
    allocate(m,n,s_,v_);
  }
  gsparse::gsparse(const gfi_array *arg) 
    : pwscmat_r(0), pwscmat_c(0), pcscmat_r(0), pcscmat_c(0), gfimat(arg) {
    if (gfi_array_get_class(arg) != GFI_SPARSE) THROW_INTERNAL_ERROR;
    v = gfi_array_is_complex(arg) ? COMPLEX : REAL; s = CSCMAT;
  }

  void gsparse::swap(gsparse &other) {
    std::swap(v, other.v); 
    std::swap(s, other.s); 
    std::swap(pcscmat_r, other.pcscmat_r);
    std::swap(pcscmat_c, other.pcscmat_c);
    std::swap(pwscmat_r, other.pwscmat_r);
    std::swap(pwscmat_c, other.pwscmat_c);
    std::swap(gfimat, other.gfimat);
  }

  gsparse& gsparse::destructive_assign(t_wscmat_r &M) { 
    destroy(); allocate(0,0, WSCMAT, REAL); (*pwscmat_r).swap(M); 
    return *this;
  }
  gsparse& gsparse::destructive_assign(t_wscmat_c &M) {
    destroy(); allocate(0,0, WSCMAT, COMPLEX); (*pwscmat_c).swap(M); 
    return *this;
  }

  gsparse& gsparse::destructive_assign(t_cscmat_r &M) {
    destroy(); allocate(0,0, CSCMAT, REAL); (*pcscmat_r).swap(M); 
    return *this;
  }

  gsparse& gsparse::destructive_assign(t_cscmat_c &M) {
    destroy(); allocate(0,0, CSCMAT, COMPLEX); (*pcscmat_c).swap(M); 
    return *this;
  }

  void gsparse::destroy() {
    if (pwscmat_r) delete pwscmat_r;
    pwscmat_r = 0;
    if (pwscmat_c) delete pwscmat_c;
    pwscmat_c = 0;
    if (pcscmat_r) delete pcscmat_r;
    pcscmat_r = 0;
    if (pcscmat_c) delete pcscmat_c;
    pcscmat_c = 0;
  }

  void gsparse::allocate(size_type m, size_type n, storage_type s_, value_type v_) {
    v = v_; s = s_;
    if (v == REAL) {
      switch (s) {
      case WSCMAT: real_wsc(new t_wscmat_r(m,n)); break;
      case CSCMAT: real_csc_w(new t_cscmat_r(m,n)); break;
      default: THROW_INTERNAL_ERROR;
      }
    } else {
      switch (s) {
      case WSCMAT: cplx_wsc(new t_wscmat_c(m,n)); break;
      case CSCMAT: cplx_csc_w(new t_cscmat_c(m,n)); break; 
      default: THROW_INTERNAL_ERROR;
      }
    }
  }

  void gsparse::deallocate(storage_type s_, value_type v_) {
    if (v_ == REAL) {
      switch (s_) {
        case WSCMAT: delete pwscmat_r; pwscmat_r = 0; break;
        case CSCMAT: delete pcscmat_r; pcscmat_r = 0; break;
        default: THROW_INTERNAL_ERROR;
      }
    } else {
      switch (s_) {
      case WSCMAT: delete pwscmat_c; pwscmat_c = 0; break;
      case CSCMAT: delete pcscmat_c; pcscmat_c = 0; break;
      default: THROW_INTERNAL_ERROR;
      }
    }
  }

  size_type gsparse::nrows() const {
    if (pwscmat_r) return gmm::mat_nrows(*pwscmat_r);
    else if (pwscmat_c) return gmm::mat_nrows(*pwscmat_c);
    else if (pcscmat_r) return gmm::mat_nrows(*pcscmat_r);
    else if (pcscmat_c) return gmm::mat_nrows(*pcscmat_c);
    else if (gfimat) return gfi_array_get_dim(gfimat)[0];
    else return 0;
  }

  size_type gsparse::ncols() const {
    if (pwscmat_r) return gmm::mat_ncols(*pwscmat_r);
    else if (pwscmat_c) return gmm::mat_ncols(*pwscmat_c);
    else if (pcscmat_r) return gmm::mat_ncols(*pcscmat_r);
    else if (pcscmat_c) return gmm::mat_ncols(*pcscmat_c);
    else if (gfimat) return gfi_array_get_dim(gfimat)[1];
    else return 0;
  }

  size_type gsparse::nnz() const {
    switch (s) {
    case WSCMAT: 
      if (pwscmat_r) return gmm::nnz(*pwscmat_r); 
      if (pwscmat_c) return gmm::nnz(*pwscmat_c);
      break;
    case CSCMAT: 
      if (pcscmat_r) return gmm::nnz(*pcscmat_r);
      if (pcscmat_c) return gmm::nnz(*pcscmat_c);
      break;
    default: THROW_INTERNAL_ERROR;
    }
    return 0;
  }

  void gsparse::to_wsc() {
    if (is_a_native_matrix_ref()) THROW_INTERNAL_ERROR;
    switch (s) {
    case WSCMAT: break;
    case CSCMAT: {
      allocate(nrows(),ncols(),WSCMAT, v);
      if (v == REAL) gmm::copy(real_csc(), real_wsc());
      else           gmm::copy(cplx_csc(), cplx_wsc());
      deallocate(CSCMAT, v);
    } break;
    default: THROW_INTERNAL_ERROR;
    }
  }

  void gsparse::to_csc() {
    switch (s) {
    case CSCMAT: break;
    case WSCMAT: {
      allocate(nrows(),ncols(),CSCMAT, v);
      if (v == REAL) gmm::copy(real_wsc(), real_csc_w());
      else           gmm::copy(cplx_wsc(), cplx_csc_w());
      deallocate(WSCMAT, v);
    } break;
    default: THROW_INTERNAL_ERROR;
    }
  }
  
  void gsparse::to_complex() {
    if (is_complex()) return;
    allocate(nrows(), ncols(), s, COMPLEX);
    switch (s) {
      case CSCMAT: gmm::copy(real_csc(), cplx_csc_w()); break;
      case WSCMAT: gmm::copy(real_wsc(), cplx_wsc()); break;
    }
    deallocate(s, REAL);
  }


  /* common templates shared between gf_spmat_* */
  template <typename MAT> 
  void setdiags(MAT &M, const std::vector<int> &v, 
		const garray<typename MAT::value_type> &w) {
    size_type m = gmm::mat_nrows(M), n = gmm::mat_ncols(M);
    for (unsigned ii=0; ii < std::min<size_type>(v.size(), w.getn()); ++ii) {
      int d = v[ii], i,j;
      if (d < 0) { i = -d; j = 0; } else { i = 0; j = d; }
      for (; i < int(m) && j < int(n) && i < int(w.getm()); ++i,++j)
	M(i,j) = w(i, ii);
    }
  }
  
  template <typename T> static void
  gf_spmat_set_diag(gsparse &gsp,
		    mexargs_in& in, bool create_matrix, T) {
    garray<T> w = in.pop().to_garray(-1, -1, T());
    if (!create_matrix && w.getm() < std::min(gsp.nrows(), gsp.ncols())) {
      THROW_BADARG("not enough rows for the diagonals (expected at least " << 
		   std::min(gsp.nrows(), gsp.ncols()) << ")");
    }
    std::vector<int> v;
    if (in.remaining()) {
      iarray vv = in.pop().to_iarray(-1);
      for (size_type i=0; i < vv.size(); ++i)
	v.push_back(vv[i]);
    } else v.push_back(0);
    if (create_matrix) {
      size_type m = w.getm(), n = w.getm();
      if (in.remaining()) { m = n = in.pop().to_integer(1,INT_MAX); }
      if (in.remaining()) { n = in.pop().to_integer(1,INT_MAX); }      
      gsp.wsc(T(), new gmm::col_matrix<gmm::wsvector<T> >(m,n));
    }
    if (w.getn() != v.size()) 
      THROW_BADARG("cannot set diagonals: inconsistent number of diags between the data (" << 
		   w.getn() << " columns and the diag numbers (" << v.size() << " elements)");
    gsp.to_wsc();
    setdiags(gsp.wsc(T()), v, w);
  }

  void
  spmat_set_diag(gsparse &gsp, mexargs_in& in, bool create_matrix) {
    if (in.front().is_complex() || (!create_matrix && gsp.is_complex()))
      gf_spmat_set_diag(gsp, in, create_matrix, complex_type());
    else gf_spmat_set_diag(gsp, in, create_matrix, scalar_type());
  }

  void spmat_load(mexargs_in& in, mexargs_out& out, 
		  mexarg_out::output_sparse_fmt fmt) {
    std::string mt = in.pop().to_string();
    std::string fname = in.pop().to_string();
    if (cmd_strmatch(mt, "hb") || 
	cmd_strmatch(mt, "harwell-boeing")) {
      gmm::HarwellBoeing_IO h; h.open(fname.c_str());
      gsparse gsp;
      if (h.is_complex()) {
	gmm::csc_matrix<complex_type> cscH;
	h.read(cscH);
	gsp.destructive_assign(cscH);
      } else {
	gmm::csc_matrix<double> cscH;
	h.read(cscH);
	gsp.destructive_assign(cscH);
      }
      out.pop().from_sparse(gsp, fmt);
      gsp.deallocate(gsp.storage(), gsp.is_complex() ? 
		     gsparse::COMPLEX : gsparse::REAL);
    } else if (cmd_strmatch(mt, "mm") || 
	       cmd_strmatch(mt, "matrix-market")) {
      gmm::MatrixMarket_IO h; h.open(fname.c_str());
      if (h.is_complex()) {
	gf_cplx_sparse_by_col H; h.read(H); 
	out.pop().from_sparse(H, fmt);
      } else {
	gf_real_sparse_by_col H; h.read(H);
	out.pop().from_sparse(H, fmt);
      }
    } else THROW_BADARG("unknown sparse matrix file-format : " << mt);
  }

}
