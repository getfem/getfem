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

#include <getfemint_gsparse_misc.h>
#include <gmm/gmm_inoutput.h>

namespace getfemint {
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

