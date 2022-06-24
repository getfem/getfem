/*===========================================================================

 Copyright (C) 2006-2020 Julien Pommier.

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
#include <getfemint_workspace.h>
#include <getfem/getfem_assembling.h>
#include <gmm/gmm_inoutput.h>
#include <gmm/gmm_MUMPS_interface.h>

using namespace getfemint;

template <typename MAT>
void copydiags(const MAT &M, const std::vector<size_type> &v,
	       garray<typename MAT::value_type> &w) {
  size_type m = gmm::mat_nrows(M), n = gmm::mat_ncols(M);
  for (unsigned ii=0; ii < v.size(); ++ii) {
    int d = int(v[ii]), i,j;
    if (d < 0) { i = -d; j = 0; } else { i = 0; j = d; }
    cout << "m=" << m << "n=" << n << ", d=" << d << ", i=" << i
	 << ", j=" << j << "\n";
    for (; i < int(m) && j < int(n); ++i,++j)
      w(i,ii) = M(i,j);
  }
}

template <typename T> static void
gf_spmat_get_full(gsparse &gsp, getfemint::mexargs_in& in, getfemint::mexargs_out& out, T) {
  size_type n,m;
  if (!in.remaining()) {
    m = gsp.nrows(); n = gsp.ncols();
    gmm::dense_matrix<T> ww(m,n);
    switch (gsp.storage()) {
    case gsparse::CSCMAT: gmm::copy(gsp.csc(T()), ww); break;
    case gsparse::WSCMAT: gmm::copy(gsp.wsc(T()), ww); break;
    default: THROW_INTERNAL_ERROR;
    }
    std::copy(ww.begin(), ww.end(),
	      out.pop().create_array(int(m),int(n),T()).begin());
  } else {
    sub_index ii = in.pop().to_sub_index().check_range(gsp.nrows());
    sub_index jj = in.remaining() ?
      in.pop().to_sub_index().check_range(gsp.ncols())
      : ii.check_range(gsp.ncols());
    m = ii.size(); n = jj.size();
    if (m == 1 && n == 1) {
      T val = T(0);
      switch (gsp.storage()) {
      case gsparse::CSCMAT: val = gsp.csc(T())(ii.first(), jj.first()); break;
      case gsparse::WSCMAT: val = gsp.wsc(T()).col(jj.first()).r(ii.first()); break;
      default: THROW_INTERNAL_ERROR;
      }
      if (gsp.is_complex()) {
	if (out.remaining()) out.pop().from_scalar(gmm::real(val));
	if (out.remaining()) out.pop().from_scalar(gmm::imag(val));
      } else {
	if (out.remaining()) out.pop().from_scalar(gmm::real(val));
	// if (out.remaining()) out.pop().from_scalar(0);
       }
    } else {
      gmm::dense_matrix<T> ww(m,n);
      switch (gsp.storage()) {
      case gsparse::CSCMAT:
	gmm::copy(gmm::sub_matrix(gsp.csc(T()),ii,jj), ww); break;
      case gsparse::WSCMAT:
	gmm::copy(gmm::sub_matrix(gsp.wsc(T()),ii,jj), ww); break;
      default: THROW_INTERNAL_ERROR;
      }
      std::copy(ww.begin(), ww.end(),
		out.pop().create_array(int(m),int(n),T()).begin());
    }
  }
  
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
    for (unsigned i=0; i < unsigned(nz); ++i) { w[i] = M.pr[i]; }
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

/*@GFCOM
  General getfem sparse matrix inquiry function. M might also be a
  native matlab sparse matrix.
@*/




// Object for the declaration of a new sub-command.

struct sub_gf_spmat_get : virtual public dal::static_stored_object {
  int arg_in_min, arg_in_max, arg_out_min, arg_out_max;
  virtual void run(getfemint::mexargs_in& in,
		   getfemint::mexargs_out& out, gsparse &gsp) = 0;
};

typedef std::shared_ptr<sub_gf_spmat_get> psub_command;

// Function to avoid warning in macro with unused arguments.
template <typename T> static inline void dummy_func(T &) {}

#define sub_command(name, arginmin, arginmax, argoutmin, argoutmax, code) { \
    struct subc : public sub_gf_spmat_get {				\
      virtual void run(getfemint::mexargs_in& in,			\
		       getfemint::mexargs_out& out, gsparse &gsp)	\
      { dummy_func(in); dummy_func(out); dummy_func(gsp); code }	\
    };									\
    psub_command psubc = std::make_shared<subc>();			\
    psubc->arg_in_min = arginmin; psubc->arg_in_max = arginmax;		\
    psubc->arg_out_min = argoutmin; psubc->arg_out_max = argoutmax;	\
    subc_tab[cmd_normalize(name)] = psubc;				\
  }                           



void gf_spmat_get(getfemint::mexargs_in& m_in,
		  getfemint::mexargs_out& m_out) {
  typedef std::map<std::string, psub_command > SUBC_TAB;
  static SUBC_TAB subc_tab;
  
  if (subc_tab.size() == 0) {
  

    /*@GET n = ('nnz')
      Return the number of non-null values stored in the sparse matrix.@*/
    sub_command
      ("nnz", 0, 0, 0, 1,
       out.pop().from_integer(int(gsp.nnz()));
       );


    /*@GET Sm = ('full'[, @list I[, @list J]])
      Return a full (sub-)matrix.
      
      The optional arguments `I` and `J`, are the sub-intervals for the
      rows and columns that are to be extracted.@*/
    sub_command
      ("full", 0, 2, 0, 1,
       if (gsp.is_complex()) gf_spmat_get_full(gsp, in, out, complex_type());
       else gf_spmat_get_full(gsp, in, out, scalar_type());
       );


    /*@GET MV = ('mult', @vec V)
      Product of the sparse matrix `M` with a vector `V`.
      
      For matrix-matrix multiplications, see SPMAT:INIT('mult').@*/
    sub_command
      ("mult", 1, 1, 0, 1,
       if (!gsp.is_complex())
	 gf_spmat_mult_or_tmult(gsp, in, out, false, scalar_type());
       else gf_spmat_mult_or_tmult(gsp, in, out, false, complex_type());
       );


    /*@GET MtV =('tmult', @vec V)
      Product of `M` transposed (conjugated if `M` is complex) with the
      vector `V`.@*/
    sub_command
      ("tmult", 1, 1, 0, 1,
       if (!gsp.is_complex())
	 gf_spmat_mult_or_tmult(gsp, in, out, true, scalar_type());
       else gf_spmat_mult_or_tmult(gsp, in, out, true, complex_type());
       );


    /*@GET D = ('diag'[, @list E])
      Return the diagonal of `M` as a vector.
      
      If `E` is used, return the sub-diagonals whose ranks are given in E.@*/
    sub_command
      ("diag", 0, 1, 0, 1,
       if (!gsp.is_complex())
	 gf_spmat_get_diag(gsp, in, out, scalar_type());
       else gf_spmat_get_diag(gsp, in, out, complex_type());
       );


    /*@GET s = ('storage')
      Return the storage type currently used for the matrix.
      
      The storage is returned as a string, either 'CSC' or 'WSC'.@*/
    sub_command
      ("storage", 0, 0, 0, 1,
       out.pop().from_string(gsp.name());
       );


    /*@GET @CELL{ni,nj} = ('size')
      Return a vector where `ni` and `nj` are the dimensions of the matrix.@*/
    sub_command
      ("size", 0, 0, 0, 1,
       iarray sz = out.pop().create_iarray_h(2);
       sz[0] = int(gsp.nrows());
       sz[1] = int(gsp.ncols());
       );


    /*@GET b = ('is_complex')
      Return 1 if the matrix contains complex values.@*/
    sub_command
      ("is_complex", 0, 0, 0, 1,
       out.pop().from_integer(gsp.is_complex());
       );


    /*@GET @CELL{JC, IR} = ('csc_ind')
      Return the two usual index arrays of CSC storage.
      
      If `M` is not stored as a CSC matrix, it is converted into CSC.@*/
    sub_command
      ("csc_ind", 0, 0, 0, 2,
       gsp.to_csc();
       if (!gsp.is_complex())
	 gf_spmat_get_data(gsp.csc(scalar_type()),  out, 0);
       else
	 gf_spmat_get_data(gsp.csc(complex_type()), out, 0);
       );

    
    /*@GET V = ('csc_val')
      Return the array of values of all non-zero entries of `M`.
      
      If `M` is not stored as a CSC matrix, it is converted into CSC.@*/
    sub_command
      ("csc_val", 0, 0, 0, 1,
       gsp.to_csc();
       if (!gsp.is_complex())
	 gf_spmat_get_data(gsp.csc(scalar_type()),  out, 1);
       else
	 gf_spmat_get_data(gsp.csc(complex_type()), out, 1);
       );


    /*@GET @CELL{N, U0} = ('dirichlet nullspace', @vec R)
    Solve the dirichlet conditions `M.U=R`.

    A solution `U0` which has a minimum L2-norm is returned, with a
    sparse matrix `N` containing an orthogonal basis of the kernel of
    the (assembled) constraints matrix `M` (hence, the PDE linear system
    should be solved on this subspace): the initial problem

    `K.U = B` with constraints `M.U = R`

    is replaced by

    `(N'.K.N).UU = N'.B` with `U = N.UU + U0`@*/
    sub_command
      ("dirichlet nullspace", 1, 1, 2, 2,
       if (gsp.is_complex())
	 gf_spmat_get_Dirichlet_nullspace(gsp, in, out, complex_type());
       else
	 gf_spmat_get_Dirichlet_nullspace(gsp, in, out, scalar_type());
       );


    /*@GET ('save', @str format, @str filename)
      Export the sparse matrix.

      the format of the file may be 'hb' for Harwell-Boeing, or 'mm'
      for Matrix-Market.@*/
    sub_command
      ("save", 2, 2, 0, 0,
       std::string fmt = in.pop().to_string();
       int ifmt;
       if (cmd_strmatch(fmt, "hb") || cmd_strmatch(fmt, "harwell-boeing")) ifmt = 0;
       else if (cmd_strmatch(fmt, "mm") || cmd_strmatch(fmt, "matrix-market")) ifmt = 1;
       else THROW_BADARG("unknown sparse matrix file-format : " << fmt);
       std::string fname = in.pop().to_string();
       gsp.to_csc();
       if (!gsp.is_complex()) {
	 if (ifmt == 0)
	   gmm::Harwell_Boeing_save(fname.c_str(), gsp.csc(scalar_type()));
	 else
           gmm::MatrixMarket_save(fname.c_str(), gsp.csc(scalar_type()));
       } else {
	 if (ifmt == 0)
	   gmm::Harwell_Boeing_save(fname.c_str(), gsp.csc(complex_type()));
	 else
           gmm::MatrixMarket_save(fname.c_str(), gsp.csc(complex_type()));
       }
       );


    /*@GET s = ('char')
      Output a (unique) string representation of the @tspmat.

      This can be used to perform comparisons between two
      different @tspmat objects.
      This function is to be completed.
      @*/
    sub_command
      ("char", 0, 0, 0, 1,
       std::stringstream s;
       if (gsp.storage() == getfemint::gsparse::WSCMAT) {
	 if (!gsp.is_complex()) s << gsp.wsc(scalar_type());
	 else                   s << gsp.wsc(complex_type());
       } else {
	 if (!gsp.is_complex()) s << gsp.csc(scalar_type());
	 else           	s << gsp.csc(complex_type());
       }
       out.pop().from_string(s.str().c_str());
       );


    /*@GET ('display')
      displays a short summary for a @tspmat object.@*/
    sub_command
      ("display", 0, 0, 0, 0,
       size_type ncc = gsp.nrows()*gsp.ncols();
       infomsg() << gsp.nrows() << "x" << gsp.ncols() << " "
       << (gsp.is_complex() ? "COMPLEX" : "REAL") << " " << gsp.name()
       << ", NNZ=" << gsp.nnz() << " (filling="
       << 100.*double(gsp.nnz())/(double(ncc == 0 ? 1 : ncc)) << "%)";
       );


#if defined(GMM_USES_MUMPS) || defined(HAVE_DMUMPS_C_H)
    /*@GET @CELL{mantissa_r, mantissa_i, exponent} = ('determinant')
      returns the matrix determinant calculated using MUMPS.@*/
    sub_command
      ("determinant", 0, 0, 0, 3,
       gsp.to_csc();
       int exponent;
       if (gsp.is_complex()) {
         complex_type det = gmm::MUMPS_determinant(gsp.csc(complex_type()),
                                                   exponent);
         if (out.remaining()) out.pop().from_scalar(gmm::real(det));
         if (out.remaining()) out.pop().from_scalar(gmm::imag(det));
       } else {
         scalar_type det = gmm::MUMPS_determinant(gsp.csc(scalar_type()),
                                                  exponent);
         if (out.remaining()) out.pop().from_scalar(det);
         if (out.remaining()) out.pop().from_scalar(0);
       }
       if (out.remaining()) out.pop().from_integer(exponent);
       );
#endif
  }

  if (m_in.narg() < 2)  THROW_BADARG( "Wrong number of input arguments");

  std::shared_ptr<gsparse> pgsp = m_in.pop().to_sparse();
  gsparse &gsp = *pgsp;
  std::string init_cmd   = m_in.pop().to_string();
  std::string cmd        = cmd_normalize(init_cmd);

  
  SUBC_TAB::iterator it = subc_tab.find(cmd);
  if (it != subc_tab.end()) {
    check_cmd(cmd, it->first.c_str(), m_in, m_out, it->second->arg_in_min,
	      it->second->arg_in_max, it->second->arg_out_min,
	      it->second->arg_out_max);
    it->second->run(m_in, m_out, gsp);
  }
  else bad_cmd(init_cmd);

}
