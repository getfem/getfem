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
    std::shared_ptr<gsparse> src = in.pop().to_sparse();
    switch (src->storage()) {
      case gsparse::WSCMAT:
	if (do_add) gmm::add (src->wsc(T()), gmm::sub_matrix(gsp.wsc(T()), ii, jj));
	else        gmm::copy(src->wsc(T()), gmm::sub_matrix(gsp.wsc(T()), ii, jj));
	break;
      case gsparse::CSCMAT:
	if (do_add) gmm::add (src->csc(T()), gmm::sub_matrix(gsp.wsc(T()), ii, jj));
	else        gmm::copy(src->csc(T()), gmm::sub_matrix(gsp.wsc(T()), ii, jj));
	break;
      default: THROW_INTERNAL_ERROR;
    }
  } else {
    garray<T> v = in.pop().to_garray(int(m),int(n),T());
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

/*@GFDOC
   Modification of the content of a getfem sparse matrix.
 @*/



// Object for the declaration of a new sub-command.

struct sub_gf_spmat_set : virtual public dal::static_stored_object {
  int arg_in_min, arg_in_max, arg_out_min, arg_out_max;
  virtual void run(getfemint::mexargs_in& in,
		   getfemint::mexargs_out& out,
		   gsparse &gsp) = 0;
};

typedef std::shared_ptr<sub_gf_spmat_set> psub_command;

// Function to avoid warning in macro with unused arguments.
template <typename T> static inline void dummy_func(T &) {}

#define sub_command(name, arginmin, arginmax, argoutmin, argoutmax, code) { \
    struct subc : public sub_gf_spmat_set {				\
      virtual void run(getfemint::mexargs_in& in,			\
		       getfemint::mexargs_out& out,			\
		       gsparse &gsp)					\
      { dummy_func(in); dummy_func(out); code }				\
    };									\
    psub_command psubc = std::make_shared<subc>();			\
    psubc->arg_in_min = arginmin; psubc->arg_in_max = arginmax;		\
    psubc->arg_out_min = argoutmin; psubc->arg_out_max = argoutmax;	\
    subc_tab[cmd_normalize(name)] = psubc;				\
  }                           



void gf_spmat_set(getfemint::mexargs_in& m_in, getfemint::mexargs_out& m_out) {
  typedef std::map<std::string, psub_command > SUBC_TAB;
  static SUBC_TAB subc_tab;

  if (subc_tab.size() == 0) {


    /*@SET ('clear'[, @list I[, @list J]])
      Erase the non-zero entries of the matrix.

      The optional arguments `I` and `J` may be specified to clear a
      sub-matrix instead of the entire matrix.@*/
    sub_command
      ("clear", 0, 2, 0, 0,
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
       );


    /*@SET ('scale', @scalar v)
      Multiplies the matrix by a scalar value `v`.@*/
    sub_command
      ("scale", 1, 1, 0, 0,
       gsp.to_wsc();
       if (!gsp.is_complex() && in.front().is_complex()) gsp.to_complex();
       if (!gsp.is_complex()) {
	 gmm::scale(gsp.real_wsc(), in.pop().to_scalar());
       } else {
	 gmm::scale(gsp.cplx_wsc(), in.pop().to_scalar(complex_type()));
       }
       );


    /*@SET ('transpose')
      Transpose the matrix.@*/
    sub_command
      ("transpose", 0, 0, 0, 0,
       if (!gsp.is_complex())
	 gf_spmat_set_transpose(gsp, TRANSP, scalar_type());
       else gf_spmat_set_transpose(gsp, TRANSP, complex_type());
       );


    /*@SET ('conjugate')
      Conjugate each element of the matrix.@*/
    sub_command
      ("conjugate", 0, 0, 0, 0,
       if (!gsp.is_complex())
	 gf_spmat_set_transpose(gsp, CONJ, scalar_type());
       else gf_spmat_set_transpose(gsp, CONJ, complex_type());
       );


    /*@SET ('transconj')
      Transpose and conjugate the matrix.@*/
    sub_command
      ("transconj", 0, 0, 0, 0,
       if (!gsp.is_complex())
	 gf_spmat_set_transpose(gsp, TRANSCONJ, scalar_type());
       else gf_spmat_set_transpose(gsp, TRANSCONJ, complex_type());
       );


    /*@SET ('to_csc')
      Convert the matrix to CSC storage.

      CSC storage is recommended for matrix-vector multiplications.@*/
    sub_command
      ("to_csc", 0, 0, 0, 0,
       gsp.to_csc();
       );


    /*@SET ('to_wsc')
      Convert the matrix to WSC storage.
      
      Read and write operation are quite fast with WSC storage.@*/
    sub_command
      ("to_wsc", 0, 0, 0, 0,
       gsp.to_wsc();
       );


    /*@SET ('to_complex')
      Store complex numbers.@*/
    sub_command
      ("to_complex", 0, 0, 0, 0,
       gsp.to_complex();
       );


    /*@SET ('diag', @dmat D [, @ivec E])
      Change the diagonal (or sub-diagonals) of the matrix.
      
      If `E` is given, `D` might be a matrix and each column of `E` will
      contain the sub-diagonal number that will be filled with the
      corresponding column of `D`.@*/
    sub_command
      ("diag", 1, 2, 0, 0,
       spmat_set_diag(gsp, in, false);
       );


    /*@SET ('assign', @ivec I, @ivec J, @mat V)
      Copy V into the sub-matrix 'M(I,J)'.
      
      `V` might be a sparse matrix or a full matrix.@*/
    sub_command
      ("assign", 3, 3, 0, 0,
       spmat_set_or_add_sub_matrix(gsp, in, false);
       );


    /*@SET ('add', @ivec I, @ivec J, @mat V)
    Add `V` to the sub-matrix 'M(I,J)'.

    `V` might be a sparse matrix or a full matrix.@*/
    sub_command
      ("add", 3, 3, 0, 0,
       spmat_set_or_add_sub_matrix(gsp, in, true);
       );

  }

  if (m_in.narg() < 2)  THROW_BADARG( "Wrong number of input arguments");


  gsparse &gsp = *(to_spmat_object(m_in.pop()));
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
