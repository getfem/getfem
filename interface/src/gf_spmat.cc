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

/*@GFDOC
  Create a new sparse matrix in GetFEM format@MATLAB{(, i.e. sparse
  matrices which are stored in the getfem workspace, not the matlab sparse
  matrices)}. These sparse matrix can be stored as CSC (compressed column
  sparse), which is the format used by Matlab, or they can be stored as WSC
  (internal format to getfem). The CSC matrices are not writable (it would
  be very inefficient), but they are optimized for multiplication with
  vectors, and memory usage. The WSC are writable, they are very fast with
  respect to random read/write operation. However their memory overhead is
  higher than CSC matrices, and they are a little bit slower for
  matrix-vector multiplications.

  By default, all newly created matrices are build as WSC matrices. This can
  be changed later with ``SPMAT:SET('to_csc',...)``, or may be changed
  automatically by getfem (for example ``::LINSOLVE()`` converts the
  matrices to CSC).

  The matrices may store REAL or COMPLEX values.
@*/



// Object for the declaration of a new sub-command.

struct sub_gf_spmat : virtual public dal::static_stored_object {
  int arg_in_min, arg_in_max, arg_out_min, arg_out_max;
  virtual void run(getfemint::mexargs_in& in,
		   getfemint::mexargs_out& out,
		   std::shared_ptr<gsparse> &gsp) = 0;
};

typedef std::shared_ptr<sub_gf_spmat> psub_command;

// Function to avoid warning in macro with unused arguments.
template <typename T> static inline void dummy_func(T &) {}

#define sub_command(name, arginmin, arginmax, argoutmin, argoutmax, code) { \
    struct subc : public sub_gf_spmat {					\
      virtual void run(getfemint::mexargs_in& in,			\
		       getfemint::mexargs_out& out,			\
		       std::shared_ptr<gsparse> &gsp)			\
      { dummy_func(in); dummy_func(out); code }				\
    };									\
    psub_command psubc = std::make_shared<subc>();			\
    psubc->arg_in_min = arginmin; psubc->arg_in_max = arginmax;		\
    psubc->arg_out_min = argoutmin; psubc->arg_out_max = argoutmax;	\
    subc_tab[cmd_normalize(name)] = psubc;				\
  }





void gf_spmat(getfemint::mexargs_in& m_in,
	      getfemint::mexargs_out& m_out) {
  typedef std::map<std::string, psub_command > SUBC_TAB;
  static SUBC_TAB subc_tab;

  if (subc_tab.size() == 0) {


    /*@INIT SM = ('empty', @int m [, @int n])
      Create a new empty (i.e. full of zeros) sparse matrix, of dimensions
      `m x n`. If `n` is omitted, the matrix dimension is `m x m`.@*/
    sub_command
      ("empty", 1, 2, 0, 1,
       size_type m = in.pop().to_integer(1,INT_MAX); size_type  n = m;
       if (in.remaining()) n = in.pop().to_integer(1, INT_MAX);
       gsp->allocate(m, n, gsparse::WSCMAT, gsparse::REAL);
       );

    /*@INIT SM = ('copy', @mat K [, @PYTHON{@list} I [, @PYTHON{@list} J]])
      Duplicate a matrix `K` (which might be a @tsp@MATLAB{ or a native matlab
      sparse matrix}). If index `I` and/or `J` are given, the matrix will
      be a submatrix of `K`. For example::

        @MATLAB{m = SPMAT:INIT('copy', sprand(50,50,.1), 1:40, [6 7 8 3 10])}
        @SCILAB{m = SPMAT:INIT('copy', sprand(50,50,.1), 1:40, [6 7 8 3 10])}
        @PYTHON{m = SPMAT:INIT('copy', SPMAT:INIT('empty',50,50), range(40), [6, 7, 8, 3, 10])}

      will return a 40x5 matrix.@*/
    sub_command
      ("copy", 1, 3, 0, 1,
       std::shared_ptr<gsparse> A = in.pop().to_sparse();
       if (!A->is_complex()) {
	 copy_spmat(*A, *gsp, in, scalar_type());
       } else {
	 copy_spmat(*A, *gsp, in, complex_type());
       }
       );

     /*@INIT SM = ('identity', @int n)
       Create a `n x n` identity matrix.@*/
    sub_command
      ("identity", 1, 1, 0, 1,
       size_type n = in.pop().to_integer(1, INT_MAX);
       gsp->real_wsc(new gsparse::t_wscmat_r(n,n));
       gmm::copy(gmm::identity_matrix(), gsp->real_wsc());
       );


     /*@INIT SM = ('mult', @tspmat A, @tspmat B)
       Create a sparse matrix as the product of the sparse matrices `A` and
       `B`. It requires that `A` and `B` be both real or both complex, you
       may have to use ``SPMAT:SET('to_complex')`` @*/
    sub_command
      ("mult", 2, 2, 0, 1,
       std::shared_ptr<gsparse> A = in.pop().to_sparse();
       std::shared_ptr<gsparse> B = in.pop().to_sparse();
       size_type m = A->nrows(); size_type n = B->ncols();

       if (A->is_complex() != B->is_complex())
	 THROW_BADARG("cannot multiply a complex matrix with a real one, use to_complex()");
       if (!A->is_complex()) gsp->real_wsc(new gsparse::t_wscmat_r(m,n));
       else gsp->cplx_wsc(new gsparse::t_wscmat_c(m,n));
       if (A->storage() == gsparse::CSCMAT
	   && B->storage() == gsparse::CSCMAT) {
	 if (!A->is_complex())
	   gmm::mult(A->real_csc(),B->real_csc(), gsp->real_wsc());
	 else
	   gmm::mult(A->cplx_csc(),B->cplx_csc(), gsp->cplx_wsc());
       }
       else if (A->storage() == gsparse::CSCMAT
		&& B->storage() == gsparse::WSCMAT) {
	 if (!A->is_complex())
	   gmm::mult(A->real_csc(),B->real_wsc(), gsp->real_wsc());
	 else
	   gmm::mult(A->cplx_csc(),B->cplx_wsc(), gsp->cplx_wsc());
       }
       else if (A->storage() == gsparse::WSCMAT
		&& B->storage() == gsparse::CSCMAT) {
	 if (!A->is_complex())
	   gmm::mult(A->real_wsc(),B->real_csc(), gsp->real_wsc());
	 else
	   gmm::mult(A->cplx_wsc(),B->cplx_csc(), gsp->cplx_wsc());
       }
       else if (A->storage() == gsparse::WSCMAT
		&& B->storage() == gsparse::WSCMAT) {
	 if (!A->is_complex())
	   gmm::mult(A->real_wsc(),B->real_wsc(), gsp->real_wsc());
	 else
	   gmm::mult(A->cplx_wsc(),B->cplx_wsc(), gsp->cplx_wsc());
       } else THROW_INTERNAL_ERROR;
       );

     /*@INIT SM = ('add', @tspmat A, @tspmat B)
       Create a sparse matrix as the sum of the sparse matrices `A` and `B`.
       Adding a real matrix with a complex matrix is possible.@*/
    sub_command
      ("add", 2, 2, 0, 1,
       std::shared_ptr<gsparse> A = in.pop().to_sparse();
       std::shared_ptr<gsparse> B = in.pop().to_sparse();
       size_type m = A->nrows(); size_type n = A->ncols();
       if (A->is_complex() != B->is_complex()) {
	 gsp->cplx_wsc(new gsparse::t_wscmat_c(m,n));
	 if (A->is_complex())
	   gf_spmat_add(*gsp, *B, *A, scalar_type(), complex_type());
	 else gf_spmat_add(*gsp, *A, *B, scalar_type(), complex_type());
       } else if (A->is_complex()) {
	 gsp->cplx_wsc(new gsparse::t_wscmat_c(m,n));
	 gf_spmat_add(*gsp, *A, *B, complex_type(), complex_type());
       } else {
	 gsp->real_wsc(new gsparse::t_wscmat_r(m,n));
	 gf_spmat_add(*gsp, *A, *B, scalar_type(), scalar_type());
       }
       );

     /*@INIT SM = ('diag', @dmat D [, @ivec E [, @int n [,@int m]]])
       Create a diagonal matrix. If `E` is given, `D` might be a matrix and
       each column of `E` will contain the sub-diagonal number that will be
       filled with the corresponding column of `D`.@*/
    sub_command
      ("diag", 1, 4, 0, 1,
       spmat_set_diag(*gsp, in, true);
       );

    /*@INIT SM = ('load','hb'|'harwell-boeing'|'mm'|'matrix-market', @str filename)
      Read a sparse matrix from an Harwell-Boeing or a Matrix-Market file
      @MATLAB{See also ``::UTIL('load matrix')``}.@*/
    sub_command
      ("load", 2, 2, 1, 1,
       load_spmat(in, *gsp);
       );

  }


  if (m_in.narg() < 1)  THROW_BADARG( "Wrong number of input arguments");

  auto gsp = std::make_shared<gsparse>();

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

  id_type id = store_spmat_object(gsp);
  m_out.pop().from_object_id(id, SPMAT_CLASS_ID);

}

/*@PYTHONEXT
  def __getitem__(self, key):
    return getfem('spmat_get',self.id, 'full',*key)
  def __setitem__(self, key, keyval):
    getfem('spmat_set', self.id, 'assign', key[0], key[1], keyval)
  def __neg__(self):
    m=Spmat('copy',self)
    m.scale(-1)
    return m
  def __add__(self, other):
    return Spmat('add',self,other)
  def __sub__(self, other):
    return Spmat('add',self,other.__neg__())
  def __mul__(self, other):
    """Multiplication of a Spmat with another Spmat or a vector or a scalar.

       The result is another Spmat object.
    """
    if isinstance(other,numbers.Number):
      m = Spmat('copy',self)
      m.set('scale',other)
    elif (isinstance(other,list) or isinstance(other, numpy.ndarray)):
      m = self.mult(other)
    else:
      m = Spmat('mult',self,other)
    return m
  def __rmul__(self, other):
    if isinstance(other,numbers.Number):
      m=Spmat('copy',self)
      m.set('scale',other)
    elif (isinstance(other,list) or isinstance(other, numpy.ndarray)):
      m=self.tmult(other)
    else:
      m=Spmat('mult',other,self)
    return m;
  @*/
