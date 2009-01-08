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

#include <getfemint_workspace.h>
#include <getfemint_precond.h>
#include <getfemint_gsparse.h>

using namespace getfemint;

/* moved out of the gprecond class to avoid
   the stupid inlining bug of gcc-3.3 and gcc-3.4
   ( http://gcc.gnu.org/bugzilla/show_bug.cgi?id=16849 ) */

   namespace getfemint {
     template<typename T> size_type
     gprecond<T>::memsize() const {
       size_type sz = sizeof(*this);
       switch (type) {
	 case IDENTITY: break;
	 case DIAG:    sz += diagonal->memsize(); break;
	 case ILUT:    sz += ilut->memsize(); break;
	 case ILU:     sz += ilu->memsize(); break;
	 case ILDLT:   sz += ildlt->memsize(); break;
	 case ILDLTT:  sz += ildltt->memsize(); break;
	 case SUPERLU:
	   sz += size_type(superlu->memsize()); break;
	 case SPMAT:   sz += gsp->memsize(); break;
       }
       return sz;
     }
   }


template <typename T> static gprecond<T>&
precond_new(mexargs_out& out, T) {
  getfemint_precond *precond = new getfemint_precond(gsparse::REAL);
  out.pop().from_object_id(workspace().push_object(precond), PRECOND_CLASS_ID);
  return precond->precond(T());
}

template <typename T> static void
precond_diagonal(gsparse &M, mexargs_out& out, T) {
  gprecond<T> &p = precond_new(out, T());
  p.type = gprecond_base::DIAG;
  p.diagonal.reset(new gmm::diagonal_precond<typename gprecond<T>::cscmat>(M.csc(T())));
}
template <typename T> static void
precond_ildlt(gsparse &M, mexargs_out& out, T) {
  gprecond<T> &p = precond_new(out, T());
  p.type = gprecond_base::ILDLT;
  p.ildlt.reset(new gmm::ildlt_precond<typename gprecond<T>::cscmat>(M.csc(T())));
}
template <typename T> static void
precond_ilu(gsparse &M, mexargs_out& out, T) {
  gprecond<T> &p = precond_new(out, T());
  p.type = gprecond_base::ILU;
  p.ilu.reset(new gmm::ilu_precond<typename gprecond<T>::cscmat>(M.csc(T())));
}
template <typename T> static void
precond_ildltt(gsparse &M, int additional_fillin, double threshold, mexargs_out& out, T) {
  gprecond<T> &p = precond_new(out, T());
  p.type = gprecond_base::ILDLTT;
  p.ildltt.reset(new gmm::ildltt_precond<typename gprecond<T>::cscmat>(M.csc(T()), additional_fillin, threshold));
}
template <typename T> static void
precond_ilut(gsparse &M, int additional_fillin, double threshold, mexargs_out& out, T) {
  gprecond<T> &p = precond_new(out, T());
  p.type = gprecond_base::ILUT;
  p.ilut.reset(new gmm::ilut_precond<typename gprecond<T>::cscmat>(M.csc(T()), additional_fillin, threshold));
}

template <typename T> static void
precond_superlu(gsparse &M, mexargs_out& out, T) {
  gprecond<T> &p = precond_new(out, T());
  p.type = gprecond_base::SUPERLU;
  p.superlu.reset(new gmm::SuperLU_factor<T>()); p.superlu.get()->build_with(M.csc(T()));
}

static void
precond_spmat(getfemint_gsparse *gsp, mexargs_out& out) {
  getfemint_precond *precond = new getfemint_precond(gsparse::REAL);
  out.pop().from_object_id(workspace().push_object(precond), PRECOND_CLASS_ID);
  precond->precond(scalar_type()).type = gprecond_base::SPMAT;
  precond->precond(scalar_type()).gsp = gsp;
  workspace().set_dependance(precond, gsp);
}

/*MLABCOM
  FUNCTION F=gf_precond(name_of_preconditioner, args)
  @TEXT PRECOND:INIT('PRECOND_init')
MLABCOM*/

/*@TEXT PRECOND:INIT('PRECOND_init')
The preconditioners may store REAL or COMPLEX values. They accept<par>
getfem sparse matrices and Matlab sparse matrices.<Par>

* PRECOND:INIT('identity')<par>
   Create a REAL identity precondioner.<par>
* PRECOND:INIT('cidentity')<par>
   Create a COMPLEX identity precondioner.<par>
* PRECOND:INIT('diagonal', @dcvec D)<par>
   Create a diagonal precondioner.<par>
* PRECOND:INIT('ildlt', @tsp m)<par>
   Create an ILDLT (Cholesky) preconditioner for the (symmetric)<par>
   sparse matrix `m`. This preconditioner has the same sparsity<par>
   pattern than `m` (no fill-in).<par>
* PRECOND:INIT('ilu', @tsp m)<par>
   Create an ILU (Incomplete LU) preconditioner for the sparse<par>
   matrix `m`. This preconditioner has the same sparsity pattern<par>
   than `m` (no fill-in).<par>
* PRECOND:INIT('ildltt', @tsp m[, @int fillin[, @scalar threshold]])<par>
   Create an ILDLT (Cholesky with filling) preconditioner for the<par>
   (symmetric) sparse matrix `m`. The preconditioner may add at most<par>
   `fillin` additional non-zero entries on each line. The default<par>
   value for `fillin` is 10, and the default threshold is1e-7.<par>
* PRECOND:INIT('ilut', @tsp m[, @int fillin[, @scalar threshold]])<par>
   Create an ILUT (Incomplete LU with filling) preconditioner for the<par>
   sparse matrix `m`. The preconditioner may add at most `fillin`<par>
   additional non-zero entries on each line. The default value for<par>
   `fillin` is 10, and the default threshold is 1e-7.<par>
* PRECOND:INIT('superlu', @tsp m)<par>
   Uses SuperLU to build an exact factorization of the sparse matrix<par>
   `m`. This preconditioner is only available if the getfem-interface<par>
   was built with SuperLU support. Note that LU factorization is likely<par>
   to eat all your memory for 3D problems.
* PRECOND:INIT('spmat', @tsp M)<par>
   ???
  @*/

void gf_precond(getfemint::mexargs_in& in, getfemint::mexargs_out& out)
{
  if (in.narg() < 1) {
    THROW_BADARG( "Wrong number of input arguments");
  }
  std::string cmd = in.pop().to_string();

  if (check_cmd(cmd, "identity", in, out, 0, 0, 0, 1)) {
    precond_new(out, scalar_type());
  } else if (check_cmd(cmd, "cidentity", in, out, 0, 0, 0, 1)) {
    precond_new(out, complex_type());
  } else if (check_cmd(cmd, "diagonal", in, out, 1, 1, 0, 1)) {
    dal::shared_ptr<gsparse> M = in.pop().to_sparse(); M->to_csc();
    if (M->is_complex()) precond_diagonal(*M, out, complex_type());
    else                 precond_diagonal(*M, out, scalar_type());
  } else if (check_cmd(cmd, "ildlt", in, out, 1, 1, 0, 1)) {
    dal::shared_ptr<gsparse> M = in.pop().to_sparse(); M->to_csc();
    if (M->is_complex()) precond_ildlt(*M, out, complex_type());
    else                 precond_ildlt(*M, out, scalar_type());
  } else if (check_cmd(cmd, "ilu", in, out, 1, 1, 0, 1)) {
    dal::shared_ptr<gsparse> M = in.pop().to_sparse(); M->to_csc();
    if (M->is_complex()) precond_ilu(*M, out, complex_type());
    else                 precond_ilu(*M, out, scalar_type());
  } else if (check_cmd(cmd, "ildltt", in, out, 1, 3, 0, 1)) {
    dal::shared_ptr<gsparse> M = in.pop().to_sparse(); M->to_csc();
    int additional_fillin = 10; scalar_type threshold = 1e-7;
    if (in.remaining()) additional_fillin = in.pop().to_integer(0, 100000);
    if (in.remaining()) threshold = in.pop().to_scalar(0,1e30);
    if (M->is_complex()) precond_ildltt(*M, additional_fillin, threshold, out, complex_type());
    else                 precond_ildltt(*M, additional_fillin, threshold, out, scalar_type());
  } else if (check_cmd(cmd, "ilut", in, out, 1, 3, 0, 1)) {
    dal::shared_ptr<gsparse> M = in.pop().to_sparse(); M->to_csc();
    int additional_fillin = 10; scalar_type threshold = 1e-7;
    if (in.remaining()) additional_fillin = in.pop().to_integer(0, 100000);
    if (in.remaining()) threshold = in.pop().to_scalar(0,1e30);
    if (M->is_complex()) precond_ilut(*M, additional_fillin, threshold, out, complex_type());
    else                 precond_ilut(*M, additional_fillin, threshold, out, scalar_type());
  } else if (check_cmd(cmd, "superlu", in, out, 1, 1, 0, 1)) {
    dal::shared_ptr<gsparse> M = in.pop().to_sparse(); M->to_csc();
    if (M->is_complex()) precond_superlu(*M, out, complex_type());
    else                 precond_superlu(*M, out, scalar_type());
  } else if (check_cmd(cmd, "spmat", in, out, 1, 1, 0, 1)) {
    getfemint_gsparse *ggsp = 0;
    if (in.front().is_gsparse()) {
      ggsp = in.pop().to_getfemint_gsparse();
    } else {
      ggsp = new getfemint_gsparse();
      dal::shared_ptr<gsparse> src = in.pop().to_sparse();
      if (src->is_complex()) {
	ggsp->sparse().allocate(src->nrows(), src->ncols(), src->storage(), complex_type());
	gmm::copy(src->csc(complex_type()), ggsp->sparse().csc_w(complex_type()));
      } else {
	ggsp->sparse().allocate(src->nrows(), src->ncols(), src->storage(), scalar_type());
	gmm::copy(src->csc(scalar_type()), ggsp->sparse().csc_w(scalar_type()));
      }
      workspace().push_object(ggsp);
    }
    precond_spmat(ggsp, out);
  } else bad_cmd(cmd);
}
