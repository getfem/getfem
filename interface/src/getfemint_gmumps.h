/* -*- c++ -*- (enables emacs c++ mode) */
/*===========================================================================

 Copyright (C) 2025-2025 Konstantinos Poulios.

 This file is a part of GetFEM

 GetFEM is free software;  you can  redistribute it  and/or modify it under
 the  terms  of the  GNU  Lesser General Public License as published by the
 Free Software Foundation;  either version 3  of  the License,  or (at your
 option) any  later  version  along with  the GCC Runtime Library Exception
 either version 3.1 or (at your option) any later version.
 This program  is  distributed  in  the  hope  that it will be useful,  but
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 License and GCC Runtime Library Exception for more details.
 You  should  have received a copy of the GNU Lesser General Public License
 along  with  this program;  if not, write to the Free Software Foundation,
 Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.

 As a special exception, you  may use  this file  as it is a part of a free
 software  library  without  restriction.  Specifically,  if   other  files
 instantiate  templates  or  use macros or inline functions from this file,
 or  you compile this  file  and  link  it  with other files  to produce an
 executable, this file  does  not  by itself cause the resulting executable
 to be covered  by the GNU Lesser General Public License.  This   exception
 does not  however  invalidate  any  other  reasons why the executable file
 might be covered by the GNU Lesser General Public License.

===========================================================================*/


#ifndef GETFEMINT_GMUMPS_H__
#define GETFEMINT_GMUMPS_H__

#include <getfemint.h>
#include <gmm/gmm_MUMPS_interface.h>

namespace getfemint {

  /*
    This class is just for redirecting the methods of a real/complex
    gmm::mumps_context, and can also be stored in the dal memory manager.
  */
  class gmumps : virtual public dal::static_stored_object {
    std::unique_ptr<gmm::mumps_context<scalar_type> > pmumps_ctx_r;
    std::unique_ptr<gmm::mumps_context<complex_type> > pmumps_ctx_c;
    bool symmetric, complex;
  public:
    gmumps(bool is_symmetric_, bool is_complex_)
      : symmetric(is_symmetric_),  complex(is_complex_) {
      // MUMPS is initiated in the constructor of gmm::mumps_context
      if (complex)
        pmumps_ctx_c = std::make_unique<gmm::mumps_context<complex_type> >(symmetric);
      else
        pmumps_ctx_r = std::make_unique<gmm::mumps_context<scalar_type> >(symmetric);
    }
    bool is_symmetric() const { return symmetric; }
    bool is_complex() const { return complex; }
//    gmm::mumps_context<scalar_type> &context_r() {
//      if (complex) THROW_ERROR("This is not a real number context");
//      return *pmumps_ctx_r;
//    }
//    gmm::mumps_context<complex_type> &context_c() {
//      if (!complex) THROW_ERROR("This is not a complex number context");
//      return *pmumps_ctx_c;
//    }
    int nrows() const
      { return complex ? pmumps_ctx_c->nrows() : pmumps_ctx_r->nrows(); }
    int &ICNTL(int I)
      { return complex ? pmumps_ctx_c->ICNTL(I) : pmumps_ctx_r->ICNTL(I); }
    const int &ICNTL(int I) const
      { return complex ? pmumps_ctx_c->ICNTL(I) : pmumps_ctx_r->ICNTL(I); }
    scalar_type &CNTL(int I)
      { return complex ? pmumps_ctx_c->CNTL(I) : pmumps_ctx_r->CNTL(I); }
    const scalar_type &CNTL(int I) const
      { return complex ? pmumps_ctx_c->CNTL(I) : pmumps_ctx_r->CNTL(I); }
    const int &INFO(int I) const
      { return complex ? pmumps_ctx_c->INFO(I) : pmumps_ctx_r->INFO(I); }
    const int &INFOG(int I) const
      { return complex ? pmumps_ctx_c->INFOG(I) : pmumps_ctx_r->INFOG(I); }
    const scalar_type &RINFO(int I) const
      { return complex ? pmumps_ctx_c->RINFO(I) : pmumps_ctx_r->RINFO(I); }
    const scalar_type &RINFOG(int I) const
      { return complex ? pmumps_ctx_c->RINFOG(I) : pmumps_ctx_r->RINFOG(I); }

    template<typename MAT>
    void set_matrix_r(const MAT& mat, bool distributed=false) {
      if (complex) THROW_ERROR("This is not a real number context.");
      pmumps_ctx_r->set_matrix(mat, distributed);
    }
    template<typename MAT>
    void set_matrix_c(const MAT& mat, bool distributed=false) {
      if (!complex) THROW_ERROR("This is not a complex number context.");
      pmumps_ctx_c->set_matrix(mat, distributed);
    }
    template<typename VEC>
    void set_vector_r(const VEC& vec) {
      if (complex) THROW_ERROR("This is not a real number context.");
      pmumps_ctx_r->set_vector(vec);
    }
    template<typename VEC>
    void set_vector_c(const VEC& vec) {
      if (!complex) THROW_ERROR("This is not a complex number context.");
      pmumps_ctx_c->set_vector(vec);
    }
    const std::vector<scalar_type>& vector_r() const {
      if (complex) THROW_ERROR("This is not a real number context");
      return pmumps_ctx_r->vector();
    }
    const std::vector<complex_type>& vector_c() const {
      if (!complex) THROW_ERROR("This is not a complex number context");
      return pmumps_ctx_c->vector();
    }

    inline void analyze() {
      if (complex) pmumps_ctx_c->analyze();
      else         pmumps_ctx_r->analyze();
    }
    inline void factorize() {
      if (complex) pmumps_ctx_c->factorize();
      else         pmumps_ctx_r->factorize();
    }
    inline void analyze_and_factorize() {
      if (complex) pmumps_ctx_c->analyze_and_factorize();
      else         pmumps_ctx_r->analyze_and_factorize();
    }
    inline void solve() {
      if (complex) pmumps_ctx_c->solve();
      else         pmumps_ctx_r->solve();
    }
    inline void factorize_and_solve() {
      if (complex) pmumps_ctx_c->factorize_and_solve();
      else         pmumps_ctx_r->factorize_and_solve();
    }
    inline void analyze_factorize_and_solve() {
      if (complex) pmumps_ctx_c->analyze_factorize_and_solve();
      else         pmumps_ctx_r->analyze_factorize_and_solve();
    }
    inline bool error_check() {
      return complex ? pmumps_ctx_c->error_check()
                     : pmumps_ctx_r->error_check();
    }

    inline void mpi_broadcast() {
      return complex ? pmumps_ctx_c->mpi_broadcast()
                     : pmumps_ctx_r->mpi_broadcast();
    }
  };

}

#endif /* GETFEMINT_GMUMPS_H__ */
