/* -*- c++ -*- (enables emacs c++ mode) */
/*===========================================================================

 Copyright (C) 2001-2020 Julien Pommier.

 This file is a part of GetFEM++

 GetFEM++  is  free software;  you  can  redistribute  it  and/or modify it
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

 As a special exception, you  may use  this file  as it is a part of a free
 software  library  without  restriction.  Specifically,  if   other  files
 instantiate  templates  or  use macros or inline functions from this file,
 or  you compile this  file  and  link  it  with other files  to produce an
 executable, this file  does  not  by itself cause the resulting executable
 to be covered  by the GNU Lesser General Public License.  This   exception
 does not  however  invalidate  any  other  reasons why the executable file
 might be covered by the GNU Lesser General Public License.

===========================================================================*/


#ifndef GETFEMINT_GSPARSE_H__
#define GETFEMINT_GSPARSE_H__

#include <getfemint.h>

namespace getfemint
{
  class gsparse : virtual public dal::static_stored_object {
  public:
    typedef enum { REAL, COMPLEX } value_type;
    typedef enum { WSCMAT, CSCMAT } storage_type;
    value_type v; storage_type s;
    typedef gmm::col_matrix<gmm::wsvector<scalar_type> > t_wscmat_r;
    typedef gmm::col_matrix<gmm::wsvector<complex_type> > t_wscmat_c;
    typedef gmm::csc_matrix<scalar_type> t_cscmat_r;
    typedef gmm::csc_matrix<complex_type> t_cscmat_c;
    typedef gmm::csc_matrix_ref<const scalar_type*, const unsigned int *,
                                const unsigned int *> t_cscmat_ref_r;
    typedef gmm::csc_matrix_ref<const complex_type*, const unsigned int *,
                                const unsigned int *> t_cscmat_ref_c;
    t_wscmat_r *pwscmat_r;
    t_wscmat_c *pwscmat_c;
    t_cscmat_r *pcscmat_r;
    t_cscmat_c *pcscmat_c;
    const gfi_array  *gfimat;

    void swap(gsparse &other);
    gsparse& destructive_assign(t_wscmat_r &M);
    gsparse& destructive_assign(t_wscmat_c &M);
    gsparse& destructive_assign(t_cscmat_r &M);
    gsparse& destructive_assign(t_cscmat_c &M);

    t_wscmat_r& real_wsc(t_wscmat_r *p = 0)
    { if (p) { v = REAL; pwscmat_r = p; } return *pwscmat_r; }
    t_wscmat_c& cplx_wsc(t_wscmat_c *p = 0)
    { if (p) { v = COMPLEX; pwscmat_c = p; }return *pwscmat_c; }
    t_cscmat_ref_r real_csc() { 
      if (gfimat && !gfi_array_is_complex(gfimat)) 
        return t_cscmat_ref_r(gfi_sparse_get_pr(gfimat), 
                              gfi_sparse_get_ir(gfimat), 
                              gfi_sparse_get_jc(gfimat),
                              gfi_array_get_dim(gfimat)[0],
                              gfi_array_get_dim(gfimat)[1]);
      else if (pcscmat_r)
        return t_cscmat_ref_r(&pcscmat_r->pr[0], &pcscmat_r->ir[0],
                              &pcscmat_r->jc[0], 
                              gmm::mat_nrows(*pcscmat_r),
                              gmm::mat_ncols(*pcscmat_r));
      else THROW_INTERNAL_ERROR;
    }
    t_cscmat_r& real_csc_w(t_cscmat_r *p=0)
    { if (p) { v = REAL; pcscmat_r = p; } return *pcscmat_r; }
    t_cscmat_ref_c cplx_csc() { 
      if (gfimat && gfi_array_is_complex(gfimat)) {
        // unsigned nc = gfi_array_get_dim(gfimat)[1];
        // cout << "jc = [";
        // for (unsigned i = 0; i <= nc; ++i)
        //   cout << (gfi_sparse_get_jc(gfimat))[i] << ", ";
        // cout << "]" << endl;
        // cout << "ir = [";
        // for (unsigned i = 0; i < (gfi_sparse_get_jc(gfimat))[nc]; ++i)
        //   cout << (gfi_sparse_get_ir(gfimat))[i] << ", ";
        // cout << "]" << endl;
        // cout << "pr = [";
        // for (unsigned i = 0; i < 2*(gfi_sparse_get_jc(gfimat))[nc]; ++i)
        //   cout << (gfi_sparse_get_pr(gfimat))[i] << ", ";
        // cout << "]" << endl;
        return t_cscmat_ref_c((complex_type*)gfi_sparse_get_pr(gfimat), 
                              gfi_sparse_get_ir(gfimat), 
                              gfi_sparse_get_jc(gfimat),
                              gfi_array_get_dim(gfimat)[0],
                              gfi_array_get_dim(gfimat)[1]);
      }
      else if (pcscmat_c)
        return t_cscmat_ref_c(&pcscmat_c->pr[0], &pcscmat_c->ir[0],
                              &pcscmat_c->jc[0], 
                              gmm::mat_nrows(*pcscmat_c),
                              gmm::mat_ncols(*pcscmat_c)); 
      else THROW_INTERNAL_ERROR;
    }
    t_cscmat_c& cplx_csc_w(t_cscmat_c *p=0)
    { if (p) { v = COMPLEX; pcscmat_c = p; } return *pcscmat_c; }

    /* overloaded versions, useful for templates */
    t_wscmat_r& wsc(double, t_wscmat_r *p=0) { return real_wsc(p); }
    t_wscmat_c& wsc(complex_type, t_wscmat_c *p=0) { return cplx_wsc(p); }
    t_cscmat_ref_r csc(double) { return real_csc(); }
    t_cscmat_ref_c csc(complex_type) { return cplx_csc(); }
    t_cscmat_r& csc_w(double, t_cscmat_r *p=0) { return real_csc_w(p); }
    t_cscmat_c& csc_w(complex_type, t_cscmat_c *p=0) { return cplx_csc_w(p); }

    gsparse() : v(REAL), s(WSCMAT), pwscmat_r(0), pwscmat_c(0),
                pcscmat_r(0), pcscmat_c(0), gfimat(0) {}
    gsparse(size_type m, size_type n, storage_type s_, value_type v_ = REAL);
    gsparse(const gfi_array *a);
    storage_type storage() const { return s; }
    bool is_complex() const { return v == COMPLEX; }
    bool is_a_native_matrix_ref() const 
    { return (gfimat != 0); }
    void destroy();
    void allocate(size_type m, size_type n, storage_type s_,
                  value_type v_ = REAL);
    void allocate(size_type m, size_type n, storage_type s_,
                  scalar_type) { allocate(m,n,s_, REAL); }
    void allocate(size_type m, size_type n, storage_type s_,
                  complex_type) { allocate(m,n,s_, COMPLEX); }
    void deallocate(storage_type s_, value_type v_);
    void deallocate() { deallocate(s,v); }
    void to_wsc();
    void to_csc();
    void to_complex();
    size_type memsize() const { return 0; /* TODO ! */ }
    size_type ncols() const;
    size_type nrows() const;
    size_type nnz() const;
    const char *name() { if (s == WSCMAT) return "WSC"; else return "CSC"; }
    template <typename V1, typename V2> void 
    mult_or_transposed_mult(const V1 &vv, V2 &ww, bool tmult) {
      typedef typename gmm::linalg_traits<V1>::value_type T;
      switch (storage()) {
        case CSCMAT: 
          if (!tmult) gmm::mult(csc(T()), vv, ww);
          else gmm::mult(gmm::conjugated(csc(T())), vv, ww);
          break;
        case WSCMAT: 
          if (!tmult) gmm::mult(wsc(T()), vv, ww);
          else gmm::mult(gmm::conjugated(wsc(T())), vv, ww);
          break;
        default: THROW_INTERNAL_ERROR;
      }
    }
    virtual ~gsparse() { destroy(); }
  };

  void spmat_set_diag(gsparse &gsp, mexargs_in& in, bool create_matrix);
  void spmat_load(mexargs_in& in, mexargs_out& out,
                  mexarg_out::output_sparse_fmt fmt);
}  /* end of namespace getfemint.                                          */

#endif /* GETFEMINT_GSPARSE_H__                                            */
