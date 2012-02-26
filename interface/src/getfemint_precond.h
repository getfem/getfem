/* -*- c++ -*- (enables emacs c++ mode) */
/*===========================================================================
 
 Copyright (C) 2001-2012 Julien Pommier.
 
 This file is a part of GETFEM++
 
 Getfem++  is  free software;  you  can  redistribute  it  and/or modify it
 under  the  terms  of the  GNU  Lesser General Public License as published
 by  the  Free Software Foundation;  either version 2.1 of the License,  or
 (at your option) any later version.
 This program  is  distributed  in  the  hope  that it will be useful,  but
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 License for more details.
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


#ifndef GETFEMINT_PRECOND_H__
#define GETFEMINT_PRECOND_H__

#include <getfemint_std.h>
#include <getfemint_object.h>
#include <gmm/gmm_precond_diagonal.h>
#include <gmm/gmm_precond_ildlt.h>
#include <gmm/gmm_precond_ildltt.h>
#include <gmm/gmm_precond_ilu.h>
#include <gmm/gmm_precond_ilut.h>
#include <getfemint.h>
#include <getfem/getfem_superlu.h>
#include <getfemint_gsparse.h>

namespace getfemint
{
  struct gprecond_base {
    size_type nrows_, ncols_;
    enum { IDENTITY, DIAG, ILDLT, ILDLTT, ILU, ILUT, SUPERLU, SPMAT } type;
    getfemint_gsparse *gsp;
    size_type nrows(void) const { return gsp ? gsp->sparse().nrows() : nrows_; }
    size_type ncols(void) const { return gsp ? gsp->sparse().ncols() : ncols_; }
    void set_dimensions(size_type m, size_type n) { nrows_ = m; ncols_ = n; }
    gprecond_base() : nrows_(0), ncols_(0), type(IDENTITY), gsp(0) {}
    const char *name() const { 
      const char *p[] = { "IDENTITY", "DIAG", "ILDLT", "ILDLTT", "ILU", "ILUT", "SUPERLU", "GSPARSE" };
      return p[type];
    }
    virtual size_type memsize() const = 0;
    virtual ~gprecond_base() {}
  };

  template <typename T> struct gprecond : public gprecond_base {
    typedef gmm::csc_matrix_ref<const T*, const unsigned int *, const unsigned int *> cscmat;
    std::auto_ptr<gmm::diagonal_precond<cscmat> > diagonal;
    std::auto_ptr<gmm::ildlt_precond<cscmat> > ildlt;
    std::auto_ptr<gmm::ildltt_precond<cscmat> > ildltt;
    std::auto_ptr<gmm::ilu_precond<cscmat> > ilu;
    std::auto_ptr<gmm::ilut_precond<cscmat> > ilut;
    std::auto_ptr<gmm::SuperLU_factor<T> > superlu;

    virtual size_type memsize() const; 
  };


  class getfemint_precond : public getfem_object {
  public:
    //typedef enum { REAL, COMPLEX } value_type;
    gsparse::value_type v;

    bool is_complex() { if (p.get() && p->gsp) return p->gsp->sparse().is_complex();
      else return v == gsparse::COMPLEX; }
    gprecond<scalar_type> &precond(scalar_type) {
      GMM_ASSERT1(!is_complex(), "cannot use a COMPLEX preconditionner with REAL data");
      return *static_cast<gprecond<scalar_type>*>(p.get()); 
    }

    gprecond<complex_type> &precond(complex_type) { 
      GMM_ASSERT1(is_complex(), "cannot use a REAL preconditionner with COMPLEX data");
      return *static_cast<gprecond<complex_type>*>(p.get());
    }

    gprecond_base &bprecond() { return *(p.get()); }
    
    size_type nrows(void) const { return p->nrows(); }
    size_type ncols(void) const { return p->ncols(); }
    void set_dimensions(size_type m, size_type n) { p->set_dimensions(m,n); }

    getfemint_precond(gsparse::value_type v_) : v(v_) { 
      if (!is_complex()) p.reset(new gprecond<scalar_type>()); 
      else p.reset(new gprecond<complex_type>());
    }
    ~getfemint_precond() {}
    id_type class_id() const { return PRECOND_CLASS_ID; }
    /*size_type memsize() const {
      return (p.get() ? p.get()->memsize() : 0);
      }*/
  private:
    std::auto_ptr<gprecond_base> p;
  };

  inline bool object_is_precond(getfem_object *o) {
    return (o->class_id() == PRECOND_CLASS_ID);
  }

  inline getfemint_precond* object_to_precond(getfem_object *o) {
    if (object_is_precond(o)) return ((getfemint_precond*)o);
    else THROW_INTERNAL_ERROR;
  }
}  /* end of namespace getfemint.                                          */

namespace gmm {
  template<typename T>
  struct linalg_traits<getfemint::gprecond<T> > {
    typedef getfemint::gprecond<T> this_type;
    typedef linalg_const is_reference;
    typedef abstract_matrix linalg_type;
    typedef T value_type;
    typedef T& reference;
    typedef value_type origin_type;
    typedef abstract_sparse storage_type;
    typedef abstract_null_type sub_row_type;
    typedef abstract_null_type const_sub_row_type;
    typedef abstract_null_type row_iterator;
    typedef abstract_null_type const_row_iterator;
    typedef abstract_null_type sub_col_type;
    typedef abstract_null_type const_sub_col_type;
    typedef abstract_null_type const_col_iterator;
    typedef abstract_null_type col_iterator;
    typedef col_major sub_orientation;
    typedef linalg_true index_sorted;
    static size_type nrows(const this_type &m) { return m.nrows(); }
    static size_type ncols(const this_type &m) { return m.ncols(); }
  };

  template <typename T, typename V1, typename V2>
  void mult_or_transposed_mult(const getfemint::gprecond<T>& precond, const V1 &v, V2 &w, bool do_mult) {
    switch (precond.type) {
      case getfemint::gprecond_base::IDENTITY: gmm::copy(v,w); break;
      case getfemint::gprecond_base::DIAG: gmm::mult(*precond.diagonal, v, w); break;
      case getfemint::gprecond_base::ILDLT: 
        if (do_mult) gmm::mult(*precond.ildlt, v, w); 
        else gmm::transposed_mult(*precond.ildlt, v, w);
        break;
      case getfemint::gprecond_base::ILDLTT: 
        if (do_mult) gmm::mult(*precond.ildltt, v, w); 
        else gmm::transposed_mult(*precond.ildltt, v, w);
        break;
      case getfemint::gprecond_base::ILU: 
        if (do_mult) gmm::mult(*precond.ilu, v, w); 
        else gmm::transposed_mult(*precond.ilu, v, w);
        break;
      case getfemint::gprecond_base::ILUT: 
        if (do_mult) gmm::mult(*precond.ilut, v, w); 
        else gmm::transposed_mult(*precond.ilut, v, w);
        break;
      case getfemint::gprecond_base::SUPERLU:
	if (do_mult) precond.superlu->solve(w,v);
	else         precond.superlu->solve(w,v,gmm::SuperLU_factor<T>::LU_TRANSP);
	break;
      case getfemint::gprecond_base::SPMAT:
	precond.gsp->sparse().mult_or_transposed_mult(v, w, !do_mult);
        break;
    }
  }
  template <typename T, typename V1, typename V2>
  void mult(const getfemint::gprecond<T>& precond, const V1 &v, V2 &w) {
    mult_or_transposed_mult(precond,v,w,true);
  }
  template <typename T, typename V1, typename V2>
  void transposed_mult(const getfemint::gprecond<T>& precond, const V1 &v, V2 &w) {
    mult_or_transposed_mult(precond,v,w,false);
  }

}
#endif /* GETFEMINT_PRECOND_H__                                         */
