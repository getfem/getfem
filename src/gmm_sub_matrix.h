/* -*- c++ -*- (enables emacs c++ mode)                                    */
/* *********************************************************************** */
/*                                                                         */
/* Library :  Generic Matrix Methods  (gmm)                                */
/* File    :  gmm_sub_vector.h : generic sub matrices.                     */
/*     									   */
/* Date : October 13, 2002.                                                */
/* Author : Yves Renard, Yves.Renard@gmm.insa-tlse.fr                      */
/*                                                                         */
/* *********************************************************************** */
/*                                                                         */
/* Copyright (C) 2001  Yves Renard.                                        */
/*                                                                         */
/* This file is a part of GETFEM++                                         */
/*                                                                         */
/* This program is free software; you can redistribute it and/or modify    */
/* it under the terms of the GNU General Public License as published by    */
/* the Free Software Foundation; version 2 of the License.                 */
/*                                                                         */
/* This program is distributed in the hope that it will be useful,         */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of          */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           */
/* GNU General Public License for more details.                            */
/*                                                                         */
/* You should have received a copy of the GNU General Public License       */
/* along with this program; if not, write to the Free Software Foundation, */
/* Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.         */
/*                                                                         */
/* *********************************************************************** */

#ifndef __GMM_SUB_MATRIX_H
#define __GMM_SUB_MATRIX_H

namespace gmm {

  /* ********************************************************************* */
  /*		sub-matrices types                                         */
  /* ********************************************************************* */

  template <class PT, class SUBI1, class SUBI2, class STORAGE, class SUBORIENT>
  struct gen_sub_matrix {
    const SUBI1 *psi1;
    const SUBI1 *psi2;
    PT pm;
    
    typedef typename std::iterator_traits<PT>::value_type M;
    typedef typename std::iterator_traits<PT>::reference ref_M;
    typedef typename linalg_traits<M>::value_type value_type;
    typedef typename mat_ref_type<PT,  M>::access_type access_type;
    
    ref_M deref(void) const { return *pm; }
    access_type operator()(size_type i, size_type j) const
    { return (*pm)(psi1->index(i), psi2->index(j)); }
    value_type operator()(size_type i, size_type j)
    { return (*pm)(psi1->index(i), psi2->index(j)); }
    size_type nrows(void) const { return psi1->size(); }
    size_type ncols(void) const { return psi2->size(); }
    
    gen_sub_matrix(ref_M mm, const SUBI1 &si1, const SUBI2 &si2)
      : psi1(&si1), psi2(&si2), pm(&mm) {}
    gen_sub_matrix() {}
  };
  
  template <class PT, class SUBI1, class SUBI2, class STORAGE>
  struct linalg_traits<gen_sub_matrix<PT, SUBI1, SUBI2, STORAGE, row_major> > {
    typedef gen_sub_matrix<PT, SUBI1, SUBI2, STORAGE, row_major> this_type;
    typedef typename std::iterator_traits<PT>::value_type M;
    typedef linalg_true is_reference;
    typedef abstract_matrix linalg_type;
    typedef typename linalg_traits<M>::value_type value_type;
    typedef typename linalg_traits<M>::reference_type reference_type;
    typedef typename linalg_traits<M>::storage_type storage_type;
    typedef abstract_null_type sub_col_type;
    typedef abstract_null_type const_sub_col_type;
    typedef typename sub_vector_type<typename linalg_traits<M>
    ::sub_row_type *, SUBI1>::vector_type sub_row_type;
    typedef typename sub_vector_type<const typename linalg_traits<M>
    ::const_sub_row_type *, SUBI1>::vector_type const_sub_row_type;
    typedef row_major sub_orientation;
    size_type nrows(const this_type &m) { return m.nrows(); }
    size_type ncols(const this_type &m) { return m.ncols(); }
    const_sub_row_type row(const this_type &m, size_type i) { 
      return sub_vector(mat_const_row(m.deref(),m.psi2->index(i)),*(m.psi1));
    }
    sub_row_type row(this_type &m, size_type i)
    { return sub_vector(mat_row(m.deref(), m.psi2->index(i)), *(m.psi1)); }
    const void* origin(const this_type &m) { return linalg_origin(m.deref()); }
    void do_clear(this_type &m)
    { for (size_type i = 0; i < nrows(m); ++i) clear(row(m,i)); }
  };

    template <class PT, class SUBI1, class SUBI2, class STORAGE>
  struct linalg_traits<gen_sub_matrix<PT, SUBI1, SUBI2, STORAGE, col_major> > {
    typedef gen_sub_matrix<PT, SUBI1, SUBI2, STORAGE, col_major> this_type;
    typedef typename std::iterator_traits<PT>::value_type M;
    typedef linalg_true is_reference;
    typedef abstract_matrix linalg_type;
    typedef typename linalg_traits<M>::value_type value_type;
    typedef typename linalg_traits<M>::reference_type reference_type;
    typedef typename linalg_traits<M>::storage_type storage_type;
    typedef abstract_null_type sub_row_type;
    typedef abstract_null_type const_sub_row_type;
    typedef typename sub_vector_type<typename linalg_traits<M>
    ::sub_col_type *, SUBI1>::vector_type sub_col_type;
    typedef typename sub_vector_type<const typename linalg_traits<M>
    ::const_sub_col_type *, SUBI1>::vector_type const_sub_col_type;
    typedef col_major sub_orientation;
    size_type nrows(const this_type &m) { return m.nrows(); }
    size_type ncols(const this_type &m) { return m.ncols(); }
    const_sub_col_type col(const this_type &m, size_type i)
    { return sub_vector(mat_const_col(m.deref(),m.psi1->index(i)),*(m.psi2)); }
    sub_col_type col(this_type &m, size_type i)
    { return sub_vector(mat_col(m.deref(), m.psi1->index(i)), *(m.psi2)); }
    const void* origin(const this_type &m) { return linalg_origin(m.deref()); }
    void do_clear(this_type &m)
    { for (size_type i = 0; i < ncols(m); ++i) clear(col(m,i)); }
  };

  /* ******************************************************************** */
  /*		sub matrices                                              */
  /* ******************************************************************** */
  
  template <class PT, class SUBI1, class SUBI2>
  struct sub_matrix_type {
    typedef typename std::iterator_traits<PT>::value_type M;
    typedef gen_sub_matrix<PT, SUBI1, SUBI2,
      typename linalg_traits<M>::storage_type,
      typename principal_orientation_type<typename
     linalg_traits<M>::sub_orientation>::potype> matrix_type;
  };

  template <class M, class SUBI1, class SUBI2>  inline
  typename sub_matrix_type<M *, SUBI1, SUBI2>::matrix_type
  sub_matrix(M &m, const SUBI1 &si1, const SUBI2 &si2) {
    return
      typename sub_matrix_type<M *, SUBI1, SUBI2>::matrix_type(m, si1, si2);
  }

  template <class M, class SUBI1>  inline
  typename sub_matrix_type<M *, SUBI1, SUBI1>::matrix_type
  sub_matrix(M &m, const SUBI1 &si1) {
    return typename sub_matrix_type<M *, SUBI1, SUBI1>::matrix_type(m, si1, si1);
  }

   template <class M, class SUBI1, class SUBI2>  inline
   typename sub_matrix_type<const M *, SUBI1, SUBI2>::matrix_type
   sub_matrix(const M &m, const SUBI1 &si1, const SUBI2 &si2) {
     return typename sub_matrix_type<const M *, SUBI1, SUBI2>::matrix_type(m,
								   si1, si2);
   }
  
  template <class M, class SUBI1>  inline
  typename sub_matrix_type<const M *, SUBI1, SUBI1>::matrix_type
  sub_matrix(const M &m, const SUBI1 &si1) {
    return typename sub_matrix_type<const M *, SUBI1, SUBI1>::matrix_type(m, si1, si1);
  }

}

#endif //  __GMM_SUB_MATRIX_H
