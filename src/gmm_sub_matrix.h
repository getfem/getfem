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
  /*		sparse sub-matrices                                        */
  /* ********************************************************************* */

  template <class PT, class IT1, class IT2> struct sparse_row_sub_matrix {
    const reverse_index *_r_i;
    IT1 begin1, end1;
    IT2 begin2, end2;
    PT m;

    typedef typename std::iterator_traits<PT>::value_type M;
    typedef typename std::iterator_traits<PT>::reference ref_M;
    typedef typename linalg_traits<M>::value_type value_type;
    typedef typename mat_ref_type<PT,  M>::access_type access_type;

    ref_M deref(void) const { return *m; }
    const reverse_index &rindex(void) const { return *_r_i; }
    access_type operator()(size_type i, size_type j) const
    { return (*m)(begin1[i], begin2[j]); }
    value_type operator()(size_type i, size_type j)
    { return (*m)(begin1[i], begin2[j]); }
    size_type nrows(void) const { return end1 - begin1; }
    size_type ncols(void) const { return end2 - begin2; }

    sparse_row_sub_matrix(ref_M mm, const IT1 &it1, const IT1 &e1,
	         const IT1 &it2, const IT1 &e2, const reverse_index &rindex1,
			  const reverse_index &)
      : _r_i(&rindex1), begin1(it1), end1(e1), begin2(it2), end2(e2), m(&mm) {}
    sparse_row_sub_matrix() {}
  };

  template <class PT, class IT1, class IT2>
  struct linalg_traits<sparse_row_sub_matrix<PT, IT1, IT2> > {
    typedef sparse_row_sub_matrix<PT, IT1, IT2> this_type;
    typedef typename std::iterator_traits<PT>::value_type M;
    typedef linalg_true is_reference;
    typedef abstract_matrix linalg_type;
    typedef typename linalg_traits<M>::value_type value_type;
    typedef typename linalg_traits<M>::reference_type reference_type;
    typedef abstract_sparse storage_type;
    typedef sparse_sub_vector<typename linalg_traits<M>
    ::sub_row_type *, IT1> sub_row_type;
    typedef sparse_sub_vector<const typename linalg_traits<M>
    ::const_sub_row_type *, IT1> const_sub_row_type;
    typedef abstract_null_type sub_col_type;
    typedef abstract_null_type const_sub_col_type;
    typedef row_major sub_orientation;
    size_type nrows(const this_type &m) { return m.nrows(); }
    size_type ncols(const this_type &m) { return m.ncols(); }
    const_sub_row_type row(const this_type &m, size_type i) { 
      return const_sub_row_type(mat_row(const_cast<const M &>(m.deref()), (m.begin2)[i]),
				m.begin1, m.end1, m.rindex());
    }
    sub_row_type row(this_type &m, size_type i)
    { return sub_row_type(mat_row(m.deref(), (m.begin2)[i]),
			  m.begin1, m.end1, m.rindex());
    }
    const void* origin(const this_type &m) { return linalg_origin(m.deref()); }
    void do_clear(this_type &m)
      { for (size_type i = 0; i < nrows(); ++i) clear(row(m,i)); }
  };

  template <class PT, class IT1, class IT2> struct sparse_col_sub_matrix {

    const reverse_index *_r_i;
    IT1 begin1, end1;
    IT2 begin2, end2;
    PT m;

    typedef typename std::iterator_traits<PT>::value_type M;
    typedef typename std::iterator_traits<PT>::reference ref_M;
    typedef typename linalg_traits<M>::value_type value_type;
    typedef typename mat_ref_type<PT,  M>::access_type access_type;

    const reverse_index &rindex(void) const { return *_r_i; }
    ref_M deref(void) const { return *m; }
    value_type operator()(size_type i, size_type j) const 
    { return (*m)(begin1[i], begin2[j]); }
    access_type operator()(size_type i, size_type j)
    { return (*m)(begin1[i], begin2[j]); }
    size_type nrows(void) const { return end1 - begin1; }
    size_type ncols(void) const { return end2 - begin2; }

    sparse_col_sub_matrix(ref_M mm, const IT1 &it1, const IT1 &e1,
		 const IT1 &it2, const IT1 &e2, const reverse_index &,
			  const reverse_index &rindex2)
      : _r_i(&rindex2), begin1(it1), end1(e1), begin2(it2), end2(e2), m(&mm) {}
    sparse_col_sub_matrix() {}

  };

  template <class PT, class IT1, class IT2>
  struct linalg_traits<sparse_col_sub_matrix<PT, IT1, IT2> > {
    typedef sparse_col_sub_matrix<PT, IT1, IT2> this_type;
    typedef typename std::iterator_traits<PT>::value_type M;
    typedef linalg_true is_reference;
    typedef abstract_matrix linalg_type;
    typedef typename linalg_traits<M>::value_type value_type;
    typedef typename linalg_traits<M>::value_type reference_type;
    typedef abstract_sparse storage_type;
    typedef sparse_sub_vector<typename linalg_traits<M>
    ::sub_col_type *, IT2> sub_col_type;
    typedef sparse_sub_vector<const typename linalg_traits<M>
    ::const_sub_col_type *, IT2> const_sub_col_type;
    typedef abstract_null_type sub_row_type;
    typedef abstract_null_type const_sub_row_type;
    typedef col_major sub_orientation;
    size_type nrows(const this_type &m) { return m.nrows(); }
    size_type ncols(const this_type &m) { return m.ncols(); }
    const_sub_col_type col(const this_type &m, size_type i) { 
      return const_sub_col_type(mat_col(const_cast<const M&>(m.deref()), (m.begin1)[i]),
				m.begin2, m.end2, m.rindex());
    }
    sub_col_type col(this_type &m, size_type i)
    { return sub_col_type(mat_row(m.deref(), (m.begin1)[i]),
			  m.begin2, m.end2, m.rindex());
    }
    const void* origin(const this_type &m) { return linalg_origin(m.deref()); }
    void do_clear(this_type &m)
      { for (size_type i = 0; i < ncols(); ++i) clear(col(m,i)); }
  };


  /* ********************************************************************* */
  /*		plain sub-matrices                                         */
  /* ********************************************************************* */

  template <class PT, class IT1, class IT2> struct plain_row_sub_matrix {
    IT1 begin1, end1;
    IT2 begin2, end2;
    PT m;

    typedef typename std::iterator_traits<PT>::value_type M;
    typedef typename std::iterator_traits<PT>::reference ref_M;
    typedef typename linalg_traits<M>::value_type value_type;
    typedef typename mat_ref_type<PT,  M>::access_type access_type;

    ref_M deref(void) const { return *m; }
    size_type nrows(void) const { return end1 - begin1; }
    size_type ncols(void) const { return end2 - begin2; }
    value_type operator()(size_type i, size_type j) const
    { return (*m)(begin1[i], begin2[j]); }
    access_type operator()(size_type i, size_type j)
    { return (*m)(begin1[i], begin2[j]); }
    plain_row_sub_matrix(ref_M mm, const IT1 &it1, const IT1 &e1,
			 const IT1 &it2, const IT1 &e2)
      : begin1(it1), end1(e1), begin2(it2), end2(e2), m(&mm) {}
    plain_row_sub_matrix() {}
  };

  template <class PT, class IT1, class IT2>
  struct linalg_traits<plain_row_sub_matrix<PT, IT1, IT2> > {
    typedef plain_row_sub_matrix<PT, IT1, IT2> this_type;
    typedef typename std::iterator_traits<PT>::value_type M;
    typedef linalg_true is_reference;
    typedef abstract_matrix linalg_type;
    typedef typename linalg_traits<M>::value_type value_type;
    typedef typename linalg_traits<M>::reference_type reference_type;
    typedef abstract_plain storage_type;
    typedef tab_ref_index_ref_with_origin<typename linalg_traits<typename
    linalg_traits<M>::sub_row_type>::iterator, IT1> sub_row_type;
    typedef tab_ref_index_ref_with_origin<typename linalg_traits<typename
    linalg_traits<M>::const_sub_row_type>::const_iterator, IT1>
    const_sub_row_type;
    typedef abstract_null_type sub_col_type;
    typedef abstract_null_type const_sub_col_type;
    typedef row_major sub_orientation;
    size_type nrows(const this_type &m) { return m.nrows(); }
    size_type ncols(const this_type &m) { return m.ncols(); }
    const_sub_row_type row(const this_type &m, size_type i) { 
      return
	const_sub_row_type(vect_const_begin(mat_const_row(m.deref(), (m.begin2)[i])),
			   m.begin1, m.end1, linalg_origin(m.deref()));
    }
    sub_row_type row(this_type &m, size_type i)
    { return sub_row_type(vect_begin(mat_row(m.deref(), (m.begin2)[i])),
			  m.begin1, m.end1, linalg_origin(m.deref()));
    }
    const void* origin(const this_type &m) { return linalg_origin(m.deref()); }
    void do_clear(this_type &m)
      { for (size_type i = 0; i < nrows(); ++i) clear(row(m,i)); }
  };

  template <class PT, class IT1, class IT2> struct plain_col_sub_matrix {

    IT1 begin1, end1;
    IT2 begin2, end2;
    PT m;

    typedef typename std::iterator_traits<PT>::value_type M;
    typedef typename std::iterator_traits<PT>::reference ref_M;
    typedef typename linalg_traits<M>::value_type value_type;
    typedef typename mat_ref_type<PT,  M>::access_type access_type;

    ref_M deref(void) const { return *m; }
    value_type operator()(size_type i, size_type j) const 
    { return (*m)(begin1[i], begin2[j]); }
    access_type operator()(size_type i, size_type j)
    { return (*m)(begin1[i], begin2[j]); }
    size_type nrows(void) const { return end1 - begin1; }
    size_type ncols(void) const { return end2 - begin2; }

    plain_col_sub_matrix(ref_M mm, const IT1 &it1, const IT1 &e1,
			 const IT1 &it2, const IT1 &e2)
      : begin1(it1), end1(e1), begin2(it2), end2(e2), m(&mm) {}
    plain_col_sub_matrix() {}

  };

  template <class PT, class IT1, class IT2>
  struct linalg_traits<plain_col_sub_matrix<PT, IT1, IT2> > {
    typedef plain_col_sub_matrix<PT, IT1, IT2> this_type;
    typedef typename std::iterator_traits<PT>::value_type M;
    typedef linalg_true is_reference;
    typedef abstract_matrix linalg_type;
    typedef typename linalg_traits<M>::value_type value_type;
    typedef typename linalg_traits<M>::value_type reference_type;
    typedef abstract_plain storage_type;
    typedef tab_ref_index_ref_with_origin<typename linalg_traits<typename
    linalg_traits<M>::sub_col_type>::iterator, IT1> sub_col_type;
    typedef tab_ref_index_ref_with_origin<typename linalg_traits<typename
    linalg_traits<M>::const_sub_col_type>::const_iterator, IT1>
    const_sub_col_type;
    typedef abstract_null_type sub_row_type;
    typedef abstract_null_type const_sub_row_type;
    typedef col_major sub_orientation;
    size_type nrows(const this_type &m) { return m.nrows(); }
    size_type ncols(const this_type &m) { return m.ncols(); }
    const_sub_col_type col(const this_type &m, size_type i) { 
      return const_sub_col_type(vect_begin(mat_col(m.deref(), (m.begin1)[i])),
				m.begin2, m.end2, linalg_origin(m.deref()));
    }
    sub_row_type col(this_type &m, size_type i)
    { return sub_col_type(vect_begin(mat_col(m.deref(), (m.begin1)[i])),
			  m.begin2, m.end2, linalg_origin(m.deref()));
    }
    const void* origin(const this_type &m) { return linalg_origin(m.deref()); }
    void do_clear(this_type &m)
      { for (size_type i = 0; i < ncols(); ++i) clear(col(m,i)); }
  };

  /* ******************************************************************** */
  /*		sub matrices with two array of indexes.                   */
  /* ******************************************************************** */
  
  template <class PT, class IT1, class IT2, class st_type, class orient>
  struct smrt_ir;
  
  template <class PT, class IT1, class IT2>
  struct smrt_ir<PT, IT1, IT2, abstract_plain, row_major> {
    typedef plain_row_sub_matrix<PT, IT1, IT2> matrix_type;
  };

  template <class PT, class IT1, class IT2>
  struct smrt_ir<PT, IT1, IT2, abstract_plain, col_major> {
    typedef plain_col_sub_matrix<PT, IT1, IT2> matrix_type;
  };

    template <class PT, class IT1, class IT2>
  struct smrt_ir<PT, IT1, IT2, abstract_sparse, row_major> {
    typedef sparse_row_sub_matrix<PT, IT1, IT2> matrix_type;
  };

  template <class PT, class IT1, class IT2>
  struct smrt_ir<PT, IT1, IT2, abstract_sparse, col_major> {
    typedef sparse_col_sub_matrix<PT, IT1, IT2> matrix_type;
  };

  template <class PT, class IT1, class IT2>
  struct sub_matrix_type {
    typedef typename std::iterator_traits<PT>::value_type M;
    typedef typename smrt_ir<PT, IT1, IT2,
      typename linalg_traits<M>::storage_type,
      typename principal_orientation_type<typename
    linalg_traits<M>::sub_orientation>::potype>::matrix_type matrix_type;
  };

  template <class M, class IT1, class IT2>  inline
  typename sub_matrix_type<M *, IT1, IT2>::matrix_type
  sub_matrix(M &m, const IT1 &it1, const IT1 &e1,
	     const IT2 &it2, const IT2 &e2) {
    return sub_matrix_st(m, it1, e1, it2, e2,
			 typename linalg_traits<M>::storage_type());
  }

  template <class M, class IT1, class IT2>  inline
  typename sub_matrix_type<const M *, IT1, IT2>::matrix_type
  sub_matrix(const M &m, const IT1 &it1, const IT1 &e1,
	     const IT2 &it2, const IT2 &e2) {
    return sub_matrix_stc(m, it1, e1, it2, e2,
			 typename linalg_traits<M>::storage_type());
  }

  template <class M, class IT1, class IT2> inline
  typename sub_matrix_type<const M *, IT1, IT2>::matrix_type
  sub_matrix(const M &m, const IT1 &it1, const IT1 &e1,
	     const IT2 &it2, const IT2 &e2,
	     const reverse_index &rindex1, const reverse_index &rindex2) {
    return sub_matrix_stc(m, it1, e1, it2, e2, 
			  typename linalg_traits<M>::storage_type(),
			  &rindex1, &rindex2);
  }

  template <class M, class IT1, class IT2>  inline
  typename sub_matrix_type<M *, IT1, IT2>::matrix_type
  sub_matrix(M &m, const IT1 &it1, const IT1 &e1,
	     const IT2 &it2, const IT2 &e2,
	     const reverse_index &rindex1, const reverse_index &rindex2) {
    return sub_matrix_st(m, it1, e1, it2, e2, 
			 typename linalg_traits<M>::storage_type(),
			 &rindex1, &rindex2);
  }

  template <class M, class IT1>  inline
  typename sub_matrix_type<M *, IT1, IT1>::matrix_type
  sub_matrix(M &m, const IT1 &it1, const IT1 &e1) {
    return sub_matrix_st(m, it1, e1, it1, e1,
			 typename linalg_traits<M>::storage_type());
  }

  template <class M, class IT1>  inline
  typename sub_matrix_type<const M *, IT1, IT1>::matrix_type
  sub_matrix(const M &m, const IT1 &it1, const IT1 &e1) {
    return sub_matrix_stc(m, it1, e1, it1, e1,
			 typename linalg_traits<M>::storage_type());
  }

  template <class M, class IT1> inline
  typename sub_matrix_type<const M *, IT1, IT1>::matrix_type
  sub_matrix(const M &m, const IT1 &it1, const IT1 &e1,
	     const reverse_index &rindex1) {
    return sub_matrix_stc(m, it1, e1, it1, e1, 
			  typename linalg_traits<M>::storage_type(),
			  &rindex1, &rindex1);
  }

  template <class M, class IT1>  inline
  typename sub_matrix_type<M *, IT1, IT1>::matrix_type
  sub_matrix(M &m, const IT1 &it1, const IT1 &e1,
	     const reverse_index &rindex1) {
    return sub_matrix_st(m, it1, e1, it1, e1, 
			 typename linalg_traits<M>::storage_type(),
			 &rindex1, &rindex1);
  }


  template <class M, class IT1, class IT2> inline
  typename sub_matrix_type<const M *, IT1, IT2>::matrix_type
  sub_matrix_stc(const M &m, const IT1 &it1, const IT1 &e1,
		 const IT2 &it2, const IT2 &e2, abstract_plain,
		 const reverse_index * = 0, const reverse_index * = 0) {
    return typename sub_matrix_type<const M *, IT1, IT2>::matrix_type
      (m, it1, e1, it2, e2);
  }

  template <class M, class IT1, class IT2> inline
  typename sub_matrix_type<M *, IT1, IT2>::matrix_type
  sub_matrix_st(M &m, const IT1 &it1, const IT1 &e1,
		const IT2 &it2, const IT2 &e2, abstract_plain,
		const reverse_index * = 0, const reverse_index * = 0) {
    return typename sub_matrix_type<M *, IT1, IT2>::matrix_type
      (m, it1, e1, it2, e2);
  }

  template <class M, class IT1, class IT2> inline
  typename sub_matrix_type<const M *, IT1, IT2>::matrix_type
  sub_matrix_stc(const M &m, const IT1 &it1, const IT1 &e1,
		 const IT2 &it2, const IT2 &e2, abstract_sparse,
		 const reverse_index *rindex1 = 0,
		 const reverse_index *rindex2 = 0) {
    return typename sub_matrix_type<const M *, IT1, IT2>::matrix_type
      (m, it1, e1, it2, e2, *rindex1, *rindex2);
  }

  template <class M, class IT1, class IT2> inline
  typename sub_matrix_type<M *, IT1, IT2>::matrix_type
  sub_matrix_st(M &m, const IT1 &it1, const IT1 &e1,
		const IT2 &it2, const IT2 &e2, abstract_sparse,
		const reverse_index *rindex1 = 0,
		const reverse_index *rindex2 = 0) {
    return typename sub_matrix_type<M *, IT1, IT2>::matrix_type
      (m, it1, e1, it2, e2, *rindex1, *rindex2);
  }

}

#endif //  __GMM_SUB_MATRIX_H
