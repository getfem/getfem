/* -*- c++ -*- (enables emacs c++ mode)                                    */
/* *********************************************************************** */
/*                                                                         */
/* Library :  Basic GEOmetric Tool  (bgeot)                                */
/* File    :  bgeot_abstract_linalg.h : generic algorithms on linear       */
/*                                      algebra                            */
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

// G.M.M. : Generic Matrix Methods
// T.M.M. : Template Matrix Methods
// A.M.A. : Abstract Matrix Algorithms
// T.M.A. : Template Matrix Algorithms

//
// A faire
//
//   . mult : optimisable dans certains cas.
//       (utilisation d'iterateurs, renvoie à add ou copie en
//        evitant les tests repetitifs)
//   . add : protection contre les écritures sur le même vecteur : contrôler
//           les origines ...
//   . faire scale et scaled sur les matrices aussi
//   . it.index() pour les vecteurs creux peut renvoyer -1 : ne faut-il pas
//       obliger l'iterateur à passer sur les -1 pour optimiser les cas ou
//       on a jamais de -1 ?
//

// Inspired from M.T.L. (http://www.osl.iu.edu/research/mtl)

#ifndef __BGEOT_ABSTRACT_LINALG_H
#define __BGEOT_ABSTRACT_LINALG_H

#include <dal_ref.h>
#include <bgeot_matrix.h>
#include <bgeot_smatrix.h>

namespace bgeot {

  /* ******************************************************************** */
  /*		Specifier types                             		  */
  /* ******************************************************************** */

  struct abstract_null_type {}; // specify an information lake.

  struct linalg_true {};
  struct linalg_false {};

  struct abstract_vector {};
  struct abstract_matrix {};
  
  struct abstract_sparse {};    // sparse matrix or vector
  struct abstract_plain {};     // plain matrix or vector
  struct abstract_indirect {};  // matrix given by the product with a vector

  struct row_major {};          // matrix with a row access.
  struct col_major {};       // matrix with a column access
  struct row_and_col {};     // both accesses but row preference
  struct col_and_row {};     // both accesses but column preference

  template <class T> struct transposed_type;
  template<> struct transposed_type<row_major>   {typedef col_major   t_type;};
  template<> struct transposed_type<col_major>   {typedef row_major   t_type;};
  template<> struct transposed_type<row_and_col> {typedef col_and_row t_type;};
  template<> struct transposed_type<col_and_row> {typedef row_and_col t_type;};

  template <class T> struct principal_orientation_type;
  template<> struct principal_orientation_type<row_major>
  {typedef row_major potype;};
  template<> struct principal_orientation_type<col_major>
  {typedef col_major potype;};
  template<> struct principal_orientation_type<row_and_col>
  {typedef row_major potype;};
  template<> struct principal_orientation_type<col_and_row>
  {typedef col_major potype;};
  
  template <class V> struct linalg_traits;

  /* ******************************************************************** */
  /*		Operations on scalars                         		  */
  /* ******************************************************************** */

  template <class T> T conj_product(T a, T b) { return a * b; }
  // à définir sur les complexes
  template <class T> T modulus(T a) { return dal::abs(a); }
  // à définir sur les complexes


  /* ******************************************************************** */
  /*		reverse indexes                             		  */
  /* ******************************************************************** */

  class reverse_index : public std::vector<size_t> {
  public :
    typedef std::vector<size_t> base_type;
    typedef base_type::value_type value_type;
    typedef base_type::iterator iterator;
    typedef base_type::const_iterator const_iterator;
    typedef base_type::reverse_iterator reverse_iterator;
    typedef base_type::const_reverse_iterator const_reverse_iterator;
    typedef base_type::pointer pointer;
    typedef base_type::const_pointer const_pointer;
    typedef base_type::reference reference;
    typedef base_type::const_reference const_reference;
    typedef base_type::size_type size_type;
    typedef base_type::difference_type difference_type;
    
    reverse_index(void) {}
    template <class IT> void init(IT it, IT ite, size_type n,
				  abstract_sparse = abstract_sparse());
    template <class IT> void init(IT, IT, size_type, abstract_plain) { }

    template <class IT, class L> reverse_index(IT it, IT ite,
					       size_type n, const L &)
    { init(it, ite, n, typename linalg_traits<L>::storage_type()); }
  };

  template <class IT>
  void reverse_index::init(IT it, IT ite, size_type n, abstract_sparse) {
      resize(n); std::fill(begin(), end(), size_type(-1));
      for (size_type i = 0; it != ite; ++it, ++i) (*this)[*it] = i;
  }

}

#include <bgeot_linalg_interface.h>

namespace bgeot {

  /* ******************************************************************** */
  /*		Identity matrix                         		  */
  /* ******************************************************************** */

  struct identity_matrix {};

  template <class V1, class V2> inline
  void mult(identity_matrix, const V1 &v1, V2 &v2) { copy(v1, v2); }
  template <class V1, class V2> inline
  void mult(identity_matrix, const V1 &v1, const V2 &v2) { copy(v1, v2); }
  template <class V1, class V2, class V3> inline
  void mult(identity_matrix, const V1 &v1, const V2 &v2, V3 &v3)
  { add(v1, v2, v3); }
  template <class V1, class V2, class V3> inline
  void mult(const identity_matrix&, const V1 &v1, const V2 &v2, const V3 &v3)
  { add(v1, v2, v3); }
  template <class M> void copy(const identity_matrix&, M &m) {
    size_type i = 0, n = std::max(mat_nrows(m), mat_ncols(m)); clear(m);
    for (; i < n; ++i) m(i,i) = typename linalg_traits<M>::value_type(1);
  }
  template <class V1, class V2> inline
  typename linalg_traits<V1>::value_type
  vect_sp(const identity_matrix &, const V1 &v1, const V2 &v2)
  { return vect_sp(v1, v2); }
  template <class M> inline void copy(const identity_matrix&, const M &m)
  { copy_ident(m, linalg_traits<M>::is_reference()); }
  template <class M> inline void copy_ident(const M &m, linalg_true)
  { copy(identity_matrix(), const_cast<M &>(m)); }
  template<class M> inline bool is_identity(const M&) { return false; }
  inline bool is_identity(const identity_matrix&) { return true; }

  /* ******************************************************************** */
  /*		                                         		  */
  /*		Generic algorithms                           		  */
  /*		                                         		  */
  /* ******************************************************************** */


  /* ******************************************************************** */
  /*		Miscellaneous                           		  */
  /* ******************************************************************** */

  template <class V> inline size_type vect_size(const V &v)
  { return linalg_traits<V>().size(v); }

  template <class MAT> inline size_type mat_nrows(const MAT &m)
  { return linalg_traits<MAT>().nrows(m); }

  template <class MAT> inline size_type mat_ncols(const MAT &m)
  { return linalg_traits<MAT>().ncols(m); }

  template <class L>
  inline const void *linalg_origin(const L &l)
  { return linalg_traits<L>().origin(l); }  

  template <class V>
  inline typename linalg_traits<V>::const_iterator vect_begin(const V &v)
  { return linalg_traits<V>().const_begin(v); }

  template <class V>
  inline typename linalg_traits<V>::const_iterator vect_end(const V &v)
  { return linalg_traits<V>().const_end(v); }

  template <class V>
  inline typename linalg_traits<V>::const_iterator vect_const_begin(const V &v)
  { return linalg_traits<V>().const_begin(v); }

  template <class V>
  inline typename linalg_traits<V>::const_iterator vect_const_end(const V &v)
  { return linalg_traits<V>().const_end(v); }

  template <class V>
  inline typename linalg_traits<V>::iterator vect_begin(V &v)
  { return linalg_traits<V>().begin(v); }

  template <class V>
  inline typename linalg_traits<V>::iterator vect_end(V &v)
  { return linalg_traits<V>().end(v); }

  template <class MAT> inline 
    typename linalg_traits<MAT>::const_sub_row_type
    mat_row(const MAT &m, size_type i)
  { return linalg_traits<MAT>().row(m, i); }

  template <class MAT> inline  
    typename linalg_traits<MAT>::const_sub_col_type
    mat_col(const MAT &m, size_type i)
  { return linalg_traits<MAT>().col(m, i); }

  template <class MAT> inline 
    typename linalg_traits<MAT>::const_sub_row_type
    mat_const_row(const MAT &m, size_type i)
  { return linalg_traits<MAT>().row(m, i); }

  template <class MAT> inline  
    typename linalg_traits<MAT>::const_sub_col_type
    mat_const_col(const MAT &m, size_type i)
  { return linalg_traits<MAT>().col(m, i); }

  template <class MAT> inline 
    typename linalg_traits<MAT>::sub_row_type
    mat_row(MAT &m, size_type i)
  { return linalg_traits<MAT>().row(m, i); }

  template <class MAT> inline  
    typename linalg_traits<MAT>::sub_col_type
    mat_col(MAT &m, size_type i)
  { return linalg_traits<MAT>().col(m, i); }

  template <class L> inline void clear(L &l)
  { return linalg_traits<L>().do_clear(l); }

  template <class L> inline
    scaled_vector_const_ref<L> scaled(const L &l,
				      typename linalg_traits<L>::value_type x)
  { return scaled_vector_const_ref<L>(l, x); }
  

  template <class L> inline transposed_ref<L> transposed(L &l)
  { return transposed_ref<L>(l); }

  template <class L> inline transposed_const_ref<L> transposed(const L &l)
  { return transposed_const_ref<L>(l); }

  inline bool _is_sparse(abstract_sparse)  { return true;  }
  inline bool _is_sparse(abstract_plain)   { return false; }
  inline bool _is_sparse(abstract_indirect)  { return false; }

  template <class L> inline bool is_sparse(const L &) 
  { return _is_sparse(linalg_traits<L>::storage_type()); }

  /* ******************************************************************** */
  /*   Selects a temporary vector type                                    */
  /*   V if V is a valid vector type,                                     */
  /*   svector if V is a reference on a sparse vector,                    */
  /*   std::vector if V is a reference on a plain vector.                 */
  /* ******************************************************************** */

  template <class R, class S, class V> struct _temporary_vector {};
  template <class V> struct _temporary_vector<linalg_true, abstract_sparse, V>
  { typedef svector<typename linalg_traits<V>::value_type> vector_type; };
  template <class V> struct _temporary_vector<linalg_true, abstract_plain, V>
  { typedef std::vector<typename linalg_traits<V>::value_type> vector_type; };
  template <class S, class V> struct _temporary_vector<linalg_false, S, V>
  { typedef V vector_type; };

  template <class V> struct temporary_vector {
    typedef typename _temporary_vector<typename linalg_traits<V>::is_reference,
      typename linalg_traits<V>::storage_type, V>::vector_type vector_type;
  };

  /* ******************************************************************** */
  /*   Selects a temporary plain vector type                              */
  /*   V if V is a valid plain vector type,                               */
  /*   std::vector if V is a reference or a sparse vector                 */
  /* ******************************************************************** */

  template <class R, class S, class V> struct _temporary_plain_vector;
  template <class S, class V> struct _temporary_plain_vector<linalg_true, S, V>
  { typedef std::vector<typename linalg_traits<V>::value_type> vector_type; };
  template <class V>
  struct _temporary_plain_vector<linalg_false, abstract_sparse, V>
  { typedef std::vector<typename linalg_traits<V>::value_type> vector_type; };
  template <class V>
  struct _temporary_plain_vector<linalg_false, abstract_plain, V>
  { typedef V vector_type; };

  template <class V> struct temporary_plain_vector {
    typedef typename _temporary_vector<typename linalg_traits<V>::is_reference,
      typename linalg_traits<V>::storage_type, V>::vector_type vector_type;
  };

  /* ******************************************************************** */
  /*   Selects a temporary sparse vector type                             */
  /*   V if V is a valid sparse vector type,                              */
  /*   std::vector if V is a reference or a plain vector                  */
  /* ******************************************************************** */

  template <class R, class S, class V> struct _temporary_sparse_vector;
  template <class S, class V>
  struct _temporary_sparse_vector<linalg_true, S, V>
  { typedef svector<typename linalg_traits<V>::value_type> vector_type; };
  template <class V>
  struct _temporary_sparse_vector<linalg_false, abstract_sparse, V>
  { typedef V vector_type; };
  template <class V>
  struct _temporary_sparse_vector<linalg_false, abstract_plain, V>
  { typedef svector<typename linalg_traits<V>::value_type> vector_type; };

  template <class V> struct temporary_sparse_vector {
    typedef typename _temporary_vector<typename linalg_traits<V>::is_reference,
      typename linalg_traits<V>::storage_type, V>::vector_type vector_type;
  };

  /* ******************************************************************** */
  /*		sub vector with an array of indexes.                      */
  /* ******************************************************************** */

  template <class V, class IT, class st_type> struct const_svrt_ir;

  template <class V, class IT>
  struct const_svrt_ir<V, IT, abstract_plain> {
    typedef tab_ref_index_ref_with_origin<typename linalg_traits<V>
    ::const_iterator, IT> vector_type;
  }; 

  template <class V, class IT>
  struct const_svrt_ir<V, IT, abstract_sparse> {
    typedef sparse_sub_vector<V, IT> vector_type;
  }; 

  template <class V, class IT, class st_type> struct svrt_ir;

  template <class V, class IT>
  struct svrt_ir<V, IT, abstract_plain> {
    typedef tab_ref_index_ref_with_origin<typename linalg_traits<V>
    ::iterator, IT> vector_type;
  }; 

  template <class V, class IT>
  struct svrt_ir<V, IT, abstract_sparse> {
    typedef const_sparse_sub_vector<V, IT> vector_type;
  };

  template <class V, class IT>
  struct sub_vector_type {
    typedef typename svrt_ir<V, IT,
      typename linalg_traits<V>::storage_type>::vector_type vector_type;
  };

  template <class V, class IT>
  struct const_sub_vector_type {
    typedef typename const_svrt_ir<V, IT,
      typename linalg_traits<V>::storage_type>::vector_type vector_type;
  };
  
  template <class V, class IT> inline
  typename const_sub_vector_type<V, IT>::vector_type
  sub_vector(const V &v, const IT &it, const IT &e) {
    return sub_vector_stc(v, it, e,
			  typename linalg_traits<V>::storage_type());
  }

  template <class V, class IT>  inline
  typename sub_vector_type<V, IT>::vector_type
  sub_vector(V &v, const IT &it, const IT &e) {
    return sub_vector_st(v, it, e,
			 typename linalg_traits<V>::storage_type());
  }

  template <class V, class IT> inline
  typename const_sub_vector_type<V, IT>::vector_type
  sub_vector(const V &v, const IT &it, const IT &e, 
	     const reverse_index &rindex) {
    return sub_vector_stc(v, it, e, 
			  typename linalg_traits<V>::storage_type(), &rindex);
  }

  template <class V, class IT>  inline
  typename sub_vector_type<V, IT>::vector_type
  sub_vector(V &v, const IT &it, const IT &e, const reverse_index &rindex) {
    return sub_vector_st(v, it, e, 
			 typename linalg_traits<V>::storage_type(), &rindex);
  }

  template <class V, class IT> inline
  typename const_sub_vector_type<V, IT>::vector_type
  sub_vector_stc(const V &v, const IT &it, const IT &e,
		 abstract_plain, const reverse_index * = 0) {
    return  typename const_sub_vector_type<V, IT>
      ::vector_type(vect_begin(v), it, e, linalg_origin(v));
  }

  template <class V, class IT> inline
  typename const_sub_vector_type<V, IT>::vector_type
  sub_vector_stc(const V &v, const IT &it, const IT &e,
		 abstract_sparse, const reverse_index *rindex = 0) {
    return typename const_sub_vector_type<V, IT>
      ::vector_type(v, it, e, *rindex);
  }

  template <class V, class IT> inline
  typename sub_vector_type<V, IT>::vector_type
  sub_vector_st(V &v, const IT &it, const IT &e,
		abstract_plain, const reverse_index * = 0) {
    return typename sub_vector_type<V, IT>
      ::vector_type(vect_begin(v), it, e, linalg_origin(v));
  }

  template <class V, class IT> inline
  typename sub_vector_type<V, IT>::vector_type
  sub_vector_st(V &v, const IT &it, const IT &e,
		abstract_sparse, const reverse_index *rindex = 0) {
    return typename sub_vector_type<V, IT>
      ::vector_type(v, it, e, *rindex);
  }

  /* ******************************************************************** */
  /*		sub matrices with two array of indexes.                   */
  /* ******************************************************************** */
  
  template <class M, class IT1, class IT2, class st_type, class orient>
  struct smrt_ir;
  
  template <class M, class IT1, class IT2>
  struct smrt_ir<M, IT1, IT2, abstract_plain, row_major> {
    typedef plain_row_sub_matrix<M, IT1, IT2> matrix_type;
  };

  template <class M, class IT1, class IT2>
  struct smrt_ir<M, IT1, IT2, abstract_plain, col_major> {
    typedef plain_col_sub_matrix<M, IT1, IT2> matrix_type;
  };

    template <class M, class IT1, class IT2>
  struct smrt_ir<M, IT1, IT2, abstract_sparse, row_major> {
    typedef sparse_row_sub_matrix<M, IT1, IT2> matrix_type;
  };

  template <class M, class IT1, class IT2>
  struct smrt_ir<M, IT1, IT2, abstract_sparse, col_major> {
    typedef sparse_col_sub_matrix<M, IT1, IT2> matrix_type;
  };

    template <class M, class IT1, class IT2, class st_type, class orient>
  struct const_smrt_ir;
  
  template <class M, class IT1, class IT2>
  struct const_smrt_ir<M, IT1, IT2, abstract_plain, row_major> {
    typedef const_plain_row_sub_matrix<M, IT1, IT2> matrix_type;
  };

  template <class M, class IT1, class IT2>
  struct const_smrt_ir<M, IT1, IT2, abstract_plain, col_major> {
    typedef const_plain_col_sub_matrix<M, IT1, IT2> matrix_type;
  };

    template <class M, class IT1, class IT2>
  struct const_smrt_ir<M, IT1, IT2, abstract_sparse, row_major> {
    typedef const_sparse_row_sub_matrix<M, IT1, IT2> matrix_type;
  };

  template <class M, class IT1, class IT2>
  struct const_smrt_ir<M, IT1, IT2, abstract_sparse, col_major> {
    typedef const_sparse_col_sub_matrix<M, IT1, IT2> matrix_type;
  };

  template <class M, class IT1, class IT2>
  struct sub_matrix_type {
    typedef typename smrt_ir<M, IT1, IT2,
      typename linalg_traits<M>::storage_type,
      typename principal_orientation_type<typename
    linalg_traits<M>::sub_orientation>::potype>::matrix_type matrix_type;
  };

  template <class M, class IT1, class IT2>
  struct const_sub_matrix_type {
    typedef typename const_smrt_ir<M, IT1, IT2,
      typename linalg_traits<M>::storage_type,
      typename principal_orientation_type<typename
    linalg_traits<M>::sub_orientation>::potype>::matrix_type matrix_type;
  };

  
  template <class M, class IT1, class IT2> inline
  typename const_sub_matrix_type<M, IT1, IT2>::matrix_type
  sub_matrix(const M &m, const IT1 &it1, const IT1 &e1,
	     const IT2 &it2, const IT2 &e2) {
    return sub_matrix_stc(m, it1, e1, it2, e2,
			  typename linalg_traits<M>::storage_type());
  }

  template <class M, class IT1, class IT2>  inline
  typename sub_matrix_type<M, IT1, IT2>::matrix_type
  sub_matrix(M &m, const IT1 &it1, const IT1 &e1,
	     const IT2 &it2, const IT2 &e2) {
    return sub_matrix_st(m, it1, e1, it2, e2,
			 typename linalg_traits<M>::storage_type());
  }

  template <class M, class IT1, class IT2> inline
  typename const_sub_matrix_type<M, IT1, IT2>::matrix_type
  sub_matrix(const M &m, const IT1 &it1, const IT1 &e1,
	     const IT2 &it2, const IT2 &e2,
	     const reverse_index &rindex1, const reverse_index &rindex2) {
    return sub_matrix_stc(m, it1, e1, it2, e2, 
			  typename linalg_traits<M>::storage_type(),
			  &rindex1, &rindex2);
  }

  template <class M, class IT1, class IT2>  inline
  typename sub_matrix_type<M, IT1, IT2>::matrix_type
  sub_matrix(M &m, const IT1 &it1, const IT1 &e1,
	     const IT2 &it2, const IT2 &e2,
	     const reverse_index &rindex1, const reverse_index &rindex2) {
    return sub_matrix_st(m, it1, e1, it2, e2, 
			 typename linalg_traits<M>::storage_type(),
			 &rindex1, &rindex2);
  }

  template <class M, class IT1, class IT2> inline
  typename const_sub_matrix_type<M, IT1, IT2>::matrix_type
  sub_matrix_stc(const M &m, const IT1 &it1, const IT1 &e1,
		 const IT2 &it2, const IT2 &e2, abstract_plain,
		 const reverse_index * = 0, const reverse_index * = 0) {
    return typename const_sub_matrix_type<M, IT1, IT2>::matrix_type
      (m, it1, e1, it2, e2);
  }

  template <class M, class IT1, class IT2> inline
  typename sub_matrix_type<M, IT1, IT2>::matrix_type
  sub_matrix_st(M &m, const IT1 &it1, const IT1 &e1,
		const IT2 &it2, const IT2 &e2, abstract_plain,
		const reverse_index * = 0, const reverse_index * = 0) {
    return typename sub_matrix_type<M, IT1, IT2>::matrix_type
      (m, it1, e1, it2, e2);
  }

  template <class M, class IT1, class IT2> inline
  typename const_sub_matrix_type<M, IT1, IT2>::matrix_type
  sub_matrix_stc(const M &m, const IT1 &it1, const IT1 &e1,
		 const IT2 &it2, const IT2 &e2, abstract_sparse,
		 const reverse_index *rindex1 = 0,
		 const reverse_index *rindex2 = 0) {
    return typename const_sub_matrix_type<M, IT1, IT2>::matrix_type
      (m, it1, e1, it2, e2, *rindex1, *rindex2);
  }

  template <class M, class IT1, class IT2> inline
  typename sub_matrix_type<M, IT1, IT2>::matrix_type
  sub_matrix_st(M &m, const IT1 &it1, const IT1 &e1,
		const IT2 &it2, const IT2 &e2, abstract_sparse,
		const reverse_index *rindex1 = 0,
		const reverse_index *rindex2 = 0) {
    return typename sub_matrix_type<M, IT1, IT2>::matrix_type
      (m, it1, e1, it2, e2, *rindex1, *rindex2);
  }

  /* ******************************************************************** */
  /*		Scalar product                             		  */
  /* ******************************************************************** */

  template <class V1, class V2>
    typename linalg_traits<V1>::value_type
    vect_sp(const V1 &v1, const V2 &v2) {
    if (vect_size(v1) != vect_size(v2))
      DAL_THROW(dimension_error,"dimensions mismatch");
    return vect_sp(v1, v2,
		   typename linalg_traits<V1>::storage_type(), 
		   typename linalg_traits<V2>::storage_type());
  }

  template <class MATSP, class V1, class V2> inline
    typename linalg_traits<V1>::value_type
    vect_sp(const MATSP &ps, const V1 &v1, const V2 &v2) {
    return vect_sp_with_mat(ps, v1, v2,
			    typename linalg_traits<MATSP>::sub_orientation());
  }

  template <class MATSP, class V1, class V2> inline
    typename linalg_traits<V1>::value_type
    vect_sp_with_mat(const MATSP &ps, const V1 &v1, const V2 &v2, row_major) {
    return vect_sp_with_matr(ps, v1, v2, 
			     typename linalg_traits<V2>::storage_type());
  }

  template <class MATSP, class V1, class V2>
    typename linalg_traits<V1>::value_type
    vect_sp_with_matr(const MATSP &ps, const V1 &v1, const V2 &v2,
		      abstract_sparse) {
    if (vect_size(v1) != mat_ncols() || vect_size(v2) != mat_nrows())
      DAL_THROW(dimension_error,"dimensions mismatch");
    size_type nr = mat_nrows(ps);
    typename linalg_traits<V2>::const_iterator
      it = vect_begin(v2), ite = vect_end(v2);
    typename linalg_traits<V1>::value_type res(0);
    for (; it != ite; ++it)
      if (it.index() != size_type(-1))
	res += conj_product(vect_sp(mat_row(ps, it.index()), v1), *it);
    return res;
  }

  template <class MATSP, class V1, class V2>
    typename linalg_traits<V1>::value_type
    vect_sp_with_matr(const MATSP &ps, const V1 &v1, const V2 &v2,
		      abstract_plain) {
    if (vect_size(v1) != mat_ncols() || vect_size(v2) != mat_nrows())
      DAL_THROW(dimension_error,"dimensions mismatch");
    size_type nr = mat_nrows(ps);
    typename linalg_traits<V2>::const_iterator
      it = vect_begin(v2), ite = vect_end(v2);
    typename linalg_traits<V1>::value_type res(0);
    for (size_type i = 0; it != ite; ++i, ++it)
      res += conj_product(vect_sp(mat_row(ps, i), v1), *it);
    return res;
  }

  template <class MATSP, class V1, class V2> inline
  typename linalg_traits<V1>::value_type
  vect_sp_with_mat(const MATSP &ps, const V1 &v1,const V2 &v2,row_and_col)
  { return vect_sp_with_mat(ps, v1, v2, row_major()); }

  template <class MATSP, class V1, class V2> inline
  typename linalg_traits<V1>::value_type
  vect_sp_with_mat(const MATSP &ps, const V1 &v1, const V2 &v2,col_major){
    return vect_sp_with_matc(ps, v1, v2,
			     typename linalg_traits<V1>::storage_type());
  }

  template <class MATSP, class V1, class V2>
    typename linalg_traits<V1>::value_type
    vect_sp_with_matc(const MATSP &ps, const V1 &v1, const V2 &v2,
		      abstract_sparse) {
    if (vect_size(v1) != mat_ncols() || vect_size(v2) != mat_nrows())
      DAL_THROW(dimension_error,"dimensions mismatch");
    typename linalg_traits<V1>::const_iterator
      it = vect_begin(v1), ite = vect_end(v1);
    typename linalg_traits<V1>::value_type res(0);
    for (; it != ite; ++it)
      if (it.index() != size_type(-1))
	res += conj_product(vect_sp(mat_col(ps, it.index()), v2), *it);
    return res;
  }

  template <class MATSP, class V1, class V2>
    typename linalg_traits<V1>::value_type
    vect_sp_with_matc(const MATSP &ps, const V1 &v1, const V2 &v2,
		      abstract_plain) {
    if (vect_size(v1) != mat_ncols() || vect_size(v2) != mat_nrows())
      DAL_THROW(dimension_error,"dimensions mismatch");
    typename linalg_traits<V1>::const_iterator
      it = vect_begin(v1), ite = vect_end(v1);
    typename linalg_traits<V1>::value_type res(0);
    for (size_type i = 0; it != ite; ++i, ++it)
      res += conj_product(vect_sp(mat_col(ps, i), v2), *it);
    return res;
  }

  template <class MATSP, class V1, class V2> inline
  typename linalg_traits<V1>::value_type
  vect_sp_with_mat(const MATSP &ps, const V1 &v1,const V2 &v2,col_and_row)
  { return vect_sp_with_mat(ps, v1, v2, column_major()); }

  template <class MATSP, class V1, class V2> inline
  typename linalg_traits<V1>::value_type
  vect_sp_with_mat(const MATSP &ps, const V1 &v1, const V2 &v2,
		   abstract_null_type) {
    V1 w(mat_nrows());
    cerr << "Warning, a temporary is used in scalar product\n";
    mult(ps, v1, w); 
    return vect_sp(w, v2);
  }

  template <class IT1, class IT2> inline
    typename std::iterator_traits<IT1>::value_type
    _vect_sp_plain(IT1 it, IT1 ite, IT2 it2) {
    typename std::iterator_traits<IT1>::value_type res(0);
    for (; it != ite; ++it, ++it2) res += conj_product(*it, *it2);
    return res;
  }

  template <class IT1, class V> inline
    typename std::iterator_traits<IT1>::value_type
    _vect_sp_sparse(IT1 it, IT1 ite, const V &v) {
    // cout << "ici aussi\n";
    typename std::iterator_traits<IT1>::value_type res(0);
    for (; it != ite; ++it) { 
      // cout << "index = " << it.index() << endl;
      if (it.index() != size_type(-1)) res += conj_product(*it, v[it.index()]); }
    // cout << "stop\n";
    return res;
  }

  template <class V1, class V2> inline
    typename linalg_traits<V1>::value_type
    vect_sp(const V1 &v1, const V2 &v2, abstract_plain, abstract_plain) {
    return _vect_sp_plain(vect_begin(v1), vect_end(v1), vect_begin(v2));
  }

  template <class V1, class V2> inline
    typename linalg_traits<V1>::value_type
    vect_sp(const V1 &v1, const V2 &v2,abstract_sparse,abstract_plain) {
    return _vect_sp_sparse(vect_begin(v1), vect_end(v1), v2);
  }

  template <class V1, class V2> inline
    typename linalg_traits<V1>::value_type
    vect_sp(const V1 &v1, const V2 &v2, abstract_plain,abstract_sparse) {
    return _vect_sp_sparse(vect_begin(v2), vect_end(v2), v1);
  }

  template <class V1, class V2> inline
    typename linalg_traits<V1>::value_type
    vect_sp(const V1 &v1, const V2 &v2,abstract_sparse,abstract_sparse) {
    return _vect_sp_sparse(vect_begin(v1), vect_end(v1), v2);
  }

  /* ******************************************************************** */
  /*		Euclidian norm                             		  */
  /* ******************************************************************** */

   template <class V>
    typename linalg_traits<V>::value_type norm2(const V &v) {
    typename linalg_traits<V>::const_iterator
      it = vect_begin(v), ite = vect_end(v);
    typename linalg_traits<V>::value_type res(0);
    for (; it != ite; ++it) res += dal::sqr(modulus(*it));
    return sqrt(res);
  }

  /* ******************************************************************** */
  /*		Inifity norm                              		  */
  /* ******************************************************************** */

  template <class V>
    typename linalg_traits<V>::value_type norminf(const V &v) {
    typename linalg_traits<V>::const_iterator
      it = vect_begin(v), ite = vect_end(v);
    typename linalg_traits<V>::value_type res(0);
    for (; it != ite; ++it) res = std::max(res, modulus(*it));
    return res;
  }
  
  /* ******************************************************************** */
  /*		norm1                                    		  */
  /* ******************************************************************** */

  template <class V>
    typename linalg_traits<V>::value_type norm1(const V &v) {
    typename linalg_traits<V>::const_iterator
      it = vect_begin(v), ite = vect_end(v);
    typename linalg_traits<V>::value_type res(0);
    for (; it != ite; ++it) res += modulus(*it);
    return res;
  }

  /* ******************************************************************** */
  /*		Copy                                    		  */
  /* ******************************************************************** */

  template <class L1, class L2> inline
  void copy(const L1& l1, L2& l2) { 
    if ((const void *)(&l1) != (const void *)(&l2)) {
      #ifdef __GETFEM_VERIFY
        if (linalg_origin(l1) == linalg_origin(l2))
	  cerr << "Warning : a conflict is possible in vector copy\n";
      #endif
      copy(l1, l2, typename linalg_traits<L1>::linalg_type(),
	   typename linalg_traits<L2>::linalg_type());
    }
  }

  template <class L1, class L2> inline
  void copy(const L1& l1, const L2& l2) {
    copy_ref(l1, l2, typename linalg_traits<L2>::is_reference());
  }

  template <class L1, class L2> inline
  void copy_ref(const L1& l1, const L2& l2, linalg_true)
  { copy(l1, const_cast<L2&>(l2)); }

  template <class L1, class L2>
  void copy(const L1& l1, L2& l2, abstract_vector, abstract_vector) {
    if (vect_size(l1) != vect_size(l2))
      DAL_THROW(dimension_error,"dimensions mismatch");
    copy_vect(l1, l2, typename linalg_traits<L1>::storage_type(),
	      typename linalg_traits<L2>::storage_type());
  }

  template <class L1, class L2>
  void copy(const L1& l1, L2& l2, abstract_matrix, abstract_matrix) {
    if (mat_ncols(l1) != mat_ncols(l2) || mat_nrows(l1) != mat_nrows(l2))
      DAL_THROW(dimension_error,"dimensions mismatch");
    copy_mat(l1, l2, typename linalg_traits<L1>::sub_orientation(),
	     typename linalg_traits<L2>::sub_orientation());
  }

  template <class V1, class V2, class C1, class C2>
  void copy_vect(const V1 &v1, const V2 &v2, C1, C2)
  { copy_vect(v1, const_cast<V2 &>(v2), C1(), C2()); }
  

  template <class L1, class L2>
  void copy_mat_by_row(const L1& l1, L2& l2)
  {
    size_type nbr = mat_nrows(l1);
    for (size_type i = 0; i < nbr; ++i)
      copy_vect(mat_row(l1, i), mat_row(l2, i),
      		typename linalg_traits<L1>::storage_type(),
		typename linalg_traits<L2>::storage_type());
  }

  template <class L1, class L2>
  void copy_mat_by_col(const L1 &l1, L2 &l2) {
    size_type nbc = mat_ncols(l1);
    for (size_type i = 0; i < nbc; ++i)
      copy_vect(mat_col(l1, i), mat_col(l2, i),
      		typename linalg_traits<L1>::storage_type(),
		typename linalg_traits<L2>::storage_type());
  }

  template <class L1, class L2> inline
  void copy_mat(const L1& l1, L2& l2, row_major, row_major)
  { copy_mat_by_row(l1, l2); }

  template <class L1, class L2> inline
  void copy_mat(const L1& l1, L2& l2, row_major, row_and_col)
  { copy_mat_by_row(l1, l2); }

  template <class L1, class L2> inline
  void copy_mat(const L1& l1, L2& l2, row_and_col, row_and_col)
  { copy_mat_by_row(l1, l2); }

  template <class L1, class L2> inline
  void copy_mat(const L1& l1, L2& l2, row_and_col, row_major)
  { copy_mat_by_row(l1, l2); }

  template <class L1, class L2> inline
  void copy_mat(const L1& l1, L2& l2, col_and_row, row_major)
  { copy_mat_by_row(l1, l2); }

  template <class L1, class L2> inline
  void copy_mat(const L1& l1, L2& l2, row_major, col_and_row)
  { copy_mat_by_row(l1, l2); }

  template <class L1, class L2> inline
  void copy_mat(const L1& l1, L2& l2, col_and_row, row_and_col)
  { copy_mat_by_row(l1, l2); }

  template <class L1, class L2> inline
  void copy_mat(const L1& l1, L2& l2, row_and_col, col_and_row)
  { copy_mat_by_row(l1, l2); }

  template <class L1, class L2> inline
  void copy_mat(const L1& l1, L2& l2, col_major, col_major)
  { copy_mat_by_col(l1, l2); }

  template <class L1, class L2> inline
  void copy_mat(const L1& l1, L2& l2, col_major, col_and_row)
  { copy_mat_by_col(l1, l2); }

  template <class L1, class L2> inline
  void copy_mat(const L1& l1, L2& l2, col_major, row_and_col)
  { copy_mat_by_col(l1, l2); }

  template <class L1, class L2> inline
  void copy_mat(const L1& l1, L2& l2, row_and_col, col_major)
  { copy_mat_by_col(l1, l2); }

  template <class L1, class L2> inline
  void copy_mat(const L1& l1, L2& l2, col_and_row, col_major)
  { copy_mat_by_col(l1, l2); }

  template <class L1, class L2> inline
  void copy_mat(const L1& l1, L2& l2, col_and_row, col_and_row)
  { copy_mat_by_col(l1, l2); }
  
  template <class L1, class L2> inline
  void copy_mat_mixed_rc(const L1& l1, L2& l2, size_type i) {
    copy_mat_mixed_rc(l1, l2, i, typename linalg_traits<L1>::storage_type());
  }

  template <class L1, class L2>
  void copy_mat_mixed_rc(const L1& l1, L2& l2, size_type i, abstract_sparse) {
    typename linalg_traits<L1>::const_iterator
      it  = vect_begin(l1), ite = vect_end(l1);
    for (; it != ite; ++it)
      if (it.index() != size_type(-1)) l2(i, it.index()) = *it;
  }

  template <class L1, class L2>
  void copy_mat_mixed_rc(const L1& l1, L2& l2, size_type i, abstract_plain) {
    typename linalg_traits<L1>::const_iterator
      it  = vect_begin(l1), ite = vect_end(l1);
    for (size_type j = 0; it != ite; ++it, ++j) l2(i, j) = *it;
  }

  template <class L1, class L2> inline
  void copy_mat_mixed_cr(const L1& l1, L2& l2, size_type i) {
    copy_mat_mixed_cr(l1, l2, i, typename linalg_traits<L1>::storage_type());
  }

  template <class L1, class L2>
  void copy_mat_mixed_cr(const L1& l1, L2& l2, size_type i, abstract_sparse) {
    typename linalg_traits<L1>::const_iterator
      it  = vect_begin(l1), ite = vect_end(l1);
    for (; it != ite; ++it)
      if (it.index() != size_type(-1)) l2(it.index(), i) = *it;
  }

  template <class L1, class L2>
  void copy_mat_mixed_cr(const L1& l1, L2& l2, size_type i, abstract_plain) {
    typename linalg_traits<L1>::const_iterator
      it  = vect_begin(l1), ite = vect_end(l1);
    for (size_type j = 0; it != ite; ++it, ++j) l2(j, i) = *it;
  }

  template <class L1, class L2>
  void copy_mat(const L1& l1, L2& l2, row_major, col_major) {
    clear(l2);
    size_type nbr = mat_nrows(l1);
    for (size_type i = 0; i < nbr; ++i)
      copy_mat_mixed_rc(mat_row(l1, i), l2, i);
  }
  
  template <class L1, class L2>
  void copy_mat(const L1& l1, L2& l2, col_major, row_major) {
    clear(l2);
    size_type nbc = mat_ncols(l1);
    for (size_type i = 0; i < nbc; ++i)
      copy_mat_mixed_cr(mat_col(l1, i), l2, i);
  }
  
  template <class L1, class L2> inline
  void copy_vect(const L1 &l1, L2 &l2, abstract_plain, abstract_plain) {
    std::copy(vect_begin(l1), vect_end(l1), vect_begin(l2));
  }
  
  template <class L1, class L2>
  void copy_vect(const L1& l1, L2& l2,
		 abstract_sparse, abstract_plain) {
    clear(l2);
    typename linalg_traits<L1>::const_iterator
      it  = vect_begin(l1), ite = vect_end(l1);
    for (; it != ite; ++it)
      if (it.index() != size_type(-1)) l2[it.index()] = *it;
  }
  
  template <class L1, class L2>
  void copy_vect(const L1& l1, L2& l2,
		 abstract_sparse, abstract_sparse) {
    clear(l2);
    typename linalg_traits<L1>::const_iterator
      it  = vect_begin(l1), ite = vect_end(l1);
    for (; it != ite; ++it)
      if (*it != (typename linalg_traits<L1>::value_type)(0)
	&& it.index() != size_type(-1)) l2[it.index()] = *it;
  }
  
  template <class L1, class L2>
  void copy_vect(const L1& l1, L2& l2,
		 abstract_plain, abstract_sparse) {
    clear(l2);
    typename linalg_traits<L1>::const_iterator
      it  = vect_begin(l1), ite = vect_end(l1);
    for (size_type i = 0; it != ite; ++it, ++i)
      if (*it != (typename linalg_traits<L1>::value_type)(0))
	l2[i] = *it;
  }

  /* ******************************************************************** */
  /*		Vector Addition                                    	  */
  /*   algorithms are build in order to avoid some conflicts whith        */
  /*   repeated arguments or with overlapping part of a same object.      */
  /*   In the latter case, conflicts are still possible.                  */
  /* ******************************************************************** */
  
  template <class L1, class L2> inline
    void add(const L1& l1, L2& l2) {
      add_spec(l1, l2, typename linalg_traits<L2>::linalg_type());
  }

  template <class L1, class L2> inline
  void add(const L1& l1, const L2& l2) {
    add_ref(l1, l2, typename linalg_traits<L2>::is_reference());
  }

  template <class L1, class L2> inline
  void add_ref(const L1& l1, const L2& l2, linalg_true)
  { add(l1, const_cast<L2 &>(l2)); }

  template <class L1, class L2>
    void add_spec(const L1& l1, L2& l2, abstract_vector) {
    if (vect_size(l1) != vect_size(l2))
      DAL_THROW(dimension_error,"dimensions mismatch");
    if (mat_ncols(l1) != 1 || mat_ncols(l2) != 1)
      DAL_THROW(to_be_done_error,"to be done.");
    add(l1, l2, typename linalg_traits<L1>::storage_type(),
	typename linalg_traits<L2>::storage_type());
  }

  template <class L1, class L2, class L3>
    void add(const L1& l1, const L2& l2, L3& l3) {
    if (vect_size(l1) != vect_size(l2) || vect_size(l1) != vect_size(l3))
      DAL_THROW(dimension_error,"dimensions mismatch"); 
    if (mat_ncols(l1) != 1 || mat_ncols(l2) != 1 || mat_ncols(l3) != 1)
      DAL_THROW(to_be_done_error,"to be done.");
    if ((const void *)(&l1) == (const void *)(&l3))
      add(l2, l3);
    else if ((const void *)(&l2) == (const void *)(&l3))
      add(l1, l3);
    else
      add(l1, l2, l3, typename linalg_traits<L1>::storage_type(),
	  typename linalg_traits<L2>::storage_type(),
	  typename linalg_traits<L3>::storage_type());
  }

  template <class L1, class L2, class L3> inline
  void add(const L1& l1, const L2& l2, const L3& l3) {
    add_ref(l1, l2, l3, typename linalg_traits<L3>::is_reference());
  }

  template <class L1, class L2, class L3> inline
  void add_ref(const L1& l1, const L2& l2, const L3& l3, linalg_true)
  { add(l1, l2, const_cast<L3 &>(l3)); }

  template <class IT1, class IT2, class IT3> inline
    void _add_full(IT1 it1, IT2 it2, IT3 it3, IT3 ite) {
    for (; it3 != ite; ++it3, ++it2, ++it1) *it3 = *it1 + *it2;
  }

  template <class IT1, class IT2, class IT3>
    void _add_almost_full(IT1 it1, IT1 ite1, IT2 it2, IT3 it3, IT3 ite3) {
    IT3 it = it3;
    for (; it != ite3; ++it, ++it2) *it = *it2;
    for (; it1 != ite1; ++it1)
      if (it.index() != size_type(-1)) *(it3 + it1.index()) += *it1;
  }

  template <class IT1, class IT2, class IT3> inline
  void _add_to_full(IT1 it1, IT1 ite1, IT2 it2, IT2 ite2,
		    IT3 it3, IT3 ite3) {
    IT3 it = it3;
    for (; it != ite3; ++it) *it = 0;
    for (; it1 != ite1; ++it1)
      if (it.index() != size_type(-1)) *(it3 + it1.index()) = *it1;
    for (; it2 != ite2; ++it2)
      if (it.index() != size_type(-1)) *(it3 + it2.index()) += *it2;    
  }
  
  template <class L1, class L2, class L3> inline
  void add(const L1& l1, const L2& l2, L3& l3,
	   abstract_plain, abstract_plain, abstract_plain) {
    _add_full(vect_begin(l1), vect_begin(l2),
	      vect_begin(l3), vect_end(l3));
  }
  
  template <class L1, class L2, class L3> inline
  void add(const L1& l1, const L2& l2, L3& l3,
	   abstract_sparse, abstract_plain, abstract_plain) {
    _add_almost_full(vect_begin(l1), vect_end(l1), vect_begin(l2),
		     vect_begin(l3), vect_end(l3));
  }
  
  template <class L1, class L2, class L3> inline
  void add(const L1& l1, const L2& l2, L3& l3,
	   abstract_plain, abstract_sparse, abstract_plain) {
    _add_almost_full(vect_begin(l2), vect_end(l2), vect_begin(l1),
		     vect_begin(l3), vect_end(l3));
  }
  
  template <class L1, class L2, class L3> inline
  void add(const L1& l1, const L2& l2, L3& l3,
	   abstract_sparse, abstract_sparse, abstract_plain) {
    _add_to_full(vect_begin(l1), vect_end(l1),
		 vect_begin(l2), vect_end(l2),
		 vect_begin(l3), vect_end(l3));
  }
  
  template <class L1, class L2, class L3> inline
  void add(const L1& l1, const L2& l2, L3& l3,
	   abstract_sparse, abstract_sparse, abstract_sparse) {
    typename linalg_traits<L1>::const_iterator
      it1 = vect_begin(l1), ite1 = vect_end(l1);
    typename linalg_traits<L2>::const_iterator
      it2 = vect_begin(l2), ite2 = vect_end(l2);
    clear(l3);
    while (it1 != ite1 && it2 != ite2) {
      while (it1.index() == size_type(-1) && it1 != ite1) ++it1;
      while (it2.index() == size_type(-1) && it2 != ite2) ++it2;
      ptrdiff_t d = it1.index() - it2.index();
      if (d < 0)
	{ l3[it1.index()] += *it1; ++it1; }
      else if (d > 0)
	{ l3[it2.index()] += *it2; ++it2; }
      else
	{ l3[it1.index()] = *it1 + *it2; ++it1; ++it2; }
    }
    for (; it1 != ite1; ++it1) l3[it1.index()] += *it1;
    for (; it2 != ite2; ++it2) l3[it2.index()] += *it2;   
  }
  
  template <class L1, class L2, class L3> inline
  void add(const L1&, const L2&, L3&,
	   abstract_plain, abstract_sparse, abstract_sparse)
  { DAL_THROW(failure_error,"Unauthorized addition"); }
  template <class L1, class L2, class L3> inline
  void add(const L1&, const L2&, L3&,
	   abstract_plain, abstract_plain, abstract_sparse)
  { DAL_THROW(failure_error,"Unauthorized addition"); }
  template <class L1, class L2, class L3> inline
  void add(const L1&, const L2&, L3&,
	   abstract_sparse, abstract_plain, abstract_sparse)
  { DAL_THROW(failure_error,"Unauthorized addition"); }
  
  template <class L1, class L2> inline
  void add(const L1& l1, L2& l2,
	   abstract_plain, abstract_plain) {
    typename linalg_traits<L1>::const_iterator
      it1 = vect_begin(l1); 
    typename linalg_traits<L2>::iterator
      it2 = vect_begin(l2), ite = vect_end(l2);
    for (; it2 != ite; ++it2, ++it1) *it2 += *it1;
  }
  
  template <class L1, class L2> inline
  void add(const L1& l1, L2& l2,
	   abstract_sparse, abstract_plain) {
    typename linalg_traits<L1>::const_iterator
      it1 = vect_begin(l1), ite1 = vect_end(l1);
    for (; it1 != ite1; ++it1) 
      if (it.index() != size_type(-1)) l2[it1.index()] += *it1;
  }
  
  template <class L1, class L2> inline
  void add(const L1& l1, L2& l2,
	   abstract_sparse, abstract_sparse) {
    typename linalg_traits<L1>::const_iterator
      it1 = vect_begin(l1), ite1 = vect_end(l1);
    for (; it1 != ite1; ++it1) 
      if (it.index() != size_type(-1)) l2[it1.index()] += *it1;
  }
  
  template <class L1, class L2> inline
  void add(const L1&, L2&, abstract_plain, abstract_sparse)
  { DAL_THROW(failure_error,"Unauthorized addition"); } 

  /* ******************************************************************** */
  /*		scale                                    	          */
  /* ******************************************************************** */

  template <class L> void scale(L& l, typename linalg_traits<L>::value_type a)
  {
    typename linalg_traits<L>::iterator it = vect_begin(l), ite = vect_end(l);
    for ( ; it != ite; ++it) *it *= a;
  }

  template <class L> inline
  void scale(const L& l, typename linalg_traits<L>::value_type a)
  { scale_const(l, a, typename linalg_traits<L>::is_reference()); }

  template <class L> inline
  void scale_const(const L& l, linalg_true)
  { scale_const(const_cast<L &>(l), a); }


  /* ******************************************************************** */
  /*		Matrix-vector mult                                    	  */
  /* ******************************************************************** */

  template <class L1, class L2, class L3>
  void mult(const L1& l1, const L2& l2, L3& l3) {
    if (mat_ncols(l1) != vect_size(l2) || mat_nrows(l1) != vect_size(l3))
      DAL_THROW(dimension_error,"dimensions mismatch");
    if (linalg_origin(l2) != linalg_origin(l3))
      mult_spec(l1, l2, l3, typename principal_orientation_type<typename
		linalg_traits<L1>::sub_orientation>::potype());
    else {
      #ifdef __GETFEM_VERIFY
        cerr << "Warning, A temporary is used for mult\n";
      #endif
      L3 temp(vect_size(l3));
      mult_spec(l1, l2, temp, typename principal_orientation_type<typename
		linalg_traits<L1>::sub_orientation>::potype());
      copy(temp, l3);
    }
  }

  template <class L1, class L2, class L3> inline
  void mult(const L1& l1, const L2& l2, const L3& l3) {
    mult_const(l1, l2, l3, typename linalg_traits<L3>::is_reference());
  }

  template <class L1, class L2, class L3> inline
  void mult_const(const L1& l1, const L2& l2, const L3& l3, linalg_true)
  { mult(l1, l2, const_cast<L3 &>(l3)); }

  template <class L1, class L2, class L3> inline
  void mult_by_row(const L1& l1, const L2& l2, L3& l3, abstract_sparse) {
    clear(l3);
    size_type nr = mat_nrows(l1);
    for (size_type i = 0; i < nr; ++i) {
      typename linalg_traits<L1>::value_type aux = vect_sp(mat_row(l1, i), l2);
      if (aux != 0) l3[i] = aux;
    }
  }

  template <class L1, class L2, class L3> inline
  void mult_by_row(const L1& l1, const L2& l2, L3& l3, abstract_plain) {
    typename linalg_traits<L3>::iterator
      it = vect_begin(l3), ite = vect_end(l3);
    for (size_type i = 0; it != ite; ++it, ++i)
      *it = vect_sp(mat_row(l1, i), l2);
  }

  template <class L1, class L2, class L3> inline
  void mult_spec(const L1& l1, const L2& l2, L3& l3, row_major)
  { mult_by_row(l1, l2, l3, typename linalg_traits<L3>::storage_type()); }

  template <class L1, class L2, class L3> inline
  void mult_spec(const L1& l1, const L2& l2, L3& l3, col_major) {
    clear(l3);
    size_type nc = mat_ncols(l1);
    for (size_type i = 0; i < nc; ++i)
      add(scaled(mat_col(l1, i), l2[i]), l3);
  }

  template <class L1, class L2, class L3> inline
  void mult_spec(const L1& l1, const L2& l2, L3& l3, abstract_null_type)
  { mult_ind(l1, l2, l3, typename linalg_traits<L1>::storage_type()); }

  template <class L1, class L2, class L3> inline
  void mult_ind(const L1& l1, const L2& l2, L3& l3, abstract_indirect) {
    DAL_THROW(failure_error,
	  "You have to define the mult(m, v1, v2) for this kind of matrix");
  }

  template <class L1, class L2, class L3, class L4>
  void mult(const L1& l1, const L2& l2, const L3& l3, L4& l4) {
    if (mat_ncols(l1) != vect_size(l2) || mat_nrows(l1) != vect_size(l3)
	|| mat_nrows(l1) != vect_size(l4))
      DAL_THROW(dimension_error,"dimensions mismatch");
    if (linalg_origin(l2) != linalg_origin(l4))
      mult_spec(l1, l2, l3, l4, typename principal_orientation_type<typename linalg_traits<L1>::sub_orientation>::potype());
    else {
      #ifdef __GETFEM_VERIFY
        cerr << "Warning, A temporary is used for mult\n";
      #endif
      L3 temp(vect_size(l3));
      mult_spec(l1,l2,l3, temp, typename principal_orientation_type<typename
		linalg_traits<L1>::sub_orientation>::potype());
      copy(temp, l4);
    }
  }
  
  template <class L1, class L2, class L3, class L4> inline
  void mult(const L1& l1, const L2& l2, const L3& l3, const L4& l4) {
    mult_const(l1, l2, l3, l4, typename linalg_traits<L3>::is_reference());
  }

  template <class L1, class L2, class L3, class L4> inline
  void mult_const(const L1& l1, const L2& l2, const L3& l3,
		  const L4& l4, linalg_true)
  { mult(l1, l2, l3, const_cast<L4 &>(l4)); }

  template <class L1, class L2, class L3, class L4> inline
  void mult_by_row(const L1& l1, const L2& l2, const L3& l3,
		   L4& l4, abstract_sparse) {
    copy(l3, l4);
    size_type nr = mat_nrows(l1);
    for (size_type i = 0; i < nr; ++i) {
      typename linalg_traits<L1>::value_type
	aux = vect_sp(mat_row(l1, i), l2);
      if (aux != 0) l4[i] += aux;
    }
  }

  template <class L1, class L2, class L3, class L4> inline
  void mult_by_row(const L1& l1, const L2& l2, const L3& l3, L4& l4,
		   abstract_plain) {
    copy(l3, l4); 
    typename linalg_traits<L4>::iterator
      it = vect_begin(l4), ite = vect_end(l4);
    for (size_type i = 0; it != ite; ++it, ++i)
      *it += vect_sp(mat_row(l1, i), l2);
  }

  template <class L1, class L2, class L3, class L4> inline
  void mult_spec(const L1& l1, const L2& l2, const L3& l3, L4& l4, row_major)
  { mult_by_row(l1, l2, l3, l4, typename linalg_traits<L4>::storage_type()); }

  template <class L1, class L2, class L3, class L4> inline
  void mult_spec(const L1& l1, const L2& l2, const L3& l3, L4& l4, col_major) {
    copy(l3, l4);
    size_type nc = mat_ncols(l1);
    for (size_type i = 0; i < nc; ++i)
      add(scaled(mat_col(l1, i), l2[i]), l4);
  }

  template <class L1, class L2, class L3, class L4> inline
  void mult_spec(const L1& l1, const L2& l2, const L3& l3,
		 L4& l4, abstract_null_type)
  { mult_ind(l1, l2, l3, l4, typename linalg_traits<L1>::storage_type()); }

  template <class L1, class L2, class L3, class L4> inline
  void mult_ind(const L1& l1, const L2& l2, const L3& l3,
		L4& l4, abstract_indirect) {
    DAL_THROW(failure_error,
	  "You have to define the mult(m, v1, v2) for this kind of matrix");
  }

}


#endif //  __BGEOT_ABSTRACT_LINALG_H
