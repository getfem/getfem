/* -*- c++ -*- (enables emacs c++ mode)                                    */
/* *********************************************************************** */
/*                                                                         */
/* Library :  Generic Matrix Methods  (gmm)                                */
/* File    :  gmm_blas.h : generic basic linear  algebra algorithms.       */
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

#ifndef __GMM_BLAS_H
#define __GMM_BLAS_H

namespace gmm {

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

  // à spécifier suivant le type de la référence

//   template <class V> inline
//   typename select_return<typename linalg_traits<V>::const_iterator,
//                          typename linalg_traits<V>::iterator,
//                          const V *>::return_type vect_begin(const V &v)
//   { return linalg_traits<V>().begin(linalg_cast(v)); }

//   template <class V> inline
//   typename select_return<typename linalg_traits<V>::const_iterator,
//                          typename linalg_traits<V>::iterator,
//                          V *>::return_type vect_begin(V &v)
//   { return linalg_traits<V>().begin(linalg_cast(v)); }



  // le reste est inutile normalement y compris les vect_const_begin ...


  template <class V>
  inline typename linalg_traits<V>::const_iterator vect_begin(const V &v)
  { return linalg_traits<V>().const_begin(v); }

  template <class V>
  inline typename linalg_traits<V>::const_iterator vect_end(const V &v)
  { return linalg_traits<V>().const_end(v); }

  template <class V>
  inline typename linalg_traits<V>::iterator vect_begin(V &v)
  { return linalg_traits<V>().begin(v); }

  template <class V>
  inline typename linalg_traits<V>::iterator vect_end(V &v)
  { return linalg_traits<V>().end(v); }


  template <class V>
  inline typename linalg_traits<V>::const_iterator vect_const_begin(const V &v)
  { return linalg_traits<V>().const_begin(v); }

  template <class V>
  inline typename linalg_traits<V>::const_iterator vect_const_end(const V &v)
  { return linalg_traits<V>().const_end(v); }
// à spécifier suivant le type de la référence
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

  template <class L> inline void clear(const L &l)
  { return linalg_traits<L>().do_clear(linalg_cast(l)); }

  template <class L> inline
    scaled_vector_const_ref<L> scaled(const L &l,
				      typename linalg_traits<L>::value_type x)
  { return scaled_vector_const_ref<L>(l, x); }
  
// à spécifier suivant le type de la référence
  template <class TYPE, class PT> struct transposed_return;
  template <class PT> struct transposed_return<abstract_vector, PT>
  { typedef vect_transposed_ref<PT> ret_type; };
  template <class PT> struct transposed_return<abstract_matrix, PT>
  { typedef transposed_ref<PT> ret_type; };

  template <class L> inline 
  transposed_ref<L *> transposed(L &l, abstract_matrix)
  { return transposed_ref<L *>(l); }

  template <class L> inline 
  transposed_ref<const L *> transposed(const L &l, abstract_matrix)
  { return transposed_ref<const L *>(l); }

  template <class L> inline
  vect_transposed_ref<L *> transposed(L &l, abstract_vector)
  { return vect_transposed_ref<L *>(l); }

  template <class L> inline 
  vect_transposed_ref<const L *> transposed(const L &l, abstract_vector)
  { return vect_transposed_ref<const L *>(l); }

  template <class L> inline
  typename transposed_return<typename linalg_traits<L>::linalg_type,
							L *>::ret_type
  transposed(L &l)
  { return transposed(l, typename linalg_traits<L>::linalg_type()); }

  template <class L> inline
  typename transposed_return<typename linalg_traits<L>::linalg_type,
							const L *>::ret_type
  transposed(const L &l)
  { return transposed(l, typename linalg_traits<L>::linalg_type()); }

  inline bool _is_sparse(abstract_sparse)  { return true;  }
  inline bool _is_sparse(abstract_plain)   { return false; }
  inline bool _is_sparse(abstract_indirect)  { return false; }

  template <class L> inline bool is_sparse(const L &) 
  { return _is_sparse(linalg_traits<L>::storage_type()); }

  /* ******************************************************************** */
  /*		Write                                   		  */
  /* ******************************************************************** */

  template <class L> void write(std::ostream &o, const L &l)
  { write(o, l, typename linalg_traits<L>::linalg_type()); }

  template <class L> inline void write(std::ostream &o, const L &l,
				       abstract_vector) {
    o << "vector(" << vect_size(l) << ") [";
    write(o, l, typename linalg_traits<L>::storage_type());
    o << " ]";
  }

  template <class L> inline void write(std::ostream &o, const L &l,
				       abstract_sparse) {
    typename linalg_traits<L>::const_iterator it = vect_begin(l),
      ite = vect_end(l);
    for (; it != ite; ++it) 
      if (it.index() != size_type(-1)
	  && *it != typename linalg_traits<L>::value_type(0))
	o << " (r" << it.index() << "," << (*it) << ")";
  }

  template <class L> inline void write(std::ostream &o, const L &l,
				       abstract_plain) {
    typename linalg_traits<L>::const_iterator it = vect_begin(l),
      ite = vect_end(l);
    if (it != ite) o << " " << (*it++);
    for (; it != ite; ++it) o << ", " << (*it);
  }

  template <class L> inline void write(std::ostream &o, const L &l,
				       abstract_matrix) {
    write(o, l, typename linalg_traits<L>::sub_orientation());
  }


  template <class L> inline void write(std::ostream &o, const L &l,
				       row_major) {
    o << "matrix(" << mat_nrows(l) << ", " << mat_ncols(l) << ")" << endl;
    for (size_type i = 0; i < mat_nrows(l); ++i) {
      o << "(";
      write(o, mat_row(l, i), typename linalg_traits<L>::storage_type());
      o << " )\n";
    }
  }

  template <class L> inline void write(std::ostream &o, const L &l,row_and_col)
  { write(o, l, row_major()); }

  template <class L> inline void write(std::ostream &o, const L &l,col_and_row)
  { write(o, l, row_major()); }

  // to be corrected
  template <class L> inline void write(std::ostream &o, const L &l,col_major) {
    o << "transposed "; write(o, transposed(l), row_major());
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
    typename temporary_vector<V1>::vector_type w(mat_nrows());
    cerr << "Warning, a temporary is used in scalar product\n";
    mult(ps, v1, w); 
    return vect_sp(w, v2);
  }

  template <class IT1, class IT2> inline
    typename std::iterator_traits<IT1>::value_type
    _vect_sp_plain(IT1 it, IT1 ite, IT2 it2) {
    typename std::iterator_traits<IT1>::value_type res(0);
    size_type n = ((ite - it) >> 3);
    for (size_type i = 0; i < n; ++i, ++it, ++it2) {
      res += conj_product(*it, *it2);
      res += conj_product(*++it, *++it2);
      res += conj_product(*++it, *++it2);
      res += conj_product(*++it, *++it2);
      res += conj_product(*++it, *++it2);
      res += conj_product(*++it, *++it2);
      res += conj_product(*++it, *++it2);
      res += conj_product(*++it, *++it2);
    }
    for (; it != ite; ++it, ++it2) res += conj_product(*it, *it2);
    return res;
  }

  template <class IT1, class V> inline
    typename std::iterator_traits<IT1>::value_type
    _vect_sp_sparse(IT1 it, IT1 ite, const V &v) {
    typename std::iterator_traits<IT1>::value_type res(0);
    for (; it != ite; ++it) 
      if (it.index() != size_type(-1)) res += conj_product(*it, v[it.index()]);
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
   typename number_traits<typename linalg_traits<V>::value_type>
   ::magnitude_type
   vect_norm2(const V &v) {
    typename linalg_traits<V>::const_iterator
      it = vect_begin(v), ite = vect_end(v);
    typename number_traits<typename linalg_traits<V>::value_type>
      ::magnitude_type res(0);
    for (; it != ite; ++it) res += dal::sqr(modulus(*it));
    return sqrt(res);
  }

  /* ******************************************************************** */
  /*		Inifity norm                              		  */
  /* ******************************************************************** */

  template <class V>
  typename number_traits<typename linalg_traits<V>::value_type>
  ::magnitude_type 
  vect_norminf(const V &v) {
    typename linalg_traits<V>::const_iterator
      it = vect_begin(v), ite = vect_end(v);
      typename number_traits<typename linalg_traits<V>::value_type>
	::magnitude_type res(0);
    for (; it != ite; ++it) res = std::max(res, modulus(*it));
    return res;
  }
  
  /* ******************************************************************** */
  /*		norm1                                    		  */
  /* ******************************************************************** */

  template <class V>
  typename number_traits<typename linalg_traits<V>::value_type>
  ::magnitude_type
  vect_norm1(const V &v) {
    typename linalg_traits<V>::const_iterator
      it = vect_begin(v), ite = vect_end(v);
    typename number_traits<typename linalg_traits<V>::value_type>
	::magnitude_type res(0);
    for (; it != ite; ++it) res += modulus(*it);
    return res;
  }

  /* ******************************************************************** */
  /*		Clean                                    		  */
  /* ******************************************************************** */

  template <class L> inline void clean(L &l, double seuil)
  { clean(l, seuil, typename linalg_traits<L>::linalg_type()); }

  template <class L> inline void clean(const L &l, double seuil)
  { clean(linalg_cast(l), seuil);}

  template <class L> inline void clean(L &l, double seuil, abstract_vector)
  { clean(l, seuil, typename linalg_traits<L>::storage_type()); }

  template <class L> void clean(L &l, double seuil, abstract_plain) {
    typename linalg_traits<L>::iterator it = vect_begin(l), ite = vect_end(l);
    for (; it != ite; ++it)
      if (modulus(*it) < seuil)
	*it = typename linalg_traits<L>::value_type(0);
  }

  template <class L> void clean(L &l, double seuil, abstract_sparse) {
    typename linalg_traits<L>::iterator it = vect_begin(l), ite = vect_end(l);
    for (; it != ite; ++it)
      if (modulus(*it) < seuil && it.index() != size_type(-1))
	l[it.index()] = typename linalg_traits<L>::value_type(0);
  }

  template <class L> inline void clean(L &l, double seuil, abstract_matrix) {
    clean(l, seuil, typename principal_orientation_type<typename
	  linalg_traits<L>::sub_orientation>::potype());
  }
  
  template <class L> void clean(L &l, double seuil, row_major) {
    for (size_type i = 0; i < mat_nrows(l); ++i)
      clean(mat_row(l, i), seuil);
  }

  template <class L> void clean(L &l, double seuil, col_major) {
    for (size_type i = 0; i < mat_ncols(l); ++i)
      clean(mat_col(l, i), seuil);
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
  void copy(const L1& l1, const L2& l2) { copy(l1, linalg_cast(l2)); }

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
      copy_vect(mat_const_row(l1, i), mat_row(l2, i),
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
  /*   algorithms are built in order to avoid some conflicts whith        */
  /*   repeated arguments or with overlapping part of a same object.      */
  /*   In the latter case, conflicts are still possible.                  */
  /* ******************************************************************** */
  
  template <class L1, class L2> inline
    void add(const L1& l1, L2& l2) {
      add_spec(l1, l2, typename linalg_traits<L2>::linalg_type());
  }

  template <class L1, class L2> inline
  void add(const L1& l1, const L2& l2) { add(l1, linalg_cast(l2)); }

  template <class L1, class L2>
    void add_spec(const L1& l1, L2& l2, abstract_vector) {
    if (vect_size(l1) != vect_size(l2))
      DAL_THROW(dimension_error,"dimensions mismatch");
    add(l1, l2, typename linalg_traits<L1>::storage_type(),
	typename linalg_traits<L2>::storage_type());
  }

  template <class L1, class L2, class L3>
    void add(const L1& l1, const L2& l2, L3& l3) {
    if (vect_size(l1) != vect_size(l2) || vect_size(l1) != vect_size(l3))
      DAL_THROW(dimension_error,"dimensions mismatch"); 
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
  void add(const L1& l1, const L2& l2, const L3& l3)
  { add(l1, l2, linalg_cast(l3)); }

  template <class IT1, class IT2, class IT3> inline
    void _add_full(IT1 it1, IT2 it2, IT3 it3, IT3 ite) {
    for (; it3 != ite; ++it3, ++it2, ++it1) *it3 = *it1 + *it2;
  }

  template <class IT1, class IT2, class IT3>
    void _add_almost_full(IT1 it1, IT1 ite1, IT2 it2, IT3 it3, IT3 ite3) {
    IT3 it = it3;
    for (; it != ite3; ++it, ++it2) *it = *it2;
    for (; it1 != ite1; ++it1)
      if (it1.index() != size_type(-1)) *(it3 + it1.index()) += *it1;
  }

  template <class IT1, class IT2, class IT3> inline
  void _add_to_full(IT1 it1, IT1 ite1, IT2 it2, IT2 ite2,
		    IT3 it3, IT3 ite3) {
    IT3 it = it3;
    for (; it != ite3; ++it) *it = 0;
    for (; it1 != ite1; ++it1)
      if (it1.index() != size_type(-1)) *(it3 + it1.index()) = *it1;
    for (; it2 != ite2; ++it2)
      if (it2.index() != size_type(-1)) *(it3 + it2.index()) += *it2;    
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
      if (it1.index() != size_type(-1)) l2[it1.index()] += *it1;
  }
  
  template <class L1, class L2> inline
  void add(const L1& l1, L2& l2,
	   abstract_sparse, abstract_sparse) {
    typename linalg_traits<L1>::const_iterator
      it1 = vect_begin(l1), ite1 = vect_end(l1);
    for (; it1 != ite1; ++it1) 
      if (it1.index() != size_type(-1)) l2[it1.index()] += *it1;
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
  { scale(linalg_cast(l), a); }


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
  void mult(const L1& l1, const L2& l2, const L3& l3)
  { mult_const(l1, l2, linalg_cast(l3)); }

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
  void mult_by_col(const L1& l1, const L2& l2, L3& l3, abstract_plain) {
    clear(l3);
    size_type nc = mat_ncols(l1);
    for (size_type i = 0; i < nc; ++i)
      add(scaled(mat_col(l1, i), l2[i]), l3);
  }

  template <class L1, class L2, class L3> inline
  void mult_by_col(const L1& l1, const L2& l2, L3& l3, abstract_sparse) {
    clear(l3);
    typename linalg_traits<L2>::const_iterator it = vect_begin(l2),
      ite = vect_end(l2);
    for (; it != ite; ++it)
      if (*it != typename linalg_traits<L2>::value_type(0))
	add(scaled(mat_col(l1, it.index()), *it), l3);
  }

  template <class L1, class L2, class L3> inline
  void mult_spec(const L1& l1, const L2& l2, L3& l3, row_major)
  { mult_by_row(l1, l2, l3, typename linalg_traits<L3>::storage_type()); }

  template <class L1, class L2, class L3> inline
  void mult_spec(const L1& l1, const L2& l2, L3& l3, col_major)
  { mult_by_col(l1, l2, l3, typename linalg_traits<L2>::storage_type()); }

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
      mult_spec(l1, l2, l3, l4, typename principal_orientation_type<typename
		linalg_traits<L1>::sub_orientation>::potype());
    else {
      #ifdef __GETFEM_VERIFY
        cerr << "Warning, A temporary is used for mult\n";
      #endif
      typename temporary_vector<L3>::vector_type temp(vect_size(l3));
      mult_spec(l1,l2,l3, temp, typename principal_orientation_type<typename
		linalg_traits<L1>::sub_orientation>::potype());
      copy(temp, l4);
    }
  }
  
  template <class L1, class L2, class L3, class L4> inline
  void mult(const L1& l1, const L2& l2, const L3& l3, const L4& l4)
  { mult_const(l1, l2, l3, linalg_cast(l4)); }

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
  void mult_by_col(const L1& l1, const L2& l2, const L3& l3, L4& l4,
		   abstract_plain) {
    copy(l3, l4);
    size_type nc = mat_ncols(l1);
    for (size_type i = 0; i < nc; ++i)
      add(scaled(mat_col(l1, i), l2[i]), l4);
  }

  template <class L1, class L2, class L3, class L4> inline
  void mult_by_col(const L1& l1, const L2& l2, const L3& l3, L4& l4,
		   abstract_sparse) {
    copy(l3, l4);
    typename linalg_traits<L2>::const_iterator it = vect_begin(l2),
      ite = vect_end(l2);
    for (; it != ite; ++it)
      if (*it != typename linalg_traits<L2>::value_type(0))
	add(scaled(mat_col(l1, it.index()), *it), l4);
  }

  template <class L1, class L2, class L3, class L4> inline
  void mult_spec(const L1& l1, const L2& l2, const L3& l3, L4& l4, row_major)
  { mult_by_row(l1, l2, l3, l4, typename linalg_traits<L4>::storage_type()); }

  template <class L1, class L2, class L3, class L4> inline
  void mult_spec(const L1& l1, const L2& l2, const L3& l3, L4& l4, col_major)
  { mult_by_col(l1, l2, l3, l4, typename linalg_traits<L2>::storage_type()); }

  template <class L1, class L2, class L3, class L4> inline
  void mult_spec(const L1& l1, const L2& l2, const L3& l3,
		 L4& l4, abstract_null_type)
  { mult_ind(l1, l2, l3, l4, typename linalg_traits<L1>::storage_type()); }

  template <class L1, class L2, class L3, class L4> inline
  void mult_ind(const L1& l1, const L2& l2, const L3& l3,
		L4& l4, abstract_indirect) {
    DAL_THROW(failure_error,
	  "You have to define gmm::mult(m, v1, v2) for this kind of matrix");
  }

}

#endif //  __GMM_BLAS_H
