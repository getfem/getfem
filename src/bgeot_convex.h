/* -*- c++ -*- (enables emacs c++ mode)                                    */
/* *********************************************************************** */
/*                                                                         */
/* Library :  Basic GEOmetric Tool  (bgeot)                                */
/* File    :  bgeot_convex.h :  convex with points                         */
/*     									   */
/* Date : December 20, 1999.                                               */
/* Author : Yves Renard, Yves.Renard@gmm.insa-tlse.fr                      */
/*                                                                         */
/* *********************************************************************** */
/*                                                                         */
/* Copyright (C) 1999-2002  Yves Renard.                                   */
/*                                                                         */
/* This file is a part of GETFEM++                                         */
/*                                                                         */
/* This program is free software; you can redistribute it and/or modify    */
/* it under the terms of the GNU Lesser General Public License as          */
/* published by the Free Software Foundation; version 2.1 of the License.  */
/*                                                                         */
/* This program is distributed in the hope that it will be useful,         */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of          */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           */
/* GNU Lesser General Public License for more details.                     */
/*                                                                         */
/* You should have received a copy of the GNU Lesser General Public        */
/* License along with this program; if not, write to the Free Software     */
/* Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307,  */
/* USA.                                                                    */
/*                                                                         */
/* *********************************************************************** */



#ifndef __BGEOT_CONVEX_H
#define __BGEOT_CONVEX_H

#include <bgeot_convex_structure.h>

namespace bgeot
{
  
  template<class PT, class PT_TAB = std::vector<PT> > class convex {
  public :
    
    typedef PT point_type;
    typedef PT_TAB point_tab_type;
    typedef typename PT_TAB::size_type size_type;
    
    typedef dal::tab_ref_index_ref< typename PT_TAB::const_iterator,
      convex_ind_ct::const_iterator> ref_convex_pt_ct;
    
    typedef dal::tab_ref_index_ref< typename PT_TAB::const_iterator,
      ref_convex_ind_ct::const_iterator> dref_convex_pt_ct;
    
  protected :
    
    pconvex_structure cvs;
    PT_TAB pts;
    
  public :
    
    ref_convex_pt_ct point_of_face(short_type i) const { 
      return ref_convex_pt_ct(pts.begin(), cvs->ind_points_of_face(i).begin(),
			      cvs->ind_points_of_face(i).end() );
    }
    
    ref_convex_pt_ct dir_points(void) const { 
      return ref_convex_pt_ct(pts.begin(), cvs->_dir_points.begin(),
			      cvs->_dir_points.end() );
    }
    
    dref_convex_pt_ct dir_points_of_face(short_type i) const {
      return dref_convex_pt_ct(pts.begin(),
			       cvs->ind_dir_points_of_face(i).begin(),
			       cvs->ind_dir_points_of_face(i).end());
    }
    pconvex_structure structure(void) const { return cvs; }
    pconvex_structure &structure(void) { return cvs; }
    const PT_TAB &points(void) const { return pts; }
    PT_TAB &points(void) { return pts; }
    short_type nb_points(void) const { return cvs->nb_points(); }
    
    void translate(const typename PT::vector_type &v);
    //template <class CONT> void base_of_orthogonal(CONT &tab);
    convex(void) { }
    convex(pconvex_structure c, const PT_TAB &t) : cvs(c), pts(t) {}
    convex(pconvex_structure c) : cvs(c) {}
  };

  template<class PT, class PT_TAB>
    void convex<PT, PT_TAB>::translate(const typename PT::vector_type &v) {
    typename PT_TAB::iterator b = pts.begin(), e = pts.end();
    for ( ; b != e ; ++b) *b += v;
  }

  /*  
  template<class PT, class PT_TAB> template<class CONT>
    void convex<PT, PT_TAB>::base_of_orthogonal(CONT &tab)
  { // programmation a revoir.
    int N = (points())[0].size();
    pconvex_structure cv = structure();
    int n = cv->dim();
    dal::dynamic_array<typename PT::vector_type> _vect;
    vsvector<double> A(N), B(N);
    ref_convex_ind_ct dptf = cv->ind_dir_points_of_face(f);
    int can_b = 0;
    
    for (int i = 0; i < n-1; i++) {
      _vect[i]  = (points())[dptf[i+1]]; _vect[i] -= (points())[dptf[0]];
      
      for (j = 0; j < i; j++)
	A[j] = vect_sp(_vect[i], _vect[j]);
      for (j = 0; j < i; j++)
	{ B = _vect[j]; B *= A[j]; _vect[i] -= B; }
      _vect[i] /= vect_norm2(_vect[i]);
    }
    
    for (int i = n; i < N; i++) {
      _vect[i] = _vect[0];
      vect_random(_vect[i]);
      for (j = 0; j < i; j++)
	A[j] = vect_sp(_vect[i], _vect[j]);
      for (j = 0; j < i; j++)
	{ B = _vect[j]; B *= A[j]; _vect[i] -= B; }
      
      if (vect_norm2(_vect[i]) < 1.0E-4 )
	i--;
      else
	_vect[i] /= vect_norm2(_vect[i]);
    }
    for (int i = n; i < N; i++) tab[i-n] = _vect[i];
  }
  */

  template<class PT, class PT_TAB>
    std::ostream &operator <<(std::ostream &o, const convex<PT, PT_TAB> &cv)
  {
    o << *(cv.structure());
    o << " points : ";
    for (size_type i = 0; i < cv.nb_points(); ++i) o << cv.points()[i] << " ";
    o << endl;
    return o;
  }

  /* ********************************************************************** */
  /* Unstabilized part.                                                     */
  /* ********************************************************************** */

  template<class PT, class PT_TAB>
    convex<PT, PT_TAB> simplex(const PT_TAB &t, int nc)
  { return convex<PT, PT_TAB>(simplex_structure(nc), t); }


  template<class PT, class PT_TAB1, class PT_TAB2>
    convex<PT> convex_product(const convex<PT, PT_TAB1> &cv1,
			      const convex<PT, PT_TAB2> &cv2)
  { // optimisable
    typename convex<PT>::point_tab_type tab;
    tab.resize(cv1.nb_points() * cv2.nb_points());
    size_type i,j,k;
    for (i = 0, k = 0; i < cv1.nb_points(); ++i)
      for (j = 0; j < cv2.nb_points(); ++j, ++k)
	{ tab[k] = (cv1.points())[i]; tab[k] += (cv2.points())[j]; }
    return convex<PT>(
	     convex_product_structure(cv1.structure(), cv2.structure()), tab);
  }

  template<class PT, class PT_TAB1, class PT_TAB2>
    convex<PT> convex_direct_product(const convex<PT, PT_TAB1> &cv1,
				     const convex<PT, PT_TAB2> &cv2)
  {
    if (cv1.nb_points() == 0 || cv2.nb_points() == 0)
      throw std::invalid_argument(
		     "convex_direct_product : null convex product");

    convex<PT> r(convex_product_structure(cv1.structure(), cv2.structure()));
    r.points().resize(r.nb_points());
    std::fill(r.points().begin(), r.points().end(), PT(r.structure()->dim()));
    dim_type dim1 = cv1.structure()->dim();
    typename PT_TAB1::const_iterator it1, it1e = cv1.points().end();
    typename PT_TAB2::const_iterator it2, it2e = cv2.points().end();
    typename convex<PT>::point_tab_type::iterator it = r.points().begin();
    for (it2 = cv2.points().begin(); it2 != it2e; ++it2)
      for (it1 = cv1.points().begin() ; it1 != it1e; ++it1, ++it)
      {
	std::copy((*it1).begin(), (*it1).end(), (*it).begin());
	std::copy((*it2).begin(), (*it2).end(), (*it).begin()+dim1);
      }
    return r;
  }

  template<class PT, class PT_TAB>
    convex<PT> convex_multiply(const convex<PT, PT_TAB> &cv, dim_type n)
  {
    if (cv.nb_points() == 0 || n == 0)
      throw std::invalid_argument(
		     "convex_multiply : null convex product");
    convex<PT> r(multiply_convex_structure(cv.structure(), n));
    r.points().resize(r.nb_points());
    std::fill(r.points().begin(), r.points().end(), PT(r.structure()->dim()));
    dim_type dim1 = cv.structure()->dim();
    typename convex<PT>::point_tab_type::iterator it = r.points().begin();
    typename PT_TAB::const_iterator it1  = cv.points().begin(), it2,
	                            it1e = cv.points().end();
    for (dim_type k = 0; k < n; ++k)
      for (it2 = it1; it2 != it1e; ++it2) *it++ = *it2;
    return r;
  }

  /* structures de reference.                                             */

//   template<class PT> convex<PT> simplex_of_reference(dim_type nc)
//   {
//     convex<PT> res(simplex_structure(nc));
//     res.points().resize(nc+1);
//     PT null(nc); null.fill(0.0);
//     std::fill(res.points().begin(), res.points().end(), null);
//     for (int i = 1; i <= nc; ++i) (res.points()[i])[i-1] = 1.0;
//     return res;
//   }

//   template<class PT> 
//     convex<PT> parallelepiped_of_reference(dim_type nc)
//   { /* optimisable. */
//     if (nc == 1) return simplex_of_reference<PT>(1);
//     else
//       return convex_direct_product<PT>(parallelepiped_of_reference<PT>(nc-1),
// 				          simplex_of_reference<PT>(1));
//   }


}  /* end of namespace bgeot.                                             */


#endif /* __BGEOT_CONVEX_H                                                */
