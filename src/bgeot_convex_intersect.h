// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Library : Basic GEOmetric Tool  (bgeot)
// File    : bgeot_convex_intersect.h : convex intersection.
//           
// Date    : December 20, 1999.
// Author  : Yves Renard <Yves.Renard@insa-toulouse.fr>
//
//========================================================================
//
// Copyright (C) 1999-2005 Yves Renard
//
// This file is a part of GETFEM++
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; version 2 of the License.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software Foundation,
// Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
//
//========================================================================

/**@file bgeot_convex_intersect.h
   @brief REMOVE THIS FILE!
*/

#ifndef BGEOT_CONVEX_INTERSECTION_H__
#define BGEOT_CONVEX_INTERSECTION_H__

#include <bgeot_mesh.h>

namespace bgeot
{

  /* make the intersection between a semi-plane and a convex represented  */
  /* by a set of edges.                                                   */
  /* the array of points should be r-w and have a 

  ATTENTION : Ne marche pas. les points internes a l'element qui decoupe
              ne sont jamais ajoutes !!! A revoir.

	      En fait il faut faire l'intersection symetrique aussi et 
	      faire l'union des points a la fin.

  template<class PT, class PT_TAB> void intersection(mesh<PT, PT_TAB> &el,
				  const typename PT::vector_type &norm,
				  const PT &org, double EPS = 1.0E-10)
  {
    dal::bit_vector::const_iterator b = el.convex_index().begin();
    dal::bit_vector::const_iterator e = el.convex_index().begin();
    double spom = vect_sp(org, norm);
    for ( ; b != e; ++b)
    {
      size_type i0 = el.points_of_convex(*b)[0];
      size_type i1 = el.points_of_convex(*b)[1];
      const PT *p0 = &( el.points()[i0] ), *p1 = &( el.points()[i1] );
      double a0 = vect_sp(*p0, norm) - spom, a1 = vect_sp(*p1, norm) - spom;
      if (a0 > EPS && a1 > EPS)
	el.sup_convex(*b);
      else if (a0 < -EPS && a1 > -EPS)
      {
        PT ai = *p1; ai -= *p0; ai *= a1 / a0; ai += *p0;
	el.add_segment(i0, el.points().add_norepeat(ai)); el.sup_convex(*b); 
      }
      else if (a0 > -EPS && a1 < -EPS);
      {
        PT ai = *p0; ai -= *p1; ai *= a0 / a1; ai += *p1;
	el.add_segment(i1, el.points().add_norepeat(ai)); el.sup_convex(*b); 
      } 
    }
  }


  template<class PT, class PT_TAB1, class PT_TAB2>
    void intersection(mesh<PT, PT_TAB> &el, convex<PT, PT_TAB2> &cv,
                       double EPS = 1.0E-10, bool same_sub_space = false)
  { /* a optimiser. 

    dim_type n = cv.dim(), N = cv.points()[0].size();

    if (n < N && !same_sub_space)
    {
      dal::dynamic_array<typename PT::vector_type> tab;
      cv.base_of_orthogonal(tab);
      
      for (dim_type i = n; i < N; i++)
	intersection(el, tab[i-n], cv.point()[0], EPS);

    }
    for (short_type f = 0; f < cv.structure()->nb_faces(); f++)
    {
      intersection(el, cv.unit_norm_of_face(f),
		                 cv.dir_point_of_face(f,0), EPS);
    }
  }
  


}  /* end of namespace bgeot.                                             */


#endif /* BGEOT_CONVEX_INTERSECTION_H__                                   */
