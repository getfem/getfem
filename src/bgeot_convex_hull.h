// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Library : Basic GEOmetric Tool  (bgeot)
// File    : bgeot_convex_hull.h : compute the convex hull.
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

/**@file bgeot_convex_hull.h
   @brief convex hull of a set of points (broken..)
*/

#ifndef BGEOT_CONVEX_HULL_H__
#define BGEOT_CONVEX_HULL_H__

#include <bgeot_mesh_structure.h>

/* a refaire ... */

namespace bgeot
{

  /**
    fonction qui a partir d'un liste de points donne l'env. convexe
    de dimension egale a celle des points. 
    Principe: toutes les arêtes possible sont mises dans le maillage,
    puis tous les simplexes possibles sont construits. Pour chaque simplexe,
    on elimine l'ensemble des arête qui croisent le simplexe.
   */
    
  template<class PT_TAB>
    void convex_hull(mesh_structure &cv, const PT_TAB &point_list, size_type N,
		                                                  double EPS)
  { /* difference entre dimension des points et du convexe pas bien geree. */

    typedef typename PT_TAB::value_type PT;
    size_type nb = point_list.begin() - point_list.end();
    size_type i,j;
    dal::dynamic_array<typename PT::vector_type> vect_list;
    std::vector<size_type> simplex(N+1);
    dal::bit_vector nne, nn;
    base_matrix A(N, N);
   
    /* ****************************************************************** */
    /*	  1 - Ajout de toutes les aretes.                                 */
    /* ****************************************************************** */

    for (i = 0; i < nb; i++)
    { for (j = 0; j < i; j++) nne.add(cv.add_segment(i, j)); }

    if (nb > N+1) {
      /* **************************************************************** */
      /*	  1 - Parcourt des simplexes.                             */
      /* **************************************************************** */

      for (i = 0; i < N+1; i++) simplex[i] = i;

      while (simplex[0] == 0) {
	for (i = 0; i < N; i++) {
	  vector_from(point_list[simplex[0]],
		      point_list[simplex[i+1]], vect_list[i]);
	  for (j = 0; j <= i; j++)
	    A(i,j) = A(j,i) = vect_sp(vect_list[i], vect_list[j]);
	}

	double s = gmm::lu_inverse(A);

	/* if the simplex is not flat .. */
	if (gmm::abs(s) > EPS) {
   
	  nn = nne;
	  while (nn.card() > 0) {
	    i = nn.take_first();
	    /* ********************************************************** */
	    /*	  1 - Parcourt des aretes.                                */
	    /* ********************************************************** */

	    size_type i1 = cv.ind_points_of_convex(i)[0],
	              i2 = cv.ind_points_of_convex(i)[1];
	    size_type ico = 0;
	    
	    for (j = 0; j < N+1; j++)
	      if (simplex[j] == i1 || simplex[j] == i2) ico++;
	    
	    // if the edge does not belong to the simplex..
	    if (ico < 2) {
	      typename PT::vector_type D0
		= vector_from(point_list[simplex[0]], point_list[i1]);
	      typename PT::vector_type DV
		= vector_from(point_list[i1], point_list[i2]);
	      base_small_vector d0(N), dd0(N), dV(N), ddV(N), v11(N);
	      for (j = 0; j < N; j++) {
		v11[j] = 1.0;
		d0[j] = vect_sp(D0, vect_list[j]);
		dV[j] = vect_sp(DV, vect_list[j]);
	      }
	      gmm::mult(A,d0,dd0);
	      gmm::mult(A,dV,ddV);
	      
	      double a = vect_sp(dd0, v11);
	      double b = vect_sp(ddV, v11);
	      double tmin, tmax;
	      bool btmin = false, btmax = false;

	      if (b > EPS)
	      { btmax = true; tmax = (1.0 - a) / b; }
	      else if (b < -EPS)
	      { btmin = true; tmin = (1.0 - a) / b; }
	      else if (a < -EPS || a > 1.0 + EPS)
	      { btmin = btmax = true; tmin = tmax = 0.0; }
	      
	      for (j = 0; j < N; j++)
	      {
		a = dd0[j]; b = ddV[j];
		if (b > EPS)
		{
		  if (!btmin) { btmin = true; tmin = (- a) / b; }
		  else { tmin = std::max(tmin, (- a) / b); }
		}
		if (b < -EPS)
		{
		  if (!btmax) { btmax = true; tmax = (- a) / b; }
		  else { tmax = std::min(tmax, (- a) / b); }
		}
		else if (a < -EPS)
		{ btmin = btmax = true; tmin = tmax = 0.0; }
	      }
	      
	      if (tmax > tmin + EPS)
	      {
		d0 = ddV; d0 *= tmin; d0 += dd0;
		dV = ddV; dV *= tmax; dV += dd0;

		a = vect_norm2(d0); b = vect_norm2(dV);

		// remove the edge if it crosses the simplex
		if (!( (a < EPS && b > 1.0 - EPS)
		    || (a > 1.0 - EPS && b > 1.0 - EPS)
		    || (a > 1.0 - EPS && b < EPS) ))
		  { nne.sup(i); cv.sup_convex(i); }

	      }
	    }
	  }
	}
	i = N; while (simplex[i] == nb-1+i-N) i--;
	simplex[i]++; while(i < N) { simplex[i+1] = simplex[i] + 1; i++; }
      }
    }
  }

  /* ******************************************************************** */
  /*	  fonction qui a partir d'un liste de points donne l'env. convexe */
  /*      de dimension egale a celle des points moins un.                 */
  /* ******************************************************************** */


  template<class PT_TAB>
    void sub_convex_hull(mesh_structure &cv, const PT_TAB &point_list,
		  const typename PT_TAB::value_type &pt0, int N, double EPS)
  {
    typedef typename PT_TAB::value_type PT;
    dal::dynamic_array<PT> sub_point_list;

    if (point_list.begin() != point_list.end())
    {

      typename PT::vector_type pn = 
	vector_from(pt0, dal::mean_value(point_list.begin(),point_list.end()));

      typename PT_TAB::const_iterator b = point_list.begin(),
	                              e = point_list.end();
      size_type i = 0;
      for ( ; b != e; ++b, ++i)
      {
	sub_point_list[i] = pn;
	sub_point_list[i] *= -1.0 * vect_sp(pn, point_list[i]);
	sub_point_list[i] += point_list[i];
      }
      sub_point_list[i] = pn;
      
      convex_hull(cv, dal::tab_ref<typename dal::dynamic_array<PT>::iterator>(
	      sub_point_list.begin(), sub_point_list.begin()+i), N-1, EPS);
    }
  }



}  /* end of namespace bgeot.                                             */


#endif /* BGEOT_CONVEX_HULL_H__                                           */
