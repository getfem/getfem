// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Library : Basic GEOmetric Tool  (bgeot)
// File    : bgeot_simplexify.h : convex simplexification.
//           
// Date    : December 20, 1999.
// Author  : Yves Renard <Yves.Renard@insa-toulouse.fr>
//
//========================================================================
//
// Copyright (C) 1999-2006 Yves Renard
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

/**@file bgeot_simplexify.h
   @brief simplexify a convex (broken..)
*/

#ifndef BGEOT_SIMPLEXIFY_H__
#define BGEOT_SIMPLEXIFY_H__

#include <bgeot_convex_hull.h>

namespace bgeot {

  /* ******************************************************************** */
  /*	  fonction qui ajoute a sl la liste des simplexes du decoupage    */
  /*         de cvt. cvt est une liste d'aretes.                          */
  /* ******************************************************************** */

  template<class PT_TAB>
    void simplexify(const mesh_structure &cvt, mesh_structure &sl,
		    const PT_TAB &point_list, size_type N, double EPS) {
    // to be done again (calling delaunay) ...
    typedef typename PT_TAB::value_type PT;
    dal::dynamic_array<mesh_structure> convex_list;
    std::vector<size_type> tab2(N+1);
    size_type nc = 0;
    convex_list[0] = cvt;
    mesh_structure *cv;

    while (nc != size_type(-1)) {
      cv = &(convex_list[nc]);
      
      size_type edge0 = (cv->convex_index()).first();
      size_type i0 = (edge0 == ST_NIL) ? ST_NIL
	: (cv->ind_points_of_convex(edge0))[0];   
      
      
      if (i0 == ST_NIL) {
	nc--;
      }
      else { 
	/* ************************************************************** */
	/* Recherche du noeud le moins attache a des aretes.              */
	/* ************************************************************** */
	
	size_type nbi0 = (cv->ind_points_to_point(i0)).size();
	
	for (size_type i = 0; i < cv->nb_max_points(); ++i) {
	  if (cv->is_point_valid(i)) {
	    size_type j = (cv->ind_points_to_point(i)).size();
	    if (nbi0 < N) { if (j > nbi0) { i0 = i; nbi0 = j; } }
	    else { if (j < nbi0 && j >= N) { i0 = i; nbi0 = j; } }
	  }
	}
	
	/* l'ideal serait de faire un test sur l'angle au solide          */
	/* minimal.                                                       */
	
	mesh_structure::ind_set tab0 = cv->ind_points_to_point(i0);
	std::vector<size_type> tab; tab.assign(tab0.begin(), tab0.end());

	if (nbi0 < N) {
	  cout << "WARNING : All points with less than " << N
	       << " neighbours\n          Verify the computation" << endl;
	  nc--;
	}
	
	if (nbi0 == N) {
	  /* ************************************************************ */
	  /* Un simplexe est detache du convexe. ajout a la liste des     */
	  /* simplexes et actualisation des listes d'aretes.              */
	  /* ************************************************************ */
	  
	  tab2[0] = i0; std::copy(tab.begin(), tab.end(), tab2.begin()+1);
	  
	  /* ensure that the simplex has direct orientation */
	  if (N>1) {
	    base_matrix M(N,N);
	    for (size_type i = 0; i < N; ++i)
	      for (size_type j = 0; j < N; ++j) 
		M(i,j) = point_list[tab2[i+1]][j] - point_list[tab2[i]][j];
	    if (gmm::lu_det(M) < 0) std::swap(tab2[0],tab2[1]);
	  }
	  
	  sl.add_simplex(N, tab2.begin());
	  
	  for (size_type i = 0; i < nbi0; ++i) {
	    cv->sup_segment(i0, tab[i]);
	    for (size_type j = i+1; j < nbi0; ++j)
	      for (size_type k = 0; k <= nc; k++)
		convex_list[k].sup_segment(tab[j], tab[i]);
	  }
	  
	  for (size_type i = 0; i < nbi0; i++)
	    for (size_type j = i+1; j < nbi0; j++)
	      for (size_type k = 0; k <= nc; k++) {
		if (convex_list[k].is_point_valid(tab[i])
		    && convex_list[k].is_point_valid(tab[j]))
		  convex_list[k].add_segment(tab[i], tab[j]);
	      }
	}
	
	if (nbi0 > N) {
	  /* ************************************************************ */
	  /* On detache le convexe consitue du point i0 et des points     */
	  /* voisins. La liste des aretes est mise a jour par             */
	  /* l'intermediaire du calcul de l'env. convexe des points       */
	  /* voisins.                                                     */
	  /* ************************************************************ */
	  
	  // cout << "decoupage d'un convexe de " << nbi0 + 1 << " points\n";
	  
	  convex_list[nc+1].clear(); 
	  for (size_type i = 0; i < nbi0; i++) {
	    cv->sup_segment(i0, tab[i]);
	    // convex_list[nc+1].add_edge(i0, tab[i]);
	  }
	  
	  std::vector<PT> sub_point_list(tab.size());
	  for (size_type jj = 0; jj < tab.size(); ++jj)
	    sub_point_list[jj] = point_list[tab[jj]];

	  sub_convex_hull(convex_list[nc+1], sub_point_list,
			  point_list[i0], N, EPS);
	
	  for (dal::bv_visitor l(convex_list[nc+1].convex_index());
	       !l.finished(); ++l) {
	    size_type i = (convex_list[nc+1].ind_points_of_convex(l))[0];
	    size_type j = (convex_list[nc+1].ind_points_of_convex(l))[1];
	    for (size_type k = 0; k <= nc; k++)
	      convex_list[k].sup_segment(i, j);
	  }
	  
	  for (dal::bv_visitor l(convex_list[nc+1].convex_index());
	       !l.finished(); ++l) {
	    size_type i = (convex_list[nc+1].ind_points_of_convex(l))[0];
	    size_type j = (convex_list[nc+1].ind_points_of_convex(l))[1];
	    
	    for (size_type k = 0; k <= nc; k++) {
	      if (convex_list[k].is_point_valid(tab[i])
		  && convex_list[k].is_point_valid(tab[j]))
		convex_list[k].add_segment(i, j);
	    }
	  }
	  
	  for (size_type i = 0; i < nbi0; i++) {
	    convex_list[nc+1].add_segment(i0, tab[i]);
	  }
	  
	  if ( cv->nb_convex() == 0 ) {
	    cout << "WARNING : convex of null volume." << endl
		 << "          Verify the computation ... " << endl;
	  }
	  else
	    nc++;
	}
      }
    }
  }
}  /* end of namespace bgeot.                                             */


#endif /* BGEOT_SIMPLEXIFY_H__                                            */
