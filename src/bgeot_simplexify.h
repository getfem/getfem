/* *********************************************************************** */
/*                                                                         */
/* Library :  Basic GEOmetric Tool  (bgeot)                                */
/* File    :  bgeot_simplexify.h : convex simplexification.                */
/*     									   */
/*                                                                         */
/* Date : December 20, 1999.                                               */
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


#ifndef __BGEOT_SIMPLEXIFY_H
#define __BGEOT_SIMPLEXIFY_H

#include <bgeot_convex_hull.h>

namespace bgeot
{

  /* ******************************************************************** */
  /*	  fonction qui ajoute a sl la liste des simplexes du decoupage    */
  /*         de cvt. cvt est une liste d'aretes.                          */
  /* ******************************************************************** */

  template<class PT_TAB>
    void simplexify(const mesh_structure &cvt, mesh_structure &sl,
		    const PT_TAB &point_list, int N, double EPS)
  { // a refaire ...
    dal::bit_vector nn;
    dal::dynamic_array<mesh_structure> convex_list;
    std::vector<size_type> tab2(N+1);
    int nc = 0;
    convex_list[0] = cvt;
    mesh_structure *cv;

    while (nc >= 0)
    {
      cv = &(convex_list[nc]);

      size_type edge0 = ((cv->convex()).index()).first();
      size_type i0 = (edge0 == ST_NIL) ? ST_NIL
                     : (cv->ind_points_of_convex(edge0))[0];   
      

      if (i0 == ST_NIL)
      {
	nc--;
      }
      else
      { 
	/* ************************************************************** */
	/* Recherche du noeud le moins attache a des aretes.              */
	/* ************************************************************** */

	size_type nbi0 = (cv->ind_points_to_point(i0)).size();
	mesh_point_st_ct::const_iterator b = (cv->points()).begin(),
	                                 e = (cv->points()).end();
	
	for (size_type i = 0; b != e; ++b, ++i)
	{
	  if ((*b).is_valid())
	  {
	    size_type j = (cv->ind_points_to_point(i)).size();
	    if (nbi0 < N)
	    { if (j > nbi0) { i0 = i; nbi0 = j; } }
	    else
	    { if (j < nbi0 && j >= N) { i0 = i; nbi0 = j; } }
	  }
	}
	
	/* l'ideal serait de faire un test sur l'angle au solide          */
	/* minimal.                                                       */

	mesh_point_search_ind_ct tab = cv->ind_points_to_point(i0);
	
	if (nbi0 < N)
	{
	  cout << "WARNING : All points with less than " << N
	       << " neighbours\n          Verify the computation" << endl;
	  nc--;
	}
	
	if (nbi0 == N)
	{
	  /* ************************************************************ */
	  /* Un simplexe est detache du convexe. ajout a la liste des     */
	  /* simplexes et actualisation des listes d'aretes.              */
	  /* ************************************************************ */

	  tab2[0] = i0; std::copy(tab.begin(), tab.end(), tab2.begin()+1);
	  sl.add_simplex(N, tab2.begin());
	  
	  for (size_type i = 0; i < nbi0; ++i)
	  {
	    cv->sup_segment(i0, tab[i]);
	    for (size_type j = i+1; j < nbi0; ++j)
	      for (size_type k = 0; k <= nc; k++)
		convex_list[k].sup_segment(tab[j], tab[i]);
	  }

	  for (size_type i = 0; i < nbi0; i++)
	    for (size_type j = i+1; j < nbi0; j++)
	      for (size_type k = 0; k <= nc; k++)
	      {
		if (((convex_list[k].points())[tab[i]].is_valid())
		    && ((convex_list[k].points())[tab[j]].is_valid()))
		  convex_list[k].add_segment(tab[i], tab[j]);
	      }
	}
	
	if (nbi0 > N)
	{
	  /* ************************************************************ */
	  /* On detache le convexe consitue du point i0 et des points     */
	  /* voisins. La liste des aretes est mise a jour par             */
	  /* l'intermediaire du calcul de l'env. convexe des points       */
	  /* voisins.                                                     */
	  /* ************************************************************ */
	  
	  // cout << "decoupage d'un convexe de " << nbi0 + 1 << " points\n";

	  convex_list[nc+1].clear(); nn.clear(); 
	  for (size_type i = 0; i < nbi0; i++)
	  {
	    cv->sup_segment(i0, tab[i]);
	    nn.add(tab[i]);
	    // convex_list[nc+1].add_edge(i0, tab[i]);
	  }
	  
	  sub_convex_hull(convex_list[nc+1],
		       dal::tab_ref_index_ref<typename PT_TAB::const_iterator,
			          mesh_point_search_ind_ct::iterator>
		  (point_list.begin(), tab.begin(), tab.begin() + tab.size()),
		  point_list[i0], N, EPS);

	  size_type l;
	  nn = convex_list[nc+1].convex_index();
	  for (l << nn; l != size_type(-1); l << nn)
	  {
	    size_type i = (convex_list[nc+1].ind_points_of_convex(l))[0];
	    size_type j = (convex_list[nc+1].ind_points_of_convex(l))[1];
	    for (size_type k = 0; k <= nc; k++)
	      convex_list[k].sup_segment(i, j);
	  }

	  nn = convex_list[nc+1].convex_index();
	  for (l << nn; l != size_type(-1); l << nn)
	  {
	    size_type i = (convex_list[nc+1].ind_points_of_convex(l))[0];
	    size_type j = (convex_list[nc+1].ind_points_of_convex(l))[1];

	    for (size_type k = 0; k <= nc; k++)
	    {
	      if (((convex_list[k].points())[tab[i]].is_valid())
		  && ((convex_list[k].points())[tab[j]].is_valid()))
		convex_list[k].add_segment(i, j);
	    }
	  }

	  for (size_type i = 0; i < nbi0; i++)
	  {
	    convex_list[nc+1].add_segment(i0, tab[i]);
	  }

	  if ( cv->nb_convex() == 0 )
	  {
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


#endif /* __BGEOT_SIMPLEXIFY_H                                            */
