/* *********************************************************************** */
/*                                                                         */
/* Library :  GEneric Tool for Finite Elements Methods (getfem)            */
/* File    :  getfem_regular_meshes.C : basic refinement of meshes.        */
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


#include <getfem_regular_meshes.h>

namespace getfem
{
  void _parallelepiped_regular_simplex_mesh(getfem_mesh &me, dim_type N,
    const base_node &org, const base_vector *ivect, const size_type *iref)
  {
    bgeot::mesh_structure cvt, sl;
    bgeot::convex<base_node>
      pararef = bgeot::parallelepiped_of_reference<base_node>(N);
    base_node a = org;
    size_type i, nbpt = pararef.nb_points();

    // cout << pararef << endl;

    for (i = 0; i < nbpt; ++i)
    {
      a.fill(0.0);
      for (dim_type n = 0; n < N; ++n)
	a.addmul(pararef.points()[i][n], ivect[n]);
      pararef.points()[i] = a;
    }

    // cout << pararef << endl;
    
    cvt.add_convex(pararef.structure(), dal::sequence_iterator<size_type>(0));
    cvt.to_edges();
    // cout << "cvt " << endl;
    // cvt.write_to_file(cout); getchar();
    bgeot::simplexify(cvt, sl, pararef.points(), N, me.eps());
    // cout << "sl " << endl;
    // sl.write_to_file(cout); getchar();

    size_type nbs = sl.nb_convex();
    std::vector<size_type> tab1(N+1), tab(N), tab3(nbpt);
    size_type total = 0;
    std::fill(tab.begin(), tab.end(), 0);
    while (tab[N-1] != iref[N-1])
    {
      for (a = org, i = 0; i < N; i++) 
	a.addmul(scalar_type(tab[i]), ivect[i]);

      for (i = 0; i < nbpt; i++)
	tab3[i] = me.add_point(a + pararef.points()[i]);
	
      for (i = 0; i < nbs; i++)
      {
	bgeot::ref_mesh_point_ind_ct tab2 = sl.ind_points_of_convex(i);
	for (dim_type l = 0; l < N+1; l++)
	  tab1[l] = tab3[(tab2[l] + ((total & 1) ? (nbpt/2) : 0)) % nbpt];
	me.add_simplex(N, tab1.begin());
      }

      for (dim_type l = 0; l < N; l++)
      {
	tab[l]++; total++;
	if (l < N-1 && tab[l] >= iref[l]) { total -= tab[l]; tab[l] = 0; }
	else break;
      }
    }
  }


  void _parallelepiped_regular_prism_mesh(getfem_mesh &me, dim_type N,
    const base_node &org, const base_vector *ivect, const size_type *iref)
  {
    getfem_mesh aux;
    _parallelepiped_regular_simplex_mesh(aux, N-1, org, ivect, iref);
    dal::bit_vector nn = aux.convex_index();
    size_type cv;
    std::vector<base_node> ptab(2 * N);
    

    for (cv << nn; cv != size_type(-1); cv << nn) {
      std::copy(aux.points_of_convex(cv).begin(),
		aux.points_of_convex(cv).end(), ptab.begin());

      for (size_type k = 0; k < iref[N-1]; ++k) {
	
	for (dim_type j = 0; j < N; ++j) ptab[j+N] = ptab[j] + ivect[N-1];
	me.add_prism_by_points(N, ptab.begin());
	
	std::copy(ptab.begin()+N, ptab.end(), ptab.begin());
      }
    }
  }



  void _parallelepiped_regular_mesh(getfem_mesh &me, dim_type N,
    const base_node &org, const base_vector *ivect, const size_type *iref)
  {
    bgeot::convex<base_node>
      pararef = bgeot::parallelepiped_of_reference<base_node>(N);
    base_node a = org;
    size_type i, nbpt = pararef.nb_points();

    for (i = 0; i < nbpt; ++i)
    {
      a.fill(0.0);
      for (dim_type n = 0; n < N; ++n)
	a.addmul(pararef.points()[i][n], ivect[n]);
      pararef.points()[i] = a;
    }

    std::vector<size_type> tab1(N+1), tab(N), tab3(nbpt);
    size_type total = 0;
    std::fill(tab.begin(), tab.end(), 0);
    while (tab[N-1] != iref[N-1])
    {
      for (a = org, i = 0; i < N; i++) 
	a.addmul(scalar_type(tab[i]), ivect[i]);

      for (i = 0; i < nbpt; i++)
	tab3[i] = me.add_point(a + pararef.points()[i]);

      me.add_convex(bgeot::parallelepiped_trans(N, 1), tab3.begin());
     
      for (dim_type l = 0; l < N; l++)
      {
	tab[l]++; total++;
	if (l < N-1 && tab[l] >= iref[l]) { total -= tab[l]; tab[l] = 0; }
	else break;
      }
    }
  }

}  /* end of namespace getfem.                                             */
