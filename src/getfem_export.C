/* *********************************************************************** */
/*                                                                         */
/* Library :  GEneric Tool for Finite Element Methods (getfem)             */
/* File    :  getfem_export.C : export and import data.                    */
/*     									   */
/* Date : October 15, 2001.                                                */
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


#include <getfem_export.h>

namespace getfem
{

  /* ********************************************************************* */
  /*                                                                       */
  /*  Gives the list of edges of a convex face.                            */
  /*                                                                       */
  /* ********************************************************************* */

  void mesh_edge_list_convex(const getfem_mesh &m, size_type i, 
			     edge_list &el, bool merge_convex)
  { // a tester ... optimisable.

    bgeot::pconvex_structure cvs = m.structure_of_convex(i);
    size_type n = cvs->dim();
    size_type nbp = cvs->nb_points();
    size_type ncv = merge_convex ? 0 : i;

    if (nbp == n+1 && cvs == bgeot::simplex_structure(n)) {
      for (dim_type k = 0; k < n; ++k)
	for (dim_type l = k+1; l <= n; ++l)
	  el.add(edge_list_elt(m.ind_points_of_convex(i)[k],
			       m.ind_points_of_convex(i)[l], ncv));
    }
    else if (nbp == (size_type(1) << n) 
	     && cvs == bgeot::parallelepiped_structure(n)) {
      for (size_type k = 0; k < (size_type(1) << n); ++k)
	for (dim_type j = 0; j < n; ++j)
	  if ((k & (1 << j)) == 0)
	    el.add(edge_list_elt(m.ind_points_of_convex(i)[k],
			      m.ind_points_of_convex(i)[k | (1 << j)], ncv));
    }
    else if (nbp == 2 * n && cvs == bgeot::prism_structure(n)) {
      for (dim_type k = 0; k < n - 1; ++k)
	for (dim_type l = k+1; l < n; ++l) {
	  el.add(edge_list_elt(m.ind_points_of_convex(i)[k],
			       m.ind_points_of_convex(i)[l], ncv));
	  el.add(edge_list_elt(m.ind_points_of_convex(i)[k+n],
			       m.ind_points_of_convex(i)[l+n], ncv));
	}
      for (dim_type k = 0; k < n; ++k)
	el.add(edge_list_elt(m.ind_points_of_convex(i)[k],
			     m.ind_points_of_convex(i)[k+n], ncv));
    }
    else {
      dal::dynamic_array<bgeot::pconvex_structure> cvstab;
      dal::dynamic_array< std::vector<size_type> > indpttab;
      size_type ncs = 1;
      cvstab[0] = cvs;
      indpttab[0].resize(cvstab[0]->nb_points());
      std::copy(m.ind_points_of_convex(i).begin(),
		m.ind_points_of_convex(i).end(), indpttab[0].begin());
      
      while (ncs != 0) {
	ncs--;
	bgeot::pconvex_structure cvs = cvstab[ncs];
	if (cvs->dim() == 1) { // il faudrait étendre aux autres cas classiques.
	  
	  for (size_type i = 1; i < cvs->nb_points(); ++i)
	    el.add(edge_list_elt((indpttab[ncs])[i-1],(indpttab[ncs])[i],ncv));
	}
	else {
	  size_type nf = cvs->nb_faces();
	  for (size_type f = 0; f < nf; ++f) {
	    if (cvs->dim() > 2) ++f;
	    cvstab[ncs+f] = (cvs->faces_structure())[f];
	    indpttab[ncs+f].resize(cvs->nb_points_of_face(f));
	    for (size_type k = 0; k < cvs->nb_points_of_face(f); ++k)
	      (indpttab[ncs+f])[k]
		= (indpttab[ncs])[(cvs->ind_points_of_face(f))[k]];
	  }
	  cvstab[ncs] = cvstab[ncs + nf - 1];
	  indpttab[ncs] = indpttab[ncs + nf - 1];
	  ncs += nf - 1;
	}
      }
    }
  }
    

  void mesh_edges_list(const getfem_mesh &m, edge_list &el, 
		       bool merge_convex)
  {
    dal::bit_vector nn = m.convex_index();
    size_type i;
    for (i << nn; i != ST_NIL; i << nn)
      mesh_edge_list_convex(m, i, el, merge_convex);
  }

}  /* end of namespace getfem.                                             */


