/* *********************************************************************** */
/*                                                                         */
/* Library :  GEneric Tool for Finite Element Methods (getfem)             */
/* File    :  getfem_export.h : export and import data.                    */
/*     									   */
/* Date : October 2, 2001.                                                 */
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


#ifndef __GETFEM_EXPORT_H
#define __GETFEM_EXPORT_H

#include <getfem_mesh_fem.h>
#include <bgeot_geotrans_inv.h>
#include <dal_tree_sorted.h>

namespace getfem
{

  /* ********************************************************************* */
  /*                                                                       */
  /*  Save a solution in a file with a Pk interpolation.                   */
  /*                                                                       */
  /* ********************************************************************* */

  template<class VECT>
    void save_solution(const std::string &filename, mesh_fem &mf,
		       const VECT &U, dim_type P, short_type K)
  { 
    dim_type N = mf.linked_mesh().dim();
    STD_NEEDED ofstream o((filename + char(0)).data());
    assert(o != NULL);
    dal::bit_vector nn = mf.convex_index();
    size_type cv;
    base_node pt1(N), pt2, pt3(P), val(1);
    base_matrix G;
    base_vector coeff;

    o << "% GETFEM++ DATA FILE\nDATA BY ELEMENT\nN = " << int(N) << "\nP = ";
    o << int(P) << "\nK = " << K << endl;
    o.precision(14);
    for (cv << nn; cv != ST_NIL; cv << nn)
    {
      bgeot::pgeometric_trans pgt = mf.linked_mesh().trans_of_convex(cv);
      pfem pfe = classical_fem(pgt, K);
      pfem pf1 = mf.fem_of_element(cv);
      size_type nbd1 = mf.nb_dof_of_element(cv);
      size_type nbd2 = pfe->nb_dof();
      size_type nbpt = pgt->nb_points();
      coeff.resize(nbd1);
      assert(pf1->target_dim() == 1);
      assert(pf1->is_equivalent());
      for (size_type i = 0; i < nbd2; ++i)
      {
	/* point corresponding to the i th node.                           */
	const base_node *ppt = &(pfe->node_of_dof(i));
	pt1.fill(0.0);
	for (size_type j = 0; j < nbpt; ++j)
	{
	  pt1.addmul(pgt->poly_vector()[j].eval(ppt->begin()),
		   mf.linked_mesh().points_of_convex(cv)[j]);
	}
	for (size_type j = 0; j < N; ++j) o << pt1[j] << " ";
	o << "  ";

	/* interpolation of the solution.                                  */
	/* faux dans le cas des éléments non tau-equivalents ou vectoriel. */
	pt2 = pfe->node_of_dof(i);
	for (size_type k = 0; k < P; ++k) {
	  for (size_type j = 0; j < nbd1; ++j) {
	    size_type dof1 = mf.ind_dof_of_element(cv)[j];
	    coeff[j] = U[dof1*P+k];
	  }
	  pf1->interpolation(pt2, G, coeff, val);
	  pt3[k] = val[0];
	}

	for (size_type j = 0; j < P; ++j) o << pt3[j] << " ";
	o << endl;
      }
      o << endl;
    }
    o.close();
  }

  /* ********************************************************************* */
  /*                                                                       */
  /*  Gives the list of edges of a mesh.                                   */
  /*                                                                       */
  /* ********************************************************************* */

  struct edge_list_elt
  {
    size_type i, j;
    inline bool operator < (const edge_list_elt &e) const
    {
      if (i < e.i) return true; if (i > e.i) return false; 
      if (j < e.j) return true; return false;
    }
    edge_list_elt(size_type ii, size_type jj)
    { i = std::min(ii, jj); j = std::max(ii, jj); }
    edge_list_elt(void) {}
  };

  typedef dal::dynamic_tree_sorted<edge_list_elt> edge_list;
  

  void mesh_edges_list(const getfem_mesh &m, edge_list &el);


  /* ********************************************************************* */
  /*                                                                       */
  /*  interpolation of a solution on another mesh.                         */
  /*   - mf_target must be of lagrange type.                               */
  /*   - the solution must be continuous.                                  */
  /*                                                                       */
  /* ********************************************************************* */

  template<class VECT>
    void interpolation_solution(mesh_fem &mf, mesh_fem &mf_target,
				const VECT &U, VECT &V, dim_type P)
  {
    dim_type N = mf.linked_mesh().dim();
    size_type cv, nb;
    base_node pt3(P), val(1);
    bgeot::geotrans_inv gti;
    dal::dynamic_array<base_node> ptab;
    dal::dynamic_array<size_type> itab;
    base_vector coeff;
    base_matrix G;

    // cout << "entering interpolation ... P = " << int(P) << " \n";
    
    nb = mf_target.nb_dof();
    for (size_type i = 0; i < nb; ++i)
      gti.add_point(mf_target.point_of_dof(i));

    dal::bit_vector nn = mf.convex_index(), ddl_touched;
    ddl_touched.add(0, nb-1);

    for (cv << nn; cv != ST_NIL; cv << nn)
    {
      // cout << "dealing with convex " << cv << endl;
      nb = gti.points_in_convex(mf.linked_mesh().convex(cv),
			    mf.linked_mesh().trans_of_convex(cv), ptab, itab);
      // cout << "nb points in this convex " << nb << endl;
      // cout << "convex : " << mf.linked_mesh().convex(cv) << endl;
      pfem pfe = mf.fem_of_element(cv);
      size_type nbd1 = pfe->nb_dof();
      coeff.resize(nbd1);
      for (size_type i = 0; i < nb; ++i)
      {
	// cout << "dealing with ddl : " << itab[i] << "  coords : " << mf_target.point_of_dof(itab[i]) << " internal coords : " << ptab[i] << endl;
	if (ddl_touched[itab[i]])
	{ // inverser les deux boucles pour gagner du temps ?
	  // Il faut verifier que le ddl est bien de Lagrange ...
	  // Il faudrait utiliser les prefem_precomp pour accelerer ...
	  pt3.fill(0.0);
	  for (size_type k = 0; k < P; ++k) {
	    for (size_type j = 0; j < nbd1; ++j) {
	      size_type dof1 = mf.ind_dof_of_element(cv)[j];
	      coeff[j] = U[dof1*P+k];
	    }
	    pfe->interpolation(ptab[i], G, coeff, val);
	    pt3[k] = val[0];

	  }
	  
	  // cout << "Résultat : " << pt3 << endl;
	  for (size_type j = 0; j < P; ++j) V[itab[i]*P+j] = pt3[j];
	  ddl_touched.sup(itab[i]);
	}
      }
    }
    // cout << "ddl untouched : " << ddl_touched << endl;
    assert(ddl_touched.card() == 0); // lancer seulement une alerte ?
  }


  
}  /* end of namespace getfem.                                             */


#endif /* __GETFEM_EXPORT_H  */
