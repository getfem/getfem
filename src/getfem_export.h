/* -*- c++ -*- (enables emacs c++ mode)                                    */
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
    std::ofstream o((filename + char(0)).data());
    if (!o) DAL_THROW(internal_error, "impossible to open file");
    dal::bit_vector nn = mf.convex_index();
    size_type cv;
    base_node pt1(N), pt2, pt3(P), val(1);
    base_matrix G;
    base_vector coeff;

    o << "% GETFEM++ DATA FILE\nBEGIN DATA ELEMENT\nN = " << int(N)
      << "\nP = " << int(P) << "\nK = " << K << endl << endl;
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

      o << "DIM = " << pgt->dim() << endl;

      if (pf1->target_dim() != 1 || !(pf1->is_equivalent()))
	DAL_THROW(to_be_done_error, "to be done ... ");
      
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
	  // il faudrait utiliser les fem_precomp pour accelerer.
	  pf1->interpolation(pt2, G, coeff, val);
	  pt3[k] = val[0];
	}

	for (size_type j = 0; j < P; ++j) o << pt3[j] << " ";
	o << endl;
      }
      o << endl;
    }
    o << "END DATA ELEMENT" << endl;
    o.close();
  }

  template<class VECT>
    void load_solution(const std::string &fi, getfem_mesh &mesh,
		       mesh_fem &mef, VECT &U, dim_type &P, short_type &K) {
    std::ifstream ist((fi + char(0)).data());
    size_type i, token, DIM, Np, N;
    char tmp[100], c;
    bool error = false, began = false;
    dal::dynamic_array<getfem::base_node> ptab;
    dal::dynamic_array<getfem::scalar_type> vtab;
    size_type nbvtab = 0;
    
    if (!ist) DAL_THROW(file_not_found_error, "File " << fi " not found");

    N = 1; P = 1; K = 1; DIM = 1; // Default parameters
    
    if (!(ftool::read_untill(ist, "BEGIN DATA ELEMENT"))) error = true;

    while (!error) {
      if (ist.eof()) { error = true; break; }
      ftool::get_token(ist, tmp, 99); token = 0;
      if (!strcmp(tmp, "N")) token = 1;
      if (!strcmp(tmp, "P")) token = 2;
      if (!strcmp(tmp, "K")) token = 3;
      if (!strcmp(tmp, "DIM")) token = 4;
      if (!strcmp(tmp, "END")) break;
      if (token) {
	if (token >= 1 && token <= 2 && began)
	  cerr << "WARNING : Ignoring principal dimensions redefinition"<<endl;
	else {
	  ftool::get_token(ist, tmp, 99);
	  if (strcmp(tmp, "=")) error = true;
	  ftool::get_token(ist, tmp, 99);
	  int k = atoi(tmp); if (k < 0 || k > 254) error = true;
	  swhitch (token) {
	  case 1 : DIM = N = i; if (k == 0) error = true; break;
	  case 2 : P = i; if (k == 0) error = true; break;
	  case 3 : K = i; break;
	  case 4 : DIM = i; break;
	  }
	}
      }
      else {
	began = true;
	for(int k = strlen(tmp)-1; k >= 0; --k) ist.putback(tmp[k]);
	size_type nbpt = 0;
	while(nbpt < 1000) {
	  if (ptab[nbpt].size() != N) ptab[nbpt].resize(N);
	  for (i = 0; i < N; ++i) ist >> (ptab[nbpt])[i];
	  for (i = 0; i < P; ++i) ist >> vtab[nbvtab * P + i];
	  nbpt++; nbvtab++;
	  while(ist.get(c))
	    { if (!(isspace(c)) || c == 10) { ist.putback(c); break; } }
	  ist.get(c); if (c != 10 && !(ist.eof())) { error = true; break; }
	  while(ist.get(c))
	    { if (!(isspace(c)) || c == 10) { ist.putback(c); break; } }
	  ist.get(c); if (c == 10 || ist.eof()) break; else ist.putback(c);
	}
	if (nbpt == 100) { error = true; break; }
	if (nbpt == 0) break;

	Np = 1; for (i = 0; i < DIM; ++i) Np *= K + 1;
	
	if (nbpt == bgeot::alpha(DIM, K)) {
	  std::vector<getfem::base_node> pnode(DIM+1);
	  for (i = 0; i <= DIM; ++i) pnode[i] = ptab[bgeot::alpha(i, K) - 1];
	  i = mesh.add_simplex_by_points(DIM, pnode.begin());
	  mef.set_finite_element(i, getfem::PK_fem(DIM, K),
				 bgeot::simplex_poly_integration(DIM));
	}
	else if (nbpt == Np) {
	  std::vector<getfem::base_node> pnode(1 << DIM);
	  size_type j;
	  for (i = 0, j = 0; i <= (1 << DIM); ++i) {
	    pnode[i] = ptab[j];
	    size_type k = i + 1, l = K-1;
	    while (!(k & 1)) { l *= K+1; k >>= 1; }
	    j += l + 1;
	  }
	  i = mesh.add_parallelepiped_by_points(DIM, pnode.begin());
	  mef.set_finite_element(i, getfem::QK_fem(DIM, K),
				 bgeot::parallelepiped_poly_integration(DIM));
	}
	else if (nbpt == bgeot::alpha(DIM - 1, K) * bgeot::alpha(1, K)) {
	  std::vector<getfem::base_node> pnode(2*DIM);
	  for (i = 0; i < DIM; ++i)
	    pnode[i] = ptab[bgeot::alpha(i, K) - 1];
	  for (i = 0; i < DIM; ++i)
	    pnode[i+DIM] = ptab[bgeot::alpha(i, K)-1 + bgeot::alpha(DIM-1, K)];
	  i = mesh.add_prism_by_points(DIM, pnode.begin());
	  mef.set_finite_element(i, getfem::PK_prism_fem(DIM, K),
				 bgeot::prism_poly_integration(DIM));
	}
	else 
	  DAL_THROW(failure_error, "Unknown element in file " << fi);
      }
    }
    
    ist.close();
    if (error) DAL_THROW(failure_error, "Format error in file " << fi);
    
    // Data repartition
    
    dal::bit_vector nn = mesh.convex_index();
    U.resize(mef.nb_dof() * P);
    std::fill(U.begin(), U.end(), 0.0);
    std::vector<int> cp(mef.nb_dof());
    std::fill(cp.begin(), cp.end(), 0);
    size_type l = 0;
    for (i << nn; i != size_type(-1); i << nn) {
      size_type nbd = mef.nb_dof_of_element(i);
      for (int j = 0; j < nbd; ++j, ++l) {
	size_type dof = mef.ind_dof_of_element(i)[j];
	for (int k = 0; k < P; ++k) U[dof*P +k] += vtab[l*P+k];
	(cp[dof])++;
      }
    }
    for (i = 0; i < mef.nb_dof(); ++i) {
      if (cp[i] == 0) DAL_THROW(internal_error, "Internal error");
      for (int k = 0; k < P; ++k) U[i*P +k] /= getfem::scalar_type(cp[i]);
    }
    
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
    if (ddl_touched.card() != 0)
      cerr << "WARNING : in interpolation_solution,"
	   << " all points have not been touched" << endl;
  }


  
}  /* end of namespace getfem.                                             */


#endif /* __GETFEM_EXPORT_H  */
