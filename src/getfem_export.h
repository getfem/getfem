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
/* Copyright (C) 2001-2002  Yves Renard.                                   */
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

#ifndef __GETFEM_EXPORT_H
#define __GETFEM_EXPORT_H

#include <getfem_mesh_fem.h>
#include <getfem_mat_elem.h>
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
		       const VECT &U, short_type K)
  { // a corriger
    dim_type N = mf.linked_mesh().dim();
    dim_type P = mf.get_qdim();
    std::ofstream o(filename.c_str());
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
      coeff.resize(nbd1);

      o << "DIM = " << int(pgt->dim()) << endl;

      if (!(pf1->is_equivalent())) 
	transfert_to_G(G, mf.linked_mesh().points_of_convex(cv));

      if (pf1->target_dim() != 1)
	DAL_THROW(to_be_done_error, "to be done ... ");
      
      for (size_type i = 0; i < nbd2; ++i)
      {
	/* point corresponding to the i th node.                           */
	pt1 = pgt->transform(pfe->node_of_dof(i),
			     mf.linked_mesh().points_of_convex(cv));
	for (size_type j = 0; j < N; ++j) o << pt1[j] << " ";
	o << "  ";

	/* interpolation of the solution.                                  */
	/* faux dans le cas des éléments vectoriel.                        */
	pt2 = pfe->node_of_dof(i);

	for (size_type k = 0; k < P; ++k) {
	  for (size_type j = 0; j < nbd1 / P; ++j)
	    coeff[j] = U[mf.ind_dof_of_element(cv)[j*P+k]];
	  // il faudrait utiliser les fem_precomp pour accelerer.
	  pf1->interpolation(pt2, G, pgt, coeff, val);
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
		       mesh_fem &mef, VECT &U, short_type &K) {
    std::ifstream ist(fi.c_str());
    size_type token, DIM, Np, N;
    char tmp[100], c;
    bool error = false, began = false;
    dal::dynamic_array<base_node> ptab;
    dal::dynamic_array<scalar_type> vtab;
    size_type nbvtab = 0;
    dim_type &P;
    
    if (!ist) DAL_THROW(file_not_found_error, "File " << fi << " not found");

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
	  // cout << "token " << token << " val = " << k << endl;
	  switch (token) {
	  case 1 : DIM = N = k; if (k == 0) error = true; break;
	  case 2 : P = k; mef.set_qdim(P); if (k == 0) error = true; break;
	  case 3 : K = k; break;
	  case 4 : DIM = k; break;
	  }
	}
      }
      else {
	began = true;
	// cout << "tmp =" << tmp << endl;
	ist.putback(' ');
	for(int k = strlen(tmp)-1; k >= 0; --k) ist.putback(tmp[k]);
	size_type nbpt = 0;
	// cout << "N = " << N << " P = " << int(P) << " K = "
	//      << K << " DIM = " << DIM << endl;
	while(nbpt < 1000) {
	  if (ptab[nbpt].size() != N) ptab[nbpt].resize(N);
	  for (dim_type i = 0; i < N; ++i) ist >> (ptab[nbpt])[i];
	  // cout << "Ptab[" << nbpt << "] = " << ptab[nbpt] << endl;
	  for (short_type i = 0; i < P; ++i) ist >> vtab[nbvtab * P + i];
	  nbpt++; nbvtab++;
	  while(ist.get(c))
	    { if (!(isspace(c)) || c == 10) { ist.putback(c); break; } }
	  ist.get(c); if (c != 10 && !(ist.eof())) { error = true; break; }
	  while(ist.get(c))
	    { if (!(isspace(c)) || c == 10) { ist.putback(c); break; } }
	  ist.get(c); if (c == 10 || ist.eof()) break; else ist.putback(c);
	}
	if (nbpt >= 1000) { error = true; break; }
	if (nbpt == 0) break;

	Np = 1; for (dim_type i = 0; i < DIM; ++i) Np *= K + 1;
	
	if (nbpt == bgeot::alpha(DIM, K)) {
	  std::vector<getfem::base_node> pnode(DIM+1);
	  for (dim_type i = 0; i <= DIM; ++i) 
	    pnode[i] = ptab[bgeot::alpha(i, K) - 1];
	  size_type i = mesh.add_simplex_by_points(DIM, pnode.begin());
	  mef.set_finite_element(i, getfem::PK_fem(DIM, K),
				 getfem::exact_simplex_im(DIM));
	}
	else if (nbpt == Np) {
	  std::vector<getfem::base_node> pnode(1 << DIM);
	  size_type j = 0;
	  for (size_type i = 0; i <= (size_type(1) << DIM); ++i) {
	    pnode[i] = ptab[j];
	    size_type k = i + 1, l = K-1;
	    while (!(k & 1)) { l *= K+1; k >>= 1; }
	    j += l + 1;
	  }
	  j = mesh.add_parallelepiped_by_points(DIM, pnode.begin());
	  mef.set_finite_element(j, getfem::QK_fem(DIM, K),
				 getfem::exact_parallelepiped_im(DIM));
	}
	else if (nbpt == bgeot::alpha(DIM - 1, K) * bgeot::alpha(1, K)) {
	  std::vector<getfem::base_node> pnode(2*DIM);
	  for (dim_type i = 0; i < DIM; ++i)
	    pnode[i] = ptab[bgeot::alpha(i, K) - 1];
	  for (dim_type i = 0; i < DIM; ++i)
	    pnode[i+DIM] = ptab[bgeot::alpha(i, K)-1 + bgeot::alpha(DIM-1, K)];
	  int i = mesh.add_prism_by_points(DIM, pnode.begin());
	  mef.set_finite_element(i, getfem::PK_prism_fem(DIM, K),
				 getfem::exact_prism_im(DIM));
	}
	else 
	  DAL_THROW(failure_error, "Unknown element in file " << fi);
      }
    }
    
    ist.close();
    if (error) DAL_THROW(failure_error, "Format error in file " << fi);
    
    // Data repartition
    // à corriger
    
    dal::bit_vector nn = mesh.convex_index();
    U.resize(mef.nb_dof());
    std::fill(U.begin(), U.end(), 0.0);
    std::vector<int> cp(mef.nb_dof() / P);
    std::fill(cp.begin(), cp.end(), 0);
    size_type l = 0, i;
    for (i << nn; i != size_type(-1); i << nn) {
      size_type nbd = mef.nb_dof_of_element(i);
      for (size_type j = 0; j < nbd / P; ++j, ++l) {
	size_type dof = mef.ind_dof_of_element(i)[j*P];
	for (short_type k = 0; k < P; ++k) U[dof +k] += vtab[l*P+k];
	(cp[dof])++;
      }
    }
    for (i = 0; i < mef.nb_dof() / P; ++i) {
      if (cp[i] == 0) DAL_THROW(internal_error, "Internal error");
      for (short_type k = 0; k < P; ++k)
	U[i*P +k] /= getfem::scalar_type(cp[i]);
    }
    
  }

  /* ********************************************************************* */
  /*                                                                       */
  /*  interpolation of a solution on same mesh.                            */
  /*   - mf_target must be of lagrange type.                               */
  /*                                                                       */
  /* ********************************************************************* */

  template<class VECT>
    void interpolation_solution_same_mesh(mesh_fem &mf_source, mesh_fem &mf_target,
					  const VECT &U, VECT &V)
  {
    dal::bit_vector nn = mf_source.convex_index();
    size_type cv;
    base_node val(1);
    base_matrix G;
    base_vector coeff;
    size_type qdim = mf_source.get_qdim();
    if ( &(mf_source.linked_mesh()) != &(mf_target.linked_mesh()))
      DAL_THROW(failure_error, "Meshes should be the same in this function.");

    if (qdim != mf_target.get_qdim())
      DAL_THROW(failure_error, "Attempt to interpolate a field of dimension " << 
		qdim << " on a mesh_fem whose Qdim is " << 
		mf_target.get_qdim());

    for (cv << nn; cv != ST_NIL; cv << nn) {
      bgeot::pgeometric_trans pgt = mf_source.linked_mesh().trans_of_convex(cv);
      pfem pf_s = mf_source.fem_of_element(cv);
      pfem pf_t = mf_target.fem_of_element(cv);
      size_type nbd_s = pf_s->nb_dof();
      size_type nbd_t = pf_t->nb_dof();
      coeff.resize(nbd_s);

      if (!(pf_s->is_equivalent())) 
	transfert_to_G(G, mf_source.linked_mesh().points_of_convex(cv));

      if (pf_s->target_dim() != 1 || pf_t->target_dim() != 1)
	DAL_THROW(to_be_done_error, "vector FEM interpolation still to be done ... ");
      
      for (size_type i = 0; i < nbd_t; ++i) {
	size_type dof_t = mf_target.ind_dof_of_element(cv)[i*qdim];
	/* interpolation of the solution.                                  */
	/* faux dans le cas des éléments vectoriel.                        */
	for (size_type k = 0; k < qdim; ++k) {
	  for (size_type j = 0; j < nbd_s; ++j) {
	    size_type dof_s = mf_source.ind_dof_of_element(cv)[j*qdim+k];
	    coeff[j] = U[dof_s];
	  }
	  // il faudrait utiliser les fem_precomp pour accelerer.
	  pf_s->interpolation(pf_t->node_of_dof(i), G, pgt, coeff, val);
	  V[dof_t + k] = val[0];
	}
      }
    }
  }

  /* ********************************************************************* */
  /*                                                                       */
  /*  interpolation of a solution on another mesh.                         */
  /*   - mf_target must be of lagrange type.                               */
  /*   - the solution must be continuous.                                  */
  /*                                                                       */
  /* ********************************************************************* */



  template<class VECT>
    void interpolation_solution(mesh_fem &mf_source, mesh_fem &mf_target,
				const VECT &U, VECT &V)
  {
    size_type cv;
    base_node val(1);
    bgeot::geotrans_inv gti;
    dal::dynamic_array<base_node> ptab;
    dal::dynamic_array<size_type> itab;
    base_vector coeff;
    base_matrix G;

    
    if (&mf_source.linked_mesh() == &mf_target.linked_mesh()) {
      interpolation_solution_same_mesh(mf_source, mf_target, U, V);
      return;
    }
    size_type qdim = mf_source.get_qdim();
    if (qdim != mf_target.get_qdim())
      DAL_THROW(failure_error, "Attempt to interpolate a field of dimension " << 
		qdim << " on a mesh_fem whose Qdim is " << 
		mf_target.get_qdim());

    dal::bit_vector tdof_added; 

    /* stocke des paires <indice, repetition> */
    std::deque<std::pair<size_type,size_type> > gti_pt_2_target_dof;

    /* stockage de tous les points dans le gti,
       et remplissage du tableau gti_pt_2_target_dof pour
       pouvoir retrouver le numéro du ddl à partir de son indice
       dans la liste de points de gti 
    */
    for (dal::bit_vector::const_iterator it = mf_target.convex_index().begin();
	 it != mf_target.convex_index().end(); ++it) {
      if (*it) {
	pfem pf_t = mf_target.fem_of_element(it.index());
	size_type qmult = mf_target.get_qdim() / pf_t->target_dim();
	
	if (pf_t->target_dim() != 1) DAL_THROW(failure_error, "still some work to do on vector FEMs!");

	for (size_type j=0; j < pf_t->nb_dof(); ++j) {
	  size_type dof_t = mf_target.ind_dof_of_element(it.index())[j*qmult];
	  if (!tdof_added[dof_t]) {
	    gti.add_point(mf_target.point_of_dof(dof_t));
	    gti_pt_2_target_dof.push_back(std::pair<size_type,size_type>(dof_t,qmult));
	    tdof_added.add(dof_t);
	  }
	}
      }
    }
    // il faudrait controler que tous les ddl de mf_target sont de
    // type lagrange
    dal::bit_vector nn = mf_source.convex_index(), ddl_touched;
    ddl_touched.add(0, mf_target.nb_dof()-1);

    for (cv << nn; cv != ST_NIL; cv << nn)
    {
      bgeot::pgeometric_trans pgt = mf_source.linked_mesh().trans_of_convex(cv);
      //cout << "dealing with convex " << cv << " remaining " << nn.card() << endl;
      size_type nb = gti.points_in_convex(mf_source.linked_mesh().convex(cv),
					  pgt, ptab, itab);
      //cout << "nb points in this convex " << nb << endl;
      // cout << "convex : " << mf_source.linked_mesh().convex(cv) << endl;
      pfem pf_s = mf_source.fem_of_element(cv);
      if (!(pf_s->is_equivalent())) 
	transfert_to_G(G, mf_source.linked_mesh().points_of_convex(cv));
      size_type nbd_s = pf_s->nb_dof();
      coeff.resize(nbd_s);
      for (size_type i = 0; i < nb; ++i)
      {
	size_type dof_t = gti_pt_2_target_dof[itab[i]].first;
	size_type nrep  = gti_pt_2_target_dof[itab[i]].second;
	assert(nrep>0);
	// cout << "dealing with ddl : " << itab[i] << "  coords : " << mf_target.point_of_dof(itab[i]) << " internal coords : " << ptab[i] << endl;
	if (ddl_touched[dof_t])
	{ // inverser les deux boucles pour gagner du temps ?
	  // Il faut verifier que le ddl est bien de Lagrange ...
	  for (size_type k = 0; k < nrep; ++k) {
	    for (size_type j = 0; j < nbd_s; ++j) {
	      size_type dof_s = mf_source.ind_dof_of_element(cv)[j*nrep+k];
	      coeff[j] = U[dof_s];
	      //cerr << "dof_s=" << dof_s << " dof_t=" << dof_t << " k=" << k << " i=" << i << endl;
	    }
	    pf_s->interpolation(ptab[i], G, pgt, coeff, val);
	    V[dof_t + k] = val[0];
	    ddl_touched.sup(dof_t+k);
	  }
	}
      }
    }
    // cout << "ddl untouched : " << ddl_touched << endl;
    if (ddl_touched.card() != 0)
      cerr << "WARNING : in interpolation_solution,"
	   << " all points have not been touched" << endl;
  }



  /*
    OLD INTERPOLATION 
  */
  template<class VECT>
    void interpolation_solution_same_mesh(mesh_fem &mf, mesh_fem &mf_target,
					  const VECT &U, VECT &V, dim_type P)
  {
    dal::bit_vector nn = mf.convex_index();
    size_type cv;
    base_node pt2, val(1);
    base_matrix G;
    base_vector coeff;

    if ( &(mf.linked_mesh()) != &(mf_target.linked_mesh()))
      DAL_THROW(failure_error, "Meshes should be the same in this function.");

    for (cv << nn; cv != ST_NIL; cv << nn)
    {
      bgeot::pgeometric_trans pgt = mf.linked_mesh().trans_of_convex(cv);
      pfem pfe = mf_target.fem_of_element(cv);
      pfem pf1 = mf.fem_of_element(cv);
      size_type nbd1 = mf.nb_dof_of_element(cv);
      size_type nbd2 = pfe->nb_dof();
      coeff.resize(nbd1);

      if (!(pf1->is_equivalent())) 
	transfert_to_G(G, mf.linked_mesh().points_of_convex(cv));

      if (pf1->target_dim() != 1)
	DAL_THROW(to_be_done_error, "to be done ... ");
      
      for (size_type i = 0; i < nbd2; ++i) {
	size_type dof2 = mf_target.ind_dof_of_element(cv)[i];
	/* interpolation of the solution.                                  */
	/* faux dans le cas des éléments vectoriel.                        */
	pt2 = pfe->node_of_dof(i);
	for (size_type k = 0; k < P; ++k) {
	  for (size_type j = 0; j < nbd1; ++j) {
	    size_type dof1 = mf.ind_dof_of_element(cv)[j];
	    coeff[j] = U[dof1*P+k];
	  }
	  // il faudrait utiliser les fem_precomp pour accelerer.
	  pf1->interpolation(pt2, G, pgt, coeff, val);
	  V[dof2*P + k] = val[0];
	}
      }
    }
  }

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
    size_type cv, nb;
    base_node pt3(P), val(1);
    bgeot::geotrans_inv gti;
    dal::dynamic_array<base_node> ptab;
    dal::dynamic_array<size_type> itab;
    base_vector coeff;
    base_matrix G;

    if (&mf.linked_mesh() == &mf_target.linked_mesh()) {
      interpolation_solution_same_mesh(mf, mf_target, U, V, P);
      return;
    }
    // cout << "entering interpolation ... P = " << int(P) << " \n";
    
    nb = mf_target.nb_dof();
    for (size_type i = 0; i < nb; ++i)
      gti.add_point(mf_target.point_of_dof(i));
    // il faudrait controler que tous les ddl de mf_target sont de
    // type lagrange

    dal::bit_vector nn = mf.convex_index(), ddl_touched;
    ddl_touched.add(0, nb-1);

    for (cv << nn; cv != ST_NIL; cv << nn)
    {
      bgeot::pgeometric_trans pgt = mf.linked_mesh().trans_of_convex(cv);
      //cout << "dealing with convex " << cv << " remaining " << nn.card() << endl;
      nb = gti.points_in_convex(mf.linked_mesh().convex(cv),
			    mf.linked_mesh().trans_of_convex(cv), ptab, itab);
      //cout << "nb points in this convex " << nb << endl;
      // cout << "convex : " << mf.linked_mesh().convex(cv) << endl;
      pfem pfe = mf.fem_of_element(cv);
      if (!(pfe->is_equivalent())) 
	transfert_to_G(G, mf.linked_mesh().points_of_convex(cv));
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
	    pfe->interpolation(ptab[i], G, pgt, coeff, val);
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


  void classical_mesh_fem(mesh_fem& mf, short_type K);
}  /* end of namespace getfem.                                             */


#endif /* __GETFEM_EXPORT_H  */
