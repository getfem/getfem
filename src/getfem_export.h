/* -*- c++ -*- (enables emacs c++ mode)                                    */
/* *********************************************************************** */
/*                                                                         */
/* Library :  GEneric Tool for Finite Element Methods (getfem)             */
/* File    :  getfem_export.h : export and import data.                    */
/*     									   */
/* Date : 2004/01/12.                                                      */
/* Author : Yves Renard, Yves.Renard@gmm.insa-tlse.fr                      */
/*          Julien Pommier, pommier@gmm.insa-tlse.fr                       */
/* *********************************************************************** */
/*                                                                         */
/* Copyright (C) 2001-2004  Yves Renard, Julien Pommier.                   */
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


#ifndef GETFEM_EXPORT_H__
#define GETFEM_EXPORT_H__

#include <getfem_mesh_fem.h>
//#include <getfem_mat_elem.h>
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
    void save_solution(const std::string &filename, const mesh_fem &mf,
		       const VECT &U, short_type K)
  { // a corriger
    dim_type N = mf.linked_mesh().dim();
    dim_type P = mf.get_qdim();
    std::ofstream o(filename.c_str());
    if (!o) DAL_THROW(internal_error, "impossible to open file");
    base_node pt1(N), pt2, pt3(P);
    base_matrix G;
    base_vector coeff, val(1);

    o << "% GETFEM++ DATA FILE\nBEGIN DATA ELEMENT\nN = " << int(N)
      << "\nP = " << int(P) << "\nK = " << K << endl << endl;
    o.precision(14);
    for (dal::bv_visitor cv(mf.convex_index()); !cv.finished(); ++cv) {
      bgeot::pgeometric_trans pgt = mf.linked_mesh().trans_of_convex(cv);
      pfem pfe = classical_fem(pgt, K);
      pfem pf1 = mf.fem_of_element(cv);
      size_type nbd1 = mf.nb_dof_of_element(cv);
      size_type nbd2 = pfe->nb_dof();
      coeff.resize(nbd1);

      o << "DIM = " << int(pgt->dim()) << endl;

      if (pf1->need_G()) 
	bgeot::vectors_to_base_matrix(G, mf.linked_mesh().points_of_convex(cv));

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
    dim_type P;
    
    if (!ist) DAL_THROW(file_not_found_error, "File " << fi << " not found");

    N = 1; P = 1; K = 1; DIM = 1; // Default parameters
    
    if (!(ftool::read_until(ist, "BEGIN DATA ELEMENT"))) error = true;

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
    
    U.resize(mef.nb_dof());
    std::fill(U.begin(), U.end(), 0.0);
    std::vector<int> cp(mef.nb_dof() / P);
    std::fill(cp.begin(), cp.end(), 0);
    size_type l = 0;
    for (dal::bv_visitor i(mesh.convex_index()); !i.finished(); ++i) {
      size_type nbd = mef.nb_dof_of_element(i);
      for (size_type j = 0; j < nbd / P; ++j, ++l) {
	size_type dof = mef.ind_dof_of_element(i)[j*P];
	for (short_type k = 0; k < P; ++k) U[dof +k] += vtab[l*P+k];
	(cp[dof])++;
      }
    }
    for (size_type i = 0; i < mef.nb_dof() / P; ++i) {
      if (cp[i] == 0) DAL_THROW(internal_error, "Internal error");
      for (short_type k = 0; k < P; ++k)
	U[i*P +k] /= getfem::scalar_type(cp[i]);
    }
    
  }

  /**
     interpolation of a solution on same mesh.
     - &mf_target.linked_mesh() == &mf_source.linked_mesh()
     - mf_target must be of lagrange type.
     - mf_target's qdim should be equal to mf_source qdim, or equal to 1
     - U.size() >= mf_source.get_qdim()
     - V.size() >= (mf_target.nb_dof() / mf_target.get_qdim()) * mf_source.get_qdim()
  */

  template<class VECT>
    void interpolation_solution_same_mesh(const mesh_fem &mf_source,
					  const mesh_fem &mf_target,
					  const VECT &U, VECT &V)
  {
    base_matrix G;
    size_type qdim = mf_source.get_qdim();
    base_vector coeff, val(qdim);
    if ( &(mf_source.linked_mesh()) != &(mf_target.linked_mesh()))
      DAL_THROW(failure_error, "Meshes should be the same in this function.");

    if (qdim != mf_target.get_qdim() && mf_target.get_qdim() != 1)
      DAL_THROW(failure_error, "Attempt to interpolate a field of dimension " << 
		qdim << " on a mesh_fem whose Qdim is " << 
		int(mf_target.get_qdim()));
    size_type qmult = mf_source.get_qdim()/mf_target.get_qdim();
    fem_precomp_pool fppool;
    /* we should sort convex by their fem */
    for (dal::bv_visitor cv(mf_source.convex_index()); !cv.finished(); ++cv) {
      bgeot::pgeometric_trans pgt = mf_source.linked_mesh().trans_of_convex(cv);
      pfem pf_s = mf_source.fem_of_element(cv);
      pfem pf_t = mf_target.fem_of_element(cv);
      size_type nbd_s = pf_s->nb_dof();
      size_type nbd_t = pf_t->nb_dof();
      coeff.resize(nbd_s*qdim);
      ref_mesh_dof_ind_ct::iterator itdof;
      itdof = mf_source.ind_dof_of_element(cv).begin();
      for (size_type k = 0; k < mf_source.nb_dof_of_element(cv); ++k, ++itdof)
	coeff[k] = U[*itdof];
      if (pf_s->need_G()) 
	bgeot::vectors_to_base_matrix(G, mf_source.linked_mesh().points_of_convex(cv));

      if (pf_s->target_dim() != 1 || pf_t->target_dim() != 1)
	DAL_THROW(to_be_done_error, "vector FEM interpolation still to be done ... ");
      pfem_precomp pfp = fppool(pf_s, pf_t->node_tab());

      itdof = mf_target.ind_dof_of_element(cv).begin();
      for (size_type i = 0; i < nbd_t; ++i, itdof+=mf_target.get_qdim()) {
	size_type dof_t = *itdof*qmult;
	/* faux dans le cas des éléments vectoriel.                        */
	pf_s->interpolation(pfp, i, G, pgt, coeff, val, qdim);
	for (size_type k=0; k < qdim; ++k) 
	  V[dof_t + k] = val[k];
      }
    }
  }


  /**
     interpolation of a solution on another mesh.
     - mf_target must be of lagrange type.
     - the solution should be continuous..
   */
  template<class VECT>
    void interpolation_solution(const mesh_fem &mf_source, const mesh_fem &mf_target,
				const VECT &U, VECT &V) {
    bgeot::geotrans_inv gti;
    dal::dynamic_array<base_node> ptab;
    dal::dynamic_array<size_type> itab;    
    std::vector<base_vector> coeff;
    base_matrix G;
    base_vector val(1);
    
    if (&mf_source.linked_mesh() == &mf_target.linked_mesh()) {
      interpolation_solution_same_mesh(mf_source, mf_target, U, V);
      return;
    }
    size_type qdim_s = mf_source.get_qdim();
    size_type qdim_t = mf_target.get_qdim();
    if (qdim_s != qdim_t && qdim_t != 1)
      DAL_THROW(failure_error, "Attempt to interpolate a field of dimension " << 
		qdim_s << " on a mesh_fem whose Qdim is " << qdim_t);
    /* initialisation of the geotrans_inv */
    {
      dal::bit_vector dof_done; dof_done.sup(0,mf_target.nb_dof());
      /* store all dof nodes into the geotrans_inv */
      for (dal::bv_visitor cv(mf_target.convex_index()); !cv.finished(); ++cv) {
	pfem pf_t = mf_target.fem_of_element(cv);
	if (pf_t->target_dim() != 1) DAL_THROW(failure_error, "still some work to do on vector FEMs!");  
	for (size_type j=0; j < pf_t->nb_dof(); ++j) {
	  size_type dof_t = mf_target.ind_dof_of_element(cv)[j*qdim_t];
	  if (!dof_done[dof_t]) {
	    // TODO: add a function in getfem_fem for inquiry about the lagrangitude of a dof ..
	    //if (dof_is_lagrange(pf_t->dof_types()[j])) {
	    gti.add_point_with_id(mf_target.point_of_dof(dof_t), dof_t);
	    dof_done.add(dof_t);
	    //} else DAL_WARNING(1, "ignoring non lagrange dof\n");
	  }
	}
      }
    }

    /* interpolation */
    dal::bit_vector ddl_touched; ddl_touched.add(0, mf_target.nb_dof());

    for (dal::bv_visitor cv(mf_source.convex_index()); !cv.finished(); ++cv) {
      bgeot::pgeometric_trans pgt=mf_source.linked_mesh().trans_of_convex(cv);
      size_type nb = gti.points_in_convex(mf_source.linked_mesh().convex(cv),
					  pgt, ptab, itab);
      if (!nb) continue;
      pfem pf_s = mf_source.fem_of_element(cv);
      size_type nbd_s = pf_s->nb_dof();
      size_type mult_s = qdim_s / pf_s->target_dim();
      val.resize(pf_s->target_dim());

      /* prepare coefficients for interpolation */
      coeff.resize(mult_s);
      for (size_type k=0; k < mult_s; ++k) coeff[k].resize(nbd_s);
      ref_mesh_dof_ind_ct::iterator itdof = mf_source.ind_dof_of_element(cv).begin();
      for (size_type j=0; j < nbd_s; ++j) 
	for (size_type k=0; k < mult_s; ++k) { coeff[k][j] = U[*itdof]; ++itdof; }
      if (pf_s->need_G()) 
	bgeot::vectors_to_base_matrix(G, mf_source.linked_mesh().points_of_convex(cv));

      for (size_type i = 0; i < nb; ++i) {
	size_type dof_t = itab[i];
	if (ddl_touched[dof_t]) {
	  ddl_touched.sup(dof_t, qdim_t);
	  for (size_type k = 0, pos=(dof_t*qdim_s)/qdim_t; k < mult_s; ++k) {
	    pf_s->interpolation(ptab[i], G, pgt, coeff[k], val);
	    for (size_type l = 0; l < pf_s->target_dim(); ++l) {
	      V[pos++] = val[l];
	    }
	  }
	}
      }
    }
    if (ddl_touched.card() != 0) {
      cerr << "WARNING : in interpolation_solution (different meshes),"
	   << ddl_touched.card() << " dof of the target mesh_fem have not been touched\n";
      for (dal::bv_visitor d(ddl_touched); !d.finished(); ++d) {
        cerr << "ddl_touched[" << d << "]=" << mf_target.point_of_dof(d) << "\n";
      }
    }
  }

  /** 
      interpolation of a solution on a set of sparse points filled in
      the provided geotrans_inv object.  The gradient is also
      interpolated of PVGRAD is non-null.
   */
  template<class VECT>
    void interpolation_solution(const mesh_fem &mf_source, bgeot::geotrans_inv &gti,
				const VECT &U, VECT &V, VECT* PVGRAD = 0) {
    size_type mdim = mf_source.linked_mesh().dim();
    size_type qdim = mf_source.get_qdim();

    base_vector val(1);
    base_matrix valg(1,mdim);
    dal::dynamic_array<base_node> ptab;
    dal::dynamic_array<size_type> itab;
    base_vector coeff;
    base_matrix G;

    dal::bit_vector ddl_touched; ddl_touched.add(0, gti.nb_points());

    for (dal::bv_visitor cv(mf_source.convex_index()); !cv.finished(); ++cv) {
      bgeot::pgeometric_trans pgt=mf_source.linked_mesh().trans_of_convex(cv);
      size_type nb = gti.points_in_convex(mf_source.linked_mesh().convex(cv),
					  pgt, ptab, itab);
      if (nb == 0) continue;
      pfem pf_s = mf_source.fem_of_element(cv);

      if (pf_s->need_G())
	bgeot::vectors_to_base_matrix(G, mf_source.linked_mesh().points_of_convex(cv));
      // cerr << "is_equiv:" << pf_s->is_equivalent() << ",inerp: G=" << G << ",nrow=" << G.nrows() << ", ncols=" << G.ncols() << endl;
      size_type nbd_s = pf_s->nb_dof();
      coeff.resize(nbd_s);
      for (size_type i = 0; i < nb; ++i)
      {
	size_type dof_t = itab[i];
	size_type nrep  = qdim / pf_s->target_dim();
	if (ddl_touched[dof_t])
	{ // inverser les deux boucles pour gagner du temps ?
	  // Il faut verifier que le ddl est bien de Lagrange ...
	  for (size_type k = 0; k < nrep; ++k) {
	    for (size_type j = 0; j < nbd_s; ++j) {
	      size_type dof_s = mf_source.ind_dof_of_element(cv)[j*nrep+k];
	      coeff[j] = U[dof_s];
	    }
	    // cerr << "cv=" << cv << ", ptab[" << i << "]=" << ptab[i] << ", coeff=" << coeff << endl;
	    pf_s->interpolation(ptab[i], G, pgt, coeff, val);
	    V[dof_t*qdim + k] = val[0];
	    
	    if (PVGRAD) {
	      pf_s->interpolation_grad(ptab[i], G, pgt, coeff, valg);
	      dal::copy_n(valg.begin(), mdim,
			  &(*PVGRAD)[dof_t*qdim*mdim]);
	    }
	    ddl_touched.sup(dof_t);
	  }
	}
      }
    }
    if (ddl_touched.card() != 0)
      cerr << "WARNING : in interpolation_solution (set of sparse points),"
	   << ddl_touched.card() << " points have not been touched\n";
  }



  /* ugly remplacement for "pre-qdim era" interpolation ( contrib/compare_solutions uses that.. ) */
  template<class VECT>
    void interpolation_solution(mesh_fem &mf, const mesh_fem &mf_target,
				const VECT &U, VECT &V, dim_type P) {
    if (mf.get_qdim() == 1) { mf.set_qdim(P); interpolation_solution(mf,mf_target,U,V); mf.set_qdim(1); }
    else interpolation_solution(mf,mf_target,U,V);
  }

  /* this function has nothing to do here .. */
  void classical_mesh_fem(mesh_fem& mf, short_type K);
}  /* end of namespace getfem.                                             */


#endif /* GETFEM_EXPORT_H__  */
