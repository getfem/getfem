/* -*- c++ -*- (enables emacs c++ mode)                                    */
/* *********************************************************************** */
/*                                                                         */
/* Library :  GEneric Tool for Finite Element Methods (getfem)             */
/* File    :  getfem_assembling.h : assemble linear system for fem.        */
/*     									   */
/* Date : November 17, 2000.                                               */
/* Authors : Yves Renard, Yves.Renard@gmm.insa-tlse.fr                     */
/*           Julien Pommier, pommier@gmm.insa-tlse.fr                      */
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


#ifndef __GETFEM_ASSEMBLING_H
#define __GETFEM_ASSEMBLING_H

#include <getfem_mesh_fem.h>
#include <bgeot_abstract_linalg.h>
namespace getfem
{

  /* ********************************************************************* */
  /*                                                                       */
  /*  Algorithmes d'assemblage pour quelques problemes elliptiques.        */
  /*                                                                       */
  /* ********************************************************************* */

  /* ********************************************************************* */
  /*	Rigidity matrix for laplacian.                                     */
  /* ********************************************************************* */

  template<class MAT, class VECT>
    void assembling_rigidity_matrix_for_laplacian(MAT &RM, const mesh_fem &mf,
						  const mesh_fem &mfdata, const VECT &A)
  { // optimisable
    size_type cv, nbd1, nbd2, N = mf.linked_mesh().dim();
    dal::bit_vector nn = mf.convex_index();
    base_tensor t;
    pfem pf1, pf2, pf1prec = 0, pf2prec = 0;
    pintegration_method pim, pimprec = 0;
    bgeot::pgeometric_trans pgt, pgtprec = 0;
    pmat_elem_computation pmec = 0;

    if (&(mf.linked_mesh()) != &(mfdata.linked_mesh()))
      DAL_THROW(std::invalid_argument,
		"This assembling procedure only works on a single mesh");

    for (cv << nn; cv != ST_NIL; cv << nn)
    {
      pf1 =     mf.fem_of_element(cv); nbd1 = pf1->nb_dof();
      pf2 = mfdata.fem_of_element(cv); nbd2 = pf2->nb_dof();
      pgt = mf.linked_mesh().trans_of_convex(cv);
      pim = mf.int_method_of_element(cv);
      if (pf1prec != pf1 || pf2prec != pf2 || pgtprec != pgt || pimprec != pim)
      {
	pmat_elem_type pme = mat_elem_product(mat_elem_product(
	    mat_elem_grad(pf1), mat_elem_grad(pf1)), mat_elem_base(pf2));
	pmec = mat_elem(pme, pim, pgt);
	pf1prec = pf1; pf2prec = pf2; pgtprec = pgt; pimprec = pim;
      }
      pmec->gen_compute(t, mf.linked_mesh().points_of_convex(cv));
      // cout << "elem matrix " << t << endl;
      base_tensor::iterator p = t.begin();
      for (size_type r = 0; r < nbd2; r++) {
	size_type dof3 = mfdata.ind_dof_of_element(cv)[r];
	for (size_type l = 0; l < N; l++) {
	  for (size_type i = 0; i < nbd1; i++) {
	    size_type dof2 = mf.ind_dof_of_element(cv)[i];
	    p += l * nbd1;
	    for (size_type j = 0; j < nbd1; j++, ++p) {
	      size_type dof1 = mf.ind_dof_of_element(cv)[j];
	      if (dof1 >= dof2) { 
		RM(dof1, dof2) += A[dof3]*(*p);
		RM(dof2, dof1) = RM(dof1, dof2);
	      }
	    }
	    p += (N-l-1) * nbd1;
	  }
	}
      }
      if (p != t.end()) DAL_THROW(dal::internal_error, "internal error"); 
    }
  }

  /* ********************************************************************* */
  /*	assembling of the term p.div(v) in Stockes mixed FEM.              */
  /* ********************************************************************* */

  template<class MAT, class VECT>
    void assembling_mixed_pressure_term(MAT &B,
					const mesh_fem &mf_u,
					const mesh_fem &mf_p,
					const mesh_fem &mf_d,
					const VECT &DATA) {
    size_type cv;
    dal::bit_vector nn = mf_u.convex_index();

    base_tensor t;

    pmat_elem_computation pmec = 0;

    pfem pf_u, pf_p, pf_d;
    pfem pf_u_prec = NULL, pf_p_prec = NULL, pf_d_prec = NULL;
    pintegration_method pim, pimprec = 0;
    bgeot::pgeometric_trans pgt, pgtprec = NULL;
    size_type nbdof_u, nbdof_p, nbdof_d;
    size_type N = mf_u.linked_mesh().dim();

    if (&(mf_u.linked_mesh()) != &(mf_p.linked_mesh())
	|| &(mf_u.linked_mesh()) != &(mf_d.linked_mesh()))
      DAL_THROW(std::invalid_argument,
		"This assembling procedure only works on a single mesh");

    /* loop over all convexes */
    for (cv << nn; cv != ST_NIL; cv << nn) {

      pf_u = mf_u.fem_of_element(cv); nbdof_u = pf_u->nb_dof();
      pf_p = mf_p.fem_of_element(cv); nbdof_p = pf_p->nb_dof();
      pf_d = mf_d.fem_of_element(cv); nbdof_d = pf_d->nb_dof();
      pgt = mf_u.linked_mesh().trans_of_convex(cv);
      pim = mf_u.int_method_of_element(cv);

      /* avoids recomputation of already known pmat_elem_computation */
      if (pf_u_prec != pf_u || pf_p_prec != pf_p || pf_d_prec != pf_d 
	  || pgtprec != pgt || pimprec != pim) {
	pmec = mat_elem(mat_elem_product(mat_elem_product(mat_elem_grad(pf_u),
						mat_elem_base(pf_p)),
			       mat_elem_base(pf_d)), pim, pgt);
	pf_u_prec = pf_u;
	pf_p_prec = pf_p;
	pf_d_prec = pf_d; pgtprec = pgt; pimprec = pim;
      }
      pmec->gen_compute(t, mf_u.linked_mesh().points_of_convex(cv));
      
      base_tensor::iterator p = t.begin();
      for (size_type i = 0; i < nbdof_d; i++) {
	size_type dof_d = mf_d.ind_dof_of_element(cv)[i];
	for (size_type j = 0; j < nbdof_p; j++) {
	  size_type dof_p = mf_p.ind_dof_of_element(cv)[j];
	  for (size_type l = 0; l < N; l++) {
	    // loop over derivation directions (d/dx, d/dy ..)
	    //	    for (size_type m = 0; m < N; m++) {
	      // loop over vector base function components (phi_x, phi_y ...)
	      for (size_type k = 0; k < nbdof_u; k++) {
		//		if (m == l) {
		/*
		  ssert(finite(DATA[dof_d])); 
		  
		  ssert(p < t.end());
		  
		  ssert(finite(*p));
		*/
		  size_type dof_u = mf_u.ind_dof_of_element(cv)[k];
		  B(dof_u*N+l, dof_p) += DATA[dof_d]*(*p);
		  //		}
		p++;
	      }
	      //	    }
	  } 
	}
      }
      if (p != t.end()) DAL_THROW(dal::internal_error, "internal error"); 
    }
  }


  /* ********************************************************************* */
  /*	Rigidity matrix for linear elasticity.                             */
  /* ********************************************************************* */

  template<class MAT, class VECT>
    void assembling_rigidity_matrix_for_linear_elasticity(MAT &RM,
							  const mesh_fem &mf, 
							  const mesh_fem &mfdata, 
							  const VECT &LAMBDA, const VECT &MU)
  { // à verifier

    size_type cv, nbd1, nbd2, N = mf.linked_mesh().dim();
    dal::bit_vector nn = mf.convex_index();
    base_tensor t;
    pfem pf1, pf2, pf1prec = NULL, pf2prec = NULL;
    pintegration_method pim, pimprec = 0;
    bgeot::pgeometric_trans pgt, pgtprec = NULL;
    pmat_elem_type pme; pmat_elem_computation pmec = 0;

    if (&(mf.linked_mesh()) != &(mfdata.linked_mesh()))
      DAL_THROW(std::invalid_argument,
		"This assembling procedure only works on a single mesh");
  
    for (cv << nn; cv != ST_NIL; cv << nn)
    {
      pf1 =     mf.fem_of_element(cv); nbd1 = pf1->nb_dof();
      pf2 = mfdata.fem_of_element(cv); nbd2 = pf2->nb_dof();
      pgt = mf.linked_mesh().trans_of_convex(cv);
      pim = mf.int_method_of_element(cv);
      if (pf1prec != pf1 || pf2prec != pf2 || pgtprec != pgt || pimprec != pim)
      {
	pme = mat_elem_product(mat_elem_product(mat_elem_grad(pf1),
						mat_elem_grad(pf1)), 
			       mat_elem_base(pf2));
	pmec = mat_elem(pme, pim, pgt);
	pf1prec = pf1; pf2prec = pf2; pgtprec = pgt; pimprec = pim;
      }
      pmec->gen_compute(t, mf.linked_mesh().points_of_convex(cv));
      base_tensor::iterator p = t.begin();
      
      size_type nbd = mf.nb_dof_of_element(cv);

      for (size_type r = 0; r < nbd2; r++)
      {
	size_type dof3 = mfdata.ind_dof_of_element(cv)[r];
	for (dim_type l = 0; l < N; l++)
	  for (size_type j = 0; j < nbd; j++)
	  {
	    size_type dof2 = mf.ind_dof_of_element(cv)[j];
	    
	    for (dim_type k = 0; k < N; k++)
	      for (size_type i = 0; i < nbd; i++, ++p)
	      {
		size_type dof1 = mf.ind_dof_of_element(cv)[i];
		
		if (dof1*N + k >= dof2*N + l)
		{
		  RM(dof1*N + k, dof2*N + l) += LAMBDA[dof3] * (*p);
		  RM(dof2*N + l, dof1*N + k) = RM(dof1*N + k, dof2*N + l);
		}
		
		if (dof1*N + l >= dof2*N + k)
		{
		  RM(dof1*N + l, dof2*N + k) += MU[dof3] * (*p);
		  RM(dof2*N + k, dof1*N + l) = RM(dof1*N + l, dof2*N + k);
		}

	      // cout << "matr elem : " << int(l) << " " << int(j) << " " << int(k) << " " << int(i) << " : " << *p << endl; getchar();
	      
		if (l == k && dof1 >= dof2)
		  for (size_type n = 0; n < N; ++n)
		  {
		    RM(dof1*N + n, dof2*N + n) += MU[dof3] * (*p);
		    RM(dof2*N + n, dof1*N + n) = RM(dof1*N + n, dof2*N + n);
		  }
		
	      }
	  }
      }
      if (p != t.end()) DAL_THROW(dal::internal_error, "internal error"); 
    }
  }

  /*
    assembles $\int_boundary{qu.v}$

    (if $u$ is a vector field of size $N$, $q$ is a square matrix $N\timesN$
    used by assem_general_boundary_conditions

    convention: Q is of the form 
       Q1_11 Q2_11 ..... Qn_11
       Q1_21 Q2_21 ..... Qn_21
       Q1_12 Q2_12 ..... Qn_12
       Q1_22 Q2_22 ..... Qn_22
    if  N = 2, and mf_d has n/N degree of freedom

    Q is a vector, so the matrice is assumed to be stored by columns
    (fortran style)
   */
  template<class MAT, class VECT>
  void assembling_boundary_qu_term(MAT &M, 
				   const mesh_fem &mf_u, size_type boundary, 
				   const mesh_fem &mf_d, const VECT &Q, dim_type N)
  {
    size_type cv;
    dal::bit_vector nn = mf_u.convex_index(), nf;
    base_tensor t;
    pfem pf_u, pf_d, pf_u_prec = NULL, pf_d_prec = NULL;
    pintegration_method pim, pimprec = 0;
    bgeot::pgeometric_trans pgt, pgtprec = NULL;
    pmat_elem_type pme; pmat_elem_computation pmec = 0;

    if (&(mf_u.linked_mesh()) != &(mf_d.linked_mesh()))
      DAL_THROW(std::invalid_argument,
		"This assembling procedure only works on a single mesh");

    for (cv << nn; cv != ST_NIL; cv << nn)
    {
      nf = mf_u.faces_of_convex_on_boundary(cv, boundary);
      if (nf.card() > 0)
      {
	size_type f, nbdof_u, nbdof_d;

	pf_u = mf_u.fem_of_element(cv); nbdof_u = pf_u->nb_dof();
	pf_d = mf_d.fem_of_element(cv); nbdof_d = pf_d->nb_dof();
	pgt = mf_u.linked_mesh().trans_of_convex(cv);
	pim = mf_u.int_method_of_element(cv);
	if (pf_u_prec != pf_u || pf_d_prec != pf_d || pgtprec!=pgt 
	    || pimprec != pim)
	{
	  pme = mat_elem_product(mat_elem_base(pf_d), 
				 mat_elem_product(mat_elem_base(pf_u),
						  mat_elem_base(pf_u)));
	  pmec = mat_elem(pme, pim, pgt);
	  pf_u_prec = pf_u; pf_d_prec = pf_d; pgtprec = pgt; pimprec = pim;
	}
	for (f << nf; f != ST_NIL; f << nf)
	{
	  pmec->gen_compute_on_face(t,mf_u.linked_mesh().points_of_convex(cv),
				    f);
	  base_tensor::iterator p = t.begin();
	  scalar_type vmax = bgeot::vect_norminf(t);

	  for (size_type j = 0; j < nbdof_u; j++) {
	    size_type dof_j = mf_u.ind_dof_of_element(cv)[j];
	    for (size_type i = 0; i < nbdof_u; i++) {
	      size_type dof_i = mf_u.ind_dof_of_element(cv)[i];
	      for (size_type id = 0; id < nbdof_d; id++) {

		size_type dof_d = mf_d.ind_dof_of_element(cv)[id];

		/* for every element of the matrix Q */
		for (int ii=0; ii < N; ii++) {
		  for (int jj=0; jj < N; jj++) {
		    /* get Q[ii][jj] for the degree of freedom 'dof_d' */
		    scalar_type data = Q[(jj*N+ii) + N*N*(dof_d)];

		    /* we filter out noise since this matrix can be used 
		       as a constraints matrix for dirichlet conditions,
		       noise may lead to 'fictive' dirichlet condition
		       (this is the case for ex. with laplace/PK(1,4)) */
		    if (data != 0. && vmax != .0 && (*p)/vmax > 1e-5) {
		      /*
			cerr << "QU : adding " << data << "*" << (*p) << " at dof_i=" << 
			dof_i << "*" << N << "+" << ii << ", dof_j=" << dof_i << "*" << 
			N << "+" << ii << endl;
		      */
		      M(dof_i*N+ii, dof_j*N+jj) += data* (*p);
		    }
		  }
		}
		p++;
	      }
	    }
	  }
	  if (p != t.end()) DAL_THROW(dal::internal_error, "internal error"); 
	}
      }
    }
  }


  template<class MAT, class VECT>
  void assembling_dirichlet_constraints(MAT &M, VECT &B, const mesh_fem &mf_u,
					size_type boundary,
					const mesh_fem &mf_rh,
				        const VECT &H, const VECT &R,
					dim_type N) {
    size_type cv;
    dal::bit_vector nn = mf_u.convex_index(), nf;
    pfem pf_u, pf_rh;
    
    assembling_boundary_qu_term(M, mf_u, boundary, mf_rh, H, N);
    assembling_Neumann_condition(B, mf_u, boundary, mf_rh, R, N);

    for (cv << nn; cv != ST_NIL; cv << nn) {
      nf = mf_u.faces_of_convex_on_boundary(cv, boundary);
      if (nf.card() > 0) {
	pf_u = mf_u.fem_of_element(cv); 
	pf_rh = mf_rh.fem_of_element(cv); 
	size_type f;
	for (f << nf; f != ST_NIL; f << nf) {
	  bgeot::pconvex_structure cvs_u = pf_u->structure();
	  bgeot::pconvex_structure cvs_rh = pf_rh->structure();
	  for (size_type i = 0; i < cvs_u->nb_points_of_face(f); ++i) {
	    size_type ind_u = cvs_u->ind_points_of_face(f)[i];
	    pdof_description tdof_u = pf_u->dof_types()[ind_u];
	    size_type dof_u = mf_u.ind_dof_of_element(cv)[ind_u];
	    for (size_type j = 0; j < cvs_rh->nb_points_of_face(f); ++j) {
	      size_type ind_rh = cvs_rh->ind_points_of_face(f)[j];
	      pdof_description tdof_rh = pf_rh->dof_types()[ind_rh];
	      if (tdof_u == tdof_rh) {
		if (bgeot::vect_dist2(pf_u->node_convex().points()[ind_u], 
				      pf_rh->node_convex().points()[ind_rh])
		    < 1.0E-7) {
		  // to be optimized (square root ..!)
		  for (size_type k = 0; k < pf_u->nb_dof(); ++k) {
		    size_type dof_k = mf_u.ind_dof_of_element(cv)[k];
		    for (int ii=0; ii < N; ii++)
		      for (int jj=0; jj < N; jj++)
			M(dof_u*N+ii, dof_k*N+jj) = 0.0;
		  }

		  size_type dof_rh = mf_rh.ind_dof_of_element(cv)[ind_rh];
		  for (int ii=0; ii < N; ii++) {
		    for (int jj=0; jj < N; jj++)
		      M(dof_u*N+ii, dof_u*N+jj) = H[(jj*N+ii) + N*N*(dof_rh)];
		    B[dof_u*N+ii] = R[dof_rh*N+ii];
		  }
		}
	      }      
	    }
	  }
	}
      }
    }
  }


  /* ********************************************************************* */
  /*	Mass matrix.                                                       */
  /* ********************************************************************* */

  template<class MATRM, class MESH_FEM>
    void mass_matrix(MATRM &M, const MESH_FEM &mf1, const MESH_FEM &mf2, dim_type N)
  {
    size_type cv, nbd1, nbd2;
    dal::bit_vector nn = mf1.convex_index();
    base_tensor t;
    pfem pf1, pf1prec = 0, pf2, pf2prec = 0;
    pintegration_method pim, pimprec = 0;
    bgeot::pgeometric_trans pgt, pgtprec = NULL;
    pmat_elem_type pme; pmat_elem_computation pmec = 0;
    // M(0,0) = 1.0;  ??

    if (&(mf1.linked_mesh()) != &(mf2.linked_mesh()))
      DAL_THROW(std::invalid_argument,
		"This assembling procedure only works on a single mesh");

    for (cv << nn; cv != ST_NIL; cv << nn)
    {
      pf1 = mf1.fem_of_element(cv); nbd1 = pf1->nb_dof();
      pf2 = mf2.fem_of_element(cv); nbd2 = pf2->nb_dof();
      pgt = mf1.linked_mesh().trans_of_convex(cv);
      pim = mf1.int_method_of_element(cv);
      if (pf1prec != pf1 || pf2prec != pf2 || pgtprec != pgt || pimprec != pim)
      {
	pme = mat_elem_product(mat_elem_base(pf1), mat_elem_base(pf2));
	pmec = mat_elem(pme, pim, pgt);
	pf1prec = pf1; pf2prec = pf2; pgtprec = pgt; pimprec = pim;
      }
      pmec->gen_compute(t, mf1.linked_mesh().points_of_convex(cv));

      // cout << "t = " << t << endl;
      
      base_tensor::iterator p = t.begin();
      for (size_type i = 0; i < nbd2; i++)
      {
	size_type dof2 = mf2.ind_dof_of_element(cv)[i];
	// cout << "cv = " << cv << " dof2 = " << dof2 << endl;
	for (size_type j = 0; j < nbd1; j++, ++p)
	{
	  size_type dof1 = mf1.ind_dof_of_element(cv)[j];
	  // cout << "dof1 = " << dof1 << " dof2 = " << dof2 << endl;
	  for (size_type k = 0; k < N; k++)
	    M(dof1*N + k, dof2*N + k) += (*p);
	}
      }
      if (p != t.end()) DAL_THROW(dal::internal_error, "internal error"); 
    }
  }

  template<class MATRM, class MESH_FEM>
    inline void mass_matrix(MATRM &M, const MESH_FEM &mf, dim_type N)
    { mass_matrix(M, mf, mf, N); }


  /* ********************************************************************* */
  /*	volumic source term.                                               */
  /* ********************************************************************* */

  template<class VECT1, class VECT2>
    void assembling_volumic_source_term(VECT1 &B, const mesh_fem &mf,
					const mesh_fem &mfdata, const VECT2 &F, dim_type N)
  {
    size_type cv, nbd1, nbd2;
    dal::bit_vector nn = mf.convex_index();
    base_tensor t;
    pfem pf1, pf2, pf1prec = NULL, pf2prec = NULL;
    pintegration_method pim, pimprec = 0;
    bgeot::pgeometric_trans pgt, pgtprec = NULL;
    pmat_elem_type pme; pmat_elem_computation pmec = 0;

    if (&(mf.linked_mesh()) != &(mfdata.linked_mesh()))
      DAL_THROW(std::invalid_argument,
		"This assembling procedure only works on a single mesh");

    for (cv << nn; cv != ST_NIL; cv << nn)
    {
      pf1 =     mf.fem_of_element(cv); nbd1 = pf1->nb_dof();
      pf2 = mfdata.fem_of_element(cv); nbd2 = pf2->nb_dof();
      pgt = mf.linked_mesh().trans_of_convex(cv);
      pim = mf.int_method_of_element(cv);
      if (pf1prec != pf1 || pf2prec != pf2 || pgtprec != pgt || pimprec != pim)
      {
	pme = mat_elem_product(mat_elem_base(pf1), mat_elem_base(pf2));
	pmec = mat_elem(pme, pim, pgt);
	pf1prec = pf1; pf2prec = pf2; pgtprec = pgt; pimprec = pim;
      }
      pmec->gen_compute(t, mf.linked_mesh().points_of_convex(cv));
      base_tensor::iterator p = t.begin();
      for (size_type i = 0; i < nbd2; i++)
      {
	size_type dof2 = mfdata.ind_dof_of_element(cv)[i];
	for (size_type j = 0; j < nbd1; j++, ++p)
	{
	  size_type dof1 = mf.ind_dof_of_element(cv)[j];
	  for (size_type k = 0; k < N; k++) B[dof1*N + k] += F[dof2*N+k]*(*p);
	}
      }
      if (p != t.end()) DAL_THROW(dal::internal_error, "internal error"); 
    }
  }

  /* ********************************************************************* */
  /*	Dirichlet Condition.                                               */
  /* ********************************************************************* */

  template<class MATRM, class VECT1, class VECT2>
    void assembling_Dirichlet_condition(MATRM &RM, VECT1 &B,
					const mesh_fem &mf,
					size_type boundary,
					const VECT2 &F, dim_type N)
  { /* Y-a-il un moyen plus performant ? */
    // Marche uniquement pour des ddl de lagrange.
    size_type cv;
    dal::bit_vector nn = mf.convex_index();
    dal::bit_vector nndof = mf.dof_on_boundary(boundary);
    pfem pf1;

    for (cv << nn; cv != ST_NIL; cv << nn)
    {
      pf1 = mf.fem_of_element(cv);
      pdof_description ldof = lagrange_dof(pf1->dim());
      size_type nbd = pf1->nb_dof();
      for (size_type i = 0; i < nbd; i++)
      {
	size_type dof1 = mf.ind_dof_of_element(cv)[i];
	if (nndof.is_in(dof1) && pf1->dof_types()[i] == ldof)
	{
	  // cout << "dof : " << i << endl;
	  for (size_type j = 0; j < nbd; j++)
	  {
	    size_type dof2 = mf.ind_dof_of_element(cv)[j];
	    for (size_type k = 0; k < N; ++k)
	      for (size_type l = 0; l < N; ++l)
	      {
		if (!(nndof.is_in(dof2)) &&
		    dof_compatibility(pf1->dof_types()[j],
				      lagrange_dof(pf1->dim())))
		  B[dof2*N+k] -= RM(dof2*N+k, dof1*N+l) * F[dof1*N+l];
		RM(dof2*N+k, dof1*N+l) = RM(dof1*N+l, dof2*N+k) = 0;
	      }
	  }
	  for (size_type k = 0; k < N; ++k)
	  { RM(dof1*N+k, dof1*N+k) = 1; B[dof1*N+k] = F[dof1*N+k]; }
	}
      }
    }
  }

//   template<class MATD, class MATG, class VECT>
//     size_type treat_Dirichlet_condition_new(const MATD &D, MATG &G,
// 					    const VECT &UD, VECT &UDD) {
//     typedef typename bgeot::VECT_TO_BE_USED<bgeot::linalg_traits<MATD>::storage_type
//       TEMP_VECT;
//     size_type nbd = D.ncols(), nbase = 0, nbr = D.nrows();
//     TEMP_VECT aux(nbr), e(nbd), f(nbd);
//     dal::dynamic_array<TEMP_VECT> base_img;
//     dal::dynamic_array<TEMP_VECT> base_img_inv;
//     size_type nb_bimg = 0;

//     // to be continued

//     // First, detection of null columns of D, and already orthogonals 
//     // vectors of the image of D.
//     dal::bit_vector nn;
//     for (size_type i = 0; i < nbd; ++i) {
//       bgeot::clear(e); e[i] = 1.0; bgeot::copy(e, f);
//       bgeot::mult(D, e, aux);
//       if (bgeot::norm2(aux) < 1.0E-8) { //à scaler sur l'ensemble de D ...
// 	G(nbase++, i) = 1.0; nn[i] = true;
//       }
//       else {
// 	bool good = true;
// 	for (size_type j = 0; j < nb_bimg; ++j)
// 	  if (dal::abs(bgeot::vect_sp(aux, base_img[j])) > 1.0E-16)
// 	    { good = false; break; }
// 	if (good) {
// 	  scalar_type n = bgeot::norm2(aux);
// 	  f /= n; aux /= n;
// 	  base_img_inv[nb_bimg] = TEMP_VECT();
// 	  base_img_inv[nb_bimg] = f;
// 	  //	  cerr << "ajout de " << aux << "\n";
// 	  aux.clean(1.0E-18);
// 	  base_img[nb_bimg++] = aux; nn[i] = true;
// 	}
//       }
//     }
//     size_type nb_triv_base = nbase;

//     for (size_type i = 0; i < nbd; ++i)
//       if (!(nn[i])) {
// 	e.clear(); e[i] = 1.0; f.clear(); f[i] = 1.0;
// 	aux = D*e;
// 	for (size_type j = 0; j < nb_bimg; ++j) { 
// 	  scalar_type c = bgeot::vect_sp(aux, base_img[j]);
// 	  //	  if (dal::abs(c > 1.0E-6) { // à scaler sur l'ensemble de D ...
// 	  if (c != 0.) {
// 	    aux -= base_img[j] * c;
// 	    f -= base_img_inv[j] * c;
// 	  }
// 	}
// 	//	cerr << "norm2(aux)= " << bgeot::vect_norm2(aux) << "\n";
// 	if (bgeot::vect_norm2(aux) < 1.0E-8) { // à scaler sur l'ensemble de D ...
// 	  G.row(nbase++) = f;
// 	}
// 	else {
// 	  scalar_type n = bgeot::vect_norm2(aux);
// 	  f /= n; aux /= n;
// 	  base_img_inv[nb_bimg] = f;
// 	  base_img[nb_bimg++] = aux;
// 	  f.clean(1.0E-18); aux.clean(1.0E-18);
// 	  //	  cerr << "ajout de " << aux << "\n";
// 	}
// 	e[i] = 0.0;
//       }

//     // Compute a solution in UDD
//     UDD.fill(0.0);
//     for (size_type i = 0; i < nb_bimg; ++i) {
//       scalar_type c = bgeot::vect_sp(base_img[i], UD);
//       UDD += base_img_inv[i].full() * c;
//     }

//     // Orthogonalisation of the basis of the kernel of D.
//     for (size_type i = nb_triv_base + 1; i < nbase; ++i) {
//       for (size_type j = nb_triv_base; j < i; ++j) {
// 	scalar_type c = bgeot::vect_sp(G.row(i), G.row(j));
// 	if (c != 0.)
// 	  G.row(i) -= G.row(j) * c;
//       }
//     }

//     // projection of UDD on the orthogonal to the kernel.
//     for (size_type j = nb_triv_base; j < nbase; ++j) {
//       scalar_type c = bgeot::vect_sp(G.row(j), UDD);
//       if (c != 0.)
// 	UDD -= G.row(j).full() * c;
//     }

//     // Test ...
//     if (bgeot::vect_norm2(D * UDD - UD) > 1.0E-12)
//       cerr << "Dirichlet condition not well inverted\n";

//     return nbase;

//   }
  

  template<class MATD, class VECT>
    size_type treat_Dirichlet_condition(const MATD &D, MATD &G,
				   const VECT &UD, VECT &UDD) {
    // To be finalized.
    //  . In order to be used with any sparse matrix type
    //  . transpose the result and give the effective dimension of
    //    the kernel
    //  . Compute the ctes / D.
    //  . Optimization (suppress temporary ...). 
    //  . Verify sizes of data

    // Build an orthogonal basis of the kernel of D in G, gives
    // the solution of minimal norm of D*U = UD in UDD and
    // return the dimension of the kernel. The function is based
    // on a Gramm-Schmidt algorithm.

    size_type nbd = D.ncols(), nbase = 0;
    bgeot::svector<scalar_type> aux(D.nrows()), e(nbd);
    bgeot::svector<scalar_type> f(nbd);
    dal::dynamic_array<bgeot::svector<scalar_type> > base_img;
    dal::dynamic_array<bgeot::svector<scalar_type> > base_img_inv;
    size_type nb_bimg = 0;

    // First, detection of null columns of D, and already orthogonals 
    // vectors of the image of D.
    dal::bit_vector nn;
    for (size_type i = 0; i < nbd; ++i) {
      e.clear(); e[i] = 1.0; f.clear(); f[i] = 1.0;
      aux = D*e;
      if (bgeot::vect_norm2(aux) < 1.0E-8) { //à scaler sur l'ensemble de D ...
	G(nbase++, i) = 1.0; nn[i] = true;
      }
      else {
	bool good = true;
	for (size_type j = 0; j < nb_bimg; ++j)
	  if (dal::abs(bgeot::vect_sp(aux, base_img[j])) > 1.0E-16)
	    { good = false; break; }
	if (good) {
	  scalar_type n = bgeot::vect_norm2(aux);
	  f /= n; aux /= n;
	  base_img_inv[nb_bimg] = f;
	  //	  cerr << "ajout de " << aux << "\n";
	  aux.clean(1.0E-18);
	  base_img[nb_bimg++] = aux; nn[i] = true;
	}
      }
    }
    size_type nb_triv_base = nbase;

    for (size_type i = 0; i < nbd; ++i)
      if (!(nn[i])) {
	e.clear(); e[i] = 1.0; f.clear(); f[i] = 1.0;
	aux = D*e;
	for (size_type j = 0; j < nb_bimg; ++j) { 
	  scalar_type c = bgeot::vect_sp(aux, base_img[j]);
	  //	  if (dal::abs(c > 1.0E-6) { // à scaler sur l'ensemble de D ...
	  if (c != 0.) {
	    aux -= base_img[j] * c;
	    f -= base_img_inv[j] * c;
	  }
	}
	//	cerr << "norm2(aux)= " << bgeot::vect_norm2(aux) << "\n";
	if (bgeot::vect_norm2(aux) < 1.0E-8) { // à scaler sur l'ensemble de D ...
	  G.row(nbase++) = f;
	}
	else {
	  scalar_type n = bgeot::vect_norm2(aux);
	  f /= n; aux /= n;
	  base_img_inv[nb_bimg] = f;
	  base_img[nb_bimg++] = aux;
	  f.clean(1.0E-18); aux.clean(1.0E-18);
	  //	  cerr << "ajout de " << aux << "\n";
	}
	e[i] = 0.0;
      }

    // Compute a solution in UDD
    UDD.fill(0.0);
    for (size_type i = 0; i < nb_bimg; ++i) {
      scalar_type c = bgeot::vect_sp(base_img[i], UD);
      UDD += base_img_inv[i].full() * c;
    }

    // Orthogonalisation of the basis of the kernel of D.
    for (size_type i = nb_triv_base + 1; i < nbase; ++i) {
      for (size_type j = nb_triv_base; j < i; ++j) {
	scalar_type c = bgeot::vect_sp(G.row(i), G.row(j));
	if (c != 0.)
	  G.row(i) -= G.row(j) * c;
      }
    }

    // projection of UDD on the orthogonal to the kernel.
    for (size_type j = nb_triv_base; j < nbase; ++j) {
      scalar_type c = bgeot::vect_sp(G.row(j), UDD);
      if (c != 0.)
	UDD -= G.row(j).full() * c;
    }

    // Test ...
    if (bgeot::vect_norm2(D * UDD - UD) > 1.0E-12)
      cerr << "Dirichlet condition not well inverted\n";

    return nbase;
  }
  

  
  /* ********************************************************************* */
  /*	Neumann Condition.                                                 */
  /* ********************************************************************* */

  template<class VECT1, class VECT2>
    void assembling_Neumann_condition(VECT1 &B, const mesh_fem &mf,
				      size_type boundary, const mesh_fem &mfdata, const VECT2 &F, dim_type N)
  {
    size_type cv, nbd1, nbd2, f;
    dal::bit_vector nn = mf.convex_index(), nf;
    base_tensor t;
    pfem pf1, pf2, pf1prec = NULL, pf2prec = NULL;
    pintegration_method pim, pimprec = 0;
    bgeot::pgeometric_trans pgt, pgtprec = NULL;
    pmat_elem_type pme; pmat_elem_computation pmec = 0;

    if (&(mf.linked_mesh()) != &(mfdata.linked_mesh()))
      DAL_THROW(std::invalid_argument,
		"This assembling procedure only works on a single mesh");
  
    for (cv << nn; cv != ST_NIL; cv << nn)
    {
      // for (int h = 0; h < mf.linked_mesh().nb_points_of_convex(cv); ++h)
      // cout << "Point " << h << " of cv " << cv << " : " 
      //     << mf.linked_mesh().points_of_convex(cv)[h] << endl;

      // for (f = 0; f < mf.linked_mesh().structure_of_convex(cv)->nb_faces(); ++f)
      // for (int h = 0; h < mf.linked_mesh().structure_of_convex(cv)->nb_points_of_face(f); ++h)
      //  cout << "Point " << h << " of cv " << cv << " on face " << f << " (" << mf.linked_mesh().structure_of_convex(cv)->ind_points_of_face(f)[h] << ") : " 
      //     << mf.linked_mesh().points_of_convex(cv)[mf.linked_mesh().structure_of_convex(cv)->ind_points_of_face(f)[h]] << endl;

      nf = mf.faces_of_convex_on_boundary(cv, boundary);
      if (nf.card() > 0)
      {
	pf1 =     mf.fem_of_element(cv); nbd1 = pf1->nb_dof();
	pf2 = mfdata.fem_of_element(cv); nbd2 = pf2->nb_dof();
	pgt = mf.linked_mesh().trans_of_convex(cv);
	pim = mf.int_method_of_element(cv);
	if (pf1prec != pf1 || pf2prec != pf2 || pgtprec!=pgt || pimprec != pim)
	{
	  pme = mat_elem_product(mat_elem_base(pf1), mat_elem_base(pf2));
	  pmec = mat_elem(pme, pim, pgt);
	  pf1prec = pf1; pf2prec = pf2; pgtprec = pgt; pimprec = pim;
	}
	for (f << nf; f != ST_NIL; f << nf)
	{
	  // cout << "cv = " << cv << " f = " << f << endl;

	  pmec->gen_compute_on_face(t,mf.linked_mesh().points_of_convex(cv),f);
	  // cout << "t = " << t << endl;
	  base_tensor::iterator p = t.begin();
	  for (size_type i = 0; i < nbd2; i++)
	  {
	    size_type dof2 = mfdata.ind_dof_of_element(cv)[i];
	    for (size_type j = 0; j < nbd1; j++, ++p)
	    {
	      size_type dof1 = mf.ind_dof_of_element(cv)[j];
	      for (size_type k = 0; k < N; k++) {
		B[dof1*N + k] += F[dof2*N+k]*(*p);
	      }
	    }
	  }
	  if (p != t.end()) DAL_THROW(dal::internal_error, "internal error"); 
	}
      }
    }
  }

  
}  /* end of namespace getfem.                                             */


#endif /* __GETFEM_ASSEMBLING_H  */
