/* -*- c++ -*- (enables emacs c++ mode)                                    */
/* *********************************************************************** */
/*                                                                         */
/* Library :  GEneric Tool for Finite Element Methods (getfem)             */
/* File    :  getfem_derivatives.h : .                                     */
/*     									   */
/* Date : June 17, 2002.                                                   */
/* Authors : Yves Renard, Yves.Renard@gmm.insa-tlse.fr                     */
/*           Julien Pommier, pommier@gmm.insa-tlse.fr                      */
/*                                                                         */
/* *********************************************************************** */
/*                                                                         */
/* Copyright (C) 2002  Yves Renard.                                        */
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

#ifndef __GETFEM_DERIVATIVES_H
#define __GETFEM_DERIVATIVES_H

#include <getfem_mesh_fem.h>

namespace getfem
{
  /*
    mf_target should be a lagrange discontinous element
    does not work with vectorial elements. ... to be done ...

    mf_target may have the same Qdim than mf, or it may
    be a scalar mesh_fem, in which case the derivatives are stored in
    the order: DxUx,DyUx,DzUx,DxUy,DyUy,...

    in any case, the size of V should be (mf_target.nbdof/mf_target.qdim) *
    (mf.qdim)*N elements (this is not checked by the function!)
  */
  template<class VECT>
  void compute_gradient(mesh_fem &mf, mesh_fem &mf_target,
			const VECT &U, VECT &V)
  {
    size_type cv;
    size_type N = mf.linked_mesh().dim();
    size_type qdim = mf.get_qdim();
    size_type target_qdim = mf_target.get_qdim();


    if (&mf.linked_mesh() != &mf_target.linked_mesh())
      DAL_THROW(std::invalid_argument, "meshes are different.");

    if (target_qdim != qdim && target_qdim != 1) {
      DAL_THROW(std::invalid_argument, "invalid Qdim for gradient mesh_fem");
    }

    base_matrix G, val;
    base_vector coeff;
 
    dal::bit_vector nn = mf.convex_index();
      
    pgeotrans_precomp pgp = NULL;
    pfem_precomp pfp = NULL;
    pfem pf, pf_target, pf_old = NULL, pf_targetold = NULL;
    bgeot::pgeometric_trans pgt;

    for (cv << nn; cv != ST_NIL; cv << nn) {
      pf = mf.fem_of_element(cv);
      pf_target = mf_target.fem_of_element(cv);
      if (pf_target->need_G() || !(pf_target->is_lagrange()))
	DAL_THROW(std::invalid_argument, 
		  "finite element target not convenient");
      
      transfert_to_G(G, mf.linked_mesh().points_of_convex(cv));

      pgt = mf.linked_mesh().trans_of_convex(cv);
      if (pf_targetold != pf_target) {
        pgp = geotrans_precomp(pgt, pf_target->node_tab());
      }
      pf_targetold = pf_target;

      if (pf_old != pf) {
	pfp = fem_precomp(pf, pf_target->node_tab());
      }
      pf_old = pf;

      size_type P = pgt->structure()->dim(); /* dimension of the convex.*/
      // base_matrix a(N, pgt->nb_points());
      base_matrix grad(N, P), TMP1(P,P), B0(P,N), B1(1, N), CS(P,P);
      base_tensor t;
      
      coeff.resize(pf->nb_dof());
      val.resize(pf->target_dim(), P);
      B1.resize(pf->target_dim(), N);
      
      for (size_type j = 0; j < pf_target->nb_dof(); ++j) {
	if (!pgt->is_linear() || j == 0) {
	  // computation of the pseudo inverse
	  bgeot::mat_product_tt(pgp->grad(j), G, grad);
	  if (P != N) {
	    bgeot::mat_product_nt(grad, grad, CS);
	    bgeot::mat_inv_cholesky(CS, TMP1);
	    bgeot::mat_product_tn(grad, CS, B0);
	  }
	  else {
	    bgeot::mat_gauss_inverse(grad, TMP1); B0 = grad;
	  }
	}

	if (pf_target->target_dim() != 1 || pf->target_dim() != 1)
	  DAL_THROW(to_be_done_error, "vectorial gradient, to be done ... ");

	pf->real_grad_base_value(pgp, pfp, j, G, B0, t);

	for (size_type q = 0; q < qdim; ++q) {
	  for (size_type l = 0; l < pf->nb_base(); ++l)
	    coeff[l] = U[mf.ind_dof_of_element(cv)[l*qdim] + q ];

	  base_tensor::const_iterator it = t.begin();
	  B1.fill(0.0);
	  for (size_type l = 0; l < N; ++l)
	    for (size_type i = 0; i < pf->nb_base(); ++i, ++it)
	      B1(0, l) += *it * coeff[i];

	  if (it != t.end()) DAL_THROW(internal_error, "internal_error");

	  if (target_qdim != 1) {
	    for (size_type l = 0; l < N; ++l)
	      V[mf_target.ind_dof_of_element(cv)[j*qdim+q]*N+l] = B1(0, l);
	  } else {
	    for (size_type l = 0; l < N; ++l)
	      V[(mf_target.ind_dof_of_element(cv)[j]*qdim + q)*N+l] = B1(0, l);
	  }
	  
	}
      }
    }
  }

//   template<class VECT>
//   void compute_gradient_old(mesh_fem &mf, mesh_fem &mf_target,
// 			    const VECT &U, VECT &V, dim_type Q)
//   {
//     size_type cv;
//     size_type N = mf.linked_mesh().dim();
    
//     if (&mf.linked_mesh() != &mf_target.linked_mesh())
//       DAL_THROW(std::invalid_argument, "meshes are different.");

//     base_matrix G, val;
//     base_vector coeff;
 
//     dal::bit_vector nn = mf.convex_index();
      
//     pgeotrans_precomp pgp = NULL;
//     pfem_precomp pfp = NULL;
//     pfem pf, pf_target, pf_old = NULL, pf_targetold = NULL;
//     bgeot::pgeometric_trans pgt;

//     for (cv << nn; cv != ST_NIL; cv << nn) {
//       pf = mf.fem_of_element(cv);
//       pf_target = mf_target.fem_of_element(cv);
//       if (!(pf_target->is_equivalent()) || !(pf_target->is_lagrange()))
// 	DAL_THROW(std::invalid_argument, 
// 		  "finite element target not convenient");
//       if (!(pf->is_equivalent())) 
// 	transfert_to_G(G, mf.linked_mesh().points_of_convex(cv));

//       pgt = mf.linked_mesh().trans_of_convex(cv);
//       if (pf_targetold != pf_target) {
//         pgp = geotrans_precomp(pgt, pf_target->node_tab());
//       }
//       pf_targetold = pf_target;

//       if (pf_old != pf) {
// 	pfp = fem_precomp(pf, pf_target->node_tab());
//       }
//       pf_old = pf;

//       size_type P = pgt->structure()->dim(); /* dimension of the convex.*/
//       base_matrix a(N, pgt->nb_points());
//       base_matrix grad(N, P), TMP1(P,P), B0(P,N), B1(1, N), CS(P,P);
      
//       /* TODO: prendre des iterateurs pour faire la copie */
//       // utiliser transfert_to_G ?
//       for (size_type j = 0; j < pgt->nb_points(); ++j) // à optimiser !!
// 	for (size_type i = 0; i < N; ++i)
// 	  a(i,j) = mf.linked_mesh().points_of_convex(cv)[j][i];
      
//       coeff.resize(pf->nb_dof());
//       val.resize(pf->target_dim(), P);
//       B1.resize(pf->target_dim(), N);
      
//       for (size_type j = 0; j < pf_target->nb_dof(); ++j) {
// 	if (!pgt->is_linear() || j == 0) {
// 	  // computation of the pseudo inverse
// 	  bgeot::mat_product(a, pgp->grad(j), grad);
// 	  if (P != N) {
// 	    bgeot::mat_product_tn(grad, grad, CS);
// 	    bgeot::mat_inv_cholesky(CS, TMP1);
// 	    bgeot::mat_product_tt(CS, grad, B0);
// 	  }
// 	  else {
// 	    bgeot::mat_gauss_inverse(grad, TMP1); B0 = grad;
// 	  }
// 	}

// 	if (pf_target->target_dim() != 1)
// 	  DAL_THROW(to_be_done_error, "vectorial gradient, to be done ... ");

// 	for (size_type q = 0; q < Q; ++q) {
// 	  for (size_type l = 0; l < pf->nb_dof(); ++l)
// 	    coeff[l] = U[mf.ind_dof_of_element(cv)[l] * Q + q ];
// 	  pf->interpolation_grad(pfp, j, G, pgt, coeff, val);
// 	  bgeot::mat_product(val, B0, B1);

// 	  for (size_type l = 0; l < N; ++l)
// 	    V[mf_target.ind_dof_of_element(cv)[j]*Q*N+q*N+l] = B1(0, l);
// 	}
//       }
	
//     }
//   }
}

#endif
