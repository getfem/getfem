/* -*- c++ -*- (enables emacs c++ mode)                                    */
/* *********************************************************************** */
/*                                                                         */
/* Library :  GEneric Tool for Finite Element Methods (getfem)             */
/* File    :  getfem_assembling.h : assemble linear system for fem.        */
/*     									   */
/* Date : June 17, 2002.                                                   */
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

#ifndef __GETFEM_DERIVATIVES_H
#define __GETFEM_DERIVATIVES_H

//#include <bgeot_geometric_trans.h>
#include <getfem_mesh_fem.h>

// module à tester ...

namespace getfem
{
  /*
    mf_target should be a lagrange discontinous element
    does not work with vectorial elements. ... to be done ...
  */

  template<class VECT>
  void compute_gradient(mesh_fem &mf, mesh_fem &mf_target,
			const VECT &U, VECT &V, dim_type Q)
  {
    size_type cv;
    size_type N = mf.linked_mesh().dim();

    assert(&mf.linked_mesh() == &mf_target.linked_mesh());

    base_matrix G, val;
    base_vector coeff;
 
    dal::bit_vector nn = mf.convex_index();
      
    pgeotrans_precomp pgp;
    pfem_precomp pfp;
    pfem pf, pf_target, pf_old = NULL, pf_targetold = NULL;
    bgeot::pgeometric_trans pgt;

    for (cv << nn; cv != ST_NIL; cv << nn) {
      pf = mf.fem_of_element(cv);
      pf_target = mf_target.fem_of_element(cv);
      assert(pf_target->is_equivalent());
      assert(pf_target->is_lagrange());

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
      base_matrix a(N, pgt->nb_points());
      base_matrix pc(pgt->nb_points() , P);
      base_matrix grad(N, P), TMP1(P,P), B0(P,N), B1(1, N), CS(P,P);
      size_type nbpt = 0, i;
      
      /* TODO: prendre des iterateurs pour faire la copie */
      for (size_type j = 0; j < pgt->nb_points(); ++j) // à optimiser !!
	for (size_type i = 0; i < N; ++i)
	  a(i,j) = mf.linked_mesh().points_of_convex(cv)[j][i];
      
      coeff.resize(pf->nb_dof());
      val.resize(pf->target_dim(), P);
      B1.resize(pf->target_dim(), N);
      
      for (size_type j = 0; j < pf_target->nb_dof(); ++j) {
	if (!pgt->is_linear() || j == 0) {
	  // computation of the pseudo inverse
	  bgeot::mat_product(a, pgp->grad(j), grad);
	  if (P != N) {
	    bgeot::mat_product_tn(grad, grad, CS);
	    bgeot::mat_inv_cholesky(CS, TMP1);
	    bgeot::mat_product_tt(CS, grad, B0);
	  }
	  else {
	    bgeot::mat_gauss_inverse(grad, TMP1); B0 = grad;
	  }
	}

	assert(pf_target->target_dim() == 1); // !!
	for (size_type q = 0; q < Q; ++q) {
	  for (size_type l = 0; l < pf->nb_dof(); ++l)
	    coeff[l] = U[mf.ind_dof_of_element(cv)[l] * Q + q ];
	  pf->interpolation_grad(pfp, j, G, coeff, val);
	  bgeot::mat_product(val, B0, B1);

	  for (size_type l = 0; l < N; ++l)
	    V[mf_target.ind_dof_of_element(cv)[j]*Q*N+q*N+l] = B1(0, l);
	}
      }
	
    }
  }







}

#endif
