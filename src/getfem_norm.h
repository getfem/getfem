/* -*- c++ -*- (enables emacs c++ mode)                                    */
/* *********************************************************************** */
/*                                                                         */
/* Library :  GEneric Tool for Finite Element Methods (getfem)             */
/* File    :  getfem_norm.h: computes various norms on pde solutions.      */
/*     									   */
/* Date : November 17, 2000.                                               */
/* Authors : Yves Renard, Yves.Renard@gmm.insa-tlse.fr                     */
/*           Julien Pommier, pommier@gmm.insa-tlse.fr                      */
/*                                                                         */
/* *********************************************************************** */
/*                                                                         */
/* Copyright (C) 2000-2002  Yves Renard.                                   */
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
 

#ifndef GETFEM_NORM_H__
#define GETFEM_NORM_H__

#include <getfem_mesh_fem.h>
#include <getfem_mat_elem.h>
namespace getfem
{



  template<class MESH_FEM, class VECT>
  scalar_type L2_norm(MESH_FEM &mf, const VECT &U, size_type N, const dal::bit_vector &cvlst)
  { /* optimisable */
    scalar_type no = 0.0;
    dal::dynamic_array<base_vector, 2> vval;
    base_tensor t;
    pfem pf1, pf1prec = NULL;
    pintegration_method pim, pimprec = 0;

    bgeot::pgeometric_trans pgt, pgtprec = NULL;
    pmat_elem_type pme; pmat_elem_computation pmec = 0;
    
    DAL_WARNING(3, "obsolete function (not qdim aware) - use asm_L2_norm");
    for (dal::bv_visitor cv(cvlst); !cv.finished(); ++cv)
      {
	pf1 =     mf.fem_of_element(cv);
	pgt = mf.linked_mesh().trans_of_convex(cv);
	pim = mf.int_method_of_element(cv);
	size_type nbd = mf.nb_dof_of_element(cv);
	if (pf1prec != pf1 || pgtprec != pgt || pimprec != pim)
	  {
	    pme = mat_elem_product(mat_elem_base(pf1), mat_elem_base(pf1));
	    pmec = mat_elem(pme, pim, pgt);
	    pf1prec = pf1; pgtprec = pgt; pimprec = pim;
	  }
	pmec->gen_compute(t, mf.linked_mesh().points_of_convex(cv), cv);
	base_tensor::iterator p = t.begin();

	for (size_type i = 0; i < nbd; i++)
	  { 
	    size_type dof1 = mf.ind_dof_of_element(cv)[i];
	    if (vval[i].size() != N) vval[i] = base_vector(N); 
	    for (size_type k = 0; k < N; k++) (vval[i])[k] = U[dof1*N+k];
	  }

	for (size_type i = 0; i < nbd; i++)
	  for (size_type j = 0; j < nbd; j++, ++p)
	    no += bgeot::vect_sp(vval[i], vval[j]) * (*p);
      
      }
    return sqrt(no);
  }

  template<class MESH_FEM, class VECT>
  scalar_type L2_norm(MESH_FEM &mf, const VECT &U, size_type N)
  {
    return L2_norm<MESH_FEM,VECT>(mf, U, N, mf.convex_index());
  }

  template<class MESH_FEM, class VECT>
  scalar_type H1_semi_norm(MESH_FEM &mf, const VECT &U, size_type N, const dal::bit_vector& cvlst)
  { /* optimisable */
    size_type NN = mf.linked_mesh().dim();
    scalar_type no = 0.0;
    dal::dynamic_array<base_vector, 2> vval;
    base_tensor t;
    pfem pf1, pf1prec = NULL;
    pintegration_method pim, pimprec = 0;
    bgeot::pgeometric_trans pgt, pgtprec = NULL;
    pmat_elem_type pme; pmat_elem_computation pmec = 0;

    DAL_WARNING(3, "obsolete function (not qdim aware) - use asm_H1_semi_norm");
    for (dal::bv_visitor cv(cvlst); !cv.finished(); ++cv)
      {
	pf1 =     mf.fem_of_element(cv);
	pgt = mf.linked_mesh().trans_of_convex(cv);
	pim = mf.int_method_of_element(cv);
	size_type nbd = mf.nb_dof_of_element(cv);
	if (pf1prec != pf1 || pgtprec != pgt || pimprec != pim)
	  {
	    pme = mat_elem_product(mat_elem_grad(pf1), mat_elem_grad(pf1));
	    pmec = mat_elem(pme, pim, pgt);
	    pf1prec = pf1; pgtprec = pgt; pimprec = pim;
	  }
	pmec->gen_compute(t, mf.linked_mesh().points_of_convex(cv), cv);
	base_tensor::iterator p = t.begin();
	for (size_type i = 0; i < nbd; i++)
	  { 
	    size_type dof1 = mf.ind_dof_of_element(cv)[i];
	    if (vval[i].size() != N) vval[i] = base_vector(N); 
	    for (size_type k = 0; k < N; k++) (vval[i])[k] = U[dof1*N+k];
	  }
	for (size_type l = 0; l < NN; l++)
	  for (size_type i = 0; i < nbd; i++)
	    for (size_type k = 0; k < NN; k++)
	      for (size_type j = 0; j < nbd; j++, ++p)
		if (k == l)
		  no += (*p) * bgeot::vect_sp(vval[i], vval[j]);
      }
    return sqrt(no);
  }

  template<class MESH_FEM, class VECT>
  scalar_type H1_semi_norm(MESH_FEM &mf, const VECT &U, size_type N)
  {
    return H1_semi_norm<MESH_FEM,VECT>(mf,U,N,mf.convex_index());
  }

  template<class MESH_FEM, class VECT>
  scalar_type H1_norm(MESH_FEM &mf, const VECT &U, size_type N, const dal::bit_vector& cvlst) {
    return sqrt( dal::sqr(L2_norm(mf, U, N, cvlst)) 
		 + dal::sqr(H1_semi_norm(mf, U, N, cvlst)));
  }

  template<class MESH_FEM, class VECT>
  scalar_type H1_norm(MESH_FEM &mf, const VECT &U, size_type N) {
    return sqrt( dal::sqr(L2_norm(mf, U, N)) 
		 + dal::sqr(H1_semi_norm(mf, U, N)));
  }
}
#endif 
