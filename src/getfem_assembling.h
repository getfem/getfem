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


#ifndef __GETFEM_ASSEMBLING_H
#define __GETFEM_ASSEMBLING_H

#include <getfem_mesh_fem.h>
#include <gmm_kernel.h>
#include <getfem_assembling_tensors.h>
namespace getfem
{
  /**
     compute $\|U\|_2$
  */
  template<class VEC>
  scalar_type asm_L2_norm(const mesh_fem &mf, const VEC &U) {
    return asm_L2_norm(mf,U,mf.convex_index());
  }
  template<class VEC>
  scalar_type asm_L2_norm(const mesh_fem &mf, const VEC &U,
			  const dal::bit_vector &cvlst) {
    generic_assembly assem;    
    if (mf.get_qdim() == 1)
      assem.set("u=data(#1); V()+=u(i).u(j).comp(Base(#1).Base(#1))(i,j)");
    else
      assem.set("u=data(#1); V()+=u(i).u(j).comp(vBase(#1).vBase(#1))(i,k,j,k)");
    assem.push_mf(mf);
    assem.push_data(U);
    bgeot::vsvector<scalar_type> v(1);
    assem.push_vec(v);
    assem.volumic_assembly(cvlst);
    return sqrt(v[0]);
  }

  /**
     compute $\|\nabla U\|_2$
   */
  template<class VEC>
  scalar_type asm_H1_semi_norm(const mesh_fem &mf, const VEC &U) {
    return asm_H1_semi_norm(mf,U,mf.convex_index());
  }
  template<class VEC>
  scalar_type asm_H1_semi_norm(const mesh_fem &mf, const VEC &U,
			       const dal::bit_vector &cvlst) {
    generic_assembly assem;    
    if (mf.get_qdim() == 1)
      assem.set("u=data(#1); V()+=u(i).u(j).comp(Grad(#1).Grad(#1))(i,d,j,d)");
    else
      assem.set("u=data(#1); V()+=u(i).u(j).comp(vGrad(#1).vGrad(#1))(i,k,d,j,k,d)");
    assem.push_mf(mf);
    assem.push_data(U);
    bgeot::vsvector<scalar_type> v(1);
    assem.push_vec(v);
    assem.volumic_assembly(cvlst);
    return sqrt(v[0]);
  }

  /** 
      compute the H1 norm of U.
  */
  template<class VEC>
  scalar_type asm_H1_norm(const mesh_fem &mf, const VEC &U) {
    return asm_H1_norm(mf,U,mf.convex_index());
  }
  template<class VEC>
  scalar_type asm_H1_norm(const mesh_fem &mf, const VEC &U,
			  const dal::bit_vector &cvlst) {
    return sqrt(dal::sqr(asm_L2_norm(mf,U,cvlst))+dal::sqr(asm_H1_semi_norm(mf,U,cvlst)));
  }
  

  /** 
      generic mass matrix assembly (on the whole mesh or on the specified boundary) 
  */
  template<class MAT>
  void asm_mass_matrix(MAT &M, const mesh_fem &mf_u1, size_type boundary=size_type(-1)) {
    asm_mass_matrix(M,mf_u1,mf_u1, boundary);
  }

  /** 
      generic mass matrix assembly (on the whole mesh or on the specified boundary) 
  */
  template<class MAT>
  void asm_mass_matrix(MAT &M, const mesh_fem &mf_u1, const mesh_fem &mf_u2,
		       size_type boundary=size_type(-1))
  {
    generic_assembly assem;
    if (&mf_u1 != &mf_u2)
      if (mf_u1.get_qdim() == 1 && mf_u2.get_qdim() == 1)
	assem.set("M(#1,#2)+=comp(Base(#1).Base(#2))");
      else
	assem.set("M(#1,#2)+=comp(vBase(#1).vBase(#2))(:,i,:,i);");
    else
      if (mf_u1.get_qdim() == 1 && mf_u2.get_qdim() == 1)
	assem.set("M(#1,#2)+=sym(comp(Base(#1).Base(#2)))");
      else
	assem.set("M(#1,#2)+=sym(comp(vBase(#1).vBase(#2))(:,i,:,i));");
    assem.push_mf(mf_u1);
    assem.push_mf(mf_u2);
    assem.push_mat(M);
    if (boundary==size_type(-1)) 
      assem.volumic_assembly();
    else assem.boundary_assembly(boundary);
  }


  /**
     source term (for both volumic sources and boundary (neumann) sources  
  */
  template<class VECT1, class VECT2>
    void asm_source_term(VECT1 &B, const mesh_fem &mf,
			 const mesh_fem &mfdata, const VECT2 &F, size_type boundary=size_type(-1))
  {
    generic_assembly assem;
    if (mf.get_qdim() == 1)
      assem.set("F=data(#2); V(#1)+=comp(Base(#1).Base(#2))(:,j).F(j);");
    else
      assem.set("F=data(qdim(#1),#2); V(#1)+=comp(vBase(#1).Base(#2))(:,i,j).F(i,j);");
    assem.push_mf(mf);
    assem.push_mf(mfdata);
    assem.push_data(F);
    assem.push_vec(B);
    if (boundary == size_type(-1))
      assem.volumic_assembly();
    else assem.boundary_assembly(boundary);
  }


  /**
    assembles $\int{qu.v}$

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

    Works for both volumic assembly and boundary assembly
   */
  template<class MAT, class VECT>
    void asm_qu_term(MAT &M, 
		     const mesh_fem &mf_u, 
		     const mesh_fem &mf_d, const VECT &Q, size_type boundary=size_type(-1))
  {
    generic_assembly assem;
    if (mf_u.get_qdim() == 1)
      assem.set("Q=data$1(#2);"
		"M(#1,#1)+=comp(Base(#1).Base(#1).Base(#2))(:,:,k).Q(k);");
    else
      assem.set("Q=data$1(qdim(#1),qdim(#1),#2);"
		"M(#1,#1)+=comp(vBase(#1).vBase(#1).Base(#2))(:,i,:,j,k).Q(i,j,k);");
    assem.push_mf(mf_u);
    assem.push_mf(mf_d);
    assem.push_data(Q);
    assem.push_mat(M);
    if (boundary == size_type(-1))
      assem.volumic_assembly();
    else
      assem.boundary_assembly(boundary);
  }


  /** 
     Stiffness matrix for linear elasticity, with Lamé coefficients
  */
  template<class MAT, class VECT>
    void asm_stiffness_matrix_for_linear_elasticity(MAT &RM,
						   const mesh_fem &mf, 
						   const mesh_fem &mfdata, 
						   const VECT &LAMBDA, const VECT &MU)
  {
    if (mf.get_qdim() != mf.linked_mesh().dim()) DAL_THROW(std::logic_error, "wrong qdim for the mesh_fem");
    /* e = strain tensor,
       M = 2*mu*e(u):e(v) + lambda*tr(e(u))*tr(e(v))
    */
    generic_assembly assem("lambda=data$1(#2); mu=data$2(#2);"
			   "t=comp(vGrad(#1).vGrad(#1).Base(#2));"
			   "e=(t{:,2,3,:,5,6,:}+t{:,3,2,:,5,6,:}+t{:,2,3,:,6,5,:}+t{:,3,2,:,6,5,:})/4;"
			   "M(#1,#1)+= sym(2*e(:,i,j,:,i,j,k).mu(k) + e(:,i,i,:,j,j,k).lambda(k))");

    assem.push_mf(mf);
    assem.push_mf(mfdata);
    assem.push_data(LAMBDA);
    assem.push_data(MU);
    assem.push_mat(RM);
    assem.volumic_assembly();
  }

  /** 
     Stiffness matrix for linear elasticity, with a general Hooke  tensor 
  */
  template<class MAT, class VECT>
    void asm_stiffness_matrix_for_linear_elasticity_Hooke(MAT &RM,
							 const mesh_fem &mf, 
							 const mesh_fem &mfdata, 
							 const VECT &H)
  {
    /* e = strain tensor,
       M = a_{i,j,k,l}e_{i,j}(u)e_{k,l}(v)
    */
    generic_assembly assem("a=data$1(qdim(#1),qdim(#1),qdim(#1),qdim(#1),#2);"
			   "t=comp(vGrad(#1).vGrad(#1).Base(#2));"
			   "e=(t{:,2,3,:,5,6,:}+t{:,3,2,:,5,6,:}+t{:,2,3,:,6,5,:}+t{:,3,2,:,6,5,:})/4;"
			   "M(#1,#1)+= sym(e(:,i,j,:,k,l,p).a(i,j,k,l,p))");
    assem.push_mf(mf);
    assem.push_mf(mfdata);
    assem.push_data(H);
    assem.push_mat(RM);
    assem.volumic_assembly();
  }

  /** two-in-one assembly of stokes equation:
      linear elasticty part and p.div(v) term are assembled at the
      same time. 
  */
  template<class MAT, class VECT>
    void asm_stokes(MAT &K, MAT &B, 
		    const mesh_fem &mf_u,
		    const mesh_fem &mf_p,
		    const mesh_fem &mf_d, const VECT &viscos)
  {
    generic_assembly assem("visc=data$1(#3); "
			   "t=comp(vGrad(#1).vGrad(#1).Base(#3));"
			   "e=(t{:,2,3,:,5,6,:}+t{:,3,2,:,5,6,:}+t{:,2,3,:,6,5,:}+t{:,3,2,:,6,5,:})/4;"
			   "M$1(#1,#1) += sym(e(:,i,j,:,i,j,k).visc(k));"          // visc*D(u):D(v)
			   "M$2(#1,#2) += comp(vGrad(#1).Base(#2))(:,i,i,:);"); // p.div v
    assem.push_mf(mf_u);
    assem.push_mf(mf_p);
    assem.push_mf(mf_d);
    assem.push_data(viscos);
    assem.push_mat(K);
    assem.push_mat(B);
    assem.volumic_assembly();
  }

  /**
     assembly of $\int_\Omega a(x)\nabla u.\nabla v$ , where $a(x)$ is scalar.
  */
  template<class MAT, class VECT>
    void asm_stiffness_matrix_for_laplacian(MAT &M, const mesh_fem &mf,
					   const mesh_fem &mfdata, const VECT &A)
  {
    generic_assembly assem("a=data$1(#2); M$1(#1,#1)+=sym(comp(Grad(#1).Grad(#1).Base(#2))(:,i,:,i,j).a(j))");
    //generic_assembly assem("a=data$1(#2); M$1(#1,#1)+=comp(Grad(#1).Grad(#1).Base(#2))(:,i,:,i,j).a(j)");
    assem.push_mf(mf);
    assem.push_mf(mfdata);
    assem.push_data(A);
    assem.push_mat(M);
    assem.volumic_assembly();
  }

  /**
     assembly of $\int_\Omega A(x)\nabla u.\nabla v$, where $A(x)$ is a matrix.
  */
  template<class MAT, class VECT>
    void asm_stiffness_matrix_for_scalar_elliptic(MAT &M, const mesh_fem &mf,
						  const mesh_fem &mfdata,
						  const VECT &A)
  {
    generic_assembly assem("a=data$1(mdim(#1),mdim(#1),#2); M$1(#1,#1)+=comp(Grad(#1).Grad(#1).Base(#2))(:,i,:,j,k).a(j,i,k)");
    assem.push_mf(mf);
    assem.push_mf(mfdata);
    assem.push_data(A);
    assem.push_mat(M);
    assem.volumic_assembly();
  }


  /* ********************************************************************* */
  /*                                                                       */
  /*  Algorithmes d'assemblage pour quelques problemes elliptiques.        */
  /*                                                                       */
  /* ********************************************************************* */

  /* ********************************************************************* */
  /*	Stiffness matrix for laplacian.                                     */
  /* ********************************************************************* */

  template<class MAT, class VECT>
    void assembling_stiffness_matrix_for_laplacian(MAT &RM, const mesh_fem &mf,
						  const mesh_fem &mfdata, const VECT &A)
  { // optimisable

    DAL_WARNING(3, "obsolete function - use asm_stiffness_matrix_for_laplacian");
    size_type nbd1, nbd2, N = mf.linked_mesh().dim();
    base_tensor t;
    pfem pf1, pf2, pf1prec = 0, pf2prec = 0;
    pintegration_method pim, pimprec = 0;
    bgeot::pgeometric_trans pgt, pgtprec = 0;
    pmat_elem_computation pmec = 0;

    if (&(mf.linked_mesh()) != &(mfdata.linked_mesh()))
      DAL_THROW(std::invalid_argument,
		"This assembling procedure only works on a single mesh");

    for (dal::bv_visitor cv(mf.convex_index()); !cv.finished(); ++cv) {
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
      pmec->gen_compute(t, mf.linked_mesh().points_of_convex(cv), cv);
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
    DAL_WARNING(3, "obsolete function - use asm_stokes");

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
    for (dal::bv_visitor cv(mf_u.convex_index()); !cv.finished(); ++cv) {
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
      pmec->gen_compute(t, mf_u.linked_mesh().points_of_convex(cv), cv);
      
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
  /*	Stiffness matrix for linear elasticity.                             */
  /* ********************************************************************* */
  template<class MAT, class VECT>
    void assembling_stiffness_matrix_for_linear_elasticity(MAT &RM,
							  const mesh_fem &mf, 
							  const mesh_fem &mfdata, 
							  const VECT &LAMBDA, const VECT &MU)
  { // à verifier
    DAL_WARNING(3, "obsolete function - use asm_stiffness_matrix_for_linear_elasticity");

    size_type nbd1, nbd2, N = mf.linked_mesh().dim();
    base_tensor t;
    pfem pf1, pf2, pf1prec = NULL, pf2prec = NULL;
    pintegration_method pim, pimprec = 0;
    bgeot::pgeometric_trans pgt, pgtprec = NULL;
    pmat_elem_type pme; pmat_elem_computation pmec = 0;

    if (&(mf.linked_mesh()) != &(mfdata.linked_mesh()))
      DAL_THROW(std::invalid_argument,
		"This assembling procedure only works on a single mesh");
  
    for (dal::bv_visitor cv(mf.convex_index()); !cv.finished(); ++cv) {
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
      pmec->gen_compute(t, mf.linked_mesh().points_of_convex(cv), cv);
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
     new version, which takes into account the qdim dimension of the mesh_fem 

     size(H) = Qdim(mf_u)^2 * nb_dof(mf_rh);
     size(R) = Qdim(mf_u)   * nb_dof(mf_rh);

     this function is able to "simplify" the dirichlet constraints (see below)
   */
  template<class MAT, class VECT>
  void asm_dirichlet_constraints(MAT &M, VECT &B, const mesh_fem &mf_u,
				 const mesh_fem &mf_rh,
				 const VECT &H, const VECT &R, size_type boundary) {
    dal::bit_vector nf;
    pfem pf_u, pf_rh;
    
    if (mf_rh.get_qdim() != 1)
      DAL_THROW(std::invalid_argument, "invalid data mesh fem for dirichlet");
    asm_qu_term(M, mf_u, mf_rh, H, boundary);
    asm_source_term(B, mf_u, mf_rh, R, boundary);

    /* step 2 : simplification of simple dirichlet conditions 
     */
    for (dal::bv_visitor cv(mf_u.convex_index()); !cv.finished(); ++cv) {
      nf = mf_u.faces_of_convex_on_boundary(cv, boundary);
      /* don't try anything with vector elements */
      if (mf_u.fem_of_element(cv)->target_dim() != 1) continue;
      if (nf.card() > 0) {
	pf_u = mf_u.fem_of_element(cv); 
	pf_rh = mf_rh.fem_of_element(cv); 
	size_type f;
	for (f << nf; f != ST_NIL; f << nf) {
	  bgeot::pconvex_structure cvs_u = pf_u->structure();
	  bgeot::pconvex_structure cvs_rh = pf_rh->structure();
	  for (size_type i = 0; i < cvs_u->nb_points_of_face(f); ++i) {

	    size_type Q = mf_u.get_qdim();  ///pf_u->target_dim() (==1)

	    size_type ind_u = cvs_u->ind_points_of_face(f)[i];
	    pdof_description tdof_u = pf_u->dof_types()[ind_u];
	    
	    for (size_type j = 0; j < cvs_rh->nb_points_of_face(f); ++j) {
	      size_type ind_rh = cvs_rh->ind_points_of_face(f)[j];
	      pdof_description tdof_rh = pf_rh->dof_types()[ind_rh];
	      /*
		same kind of dof and same location of dof ? 
		=> then the previous was not useful for this dofs (introducing
		a mass matrix which is not diagonal in the constraints matrix)
		-> the constraint is simplified:
		we replace \int{(H_j.psi_j)*phi_i}=\int{R_j.psi_j}  (sum over j)
		with             H_j*phi_i = R_j     
	      */
	      if (tdof_u == tdof_rh &&
		  bgeot::vect_dist2_sqr(pf_u->node_convex().points()[ind_u], 
					pf_rh->node_convex().points()[ind_rh]) < 1.0E-14) {
		/* the dof might be "duplicated" */
		for (size_type q = 0; q < Q; ++q) {
		  size_type dof_u = mf_u.ind_dof_of_element(cv)[ind_u*Q + q];

		  /* "erase" the row */
		  for (size_type k = 0; k < mf_u.nb_dof_of_element(cv); ++k) {
		    size_type dof_k = mf_u.ind_dof_of_element(cv)[k];
		    M(dof_u, dof_k) = 0.0;
		  }

		  size_type dof_rh = mf_rh.ind_dof_of_element(cv)[ind_rh];
		  /* set the "simplified" row */
		  for (unsigned jj=0; jj < Q; jj++) {
		    size_type dof_u2 = mf_u.ind_dof_of_element(cv)[ind_u*Q + jj];		    
		    M(dof_u, dof_u2) = H[(jj*Q+q) + Q*Q*(dof_rh)];
		  }
		  B[dof_u] = R[dof_rh*Q+q];
		}
	      }      
	    }
	  }
	}
      }
    }
  }

  template<class MAT, class VECT>
  void asm_dirichlet_constraints(MAT &M, VECT &B, const mesh_fem &mf_u,
				 const mesh_fem &mf_rh,
				 const VECT &R, size_type boundary) {
    if (mf_rh.get_qdim() != 1) 
      DAL_THROW(std::invalid_argument,
		"mf_rh should be a scalar (qdim=1) mesh_fem");
    size_type N = mf_rh.nb_dof(), Q=mf_u.get_qdim();
    VECT H(dal::sqr(mf_u.get_qdim())*N); gmm::clear(H);
    
    for (size_type i=0; i < N; ++i) {
      for (size_type q=0; q < Q; ++q) {
	H[i*Q*Q+q*Q+q]=1;
      }
    }
    asm_dirichlet_constraints(M, B, mf_u, mf_rh, H, R, boundary);
  }

  /* old version, pre-Qdim */
  template<class MAT, class VECT>
  void assembling_dirichlet_constraints(MAT &M, VECT &B, const mesh_fem &mf_u,
					size_type boundary,
					const mesh_fem &mf_rh,
				        const VECT &H, const VECT &R,
					dim_type N) {
    DAL_WARNING(3, "obsolete function - use asm_dirichlet_constraints");
    dal::bit_vector nf;
    pfem pf_u, pf_rh;
    
    assembling_boundary_qu_term(M, mf_u, boundary, mf_rh, H, N);
    assembling_Neumann_condition(B, mf_u, boundary, mf_rh, R, N);
    for (dal::bv_visitor cv(mf_u.convex_index()); !cv.finished(); ++cv) {
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
    void mass_matrix(MATRM &M, const MESH_FEM &mf1,
		     const MESH_FEM &mf2, dim_type N)
  {
    DAL_WARNING(3, "obsolete function - use asm_mass_matrix");
    size_type nbd1, nbd2;
    base_tensor t;
    pfem pf1, pf1prec = 0, pf2, pf2prec = 0;
    pintegration_method pim, pimprec = 0;
    bgeot::pgeometric_trans pgt, pgtprec = NULL;
    pmat_elem_type pme; pmat_elem_computation pmec = 0;
    // M(0,0) = 1.0;  ??

    if (&(mf1.linked_mesh()) != &(mf2.linked_mesh()))
      DAL_THROW(std::invalid_argument,
		"This assembling procedure only works on a single mesh");

    for (dal::bv_visitor cv(mf1.convex_index()); !cv.finished(); ++cv) {
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
      pmec->gen_compute(t, mf1.linked_mesh().points_of_convex(cv), cv);

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

  template<class MATRM, class MESH_FEM>
  void mass_matrix_on_boundary(MATRM &M, const MESH_FEM &mf1,
		    const MESH_FEM &mf2, size_type boundary, dim_type N)
  {
    DAL_WARNING(3, "obsolete function - use asm_mass_matrix");
    size_type nbd1, nbd2, f;
    dal::bit_vector nf;
    base_tensor t;
    pfem pf1, pf1prec = 0, pf2, pf2prec = 0;
    pintegration_method pim, pimprec = 0;
    bgeot::pgeometric_trans pgt, pgtprec = NULL;
    pmat_elem_type pme; pmat_elem_computation pmec = 0;
    // M(0,0) = 1.0;  ??

    if (&(mf1.linked_mesh()) != &(mf2.linked_mesh()))
      DAL_THROW(std::invalid_argument,
		"This assembling procedure only works on a single mesh");

    for (dal::bv_visitor cv(mf1.convex_index()); !cv.finished(); ++cv) {
      nf = mf1.faces_of_convex_on_boundary(cv, boundary);
      if (nf.card() > 0) {
	pf1 = mf1.fem_of_element(cv); nbd1 = pf1->nb_dof();
	pf2 = mf2.fem_of_element(cv); nbd2 = pf2->nb_dof();
	pgt = mf1.linked_mesh().trans_of_convex(cv);
	pim = mf1.int_method_of_element(cv);
	if (pf1prec != pf1 || pf2prec != pf2 || pgtprec != pgt 
	    || pimprec != pim) {
	  pme = mat_elem_product(mat_elem_base(pf1), mat_elem_base(pf2));
	  pmec = mat_elem(pme, pim, pgt);
	  pf1prec = pf1; pf2prec = pf2; pgtprec = pgt; pimprec = pim;
	}

	for (f << nf; f != ST_NIL; f << nf) {

	  pmec->gen_compute_on_face(t, mf1.linked_mesh().points_of_convex(cv),
				    f, cv);
	  
	  base_tensor::iterator p = t.begin();
	  for (size_type i = 0; i < nbd2; i++) {
	    size_type dof2 = mf2.ind_dof_of_element(cv)[i];
	    // cout << "cv = " << cv << " dof2 = " << dof2 << endl;
	    for (size_type j = 0; j < nbd1; j++, ++p) {
	      size_type dof1 = mf1.ind_dof_of_element(cv)[j];
	      // cout << "dof1 = " << dof1 << " dof2 = " << dof2 << endl;
	      for (size_type k = 0; k < N; k++)
		M(dof1*N + k, dof2*N + k) += (*p);
	    }
	  }
	  if (p != t.end()) DAL_THROW(dal::internal_error, "internal error"); 
	}
      }
    }
  }

  template<class MATRM, class MESH_FEM>
  inline void mass_matrix_on_boundary(MATRM &M, const MESH_FEM &mf,
				      size_type boundary, dim_type N)
  { mass_matrix_on_boundary(M, mf, mf, boundary, N); }


  /* ********************************************************************* 
     volumic source term
     ********************************************************************* */
  template<class VECT1, class VECT2>
    void assembling_volumic_source_term(VECT1 &B, const mesh_fem &mf,
					const mesh_fem &mfdata, const VECT2 &F, dim_type N)
  {
    DAL_WARNING(3, "obsolete function - use asm_source_term");
    size_type nbd1, nbd2;
    base_tensor t;
    pfem pf1, pf2, pf1prec = NULL, pf2prec = NULL;
    pintegration_method pim, pimprec = 0;
    bgeot::pgeometric_trans pgt, pgtprec = NULL;
    pmat_elem_type pme; pmat_elem_computation pmec = 0;

    if (&(mf.linked_mesh()) != &(mfdata.linked_mesh()))
      DAL_THROW(std::invalid_argument,
		"This assembling procedure only works on a single mesh");

    for (dal::bv_visitor cv(mf.convex_index()); !cv.finished(); ++cv) {
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
      pmec->gen_compute(t, mf.linked_mesh().points_of_convex(cv), cv);
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
					const VECT2 &F)
  { // Marche uniquement pour des ddl de lagrange.
    size_type Q=mf.get_qdim();
    dal::bit_vector nndof = mf.dof_on_boundary(boundary);
    pfem pf1;

    for (dal::bv_visitor cv(mf.convex_index()); !cv.finished(); ++cv) {
      pf1 = mf.fem_of_element(cv);
      pdof_description ldof = lagrange_dof(pf1->dim());
      size_type nbd = pf1->nb_dof();
      for (size_type i = 0; i < nbd; i++)
      {
	size_type dof1 = mf.ind_dof_of_element(cv)[i*Q];
	if (nndof.is_in(dof1) && pf1->dof_types()[i] == ldof)
	{
	  // cout << "dof : " << i << endl;
	  for (size_type j = 0; j < nbd; j++)
	  {
	    size_type dof2 = mf.ind_dof_of_element(cv)[j*Q];
	    for (size_type k = 0; k < Q; ++k)
	      for (size_type l = 0; l < Q; ++l)
	      {
		if (!(nndof.is_in(dof2)) &&
		    dof_compatibility(pf1->dof_types()[j],
				      lagrange_dof(pf1->dim())))
		  B[dof2+k] -= RM(dof2+k, dof1+l) * F[dof1+l];
		RM(dof2+k, dof1+l) = RM(dof1+l, dof2+k) = 0;
	      }
	  }
	  for (size_type k = 0; k < Q; ++k)
	  { RM(dof1+k, dof1+k) = 1; B[dof1+k] = F[dof1+k]; }
	}
      }
    }
  }

  /* add a dirichlet condition on a single dof, modifiying the matrix
     RM and the rhs B.  (keeping the symmetry properties) */
  template<class MATRM, class VECT1>
    void add_Dirichlet_dof(MATRM &RM, VECT1 &B,
			   const mesh_fem &mf,
			   size_type dof, scalar_type dof_val) {     
    size_type Q=mf.get_qdim();
    bgeot::mesh_convex_ind_ct dofcv = mf.convex_to_dof(dof);
    pfem pf1;

    for (bgeot::mesh_convex_ind_ct::const_iterator it = dofcv.begin();
	 it != dofcv.end(); ++it) {
      pf1 = mf.fem_of_element(*it);
      if (pf1->target_dim() != 1)
	DAL_THROW(to_be_done_error, "sorry, to be done ... ");
      size_type nbd = pf1->nb_dof();
      for (size_type i = 0; i < nbd * Q; i++) {
	size_type dof1 = mf.ind_dof_of_element(*it)[i];
	if (dof == dof1) {
	  for (size_type j = 0; j < nbd * Q; j++) {
	    size_type dof2 = mf.ind_dof_of_element(*it)[j];
	    if (!(dof == dof2)) {
	      B[dof2] -= RM(dof2, dof1) * dof_val;
	      RM(dof2, dof1) = RM(dof1, dof2) = 0;
	    }
	  }
	  RM(dof1, dof1) = 1; B[dof1] = dof_val;
	}
      }
    }
  }


  template<class MATRM, class VECT1, class VECT2>
    void old_assembling_Dirichlet_condition(MATRM &RM, VECT1 &B,
					const mesh_fem &mf,
					size_type boundary,
					const VECT2 &F, dim_type N)
  { /* Y-a-il un moyen plus performant ? */
    // Marche uniquement pour des ddl de lagrange.
    dal::bit_vector nndof = mf.dof_on_boundary(boundary);
    pfem pf1;

    for (dal::bv_visitor cv(mf.convex_index()); !cv.finished(); ++cv) {
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


  template<class MATD, class MATG, class VECT>
  size_type Dirichlet_nullspace(const MATD &D, MATG &G,
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
    typedef typename gmm::linalg_traits<MATD>::value_type T;
    typedef typename gmm::number_traits<T>::magnitude_type R;
    typedef typename gmm::temporary_vector<MATD>::vector_type TEMP_VECT;
    R tol = gmm::default_tol(R());
    R norminfD = gmm::mat_norminf(D);
    size_type nbd = gmm::mat_ncols(D), nbase = 0, nbr = gmm::mat_nrows(D);
    TEMP_VECT aux(nbr), e(nbd), f(nbd);
    dal::dynamic_array<TEMP_VECT> base_img;
    dal::dynamic_array<TEMP_VECT> base_img_inv;
    size_type nb_bimg = 0;

    if (!(gmm::is_col_matrix(D)))
      DAL_WARNING(3,
		  "Dirichlet_nullspace is inefficient when D is a row matrix");
    // First, detection of null columns of D, and already orthogonals 
    // vectors of the image of D.
    dal::bit_vector nn;
    for (size_type i = 0; i < nbd; ++i) {
      gmm::clear(e); e[i] = T(1);
      gmm::mult(D, e, aux);
      R n = gmm::vect_norm2(aux);

      if (n < norminfD*tol) {
	G(i, nbase++) = T(1); nn[i] = true;
      }
      else {
	bool good = true;
	for (size_type j = 0; j < nb_bimg; ++j)
	  if (dal::abs(gmm::vect_sp(aux, base_img[j])) > T(0))
	    { good = false; break; }
	if (good) {
	  gmm::copy(e, f);
	  gmm::scale(f, T(R(1)/n)); gmm::scale(aux, T(R(1)/n));
	  base_img_inv[nb_bimg] = TEMP_VECT(nbd);
	  gmm::copy(f, base_img_inv[nb_bimg]);
	  gmm::clean(aux, gmm::vect_norminf(aux)*tol);
	  base_img[nb_bimg] = TEMP_VECT(nbr);
	  gmm::copy(aux, base_img[nb_bimg++]);
	  nn[i] = true;
	}
      }
    }
    size_type nb_triv_base = nbase;

    for (size_type i = 0; i < nbd; ++i)
      if (!(nn[i])) {
	gmm::clear(e); e[i] = T(1); gmm::clear(f); f[i] = T(1);
	gmm::mult(D, e, aux);
	for (size_type j = 0; j < nb_bimg; ++j) { 
	  T c = gmm::vect_sp(aux, base_img[j]);
	  //	  if (dal::abs(c > 1.0E-6) { // à scaler sur l'ensemble de D ...
	  if (c != T(0)) {
	    gmm::add(gmm::scaled(base_img[j], -c), aux);
	    gmm::add(gmm::scaled(base_img_inv[j], -c), f);
	  }
	}
	if (gmm::vect_norm2(aux) < norminfD*tol*R(10000)) {
	  gmm::copy(f, gmm::mat_col(G, nbase++));
	}
	else {
	  R n = gmm::vect_norm2(aux);
	  gmm::scale(f, T(R(1)/n)); gmm::scale(aux, T(R(1)/n));
	  gmm::clean(f, tol*gmm::vect_norminf(f)); gmm::clean(aux, tol*gmm::vect_norminf(aux));
	  base_img_inv[nb_bimg] = TEMP_VECT(nbd);
	  gmm::copy(f, base_img_inv[nb_bimg]);
	  base_img[nb_bimg] = TEMP_VECT(nbr);
	  gmm::copy(aux, base_img[nb_bimg++]);
	}
      }
    // Compute a solution in UDD
    gmm::clear(UDD);
    for (size_type i = 0; i < nb_bimg; ++i) {
      scalar_type c = gmm::vect_sp(base_img[i], UD);
      gmm::add(gmm::scaled(base_img_inv[i], c), UDD);
    }
    // Orthogonalisation of the basis of the kernel of D.
    for (size_type i = nb_triv_base + 1; i < nbase; ++i) {
      for (size_type j = nb_triv_base; j < i; ++j) {
	T c = gmm::vect_sp(gmm::mat_col(G,i), gmm::mat_col(G,j));
	if (c != T(0))
	  gmm::add(gmm::scaled(gmm::mat_col(G,j), -c), gmm::mat_col(G,i));
      }
    }
    // projection of UDD on the orthogonal to the kernel.
    for (size_type j = nb_triv_base; j < nbase; ++j) {
      T c = gmm::vect_sp(gmm::mat_col(G,j), UDD);
      if (c != T(0))
	gmm::add(gmm::scaled(gmm::mat_col(G,j), -c), UDD);
    }
    // Test ...
    gmm::mult(D, UDD, gmm::scaled(UD, T(-1)), aux);
    if (gmm::vect_norm2(aux) > gmm::vect_norm2(UDD)*tol*R(10000))
      DAL_WARNING(2, "Dirichlet condition not well inverted: residu="
		  << gmm::vect_norm2(aux));
    return nbase;
  }
  

//   template<class MATD, class VECT>
//     size_type treat_Dirichlet_condition(const MATD &D, MATD &G,
// 				   const VECT &UD, VECT &UDD) {

//     size_type nbd = D.ncols(), nbase = 0;
//     gmm::wsvector<scalar_type> aux(D.nrows()), e(nbd);
//     gmm::wsvector<scalar_type> f(nbd);
//     dal::dynamic_array<gmm::wsvector<scalar_type> > base_img;
//     dal::dynamic_array<gmm::wsvector<scalar_type> > base_img_inv;
//     size_type nb_bimg = 0;

//     // First, detection of null columns of D, and already orthogonals 
//     // vectors of the image of D.
//     dal::bit_vector nn;
//     for (size_type i = 0; i < nbd; ++i) {
//       e.clear(); e[i] = 1.0; f.clear(); f[i] = 1.0;
//       aux = D*e;
//       if (gmm::vect_norm2(aux) < 1.0E-8) { //à scaler sur l'ensemble de D ...
// 	G(nbase++, i) = 1.0; nn[i] = true;
//       }
//       else {
// 	bool good = true;
// 	for (size_type j = 0; j < nb_bimg; ++j)
// 	  if (dal::abs(gmm::vect_sp(aux, base_img[j])) > 1.0E-16)
// 	    { good = false; break; }
// 	if (good) {
// 	  scalar_type n = gmm::vect_norm2(aux);
// 	  f /= n; aux /= n;
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
// 	  scalar_type c = gmm::vect_sp(aux, base_img[j]);
// 	  //	  if (dal::abs(c > 1.0E-6) { // à scaler sur l'ensemble de D ...
// 	  if (c != 0.) {
// 	    aux -= base_img[j] * c;
// 	    f -= base_img_inv[j] * c;
// 	  }
// 	}
// 	//	cerr << "norm2(aux)= " << gmm::vect_norm2(aux) << "\n";
// 	if (gmm::vect_norm2(aux) < 1.0E-8) { // à scaler sur l'ensemble de D ...
// 	  G.row(nbase++) = f;
// 	}
// 	else {
// 	  scalar_type n = gmm::vect_norm2(aux);
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
//       scalar_type c = gmm::vect_sp(base_img[i], UD);
//       UDD += base_img_inv[i].full() * c;
//     }

//     // Orthogonalisation of the basis of the kernel of D.
//     for (size_type i = nb_triv_base + 1; i < nbase; ++i) {
//       for (size_type j = nb_triv_base; j < i; ++j) {
// 	scalar_type c = gmm::vect_sp(G.row(i), G.row(j));
// 	if (c != 0.)
// 	  G.row(i) -= G.row(j) * c;
//       }
//     }

//     // projection of UDD on the orthogonal to the kernel.
//     for (size_type j = nb_triv_base; j < nbase; ++j) {
//       scalar_type c = gmm::vect_sp(G.row(j), UDD);
//       if (c != 0.)
// 	UDD -= G.row(j).full() * c;
//     }

//     // Test ...
//     if (gmm::vect_norm2(D * UDD - UD) > 1.0E-12)
//       cerr << "Dirichlet condition not well inverted\n";

//     return nbase;
//   }
  

  
  /* ********************************************************************* */
  /*	Neumann Condition.                                                 */
  /* ********************************************************************* */

  template<class VECT1, class VECT2>
    void assembling_Neumann_condition(VECT1 &B, const mesh_fem &mf,
				      size_type boundary, const mesh_fem &mfdata, const VECT2 &F, dim_type N)
  {
    DAL_WARNING(3, "obsolete function - use asm_source_term");
    size_type nbd1, nbd2, f;
    dal::bit_vector nf;
    base_tensor t;
    pfem pf1, pf2, pf1prec = NULL, pf2prec = NULL;
    pintegration_method pim, pimprec = 0;
    bgeot::pgeometric_trans pgt, pgtprec = NULL;
    pmat_elem_type pme; pmat_elem_computation pmec = 0;

    if (&(mf.linked_mesh()) != &(mfdata.linked_mesh()))
      DAL_THROW(std::invalid_argument,
		"This assembling procedure only works on a single mesh");
  
    for (dal::bv_visitor cv(mf.convex_index()); !cv.finished(); ++cv) {
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

	  pmec->gen_compute_on_face(t,mf.linked_mesh().points_of_convex(cv), 
				    f, cv);
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
