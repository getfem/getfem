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


#ifndef GETFEM_ASSEMBLING_H__
#define GETFEM_ASSEMBLING_H__

#include <getfem_assembling_tensors.h>

namespace getfem
{
  /**
   * compute $\|U\|2_$, U might be real or complex
   */
  template<typename VEC>
  scalar_type asm_L2_norm(const mesh_fem &mf, const VEC &U) {
    return asm_L2_norm(mf,U,mf.convex_index());
  }

  template<typename VEC>
  scalar_type asm_L2_norm(const mesh_fem &mf, const VEC &U,
			  const dal::bit_vector &cvlst) {
    return sqrt(asm_L2_norm_sqr(mf,U,cvlst,
			       typename gmm::linalg_traits<VEC>::value_type()));
  }

  template<typename VEC, typename T>
  scalar_type asm_L2_norm_sqr(const mesh_fem &mf, const VEC &U,
			      const dal::bit_vector &cvlst, T) {
    generic_assembly assem;    
    if (mf.get_qdim() == 1)
      assem.set("u=data(#1); V()+=u(i).u(j).comp(Base(#1).Base(#1))(i,j)");
    else
      assem.set("u=data(#1);"
		"V()+=u(i).u(j).comp(vBase(#1).vBase(#1))(i,k,j,k)");
    assem.push_mf(mf);
    assem.push_data(U);
    bgeot::vsvector<scalar_type> v(1);
    assem.push_vec(v);
    assem.volumic_assembly(cvlst);
    return v[0];
  }

  template<typename VEC, typename T>
  scalar_type asm_L2_norm_sqr(const mesh_fem &mf, const VEC &U,
			      const dal::bit_vector &cvlst, std::complex<T>) {
    return asm_L2_norm_sqr(mf,gmm::real_part(U),cvlst,T()) + 
      asm_L2_norm_sqr(mf,gmm::imag_part(U),cvlst,T());
  }

  
  /**
   * compute $\|\nabla U\|2_$, U might be real or complex
   */
  template<typename VEC>
  scalar_type asm_H1_semi_norm(const mesh_fem &mf, const VEC &U) {
    return asm_H1_semi_norm(mf,U,mf.convex_index());
  }

  template<typename VEC>
  scalar_type asm_H1_semi_norm(const mesh_fem &mf, const VEC &U,
			       const dal::bit_vector &cvlst) {
    return sqrt(asm_H1_semi_norm_sqr(mf,U,cvlst,typename gmm::linalg_traits<VEC>::value_type()));
  }

  template<typename VEC, typename T>
  scalar_type asm_H1_semi_norm_sqr(const mesh_fem &mf, const VEC &U,
				   const dal::bit_vector &cvlst, T) {
    generic_assembly assem;    
    if (mf.get_qdim() == 1)
      assem.set("u=data(#1); V()+=u(i).u(j).comp(Grad(#1).Grad(#1))(i,d,j,d)");
    else
      assem.set("u=data(#1);"
		"V()+=u(i).u(j).comp(vGrad(#1).vGrad(#1))(i,k,d,j,k,d)");
    assem.push_mf(mf);
    assem.push_data(U);
    bgeot::vsvector<scalar_type> v(1);
    assem.push_vec(v);
    assem.volumic_assembly(cvlst);
    return v[0];
  }

  template<typename VEC, typename T>
  scalar_type asm_H1_semi_norm_sqr(const mesh_fem &mf, const VEC &U,
			       const dal::bit_vector &cvlst, std::complex<T>) {
    return asm_H1_semi_norm_sqr(mf, gmm::real_part(U), cvlst, T()) + 
      asm_H1_semi_norm_sqr(mf, gmm::imag_part(U), cvlst, T());
  }


  /** 
   *   compute the H1 norm of U.
   */
  template<typename VEC>
  scalar_type asm_H1_norm(const mesh_fem &mf, const VEC &U) {
    return asm_H1_norm(mf,U,mf.convex_index());
  }

  template<typename VEC>
  scalar_type asm_H1_norm(const mesh_fem &mf, const VEC &U,
			  const dal::bit_vector &cvlst) {
    return sqrt(dal::sqr(asm_L2_norm(mf,U,cvlst))
		+dal::sqr(asm_H1_semi_norm(mf,U,cvlst)));
  }
  
  /** 
   *  generic mass matrix assembly (on the whole mesh or on the specified
   *  boundary) 
   */
  template<typename MAT>
  void asm_mass_matrix(const MAT &M, const mesh_fem &mf_u1,
		       size_type boundary=size_type(-1)) {
    asm_mass_matrix(const_cast<MAT &>(M),mf_u1,mf_u1, boundary);
  }

  /** 
   *  generic mass matrix assembly (on the whole mesh or on the specified
   *  boundary) 
   */
  template<typename MAT>
  void asm_mass_matrix(MAT &M, const mesh_fem &mf_u1, const mesh_fem &mf_u2,
		       size_type boundary=size_type(-1)) {
    generic_assembly assem;
    if (&mf_u1 != &mf_u2)
      if (mf_u1.get_qdim() == 1 && mf_u2.get_qdim() == 1)
	assem.set("M(#1,#2)+=comp(Base(#1).Base(#2))");
      else
	assem.set("M(#1,#2)+=comp(vBase(#1).vBase(#2))(:,i,:,i);");
    else
      if (mf_u1.get_qdim() == 1)
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
   *  generic mass matrix assembly with an additional parameter
   *  (on the whole mesh or on the specified boundary) 
   */
  template<typename MAT, typename VECT>
  void asm_mass_matrix_param(MAT &M, const mesh_fem &mf_u, const mesh_fem &mfdata,
			     const VECT &F, size_type boundary=size_type(-1)) {
    generic_assembly assem;
    if (mf_u.get_qdim() == 1)
      assem.set("F=data(#2);"
		"M(#1,#1)+=sym(comp(Base(#1).Base(#1).Base(#2))(:,:,i).F(i))");
    else
      assem.set("F=data(#2);"
	 "M(#1,#1)+=sym(comp(vBase(#1).vBase(#1).Base(#2))(:,i,:,i,j).F(j));");
    assem.push_mf(mf_u);
    assem.push_mf(mfdata);
    assem.push_data(F);
    assem.push_mat(M);
    if (boundary==size_type(-1)) 
      assem.volumic_assembly();
    else assem.boundary_assembly(boundary);
  }

  /*  source term (for both volumic sources and boundary (neumann) sources.
   *  real version.
   */
  template<typename VECT1, typename VECT2, typename T>
    void asm_source_term(const VECT1 &B, const mesh_fem &mf,
			 const mesh_fem &mfdata, const VECT2 &F,
			 size_type boundary, T) {
    generic_assembly assem;
    if (mf.get_qdim() == 1)
      assem.set("F=data(#2); V(#1)+=comp(Base(#1).Base(#2))(:,j).F(j);");
    else
      assem.set("F=data(qdim(#1),#2);"
		"V(#1)+=comp(vBase(#1).Base(#2))(:,i,j).F(i,j);");
    assem.push_mf(mf);
    assem.push_mf(mfdata);
    assem.push_data(F);
    assem.push_vec(const_cast<VECT1&>(B));
    (boundary == size_type(-1)) ?
      assem.volumic_assembly() : assem.boundary_assembly(boundary);
  }

  /*  source term (for both volumic sources and boundary (neumann) sources.
   *  complex version.
   */
  template<typename VECT1, typename VECT2, typename T>
    void asm_source_term(VECT1 &B, const mesh_fem &mf,
			 const mesh_fem &mfdata, const VECT2 &F,
			 size_type boundary, std::complex<T>) {
    /*    generic_assembly assem;
    if (mf.get_qdim() == 1)
      assem.set("Fr=data$1(#2); V$1(#1)+=comp(Base(#1).Base(#2))(:,j).Fr(j);"
		"Fi=data$2(#2); V$2(#1)+=comp(Base(#1).Base(#2))(:,j).Fi(j);");
    else
      assem.set("Fr=data$1(qdim(#1),#2); Fi=data$2(qdim(#1),#2);"
		"V$1(#1)+=comp(vBase(#1).Base(#2))(:,i,j).Fr(i,j);"
		"V$2(#1)+=comp(vBase(#1).Base(#2))(:,i,j).Fi(i,j);");
    assem.push_mf(mf);
    assem.push_mf(mfdata);
    assem.push_data(gmm::real_part(F));
    assem.push_data(gmm::imag_part(F));
    assem.push_vec(gmm::real_part(B));
    assem.push_vec(gmm::imag_part(B));
    (boundary == size_type(-1)) ?
    assem.volumic_assembly() : assem.boundary_assembly(boundary);*/
    asm_source_term(gmm::real_part(B), mf, mfdata, gmm::real_part(F), boundary, T());
    asm_source_term(gmm::imag_part(B), mf, mfdata, gmm::imag_part(F), boundary, T());
  }

  /** 
   *  source term (for both volumic sources and boundary (neumann) sources.
   */
  template<typename VECT1, typename VECT2>
  void asm_source_term(VECT1 &B, const mesh_fem &mf,
		       const mesh_fem &mfdata, const VECT2 &F,
		       size_type boundary=size_type(-1)) {
    if (mfdata.get_qdim() != 1)
      DAL_THROW(invalid_argument, "invalid data mesh fem (Qdim=1 required)");
    asm_source_term(B, mf, mfdata, F, boundary,
		    typename gmm::linalg_traits<VECT2>::value_type());
  }



  template <typename V> bool is_Q_symmetric(const V& Q, size_type q, size_type nbd) {
    /* detect the symmetricity of Q (in that case the symmetricity of
     * the final matrix will be ensured, and computations will be
     * slightly speed up */
    for (size_type k=0; k < nbd; ++k)
      for (size_type i=1; i < q; ++i)
	for (size_type j=0; j < i; ++j)
	  if (Q[k*q*q+i*q+j] != Q[k*q*q+j*q+i])
	    return false;
    return true;
  }

  template<typename MAT, typename VECT, typename T>
    void asm_qu_term_(const MAT &M, 
		     const mesh_fem &mf_u, 
		     const mesh_fem &mf_d, const VECT &Q, 
		      size_type boundary, T) {
    generic_assembly assem;    
    if (mf_u.get_qdim() == 1)
      assem.set("Q=data$1(#2);"
		"M(#1,#1)+=comp(Base(#1).Base(#1).Base(#2))(:,:,k).Q(k);");
    else {
      if (is_Q_symmetric(Q,mf_u.get_qdim(),mf_d.nb_dof()))
	assem.set("Q=data$1(qdim(#1),qdim(#1),#2);"
		  "M(#1,#1)+=sym(comp(vBase(#1).vBase(#1).Base(#2))"
		  "(:,i,:,j,k).Q(i,j,k));");
      else
	assem.set("Q=data$1(qdim(#1),qdim(#1),#2);"
		  "M(#1,#1)+=comp(vBase(#1).vBase(#1).Base(#2))"
		  "(:,i,:,j,k).Q(i,j,k);");
    }
    assem.push_mf(mf_u);
    assem.push_mf(mf_d);
    assem.push_data(Q);
    assem.push_mat(const_cast<MAT&>(M));
    if (boundary == size_type(-1))
      assem.volumic_assembly();
    else
      assem.boundary_assembly(boundary);
  }

 template<typename MAT, typename VECT, typename T>
    void asm_qu_term_(MAT &M, 
		     const mesh_fem &mf_u, 
		     const mesh_fem &mf_d, const VECT &Q, 
		      size_type boundary, std::complex<T>) {
   asm_qu_term_(gmm::real_part(M), mf_u, mf_d, gmm::real_part(Q), boundary, T());
   asm_qu_term_(gmm::imag_part(M), mf_u, mf_d, gmm::imag_part(Q), boundary, T());
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

    Q is a vector, so the matrix is assumed to be stored by columns
    (fortran style)

    Works for both volumic assembly and boundary assembly
  */
  template<typename MAT, typename VECT>
    void asm_qu_term(const MAT &M, 
		     const mesh_fem &mf_u, 
		     const mesh_fem &mf_d, const VECT &Q, 
		     size_type boundary=size_type(-1)) {
    if (mf_d.get_qdim() != 1)
      DAL_THROW(invalid_argument, "invalid data mesh fem (Qdim=1 required)");
    asm_qu_term_(const_cast<MAT&>(M), mf_u, mf_d, Q, boundary, typename gmm::linalg_traits<VECT>::value_type());
  }


  /** 
     Stiffness matrix for linear elasticity, with Lamé coefficients
  */
  template<class MAT, class VECT>
    void asm_stiffness_matrix_for_linear_elasticity(const MAT &RM_,
					   const mesh_fem &mf,
					   const mesh_fem &mfdata,
					   const VECT &LAMBDA,const VECT &MU) {
    MAT &RM = const_cast<MAT &>(RM_);
    if (mfdata.get_qdim() != 1)
      DAL_THROW(invalid_argument, "invalid data mesh fem (Qdim=1 required)");
    
    if (mf.get_qdim() != mf.linked_mesh().dim())
      DAL_THROW(std::logic_error, "wrong qdim for the mesh_fem");
    /* e = strain tensor,
       M = 2*mu*e(u):e(v) + lambda*tr(e(u))*tr(e(v))
    */
    generic_assembly assem("lambda=data$1(#2); mu=data$2(#2);"
			   "t=comp(vGrad(#1).vGrad(#1).Base(#2));"
			   //"e=(t{:,2,3,:,5,6,:}+t{:,3,2,:,5,6,:}"
			   //"+t{:,2,3,:,6,5,:}+t{:,3,2,:,6,5,:})/4;"
			   //"e=(t{:,2,3,:,5,6,:}+t{:,3,2,:,5,6,:})*0.5;"
			   /*"M(#1,#1)+= sym(2*e(:,i,j,:,i,j,k).mu(k)"
                             " + e(:,i,i,:,j,j,k).lambda(k))");*/
                           "M(#1,#1)+= sym(t(:,i,j,:,i,j,k).mu(k)+t(:,j,i,:,i,j,k).mu(k)"
                             " + t(:,i,i,:,j,j,k).lambda(k))");
    assem.push_mf(mf);
    assem.push_mf(mfdata);
    assem.push_data(LAMBDA);
    assem.push_data(MU);
    assem.push_mat(RM);
    assem.volumic_assembly();
  }

  /** 
     Stiffness matrix for linear elasticity, with a general Hooke tensor. This is more a
     demonstration of generic assembly than something useful !
  */
  template<typename MAT, typename VECT> void
  asm_stiffness_matrix_for_linear_elasticity_Hooke(MAT &RM,
						   const mesh_fem &mf, 
						   const mesh_fem &mfdata, 
						   const VECT &H) {
    if (mfdata.get_qdim() != 1)
      DAL_THROW(invalid_argument, "invalid data mesh fem (Qdim=1 required)");
    /* e = strain tensor,
       M = a_{i,j,k,l}e_{i,j}(u)e_{k,l}(v)
    */
    generic_assembly assem("a=data$1(qdim(#1),qdim(#1),qdim(#1),qdim(#1),#2);"
			   "t=comp(vGrad(#1).vGrad(#1).Base(#2));"
			   "e=(t{:,2,3,:,5,6,:}+t{:,3,2,:,5,6,:}"
			   "  +t{:,2,3,:,6,5,:}+t{:,3,2,:,6,5,:})/4;"
			   "M(#1,#1)+= sym(e(:,i,j,:,k,l,p).a(i,j,k,l,p))");
    assem.push_mf(mf);
    assem.push_mf(mfdata);
    assem.push_data(H);
    assem.push_mat(RM);
    assem.volumic_assembly();
  }

  /** two-in-one assembly of stokes equation:
   *  linear elasticty part and p.div(v) term are assembled at the
   *  same time. 
   */
  template<typename MAT, typename VECT>
    void asm_stokes(MAT &K, MAT &BT, 
		    const mesh_fem &mf_u, const mesh_fem &mf_p,
		    const mesh_fem &mf_d, const VECT &viscos) {
    if (mf_d.get_qdim() != 1)
      DAL_THROW(invalid_argument, "invalid data mesh fem (Qdim=1 required)");
    generic_assembly assem("visc=data$1(#3); "
			   "t=comp(vGrad(#1).vGrad(#1).Base(#3));"
			   "e=(t{:,2,3,:,5,6,:}+t{:,3,2,:,5,6,:}"
			   "  +t{:,2,3,:,6,5,:}+t{:,3,2,:,6,5,:})/4;"
			   // visc*D(u):D(v)
			   "M$1(#1,#1)+=sym(e(:,i,j,:,i,j,k).visc(k));"
			   // p.div v
			   "M$2(#1,#2)+=comp(vGrad(#1).Base(#2))(:,i,i,:);");
    assem.push_mf(mf_u);
    assem.push_mf(mf_p);
    assem.push_mf(mf_d);
    assem.push_data(viscos);
    assem.push_mat(K);
    assem.push_mat(BT);
    assem.volumic_assembly();
  }

  template<typename MAT>
    void asm_stokes_B(MAT &B, const mesh_fem &mf_u, const mesh_fem &mf_p) {
    if (mf_p.get_qdim() != 1)
      DAL_THROW(invalid_argument, "invalid data mesh fem (Qdim=1 required)");
    // generic_assembly assem("M$1(#1,#2)+=comp(vGrad(#1).Base(#2))(:,i,i,:);");
    generic_assembly assem("M$1(#1,#2)+=-comp(Base(#1).vGrad(#2))(:,:,i,i);");
    assem.push_mf(mf_p);
    assem.push_mf(mf_u);
    assem.push_mat(B);
    assem.volumic_assembly();
  }


  /**
     assembly of $\int_\Omega a(x)\nabla u.\nabla v$ , where $a(x)$ is scalar.
  */
  template<typename MAT>
    void asm_stiffness_matrix_for_homogeneous_laplacian(const MAT &M_,
							const mesh_fem &mf) {
    MAT &M = const_cast<MAT &>(M_);
    generic_assembly 
      assem("M$1(#1,#1)+=sym(comp(Grad(#1).Grad(#1))(:,i,:,i))");
    assem.push_mf(mf);
    assem.push_mat(M);
    assem.volumic_assembly();
  }


  /**
     assembly of $\int_\Omega a(x)\nabla u.\nabla v$ , where $a(x)$ is scalar.
  */
  template<typename MAT, typename VECT>
    void asm_stiffness_matrix_for_laplacian(MAT &M, const mesh_fem &mf,
					    const mesh_fem &mfdata,
					    const VECT &A) {
    if (mfdata.get_qdim() != 1)
      DAL_THROW(invalid_argument, "invalid data mesh fem (Qdim=1 required)");
    generic_assembly
      assem("a=data$1(#2); M$1(#1,#1)+="
	    "sym(comp(Grad(#1).Grad(#1).Base(#2))(:,i,:,i,j).a(j))");
    //generic_assembly assem("a=data$1(#2); M$1(#1,#1)"
    //        "+=comp(Grad(#1).Grad(#1).Base(#2))(:,i,:,i,j).a(j)");
    assem.push_mf(mf);
    assem.push_mf(mfdata);
    assem.push_data(A);
    assem.push_mat(M);
    assem.volumic_assembly();
  }

  /**
     assembly of $\int_\Omega A(x)\nabla u.\nabla v$, where $A(x)$
     is a NxN matrix.
     Arguments:
      - M  : a sparse matrix of dimensions mf.nb_dof() x mf.nb_dof()

      - mf : the mesh_fem that describes the solution, with
      mf.get_qdim() == N.

      - mfdata : the mesh_fem that describes the coefficients of A
      (mfdata.get_qdim() == 1).

      - A__ : a (very large) vector, which is a flattened (n x n x
      mfdata.nb_dof()) 3D array. For each dof of mfdata, it contains
      the n x n coefficients of A. As usual, the order is the
      "fortran-order", i.e. A__ = [A_11(dof1) A_21(dof1) A_31(dof1)
      A_12(dof1) A_22(dof1) ... A_33(dof) A_11(dof2)
      .... A_33(lastdof)]
  */
  template<typename MAT, typename VECT>
    void asm_stiffness_matrix_for_scalar_elliptic(MAT &M, const mesh_fem &mf,
						  const mesh_fem &mfdata,
						  const VECT &A)
  {
    if (mfdata.get_qdim() != 1)
      DAL_THROW(invalid_argument, "invalid data mesh fem (Qdim=1 required)");
    generic_assembly assem("a=data$1(mdim(#1),mdim(#1),#2); M$1(#1,#1)+=comp(Grad(#1).Grad(#1).Base(#2))(:,i,:,j,k).a(j,i,k)");
    assem.push_mf(mf);
    assem.push_mf(mfdata);
    assem.push_data(A);
    assem.push_mat(M);
    assem.volumic_assembly();
  }


  /** 
      assembly of the term $\int_\Omega Ku.v - \nabla u.\nabla v$, 
      for the helmholtz equation ($\Delta u$ + k^2u = 0$, with $K=k^2$).

      The argument K_squared may be a real or a complex-valued vector.
  */
  template<typename MAT, typename VECT>
  void asm_Helmholtz(MAT &M, const mesh_fem &mf_u, const mesh_fem &mf_data,
		     const VECT &K_squared) {
    asm_Helmholtz(M, mf_u, mf_data, K_squared,
		  typename gmm::linalg_traits<VECT>::value_type());
  }

  template<typename MAT, typename VECT, typename T>
  void asm_Helmholtz(MAT &M, const mesh_fem &mf_u, const mesh_fem &mf_data,
		     const VECT &K_squared, T) {
    asm_Helmholtz_real(M, mf_u, mf_data, K_squared);
  }

  template<typename MAT, typename VECT, typename T>
  void asm_Helmholtz(MAT &M, const mesh_fem &mf_u, const mesh_fem &mf_data,
		     const VECT &K_squared, std::complex<T>) {
    asm_Helmholtz_cplx(gmm::real_part(M), gmm::imag_part(M), mf_u, mf_data,
		       gmm::real_part(K_squared), gmm::imag_part(K_squared));
  }


  template<typename MATr, typename MATi, typename VECTr, typename VECTi>  
  void asm_Helmholtz_cplx(const MATr &Mr, const MATi &Mi, const mesh_fem &mf_u,
			  const mesh_fem &mf_data,
			  const VECTr &K_squaredr, const VECTi &K_squaredi) {
    generic_assembly assem("Kr=data$1(#2); Ki=data$2(#2);"
			   "m = comp(Base(#1).Base(#1).Base(#2)); "
			   "M$1(#1,#1)+=sym(m(:,:,i).Kr(i) - "
			   "comp(Grad(#1).Grad(#1))(:,i,:,i));"
			   "M$2(#1,#1)+=sym(m(:,:,i).Ki(i));");
    assem.push_mf(mf_u);
    assem.push_mf(mf_data);
    assem.push_data(K_squaredr); //gmm::real_part(K_squared));
    assem.push_data(K_squaredi); //gmm::imag_part(K_squared));
    assem.push_mat(const_cast<MATr&>(Mr));
    assem.push_mat(const_cast<MATi&>(Mi));
    assem.volumic_assembly();
  }

  template<typename MAT, typename VECT>  
  void asm_Helmholtz_real(const MAT &M, const mesh_fem &mf_u,
			  const mesh_fem &mf_data, const VECT &K_squared) {
    generic_assembly assem("K=data$1(#2);"
			   "m = comp(Base(#1).Base(#1).Base(#2)); "
			   "M$1(#1,#1)+=sym(m(:,:,i).K(i) - "
			   "comp(Grad(#1).Grad(#1))(:,i,:,i));");
    assem.push_mf(mf_u);
    assem.push_mf(mf_data);
    assem.push_data(K_squared);
    assem.push_mat(const_cast<MAT&>(M));
    assem.volumic_assembly();
  }

  enum { ASMDIR_BUILDH = 1, ASMDIR_BUILDR = 2, ASMDIR_SIMPLIFY = 4,
	 ASMDIR_BUILDALL = 7 };

  /**
     Assembly of Dirichlet constraints h(x)u(x) = r(x), where h is a
     QxQ matrix field (Q == mf_u.get_qdim()), outputs a
     (under-determined) linear system MU=B.

     size(h_data) = Q^2 * nb_dof(mf_rh);
     size(r_data) = Q   * nb_dof(mf_rh);

     This function tries hard to make M diagonal or mostly diagonal:
     this function is able to "simplify" the dirichlet constraints (see below)
     version = |ASMDIR_BUILDH : build H
	       |ASMDIR_BUILDR : build R
	       |ASMDIR_SIMPLIFY : simplify
	       |ASMDIR_BUILDALL : do everything.
  */

  template<typename MAT, typename VECT1, typename VECT2, typename VECT3>
  void asm_dirichlet_constraints(MAT &H, VECT1 &R, const mesh_fem &mf_u,
				 const mesh_fem &mf_rh,
				 const VECT2 &h_data, const VECT3 &r_data,
				 size_type boundary,
				 int version =  ASMDIR_BUILDALL) {
    mesh_cvf_set::face_bitset nf;
    pfem pf_u, pf_rh;
    
    if (mf_rh.get_qdim() != 1)
      DAL_THROW(invalid_argument, "invalid data mesh fem (Qdim=1 required)");
    if (version & ASMDIR_BUILDH) {
      asm_qu_term(H, mf_u, mf_rh, h_data, boundary);
      std::vector<size_type> ind(0);
      dal::bit_vector bdof = mf_u.dof_on_set(boundary);
      for (size_type i = 0; i < mf_u.nb_dof(); ++i)
	if (!(bdof[i])) ind.push_back(i);
      gmm::clear(gmm::sub_matrix(H, gmm::sub_index(ind)));
    }
    if (version & ASMDIR_BUILDR) asm_source_term(R, mf_u, mf_rh, r_data, boundary);
    if (!(version & ASMDIR_SIMPLIFY)) return;

    /* step 2 : simplification of simple dirichlet conditions */
    dal::bit_vector bv = mf_u.linked_mesh().convex_in_set(boundary);
    for (dal::bv_visitor cv(bv); !cv.finished(); ++cv) {
      nf = mf_u.linked_mesh().faces_of_convex_in_set(cv, boundary);
      size_type nbf = mf_u.linked_mesh().structure_of_convex(cv)->nb_faces();
      /* don't try anything with vector elements */
      if (mf_u.fem_of_element(cv)->target_dim() != 1) continue;
      if (nf.count() > 0) {
	pf_u = mf_u.fem_of_element(cv); 
	pf_rh = mf_rh.fem_of_element(cv);
	
	for (size_type f = 0; f < nbf; ++f)
	  if (nf[f]) {
	    bgeot::pconvex_structure cvs_u = pf_u->structure(cv);
	    bgeot::pconvex_structure cvs_rh = pf_rh->structure(cv);
	    for (size_type i = 0; i < cvs_u->nb_points_of_face(f); ++i) {
	      
	      size_type Q = mf_u.get_qdim();  // pf_u->target_dim() (==1)
	      
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
		we replace \int{(H_j.psi_j)*phi_i}=\int{R_j.psi_j} (sum over j)
		with             H_j*phi_i = R_j     
	      */
		if (tdof_u == tdof_rh &&
		    bgeot::vect_dist2_sqr((*(pf_u->node_tab(cv)))[ind_u], 
					  (*(pf_rh->node_tab(cv)))[ind_rh])
		    < 1.0E-14) {
		  /* the dof might be "duplicated" */
		  for (size_type q = 0; q < Q; ++q) {
		    size_type dof_u = mf_u.ind_dof_of_element(cv)[ind_u*Q + q];
		    
		    /* "erase" the row */
		    if (version & ASMDIR_BUILDH)
		      for (size_type k=0; k < mf_u.nb_dof_of_element(cv); ++k)
			H(dof_u, mf_u.ind_dof_of_element(cv)[k]) = 0.0;
		    
		    size_type dof_rh = mf_rh.ind_dof_of_element(cv)[ind_rh];
		    /* set the "simplified" row */
		    if (version & ASMDIR_BUILDH)
		      for (unsigned jj=0; jj < Q; jj++) {
			size_type dof_u2
			  = mf_u.ind_dof_of_element(cv)[ind_u*Q+jj];
			H(dof_u, dof_u2) = h_data[(jj*Q+q) + Q*Q*(dof_rh)];
		      }
		    if (version & ASMDIR_BUILDR) R[dof_u] = r_data[dof_rh*Q+q];
		  }
		}
	      }
	    }
	  }
      }
    }
  }
  
  template<typename MAT, typename VECT>
  void asm_dirichlet_constraints(MAT &H, VECT &R, const mesh_fem &mf_u,
				 const mesh_fem &mf_rh,
				 const VECT &r_data, size_type boundary,
				 int version = ASMDIR_BUILDALL) {
    if (mf_rh.get_qdim() != 1) 
      DAL_THROW(invalid_argument,"mf_rh should be a scalar (qdim=1) mesh_fem");
    size_type N = mf_rh.nb_dof(), Q=mf_u.get_qdim();
    VECT h_data(dal::sqr(mf_u.get_qdim())*N); gmm::clear(H);
    
    for (size_type i=0; i < N; ++i)
      for (size_type q=0; q < Q; ++q)  h_data[i*Q*Q+q*Q+q]=1;
   
    asm_dirichlet_constraints(H, R, mf_u, mf_rh, h_data, r_data, boundary, version);
  }

  /**
     Faster (and simpler) assembly of simple Dirichlet conditions (
     u(x) = F(x) on a boundary) 
     
     The input matrix RM and the right hand side B are modified to
     enforce the Dirichlet condition. The symmetry properties of RM
     are kept.
  */
  template<typename MATRM, typename VECT1, typename VECT2>
    void assembling_Dirichlet_condition(MATRM &RM, VECT1 &B,
					const mesh_fem &mf,
					size_type boundary,
					const VECT2 &F)
  { // Marche uniquement pour des ddl de lagrange.
    size_type Q=mf.get_qdim();
    dal::bit_vector nndof = mf.dof_on_set(boundary);
    pfem pf1;
    for (dal::bv_visitor cv(mf.convex_index()); !cv.finished(); ++cv) {
      pf1 = mf.fem_of_element(cv);
      pdof_description ldof = lagrange_dof(pf1->dim());
      size_type nbd = pf1->nb_dof(cv);
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
  template<typename MATRM, typename VECT1>
    void add_Dirichlet_dof(MATRM &RM, VECT1 &B,
			   const mesh_fem &mf,
			   size_type dof, 
			   typename gmm::linalg_traits<MATRM>::value_type dof_val) {
    size_type Q=mf.get_qdim();
    bgeot::mesh_convex_ind_ct dofcv = mf.convex_to_dof(dof);
    pfem pf1;

    for (bgeot::mesh_convex_ind_ct::const_iterator it = dofcv.begin();
	 it != dofcv.end(); ++it) {
      pf1 = mf.fem_of_element(*it);
      if (pf1->target_dim() != 1)
	DAL_THROW(to_be_done_error, "sorry, to be done ... ");
      size_type nbd = pf1->nb_dof(*it);
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

  /** 
      Build an orthogonal basis of the kernel of H in NS, gives the
      solution of minimal norm of H*U = R in U0 and return the
      dimension of the kernel. The function is based on a
      Gramm-Schmidt algorithm.
  */
   template<typename MAT1, typename MAT2, typename VECT>
  size_type Dirichlet_nullspace(const MAT1 &H, MAT2 &NS,
				const VECT &R, VECT &U0) {

    // To be finalized.
    //  . In order to be used with any sparse matrix type
    //  . transpose the result and give the effective dimension of the kernel
    //  . Compute the ctes / H.
    //  . Optimization (suppress temporary ...). 
    //  . Verify sizes of data

   typedef typename gmm::linalg_traits<MAT1>::value_type T;
    typedef typename gmm::number_traits<T>::magnitude_type MAGT;
    typedef typename gmm::temporary_vector<MAT1>::vector_type TEMP_VECT;
    MAGT tol = gmm::default_tol(MAGT());
    MAGT norminfH = gmm::mat_maxnorm(H);
    size_type nbd = gmm::mat_ncols(H), nbase = 0, nbr = gmm::mat_nrows(H);
    TEMP_VECT aux(nbr), e(nbd), f(nbd);
    dal::dynamic_array<TEMP_VECT> base_img;
    dal::dynamic_array<TEMP_VECT> base_img_inv;
    size_type nb_bimg = 0;

    if (!(gmm::is_col_matrix(H)))
      DAL_WARNING(2, "Dirichlet_nullspace inefficient for a row matrix H");
    // First, detection of null columns of H, and already orthogonals 
    // vectors of the image of H.
    dal::bit_vector nn;
    for (size_type i = 0; i < nbd; ++i) {
      gmm::clear(e); e[i] = T(1);
      gmm::mult(H, e, aux);
      MAGT n = gmm::vect_norm2(aux);

      if (n < 1e-8) { //norminfH*tol) {
	NS(i, nbase++) = T(1); nn[i] = true;
      }
      else {
	bool good = true;
	for (size_type j = 0; j < nb_bimg; ++j)
	  if (dal::abs(gmm::vect_sp(aux, base_img[j])) > MAGT(0))
	    { good = false; break; }
	if (good) {
	  gmm::copy(e, f);
	  gmm::scale(f, T(MAGT(1)/n)); gmm::scale(aux, T(MAGT(1)/n));
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
	gmm::mult(H, e, aux);
	for (size_type j = 0; j < nb_bimg; ++j) { 
	  T c = gmm::vect_sp(aux, base_img[j]);
	  //	  if (dal::abs(c > 1.0E-6) { // à scaler sur l'ensemble de H ...
	  if (c != T(0)) {
	    gmm::add(gmm::scaled(base_img[j], -c), aux);
	    gmm::add(gmm::scaled(base_img_inv[j], -c), f);
	  }
	}
	if (gmm::vect_norm2(aux) < norminfH*tol*MAGT(10000)) {
	  gmm::copy(f, gmm::mat_col(NS, nbase++));
	}
	else {
	  MAGT n = gmm::vect_norm2(aux);
	  gmm::scale(f, T(MAGT(1)/n)); gmm::scale(aux, T(MAGT(1)/n));
	  gmm::clean(f, tol*gmm::vect_norminf(f));
	  gmm::clean(aux, tol*gmm::vect_norminf(aux));
	  base_img_inv[nb_bimg] = TEMP_VECT(nbd);
	  gmm::copy(f, base_img_inv[nb_bimg]);
	  base_img[nb_bimg] = TEMP_VECT(nbr);
	  gmm::copy(aux, base_img[nb_bimg++]);
	}
      }
    // Compute a solution in U0
    gmm::clear(U0);
    for (size_type i = 0; i < nb_bimg; ++i) {
      T c = gmm::vect_sp(base_img[i], R);
      gmm::add(gmm::scaled(base_img_inv[i], c), U0);
    }
    // Orthogonalisation of the basis of the kernel of H.
    for (size_type i = nb_triv_base + 1; i < nbase; ++i) {
      for (size_type j = nb_triv_base; j < i; ++j) {
	T c = gmm::vect_sp(gmm::mat_col(NS,i), gmm::mat_col(NS,j));
	if (c != T(0))
	  gmm::add(gmm::scaled(gmm::mat_col(NS,j), -c), gmm::mat_col(NS,i));
      }
    }
    // projection of U0 on the orthogonal to the kernel.
    for (size_type j = nb_triv_base; j < nbase; ++j) {
      T c = gmm::vect_sp(gmm::mat_col(NS,j), U0);
      if (c != T(0))
	gmm::add(gmm::scaled(gmm::mat_col(NS,j), -c), U0);
    }
    // Test ...
    gmm::mult(H, U0, gmm::scaled(R, T(-1)), aux);
    if (gmm::vect_norm2(aux) > gmm::vect_norm2(U0)*tol*MAGT(10000))
      DAL_WARNING(2, "Dirichlet condition not well inverted: residu="
		  << gmm::vect_norm2(aux));
    return nbase;
  }

}  /* end of namespace getfem.                                             */


#endif /* GETFEM_ASSEMBLING_H__  */
