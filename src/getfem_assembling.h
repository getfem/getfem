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
     compute $\|U\|2_$
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
     compute $\|\nabla U\|2_$
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
    return sqrt(dal::sqr(asm_L2_norm(mf,U,cvlst))
		+dal::sqr(asm_H1_semi_norm(mf,U,cvlst)));
  }
  

  /** 
      generic mass matrix assembly (on the whole mesh or on the specified
      boundary) 
  */
  template<class MAT>
  void asm_mass_matrix(MAT &M, const mesh_fem &mf_u1,
		       size_type boundary=size_type(-1)) {
    asm_mass_matrix(M,mf_u1,mf_u1, boundary);
  }

  /** 
      generic mass matrix assembly (on the whole mesh or on the specified
      boundary) 
  */
  template<class MAT>
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
     source term (for both volumic sources and boundary (neumann) sources  
  */
  template<class VECT1, class VECT2>
    void asm_source_term(VECT1 &B, const mesh_fem &mf,
			 const mesh_fem &mfdata, const VECT2 &F,
			 size_type boundary=size_type(-1)) {
    if (mfdata.get_qdim() != 1)
      DAL_THROW(std::invalid_argument,
		"invalid data mesh fem (Qdim=1 required)");
    generic_assembly assem;
    if (mf.get_qdim() == 1)
      assem.set("F=data(#2); V(#1)+=comp(Base(#1).Base(#2))(:,j).F(j);");
    else
      assem.set("F=data(qdim(#1),#2);"
		"V(#1)+=comp(vBase(#1).Base(#2))(:,i,j).F(i,j);");
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

    Q is a vector, so the matrix is assumed to be stored by columns
    (fortran style)

    Works for both volumic assembly and boundary assembly
  */
  template<class MAT, class VECT>
    void asm_qu_term(MAT &M, 
		     const mesh_fem &mf_u, 
		     const mesh_fem &mf_d, const VECT &Q, 
		     size_type boundary=size_type(-1)) {
    if (mf_d.get_qdim() != 1)
      DAL_THROW(std::invalid_argument,
		"invalid data mesh fem (Qdim=1 required)");
    generic_assembly assem;
    if (mf_u.get_qdim() == 1)
      assem.set("Q=data$1(#2);"
		"M(#1,#1)+=comp(Base(#1).Base(#1).Base(#2))(:,:,k).Q(k);");
    else {
      /* detect the symmetricity of Q (in that case the symmetricity of
       * the final matrix will be ensured, and computations will be
       * slightly speed up */
      bool Q_symmetric = true;
      size_type q = mf_u.get_qdim();
      for (size_type k=0; k < mf_d.nb_dof(); ++k)
	for (size_type i=1; i < q; ++i)
	  for (size_type j=0; j < i; ++j)
	    if (Q[k*q*q+i*q+j] != Q[k*q*q+j*q+i])
	      { Q_symmetric = false; goto bye; }
    bye:
      if (Q_symmetric)
	assem.set("Q=data$1(qdim(#1),qdim(#1),#2);"
		  "M(#1,#1)+=sym(comp(vBase(#1).vBase(#1).Base(#2))(:,i,:,j,k).Q(i,j,k));");
      else
	assem.set("Q=data$1(qdim(#1),qdim(#1),#2);"
		  "M(#1,#1)+=comp(vBase(#1).vBase(#1).Base(#2))(:,i,:,j,k).Q(i,j,k);");
    }
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
    void asm_stiffness_matrix_for_linear_elasticity(const MAT &RM_,
					   const mesh_fem &mf, 
					   const mesh_fem &mfdata, 
					   const VECT &LAMBDA,const VECT &MU) {
    MAT &RM = const_cast<MAT &>(RM_);
    if (mfdata.get_qdim() != 1)
      DAL_THROW(std::invalid_argument,
		"invalid data mesh fem (Qdim=1 required)");
    
    if (mf.get_qdim() != mf.linked_mesh().dim())
      DAL_THROW(std::logic_error, "wrong qdim for the mesh_fem");
    /* e = strain tensor,
       M = 2*mu*e(u):e(v) + lambda*tr(e(u))*tr(e(v))
    */
    generic_assembly assem("lambda=data$1(#2); mu=data$2(#2);"
			   "t=comp(vGrad(#1).vGrad(#1).Base(#2));"
			   //"e=(t{:,2,3,:,5,6,:}+t{:,3,2,:,5,6,:}"
			   //"+t{:,2,3,:,6,5,:}+t{:,3,2,:,6,5,:})/4;"
			   "e=(t{:,2,3,:,5,6,:}+t{:,3,2,:,5,6,:})*0.5;"
			   "M(#1,#1)+= sym(2*e(:,i,j,:,i,j,k).mu(k)"
			   " + e(:,i,i,:,j,j,k).lambda(k))");

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
    if (mfdata.get_qdim() != 1)
      DAL_THROW(std::invalid_argument,
		"invalid data mesh fem (Qdim=1 required)");
    /* e = strain tensor,
       M = a_{i,j,k,l}e_{i,j}(u)e_{k,l}(v)
    */
    generic_assembly assem("a=data$1(qdim(#1),qdim(#1),qdim(#1),qdim(#1),#2);"
			   "t=comp(vGrad(#1).vGrad(#1).Base(#2));"
			   "e=(t{:,2,3,:,5,6,:}+t{:,3,2,:,5,6,:}"
			   "+t{:,2,3,:,6,5,:}+t{:,3,2,:,6,5,:})/4;"
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
    if (mf_d.get_qdim() != 1)
      DAL_THROW(std::invalid_argument, "invalid data mesh fem for asm_stokes (Qdim=1 required)");
    generic_assembly assem("visc=data$1(#3); "
			   "t=comp(vGrad(#1).vGrad(#1).Base(#3));"
			   "e=(t{:,2,3,:,5,6,:}+t{:,3,2,:,5,6,:}+t{:,2,3,:,6,5,:}+t{:,3,2,:,6,5,:})/4;"
			   "M$1(#1,#1) += sym(e(:,i,j,:,i,j,k).visc(k));"          // visc*D(u):D(v)
			   "M$2(#1,#2) += comp(vGrad(#1).Base(#2))(:,i,i,:);");    // p.div v
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
					    const mesh_fem &mfdata,
					    const VECT &A) {
    if (mfdata.get_qdim() != 1)
      DAL_THROW(std::invalid_argument,
		"invalid data mesh fem (Qdim=1 required)");
    generic_assembly assem("a=data$1(#2); M$1(#1,#1)+=sym(comp(Grad(#1).Grad(#1).Base(#2))(:,i,:,i,j).a(j))");
    //generic_assembly assem("a=data$1(#2); M$1(#1,#1)+=comp(Grad(#1).Grad(#1).Base(#2))(:,i,:,i,j).a(j)");
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
  template<class MAT, class VECT>
    void asm_stiffness_matrix_for_scalar_elliptic(MAT &M, const mesh_fem &mf,
						  const mesh_fem &mfdata,
						  const VECT &A)
  {
    if (mfdata.get_qdim() != 1)
      DAL_THROW(std::invalid_argument, "invalid data mesh fem (Qdim=1 required)");
    generic_assembly assem("a=data$1(mdim(#1),mdim(#1),#2); M$1(#1,#1)+=comp(Grad(#1).Grad(#1).Base(#2))(:,i,:,j,k).a(j,i,k)");
    assem.push_mf(mf);
    assem.push_mf(mfdata);
    assem.push_data(A);
    assem.push_mat(M);
    assem.volumic_assembly();
  }

  template<class MAT, class VECT>
  void asm_Helmholtz(MAT &M, const mesh_fem &mf_u, const mesh_fem &mfdata, const VECT &K_squared) {
    generic_assembly assem("Kr=data$1(#2); Ki=data$2(#2);"
			   "m = comp(Base(#1).Base(#1).Base(#2)); "
			   "M$1(#1,#1)+=sym(comp(Grad(#1).Grad(#1))(:,i,:,i) + m(:,:,i).Kr(i)); "
			   "M$2(#1,#1)+=sym(m(:,:,i).Ki(i));");
    assem.push_mf(mf_u);
    assem.push_mf(mfdata);
    assem.push_data(gmm::real_part(K_squared));
    assem.push_data(gmm::imag_part(K_squared));
    assem.push_mat(gmm::real_part(M));
    assem.push_mat(gmm::imag_part(M));
    assem.volumic_assembly();
  }


  /* 
     new version, which takes into account the qdim dimension of the mesh_fem 

     size(H) = Qdim(mf_u)^2 * nb_dof(mf_rh);
     size(R) = Qdim(mf_u)   * nb_dof(mf_rh);

     this function is able to "simplify" the dirichlet constraints (see below)
     version = +1 : build H
	       +2 : build R
	       +4 : simplify
	       (for instance version = 7 do everything).
     TODO : a version without H ?
   */
  template<class MAT, class VECT>
  void asm_dirichlet_constraints(MAT &M, VECT &B, const mesh_fem &mf_u,
				 const mesh_fem &mf_rh,
				 const VECT &H, const VECT &R,
				 size_type boundary, int version = 7) {
    dal::bit_vector nf;
    pfem pf_u, pf_rh;
    
    if (mf_rh.get_qdim() != 1)
      DAL_THROW(std::invalid_argument,
		"invalid data mesh fem for dirichlet (Qdim=1 required)");
    if (version & 1) {
      asm_qu_term(M, mf_u, mf_rh, H, boundary);
      std::vector<size_type> ind(0);
      dal::bit_vector bdof = mf_u.dof_on_boundary(boundary);
      for (size_type i = 0; i < mf_u.nb_dof(); ++i)
	if (!(bdof[i])) ind.push_back(i);
      gmm::clear(gmm::sub_matrix(M, gmm::sub_index(ind)));
    }
    if (version & 2) asm_source_term(B, mf_u, mf_rh, R, boundary);
    if (!(version & 4)) return;

    /* step 2 : simplification of simple dirichlet conditions */
    dal::bit_vector bv = mf_u.convex_on_boundary(boundary);
    for (dal::bv_visitor cv(bv); !cv.finished(); ++cv) {
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
		  bgeot::vect_dist2_sqr(pf_u->node_convex().points()[ind_u], 
					pf_rh->node_convex().points()[ind_rh])
		  < 1.0E-14) {
		/* the dof might be "duplicated" */
		for (size_type q = 0; q < Q; ++q) {
		  size_type dof_u = mf_u.ind_dof_of_element(cv)[ind_u*Q + q];
		  
		  /* "erase" the row */
		  if (version & 1)
		    for (size_type k=0; k < mf_u.nb_dof_of_element(cv); ++k)
		      M(dof_u, mf_u.ind_dof_of_element(cv)[k]) = 0.0;
		  
		  size_type dof_rh = mf_rh.ind_dof_of_element(cv)[ind_rh];
		  /* set the "simplified" row */
		  if (version & 1)
		    for (unsigned jj=0; jj < Q; jj++) {
		      size_type dof_u2
			= mf_u.ind_dof_of_element(cv)[ind_u*Q+jj];
		      M(dof_u, dof_u2) = H[(jj*Q+q) + Q*Q*(dof_rh)];
		    }
		  if (version & 2) B[dof_u] = R[dof_rh*Q+q];
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
				 const VECT &R, size_type boundary,
				 int version = 7) {
    if (mf_rh.get_qdim() != 1) 
      DAL_THROW(std::invalid_argument,
		"mf_rh should be a scalar (qdim=1) mesh_fem");
    size_type N = mf_rh.nb_dof(), Q=mf_u.get_qdim();
    VECT H(dal::sqr(mf_u.get_qdim())*N); gmm::clear(H);
    
    for (size_type i=0; i < N; ++i)
      for (size_type q=0; q < Q; ++q)  H[i*Q*Q+q*Q+q]=1;
   
    asm_dirichlet_constraints(M, B, mf_u, mf_rh, H, R, boundary, version);
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
    R norminfD = gmm::mat_maxnorm(D);
    size_type nbd = gmm::mat_ncols(D), nbase = 0, nbr = gmm::mat_nrows(D);
    TEMP_VECT aux(nbr), e(nbd), f(nbd);
    dal::dynamic_array<TEMP_VECT> base_img;
    dal::dynamic_array<TEMP_VECT> base_img_inv;
    size_type nb_bimg = 0;

    if (!(gmm::is_col_matrix(D)))
      DAL_WARNING(2,
		  "Dirichlet_nullspace is inefficient when D is a row matrix");
    // First, detection of null columns of D, and already orthogonals 
    // vectors of the image of D.
    dal::bit_vector nn;
    for (size_type i = 0; i < nbd; ++i) {
      gmm::clear(e); e[i] = T(1);
      gmm::mult(D, e, aux);
      R n = gmm::vect_norm2(aux);

      if (n < 1e-8) { //norminfD*tol) {
	G(i, nbase++) = T(1); nn[i] = true;
      }
      else {
	bool good = true;
	for (size_type j = 0; j < nb_bimg; ++j)
	  if (dal::abs(gmm::vect_sp(aux, base_img[j])) > R(0))
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
      T c = gmm::vect_sp(base_img[i], UD);
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

}  /* end of namespace getfem.                                             */


#endif /* GETFEM_ASSEMBLING_H__  */
