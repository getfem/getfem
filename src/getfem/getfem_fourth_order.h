/* -*- c++ -*- (enables emacs c++ mode) */
/*===========================================================================

 Copyright (C) 2006-2015 Yves Renard

 This file is a part of GetFEM++

 GetFEM++  is  free software;  you  can  redistribute  it  and/or modify it
 under  the  terms  of the  GNU  Lesser General Public License as published
 by  the  Free Software Foundation;  either version 3 of the License,  or
 (at your option) any later version along with the GCC Runtime Library
 Exception either version 3.1 or (at your option) any later version.
 This program  is  distributed  in  the  hope  that it will be useful,  but
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 License and GCC Runtime Library Exception for more details.
 You  should  have received a copy of the GNU Lesser General Public License
 along  with  this program;  if not, write to the Free Software Foundation,
 Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.

 As a special exception, you  may use  this file  as it is a part of a free
 software  library  without  restriction.  Specifically,  if   other  files
 instantiate  templates  or  use macros or inline functions from this file,
 or  you compile this  file  and  link  it  with other files  to produce an
 executable, this file  does  not  by itself cause the resulting executable
 to be covered  by the GNU Lesser General Public License.  This   exception
 does not  however  invalidate  any  other  reasons why the executable file
 might be covered by the GNU Lesser General Public License.

===========================================================================*/

/**@file getfem_fourth_order.h
   @author  Yves Renard <Yves.Renard@insa-lyon.fr>,
            Julien Pommier <Julien.Pommier@insa-toulouse.fr>
   @date January 6, 2006.
   @brief assembly procedures and bricks for fourth order pdes.
*/
#ifndef GETFEM_FOURTH_ORDER_H_
#define GETFEM_FOURTH_ORDER_H__

#include "getfem_models.h"
#include "getfem_assembling.h"

namespace getfem {
  
  /* ******************************************************************** */
  /*            Bilaplacian assembly routines.                            */
  /* ******************************************************************** */

  /**
     assembly of @f$\int_\Omega \Delta u \Delta v@f$.
     @ingroup asm
  */
  template<typename MAT, typename VECT>
  void asm_stiffness_matrix_for_bilaplacian
  (const MAT &M, const mesh_im &mim, const mesh_fem &mf,
   const mesh_fem &mf_data, const VECT &A,
   const mesh_region &rg = mesh_region::all_convexes()) {
    generic_assembly assem
      ("a=data$1(#2);"
       "M(#1,#1)+=sym(comp(Hess(#1).Hess(#1).Base(#2))(:,i,i,:,j,j,k).a(k))");
    assem.push_mi(mim);
    assem.push_mf(mf);
    assem.push_mf(mf_data);
    assem.push_data(A);
    assem.push_mat(const_cast<MAT &>(M));
    assem.assembly(rg);
  }

  template<typename MAT, typename VECT>
  void asm_stiffness_matrix_for_homogeneous_bilaplacian
  (const MAT &M, const mesh_im &mim, const mesh_fem &mf,
   const VECT &A, const mesh_region &rg = mesh_region::all_convexes()) {
    generic_assembly assem
      ("a=data$1(1);"
       "M(#1,#1)+=sym(comp(Hess(#1).Hess(#1))(:,i,i,:,j,j).a(1))");
    assem.push_mi(mim);
    assem.push_mf(mf);
    assem.push_data(A);
    assem.push_mat(const_cast<MAT &>(M));
    assem.assembly(rg);
  }


  template<typename MAT, typename VECT>
  void asm_stiffness_matrix_for_bilaplacian_KL
  (const MAT &M, const mesh_im &mim, const mesh_fem &mf,
   const mesh_fem &mf_data, const VECT &D_, const VECT &nu_,
   const mesh_region &rg = mesh_region::all_convexes()) {
    generic_assembly assem
      ("d=data$1(#2); n=data$2(#2);"
       "t=comp(Hess(#1).Hess(#1).Base(#2).Base(#2));"
       "M(#1,#1)+=sym(t(:,i,j,:,i,j,k,l).d(k)-t(:,i,j,:,i,j,k,l).d(k).n(l)"
       "+t(:,i,i,:,j,j,k,l).d(k).n(l))");
    assem.push_mi(mim);
    assem.push_mf(mf);
    assem.push_mf(mf_data);
    assem.push_data(D_);
    assem.push_data(nu_);
    assem.push_mat(const_cast<MAT &>(M));
    assem.assembly(rg);
  }

  template<typename MAT, typename VECT>
  void asm_stiffness_matrix_for_homogeneous_bilaplacian_KL
  (const MAT &M, const mesh_im &mim, const mesh_fem &mf,
   const VECT &D_, const VECT &nu_,
   const mesh_region &rg = mesh_region::all_convexes()) {
    generic_assembly assem
      ("d=data$1(1); n=data$2(1);"
       "t=comp(Hess(#1).Hess(#1));"
       "M(#1,#1)+=sym(t(:,i,j,:,i,j).d(1)-t(:,i,j,:,i,j).d(1).n(1)"
       "+t(:,i,i,:,j,j).d(1).n(1))");
    assem.push_mi(mim);
    assem.push_mf(mf);
    assem.push_data(D_);
    assem.push_data(nu_);
    assem.push_mat(const_cast<MAT &>(M));
    assem.assembly(rg);
  }

  /* ******************************************************************** */
  /*            Bilaplacian bricks.                                       */
  /* ******************************************************************** */

  
  /** Adds a bilaplacian brick on the variable
      `varname` and on the mesh region `region`. 
      This represent a term :math:`\Delta(D \Delta u)`. 
      where :math:`D(x)` is a coefficient determined by `dataname` which
      could be constant or described on a f.e.m. The corresponding weak form
      is :math:`\int D(x)\Delta u(x) \Delta v(x) dx`.
  */
  size_type add_bilaplacian_brick
  (model &md, const mesh_im &mim, const std::string &varname,
   const std::string &dataname, size_type region = size_type(-1));

  /** Adds a bilaplacian brick on the variable
      `varname` and on the mesh region `region`.
      This represent a term :math:`\Delta(D \Delta u)` where :math:`D(x)`
      is a the flexion modulus determined by `dataname1`. The term is
      integrated by part following a Kirchhoff-Love plate model
      with `dataname2` the poisson ratio.
  */
  size_type add_bilaplacian_brick_KL
  (model &md, const mesh_im &mim, const std::string &varname,
   const std::string &dataname1, const std::string &dataname2,
   size_type region = size_type(-1));


  /* ******************************************************************** */
  /*            Normale derivative source term assembly routines.         */
  /* ******************************************************************** */

  /**
     assembly of @f$\int_\Gamma{\partial_n u f}@f$.
     @ingroup asm
  */
  template<typename VECT1, typename VECT2>
  void asm_normal_derivative_source_term
  (VECT1 &B, const mesh_im &mim, const mesh_fem &mf, const mesh_fem &mf_data,
   const VECT2 &F, const mesh_region &rg) {
    GMM_ASSERT1(mf_data.get_qdim() == 1,
                "invalid data mesh fem (Qdim=1 required)");

    size_type Q = gmm::vect_size(F) / mf_data.nb_dof();

    const char *s;
    if (mf.get_qdim() == 1 && Q == 1)
      s = "F=data(#2);"
        "V(#1)+=comp(Grad(#1).Normal().Base(#2))(:,i,i,j).F(j);";
    else if (mf.get_qdim() == 1 && Q == gmm::sqr(mf.linked_mesh().dim()))
      s = "F=data(mdim(#1),mdim(#1),#2);"
        "V(#1)+=comp(Grad(#1).Normal().Normal().Normal().Base(#2))"
        "(:,i,i,k,l,j).F(k,l,j);";
    else if (mf.get_qdim() > size_type(1) && Q == mf.get_qdim())
      s = "F=data(qdim(#1),#2);"
        "V(#1)+=comp(vGrad(#1).Normal().Base(#2))(:,i,k,k,j).F(i,j);";
    else if (mf.get_qdim() > size_type(1) &&
             Q == size_type(mf.get_qdim()*gmm::sqr(mf.linked_mesh().dim())))
      s = "F=data(qdim(#1),mdim(#1),mdim(#1),#2);"
        "V(#1)+=comp(vGrad(#1).Normal().Normal().Normal().Base(#2))"
        "(:,i,k,k,l,m,j).F(i,l,m,j);";
    else
      GMM_ASSERT1(false, "invalid rhs vector");
    asm_real_or_complex_1_param(B, mim, mf, mf_data, F, rg, s);
  }

  template<typename VECT1, typename VECT2>
  void asm_homogeneous_normal_derivative_source_term
  (VECT1 &B, const mesh_im &mim, const mesh_fem &mf,
   const VECT2 &F, const mesh_region &rg) {
    
    size_type Q = gmm::vect_size(F);

    const char *s;
    if (mf.get_qdim() == 1 && Q == 1)
      s = "F=data(1);"
        "V(#1)+=comp(Grad(#1).Normal())(:,i,i).F(1);";
    else if (mf.get_qdim() == 1 && Q == gmm::sqr(mf.linked_mesh().dim()))
      s = "F=data(mdim(#1),mdim(#1));"
        "V(#1)+=comp(Grad(#1).Normal().Normal().Normal())"
        "(:,i,i,l,j).F(l,j);";
    else if (mf.get_qdim() > size_type(1) && Q == mf.get_qdim())
      s = "F=data(qdim(#1));"
        "V(#1)+=comp(vGrad(#1).Normal())(:,i,k,k).F(i);";
    else if (mf.get_qdim() > size_type(1) &&
             Q == size_type(mf.get_qdim()*gmm::sqr(mf.linked_mesh().dim())))
      s = "F=data(qdim(#1),mdim(#1),mdim(#1));"
        "V(#1)+=comp(vGrad(#1).Normal().Normal().Normal())"
        "(:,i,k,k,l,m).F(i,l,m);";
    else
      GMM_ASSERT1(false, "invalid rhs vector");
    asm_real_or_complex_1_param(B, mim, mf, mf, F, rg, s);
  }


  /* ******************************************************************** */
  /*                Normale derivative source term brick.                 */
  /* ******************************************************************** */


  /** Adds a normal derivative source term brick
      :math:`F = \int b.\partial_n v` on the variable `varname` and the
      mesh region `region`.
     
      Update the right hand side of the linear system.
      `dataname` represents `b` and `varname` represents `v`.
  */
  size_type add_normal_derivative_source_term_brick
  (model &md, const mesh_im &mim, const std::string &varname,
   const std::string &dataname, size_type region);


  /* ******************************************************************** */
  /*           Special boundary condition for Kirchhoff-Love model.       */
  /* ******************************************************************** */

  /*
    assembly of the special boundary condition for Kirchhoff-Love model.
    @ingroup asm
  */
  template<typename VECT1, typename VECT2>
  void asm_neumann_KL_term
  (VECT1 &B, const mesh_im &mim, const mesh_fem &mf, const mesh_fem &mf_data,
   const VECT2 &M, const VECT2 &divM, const mesh_region &rg) {
    GMM_ASSERT1(mf_data.get_qdim() == 1,
                "invalid data mesh fem (Qdim=1 required)");

    generic_assembly assem
      ("MM=data$1(mdim(#1),mdim(#1),#2);"
       "divM=data$2(mdim(#1),#2);"
       "V(#1)+=comp(Base(#1).Normal().Base(#2))(:,i,j).divM(i,j);"
       "V(#1)+=comp(Grad(#1).Normal().Base(#2))(:,i,j,k).MM(i,j,k)*(-1);"
       "V(#1)+=comp(Grad(#1).Normal().Normal().Normal().Base(#2))(:,i,i,j,k,l).MM(j,k,l);");
    
    assem.push_mi(mim);
    assem.push_mf(mf);
    assem.push_mf(mf_data);
    assem.push_data(M);
    assem.push_data(divM);
    assem.push_vec(B);
    assem.assembly(rg);
  }

  template<typename VECT1, typename VECT2>
  void asm_neumann_KL_homogeneous_term
  (VECT1 &B, const mesh_im &mim, const mesh_fem &mf,
   const VECT2 &M, const VECT2 &divM, const mesh_region &rg) {

    generic_assembly assem
      ("MM=data$1(mdim(#1),mdim(#1));"
       "divM=data$2(mdim(#1));"
       "V(#1)+=comp(Base(#1).Normal())(:,i).divM(i);"
       "V(#1)+=comp(Grad(#1).Normal())(:,i,j).MM(i,j)*(-1);"
       "V(#1)+=comp(Grad(#1).Normal().Normal().Normal())(:,i,i,j,k).MM(j,k);");
    
    assem.push_mi(mim);
    assem.push_mf(mf);
    assem.push_data(M);
    assem.push_data(divM);
    assem.push_vec(B);
    assem.assembly(rg);
  }

  /* ******************************************************************** */
  /*                Kirchhoff Love Neumann term brick.                    */
  /* ******************************************************************** */


  /** Adds a Neumann term brick for Kirchhoff-Love model
      on the variable `varname` and the mesh region `region`.
      `dataname1` represents the bending moment tensor and  `dataname2`
      its divergence.
  */
  size_type add_Kirchhoff_Love_Neumann_term_brick
  (model &md, const mesh_im &mim, const std::string &varname,
   const std::string &dataname1, const std::string &dataname2,
   size_type region);


  /* ******************************************************************** */
  /*            Normal derivative Dirichlet assembly routines.            */
  /* ******************************************************************** */

  /**
     Assembly of normal derivative Dirichlet constraints
     @f$ \partial_n u(x) = r(x) @f$ in a weak form
     @f[ \int_{\Gamma} \partial_n u(x)v(x)=\int_{\Gamma} r(x)v(x) \forall v@f],
     where @f$ v @f$ is in
     the space of multipliers corresponding to mf_mult.

     size(r_data) = Q   * nb_dof(mf_rh);

     version = |ASMDIR_BUILDH : build H
     |ASMDIR_BUILDR : build R
     |ASMDIR_BUILDALL : do everything.

     @ingroup asm
  */

  template<typename MAT, typename VECT1, typename VECT2>
  void asm_normal_derivative_dirichlet_constraints
  (MAT &H, VECT1 &R, const mesh_im &mim, const mesh_fem &mf_u,
   const mesh_fem &mf_mult, const mesh_fem &mf_r,
   const VECT2 &r_data, const mesh_region &rg, bool R_must_be_derivated, 
   int version) {
    typedef typename gmm::linalg_traits<VECT1>::value_type value_type;
    typedef typename gmm::number_traits<value_type>::magnitude_type magn_type;
    
    rg.from_mesh(mim.linked_mesh()).error_if_not_faces();
    
    if (version & ASMDIR_BUILDH) {
      const char *s;
      if (mf_u.get_qdim() == 1 && mf_mult.get_qdim() == 1)
        s = "M(#1,#2)+=comp(Base(#1).Grad(#2).Normal())(:,:,i,i)";
      else
        s = "M(#1,#2)+=comp(vBase(#1).vGrad(#2).Normal())(:,i,:,i,j,j);";
      
      generic_assembly assem(s);
      assem.push_mi(mim);
      assem.push_mf(mf_mult);
      assem.push_mf(mf_u);
      assem.push_mat(H);
      assem.assembly(rg);
      gmm::clean(H, gmm::default_tol(magn_type())
                 * gmm::mat_maxnorm(H) * magn_type(1000));
    }
    if (version & ASMDIR_BUILDR) {
      GMM_ASSERT1(mf_r.get_qdim() == 1,
                "invalid data mesh fem (Qdim=1 required)");
      if (!R_must_be_derivated) {
        asm_normal_source_term(R, mim, mf_mult, mf_r, r_data, rg);
      } else {
        asm_real_or_complex_1_param
          (R, mim, mf_mult, mf_r, r_data, rg,
           "R=data(#2); V(#1)+=comp(Base(#1).Grad(#2).Normal())(:,i,j,j).R(i)");
      }
    }
  }

  /* ******************************************************************** */
  /*            Normal derivative Dirichlet condition bricks.             */
  /* ******************************************************************** */

  /** Adds a Dirichlet condition on the normal derivative of the variable
      `varname` and on the mesh region `region` (which should be a boundary. 
      The general form is
      :math:`\int \partial_n u(x)v(x) = \int r(x)v(x) \forall v`
      where :math:`r(x)` is
      the right hand side for the Dirichlet condition (0 for
      homogeneous conditions) and :math:`v` is in a space of multipliers
      defined by the variable `multname` on the part of boundary determined
      by `region`. `dataname` is an optional parameter which represents
      the right hand side of the Dirichlet condition.
      If `R_must_be_derivated` is set to `true` then the normal
      derivative of `dataname` is considered.
  */
  size_type add_normal_derivative_Dirichlet_condition_with_multipliers
  (model &md, const mesh_im &mim, const std::string &varname,
   const std::string &multname, size_type region,
   const std::string &dataname = std::string(),
   bool R_must_be_derivated = false);
  

  /** Adds a Dirichlet condition on the normal derivative of the variable
      `varname` and on the mesh region `region` (which should be a boundary. 
      The general form is
      :math:`\int \partial_n u(x)v(x) = \int r(x)v(x) \forall v`
      where :math:`r(x)` is
      the right hand side for the Dirichlet condition (0 for
      homogeneous conditions) and :math:`v` is in a space of multipliers
      defined by the trace of mf_mult on the part of boundary determined
      by `region`. `dataname` is an optional parameter which represents
      the right hand side of the Dirichlet condition.
      If `R_must_be_derivated` is set to `true` then the normal
      derivative of `dataname` is considered.
  */
  size_type add_normal_derivative_Dirichlet_condition_with_multipliers
  (model &md, const mesh_im &mim, const std::string &varname,
   const mesh_fem &mf_mult, size_type region,
   const std::string &dataname = std::string(),
   bool R_must_be_derivated = false);

  /** Adds a Dirichlet condition on the normal derivative of the variable
      `varname` and on the mesh region `region` (which should be a boundary. 
      The general form is
      :math:`\int \partial_n u(x)v(x) = \int r(x)v(x) \forall v`
      where :math:`r(x)` is
      the right hand side for the Dirichlet condition (0 for
      homogeneous conditions) and :math:`v` is in a space of multipliers
      defined by the trace of a Lagranfe finite element method of degree
      `degree` and on the boundary determined
      by `region`. `dataname` is an optional parameter which represents
      the right hand side of the Dirichlet condition.
      If `R_must_be_derivated` is set to `true` then the normal
      derivative of `dataname` is considered.
  */
  size_type add_normal_derivative_Dirichlet_condition_with_multipliers
  (model &md, const mesh_im &mim, const std::string &varname,
   dim_type degree, size_type region,
   const std::string &dataname = std::string(),
   bool R_must_be_derivated = false);

    /** Adds a Dirichlet condition on the normal derivative of the variable
      `varname` and on the mesh region `region` (which should be a boundary. 
      The general form is
      :math:`\int \partial_n u(x)v(x) = \int r(x)v(x) \forall v`
      where :math:`r(x)` is
      the right hand side for the Dirichlet condition (0 for
      homogeneous conditions). For this brick the condition is enforced with
      a penalisation with a penanalization parameter `penalization_coeff` on
      the boundary determined by `region`.
      `dataname` is an optional parameter which represents
      the right hand side of the Dirichlet condition.
      If `R_must_be_derivated` is set to `true` then the normal
      derivative of `dataname` is considered.
      Note that is is possible to change the penalization coefficient
      using the function `getfem::change_penalization_coeff` of the standard
      Dirichlet condition.
  */
  size_type add_normal_derivative_Dirichlet_condition_with_penalization
  (model &md, const mesh_im &mim, const std::string &varname,
   scalar_type penalisation_coeff, size_type region, 
   const std::string &dataname = std::string(),
   bool R_must_be_derivated = false);
  





}  /* end of namespace getfem.                                             */


#endif /* GETFEM_FOURTH_ORDER_H__ */
