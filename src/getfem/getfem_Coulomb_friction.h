// -*- c++ -*- (enables emacs c++ mode)
//===========================================================================
//
// Copyright (C) 2004-2010 Yves Renard
//
// This file is a part of GETFEM++
//
// Getfem++  is  free software;  you  can  redistribute  it  and/or modify it
// under  the  terms  of the  GNU  Lesser General Public License as published
// by  the  Free Software Foundation;  either version 2.1 of the License,  or
// (at your option) any later version.
// This program  is  distributed  in  the  hope  that it will be useful,  but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
// License for more details.
// You  should  have received a copy of the GNU Lesser General Public License
// along  with  this program;  if not, write to the Free Software Foundation,
// Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.
//
// As a special exception, you  may use  this file  as it is a part of a free
// software  library  without  restriction.  Specifically,  if   other  files
// instantiate  templates  or  use macros or inline functions from this file,
// or  you compile this  file  and  link  it  with other files  to produce an
// executable, this file  does  not  by itself cause the resulting executable
// to be covered  by the GNU Lesser General Public License.  This   exception
// does not  however  invalidate  any  other  reasons why the executable file
// might be covered by the GNU Lesser General Public License.
//
//===========================================================================

/** @file getfem_Coulomb_friction.h
    @author Yves Renard <Yves.Renard@insa-lyon.fr>
    @author Julien Pommier <Julien.Pommier@insa-toulouse.fr>
    @date July 6, 2004.
    @brief Unilateral contact and Coulomb friction condition brick.
 */
#ifndef GETFEM_COULOMB_FRICTION_H__
#define GETFEM_COULOMB_FRICTION_H__

#include "getfem_modeling.h"
#include "getfem_models.h"

namespace getfem {

 typedef gmm::row_matrix<gmm::rsvector<scalar_type> > CONTACT_B_MATRIX;

  /** Add a frictionless contact condition to the model. If U is the vector
      of degrees of freedom on which the unilateral constraint is applied,
      the matrix `BN` has to be such that this condition is defined by
      $B_N U \le gap$. The constraint is prescribed thank to a multiplier
      `multname_n` whose dimension should be equal to the number of lines of
      `BN`. The augmentation parameter `r` should be chosen in a range of
      acceptabe values (see Getfem user documentation). `dataname_gap` is an
      optional parameter representing the initial gap. It can be a single value
      or a vector of value. `dataname_alpha` is an optional homogenization
      parameter for the augmentation parameter
      (see Getfem user documentation). The parameter `symmetrized` indicates
      that the symmetry of the tangent matrix will be kept or not. 
  */
  size_type add_basic_contact_brick
  (model &md, const std::string &varname_u, const std::string &multname_n,
   const std::string &dataname_r, CONTACT_B_MATRIX &BN,
   std::string dataname_gap = "", std::string dataname_alpha = "",
   bool symmetrized = false, bool Hughes_stabilized = false);


  /** Add a contact with friction condition to the model. If U is the vector
      of degrees of freedom on which the condition is applied,
      the matrix `BN` has to be such that the contact condition is defined
      by $B_N U \le gap$ and `BT` have to be such that the relative tangential
      displacement is $B_T U$. The matrix `BT` should have as many rows as
      `BN` multiplied b $d-1$ where $d$ is the domain dimension.
      The contact condition is prescribed thank to a multiplier
      `multname_n` whose dimension should be equal to the number of rows of
      `BN` and the friction condition by a mutliplier `multname_t` whise size
      should be the number of rows of `BT`.
      The parameter `dataname_friction_coeff` describe the friction
      coefficient. It could be a scalar or a vector describing the
      coefficient on each contact condition. 
      The augmentation parameter
      `r` should be chosen in a range of acceptabe values
      (see Getfem user documentation). `dataname_gap` is an
      optional parameter representing the initial gap. It can be a single value
      or a vector of value. `dataname_alpha` is an optional homogenization
      parameter for the augmentation parameter
      (see Getfem user documentation). The parameter `symmetrized` indicates
      that a part of the symmetry of the tangent matrix will be kept or not
      (except for the coupling bewteen contact and friction). 
  */
  size_type add_basic_contact_with_friction_brick
  (model &md, const std::string &varname_u, const std::string &multname_n,
   const std::string &multname_t, const std::string &dataname_r,
   CONTACT_B_MATRIX &BN, CONTACT_B_MATRIX &BT,
   std::string dataname_friction_coeff, 
   std::string dataname_gap, std::string dataname_alpha,
   bool symmetrized, bool Hughes_stabilized = false) ;

  /** Can be used to change the matrix BN of a basic contact/friction brick
   */
  CONTACT_B_MATRIX &contact_brick_set_BN(model &md, size_type indbrick);
  
  /** Can be used to change the matrix DN of a basic contact/friction brick
   */
  CONTACT_B_MATRIX &contact_brick_set_DN(model &md, size_type indbrick);

 /** Can be used to change the matrix DT of a basic contact/friction brick
   */
  CONTACT_B_MATRIX &contact_brick_set_DT(model &md, size_type indbrick);

  /** Can be used to change the matrix BT of a basic contact/friction brick
   */
  CONTACT_B_MATRIX &contact_brick_set_BT(model &md, size_type indbrick);

/** Add Hughes stabilized frictionless contact condition to the model. If U 
    is the vector of degrees of freedom on which the unilateral constraint is applied, 
    and Lambda the multiplier Vector of contact force.Then Hughes stabilized frictionless
    contact condition is defined by the matrix `BN` and 'DN' have to be such that this
    condition is defined by $B_N U - DN Lambda \le 0$. where 'DN' is the masse matrix 
    relative to stabilzed term.
    The augmentation parameter `r` should be chosen in a range of acceptabe values. 
    `dataname_gap` is an optional parameter representing the initial gap. It can be
    a single value or a vector of value. `dataname_alpha` is an optional homogenization 
    parameter for the augmentation parameter. The parameter `symmetrized` indicates that 
    a part of the symmetry of the tangent matrix will be kept or not
  */
   inline size_type add_Hughes_stab_basic_contact_brick
   (model &md, const std::string &varname_u, const std::string &multname_n,
    const std::string &dataname_r, CONTACT_B_MATRIX &BN, CONTACT_B_MATRIX &DN,
    std::string dataname_gap="", std::string dataname_alpha="",
    bool symmetrized=false) {

    size_type indbrick = add_basic_contact_brick
      (md, varname_u, multname_n, dataname_r, BN,
       dataname_gap, dataname_alpha, symmetrized, true);
    gmm::resize(contact_brick_set_DN(md, indbrick),
                gmm::mat_nrows(DN), gmm::mat_ncols(DN));
    gmm::copy(DN, contact_brick_set_DN(md, indbrick));
    return indbrick;
   } 

  /**  Add Hughes stabilized friction contact condition to the model. If U is the vector
      of degrees of freedom on which the condition is applied,
      the matrix `BN` have to be such that the contact condition is defined
      by $B_N U+DN Lambda \le 0$ (where 'DN' is the masse matrix 
      relative to stabilzed term) and `BT` have to be such that the relative tangential
      displacement is $B_T U$. The matrix `BT` should have as many rows as
      `BN` multiplied b $d-1$ where $d$ is the domain dimension.
      The contact condition is prescribed thank to a multiplier
      `multname_n` whose dimension should be equal to the number of rows of
      `BN` and the friction condition by a mutliplier `multname_t` whise size
      should be the number of rows of `BT`.
      The parameter `dataname_friction_coeff` describe the friction
      coefficient. It could be a scalar or a vector describing the
      coefficient on each contact condition. 
      The augmentation parameter
      `r` should be chosen in a range of acceptabe values
      (see Getfem user documentation). `dataname_gap` is an
      optional parameter representing the initial gap. It can be a single value
      or a vector of value. `dataname_alpha` is an optional homogenization
      parameter for the augmentation parameter
      (see Getfem user documentation). The parameter `symmetrized` indicates
      that a part of the symmetry of the tangent matrix will be kept or not
      (except for the coupling bewteen contact and friction). 
   **/
  inline size_type add_Hughes_stab_friction_contact_brick
  (model &md, const std::string &varname_u, const std::string &multname_n,
   const std::string &multname_t, const std::string &dataname_r,
   CONTACT_B_MATRIX &BN, CONTACT_B_MATRIX &BT, CONTACT_B_MATRIX &DN,CONTACT_B_MATRIX &DT,
   std::string dataname_friction_coeff, 
   std::string dataname_gap="", std::string dataname_alpha="",
   bool symmetrized=false){
  
    size_type indbrick =add_basic_contact_with_friction_brick
      (md, varname_u, multname_n, multname_t, dataname_r, BN, BT,
       dataname_friction_coeff, dataname_gap, dataname_alpha, symmetrized, true);
    gmm::resize(contact_brick_set_DN(md, indbrick),
		gmm::mat_nrows(DN), gmm::mat_ncols(DN));
    gmm::copy(DN, contact_brick_set_DN(md, indbrick));
    
    gmm::resize(contact_brick_set_DT(md, indbrick),
		gmm::mat_nrows(DT), gmm::mat_ncols(DT));
    gmm::copy(DT, contact_brick_set_DT(md, indbrick));
    return indbrick;
  } 



  /** Add a frictionless contact condition with a rigid obstacle
      to the model. The condition is applied on the variable `varname_u`
      on the boundary corresponding to `region`. The rigid obstacle should
      be described with the string `obstacle` being a signed distance to
      the obstacle. This string should be an expression where the coordinates
      are 'x', 'y' in 2D and 'x', 'y', 'z' in 3D. For instance, if the rigid
      obstacle correspond to $z \le 0$, the corresponding signed distance will
      be simply "z". `multname_n` should be a fixed size variable whose size is
      the number of degrees of freedom on boundary `region`. It represents the
      contact equivalent nodal forces. 
      The augmentation parameter `r` should be chosen in a
      range of acceptabe values (close to the Young modulus of the elastic
      body, see Getfem user documentation). The
      parameter `symmetrized` indicates that the symmetry of the tangent
      matrix will be kept or not. Basically, this brick computes the matrix BN
      and the vectors gap and alpha and calls the basic contact brick.
  */
  size_type add_contact_with_rigid_obstacle_brick
  (model &md, const mesh_im &mim, const std::string &varname_u,
   const std::string &multname_n, const std::string &dataname_r,
   size_type region, const std::string &obstacle, bool symmetrized = false);


  /** Add a contact with friction condition with a rigid obstacle
      to the model. The condition is applied on the variable `varname_u`
      on the boundary corresponding to `region`. The rigid obstacle should
      be described with the string `obstacle` being a signed distance to
      the obstacle. This string should be an expression where the coordinates
      are 'x', 'y' in 2D and 'x', 'y', 'z' in 3D. For instance, if the rigid
      obstacle correspond to $z \le 0$, the corresponding signed distance will
      be simply "z". `multname_n` should be a fixed size variable whose size is
      the number of degrees of freedom on boundary `region`. It represents the
      contact equivalent nodal forces. 
      `multname_t` should be a fixed size variable whose size is
      the number of degrees of freedom on boundary `region` multiplied by
      $d-1$ where $d$ is the domain dimension. It represents the
      friction equivalent nodal forces.
      The augmentation parameter `r` should be chosen in a
      range of acceptabe values (close to the Young modulus of the elastic
      body, see Getfem user documentation). `dataname_friction_coeff` is
      the friction coefficient. It could be a scalar or a vector of values
      representing the friction coefficient on each contact node.
      The parameter `symmetrized` indicates that the symmetry of the tangent
      matrix will be kept or not. Basically, this brick computes the matrix BN
      and the vectors gap and alpha and calls the basic contact brick.
  */
  size_type add_contact_with_friction_with_rigid_obstacle_brick
  (model &md, const mesh_im &mim, const std::string &varname_u,
   const std::string &multname_n, const std::string &multname_t,
   const std::string &dataname_r, const std::string &dataname_friction_coeff,
   size_type region, const std::string &obstacle, bool symmetrized); 

  /** Add a frictionless contact condition with a rigid obstacle
      to the model. This brick add a contact which is defined
      in an integral way. Is it the direct approximation of an augmented
      Lagrangian formulation (see Getfem user documentation) defined at the
      continuous level. The advantage should be a better scalability:
      the number of
      Newton iterations should be more or less independent of the mesh size.
      The condition is applied on the variable `varname_u`
      on the boundary corresponding to `region`. The rigid obstacle should
      be described with the data `dataname_obstacle` being a signed distance to
      the obstacle (interpolated on a finite element method).
      `multname_n` should be a fem variable representing the contact stress.
      An inf-sup condition beetween `multname_n` and `varname_u` is required.
      The augmentation parameter `dataname_r` should be chosen in a
      range of acceptabe values.
      The possible value for `option` is 1 for the non-symmetric
      Alart-Curnier version, 2 for the symmetric one and 3 for the
      non-symmetric Alart-Curnier with an additional augmentation.
  */
  size_type add_continuous_contact_with_rigid_obstacle_brick
  (model &md, const mesh_im &mim, const std::string &varname_u,
   const std::string &multname_n, const std::string &dataname_obs,
   const std::string &dataname_r, size_type region, int option); 

  /** Add a contact with friction condition with a rigid obstacle
      to the model. This brick add a contact which is defined
      in an integral way. Is it the direct approximation of an augmented
      Lagrangian formulation (see Getfem user documentation) defined at the
      continuous level. The advantage should be a better scalability:
      the number of
      Newton iterations should be more or less independent of the mesh size.
      The condition is applied on the variable `varname_u`
      on the boundary corresponding to `region`. The rigid obstacle should
      be described with the data `dataname_obstacle` being a signed distance to
      the obstacle (interpolated on a finite element method).
      `multname_n` should be a fem variable representing the contact stress.
      An inf-sup condition beetween `multname_n` and `varname_u` is required.
      The augmentation parameter `dataname_r` should be chosen in a
      range of acceptabe values. `dataname_friction_coeff` is the friction
      coefficient which could be constant or defined on a finite element
      method. 
      The possible value for `option` is 1 for the non-symmetric
      Alart-Curnier version, 2 for the symmetric one and 3 for the
      non-symmetric Alart-Curnier with an additional augmentation.
      'dataname_alpha' and 'dataname_wt' are optional parameters to solve
      dynamical friction problems.
  */
  size_type add_continuous_contact_with_friction_with_rigid_obstacle_brick
  (model &md, const mesh_im &mim, const std::string &varname_u,
   const std::string &multname_n, const std::string &dataname_obs,
   const std::string &dataname_r, const std::string &dataname_friction_coeff,
   size_type region, int option,
   const std::string &dataname_alpha = "",
   const std::string &dataname_wt = "");


  /** Add a frictionless contact condition between two faces of one or two
      elastic bodies. The condition is applied on the variable `varname_u` or
      the variables `varname_u1` and `varname_u2` depending if a single or
      two distinct displacement fields are given. Vectors `rg1` and `rg2`
      contain pairs of regions expected to come in contact with each other. In
      case of a single region per side, `rg1` and `rg2` can be given as normal
      integers. In the single displacement variable case the regions defined in
      both `rg1` and `rg2` refer to the variable `varname_u`. In the case of
      two displacement variables, `rg1` refers to `varname_u1` and `rg2` refers
      to `varname_u2`. `multname_n` should be a fixed size variable whose size
      is the number of degrees of freedom on those regions among the ones
      defined in `rg1` and `rg2` which are characterized as "slaves". It
      represents the contact equivalent nodal forces. The augmentation
      parameter `r` should be chosen in a range of acceptabe values (close to
      the Young modulus of the elastic body, see Getfem user documentation).
      The optional parameters `slave1` and `slave2` declare if the regions
      defined in `rg1` and `rg2` are correspondingly considered as "slaves".
      By default `slave1` is true and `slave2` is false, i.e. `rg1` contains
      the slave surfaces, while 'rg2' the master surfaces. Preferably only
      one of `slave1` and `slave2` is set to true. The parameter `symmetrized`
      indicates that the symmetry of the tangent matrix will be kept or not.
      Basically, this brick computes the matrix BN and the vectors gap and
      alpha and calls the basic contact brick.
  */
  size_type add_nonmatching_meshes_contact_brick
  (model &md, const mesh_im &mim1, const mesh_im &mim2,
   const std::string &varname_u1, const std::string &varname_u2,
   std::string &multname_n, const std::string &dataname_r,
   const std::vector<size_type> &rg1, const std::vector<size_type> &rg2,
   bool slave1=true, bool slave2=false, bool symmetrized=false);

  inline size_type add_nonmatching_meshes_contact_brick
  (model &md, const mesh_im &mim1, const mesh_im &mim2,
   const std::string &varname_u1, const std::string &varname_u2,
   std::string &multname_n, const std::string &dataname_r,
   size_type rg1, size_type rg2, bool slave1=true, bool slave2=false,
   bool symmetrized=false) {

    std::vector<size_type> vrg1(1,rg1);
    std::vector<size_type> vrg2(1,rg2);
    return add_nonmatching_meshes_contact_brick
      (md, mim1, mim2, varname_u1, varname_u2, multname_n, dataname_r,
       vrg1, vrg2, slave1, slave2, symmetrized);
  }

  inline size_type add_nonmatching_meshes_contact_brick
  (model &md, const mesh_im &mim, const std::string &varname_u,
   std::string &multname_n, const std::string &dataname_r,
   const std::vector<size_type> &rg1, const std::vector<size_type> &rg2,
   bool slave1=true, bool slave2=false, bool symmetrized=false) {

    return add_nonmatching_meshes_contact_brick
      (md, mim, mim, varname_u, varname_u, multname_n, dataname_r,
       rg1, rg2, slave1, slave2, symmetrized);
  }

  inline size_type add_nonmatching_meshes_contact_brick
  (model &md, const mesh_im &mim, const std::string &varname_u,
   std::string &multname_n, const std::string &dataname_r,
   size_type rg1, size_type rg2, bool slave1=true, bool slave2=false,
   bool symmetrized=false) {

    std::vector<size_type> vrg1(1,rg1);
    std::vector<size_type> vrg2(1,rg2);
    return add_nonmatching_meshes_contact_brick
      (md, mim, mim, varname_u, varname_u, multname_n, dataname_r,
       vrg1, vrg2, slave1, slave2, symmetrized);
  }


  /** Add a contact with friction condition between two faces of one or two
      elastic bodies. The condition is applied on the variable `varname_u` or
      the variables `varname_u1` and `varname_u2` depending if a single or
      two distinct displacement fields are given. Vectors `rg1` and `rg2`
      contain pairs of regions expected to come in contact with each other. In
      case of a single region per side, `rg1` and `rg2` can be given as normal
      integers. In the single displacement variable case the regions defined in
      both `rg1` and `rg2` refer to the variable `varname_u`. In the case of
      two displacement variables, `rg1` refers to `varname_u1` and `rg2` refers
      to `varname_u2`. `multname_n` should be a fixed size variable whose size
      is the number of degrees of freedom on those regions among the ones
      defined in `rg1` and `rg2` which are characterized as "slaves". It
      represents the contact equivalent nodal normal forces. `multname_t`
      should be a fixed size variable whose size corresponds to the size of
      `multname_n` multiplied by qdim - 1 . It represents the contact
      equivalent nodal tangent (frictional) forces. The augmentation parameter
      `r` should be chosen in a range of acceptabe values (close to the Young
      modulus of the elastic body, see Getfem user documentation). The friction
      coefficient stored in the parameter `friction_coeff` is either a single
      value or a vector of the same size as `multname_n`. The optional
      parameters `slave1` and `slave2` declare if the regions defined in `rg1`
      and `rg2` are correspondingly considered as "slaves". By default `slave1`
      is true and `slave2` is false, i.e. `rg1` contains the slave surfaces,
      while 'rg2' the master surfaces. Preferably only one of `slave1` and
      `slave2` is set to true. The parameter `symmetrized` indicates that the
      symmetry of the tangent matrix will be kept or not.
      Basically, this brick computes the matrices BN and BT as well the vectors
      gap and alpha and calls the basic contact brick.
  */
  size_type add_nonmatching_meshes_contact_with_friction_brick
  (model &md, const mesh_im &mim1, const mesh_im &mim2,
   const std::string &varname_u1, const std::string &varname_u2,
   std::string &multname_n, std::string &multname_t,
   const std::string &dataname_r, const std::string &dataname_friction_coeff,
   const std::vector<size_type> &rg1, const std::vector<size_type> &rg2,
   bool slave1=true, bool slave2=false, bool symmetrized=false);

  inline size_type add_nonmatching_meshes_contact_with_friction_brick
  (model &md, const mesh_im &mim1, const mesh_im &mim2,
   const std::string &varname_u1, const std::string &varname_u2,
   std::string &multname_n, std::string &multname_t,
   const std::string &dataname_r, const std::string &dataname_friction_coeff,
   size_type rg1, size_type rg2, bool slave1=true, bool slave2=false,
   bool symmetrized=false) {

    std::vector<size_type> vrg1(1,rg1);
    std::vector<size_type> vrg2(1,rg2);
    return add_nonmatching_meshes_contact_with_friction_brick
      (md, mim1, mim2, varname_u1, varname_u2, multname_n, multname_t,
       dataname_r, dataname_friction_coeff,
       vrg1, vrg2, slave1, slave2, symmetrized);
  }

  inline size_type add_nonmatching_meshes_contact_with_friction_brick
  (model &md, const mesh_im &mim, const std::string &varname_u,
   std::string &multname_n, std::string &multname_t,
   const std::string &dataname_r, const std::string &dataname_friction_coeff,
   const std::vector<size_type> &rg1, const std::vector<size_type> &rg2,
   bool slave1=true, bool slave2=false, bool symmetrized=false) {

    return add_nonmatching_meshes_contact_with_friction_brick
      (md, mim, mim, varname_u, varname_u, multname_n, multname_t,
       dataname_r, dataname_friction_coeff,
       rg1, rg2, slave1, slave2, symmetrized);
  }

  inline size_type add_nonmatching_meshes_contact_with_friction_brick
  (model &md, const mesh_im &mim, const std::string &varname_u,
   std::string &multname_n, std::string &multname_t,
   const std::string &dataname_r, const std::string &dataname_friction_coeff,
   size_type rg1, size_type rg2, bool slave1=true, bool slave2=false,
   bool symmetrized=false) {

    std::vector<size_type> vrg1(1,rg1);
    std::vector<size_type> vrg2(1,rg2);
    return add_nonmatching_meshes_contact_with_friction_brick
      (md, mim, mim, varname_u, varname_u, multname_n, multname_t,
       dataname_r, dataname_friction_coeff,
       vrg1, vrg2, slave1, slave2, symmetrized);
  }

  class friction_nonlinear_term : public nonlinear_elem_term {
    
  public:
    dim_type N;
    const mesh_fem &mf_u;
    const mesh_fem &mf_lambda;
    const mesh_fem &mf_obs;
    const mesh_fem *mf_coeff;
    std::vector<scalar_type> U, lambda_n, obs, friction_coeff, WT;
    bgeot::multi_index sizes_;

    base_small_vector no, aux1, lt, lambda, auxN, zt, V, coeff;
    // base_vector coeff, V;
    base_matrix grad, GP;
    scalar_type ln, un, g, r, alpha, f_coeff;
    bool contact_only;
    size_type option;

    template <class VECT> friction_nonlinear_term
    (const mesh_fem &mf_u_, const VECT &U_, 
     const mesh_fem &mf_lambda_, const VECT &lambda_n_, 
     const mesh_fem &mf_obs_, const VECT &obs_,
     scalar_type r__, size_type option_, bool contact_only_ = true,
     scalar_type alpha_ = scalar_type(-1), const mesh_fem *mf_coeff_ = 0,
     const VECT *f_coeff_ = 0, const VECT *WT_ = 0) :
      N(mf_u_.linked_mesh().dim()), mf_u(mf_u_), mf_lambda(mf_lambda_),
      mf_obs(mf_obs_), U(mf_u.nb_basic_dof()),
      lambda_n(mf_lambda_.nb_basic_dof()), obs(mf_obs_.nb_basic_dof()),
      r(r__), alpha(alpha_), contact_only(contact_only_), option(option_) {
      
      sizes_.resize(1); sizes_[0] = 1;
      switch (option) {
      case 0: case 2: case 4: case 5: case 7: case 11: case 13: case 19: case 21: 
	sizes_[0] = N; break;
      case 1: case 3: case 10: case 12: case 15: case 17: case 18: case 20:
	sizes_.resize(2); sizes_[0] = sizes_[1] = N;  break;
      }
      
      mf_u.extend_vector(U_, U);
      mf_lambda.extend_vector(lambda_n_, lambda_n);
      mf_obs.extend_vector(obs_, obs);
      
      V.resize(N); no.resize(N); aux1.resize(1); auxN.resize(N);
      lambda.resize(N); lt.resize(N); zt.resize(N);
      gmm::resize(grad, 1, N);
      gmm::resize(GP, N, N);

      if (!contact_only) {
	mf_coeff = mf_coeff_;
	if (!mf_coeff)
	  f_coeff = (*f_coeff_)[0];
	else {
	  friction_coeff.resize(mf_coeff->nb_basic_dof());
	  mf_coeff->extend_vector(*f_coeff_, friction_coeff);
	}
	if (WT_ && gmm::vect_size(*WT_)) {
	  WT.resize(mf_u.nb_basic_dof());
	  mf_u.extend_vector(*WT_, WT);
	}
	else
	  WT.resize(0);
      }
      
      GMM_ASSERT1(mf_u.get_qdim() == N, "wrong qdim for the mesh_fem");
    }
    
    const bgeot::multi_index &sizes() const { return sizes_; }
    
    virtual void compute(fem_interpolation_context&, bgeot::base_tensor &t);
    virtual void prepare(fem_interpolation_context& ctx, size_type nb);

  };


  /** Specific assembly procedure for the use of an Uzawa algorithm to solve
      contact problems.
  */
  template<typename VECT1> 
  void asm_Coulomb_friction_continuous_Uzawa_proj
  (VECT1 &R, const mesh_im &mim, const getfem::mesh_fem &mf_u,
   const VECT1 &U, const getfem::mesh_fem &mf_lambda, const VECT1 &lambda_n,
   const getfem::mesh_fem &mf_obs, const VECT1 &obs, scalar_type r,
   const mesh_region &rg = mesh_region::all_convexes()) {
    
    friction_nonlinear_term nterm1(mf_u, U, mf_lambda, lambda_n, mf_obs,
				   obs, r, 9);

    getfem::generic_assembly assem;
    assem.set("V(#2)+=comp(NonLin$1(#1,#1,#2,#3).Base(#2))(i,:); ");
    assem.push_mi(mim);
    assem.push_mf(mf_u);
    assem.push_mf(mf_lambda);
    assem.push_mf(mf_obs);
    assem.push_nonlinear_term(&nterm1);
    assem.push_vec(R);
    assem.assembly(rg); 
  }


//===========================================================================
//
//  Brick for the old brick system
//
//===========================================================================

# define MDBRICK_COULOMB_FRICTION 434245

  /**
   * Unilateral contact and Coulomb friction condition brick.
   * (for conformal small displacement problems)
   * @ingroup bricks
   */
  template<typename MODEL_STATE = standard_model_state>
  class mdbrick_Coulomb_friction : public mdbrick_abstract<MODEL_STATE>  {
    
    TYPEDEF_MODEL_STATE_TYPES;

    mdbrick_abstract<MODEL_STATE> &sub_problem;
    size_type num_fem;

    T_MATRIX BN, BT;
    typedef gmm::row_matrix<gmm::rsvector<value_type> > RT_MATRIX;
    typedef gmm::dense_matrix<bool> CH_MATRIX;
    RT_MATRIX AUG_M; // For Hansbo augmentation.
    CH_MATRIX CH_M; // For determining the Jacobian manually; only for 2D.
    VECTOR gap, threshold, WT, WN, friction_coef, RLN, RLT;
    value_type r, alpha, beta;
    size_type d, nbc;

    const mesh_fem *mf_u;
    gmm::sub_interval SUBU, SUBN, SUBT;
    
    bool Tresca_version, symmetrized, contact_only, really_stationary;

    template<typename VEC> static void ball_projection(const VEC &x,
                                                       value_type radius) {
      value_type a = gmm::vect_norm2(x);
      if (radius <= 0) gmm::clear(const_cast<VEC&>(x));
      else if (a > radius) gmm::scale(const_cast<VEC&>(x), radius/a); 
    }

    template<class VEC, class VECR>
    static void ball_projection_grad_r(const VEC &x, value_type radius,
                                       VECR &g) {
      value_type a = gmm::vect_norm2(x);
      if (radius > 0 && a >= radius)
        gmm::copy(gmm::scaled(x, value_type(1)/a), g);
      else gmm::clear(g);
    }

    template <class VEC, class MAT>
    static void ball_projection_grad(const VEC &x, double radius, MAT &g) {
      if (radius <= value_type(0)) { gmm::clear(g); return; }
      gmm::copy(gmm::identity_matrix(), g);
      value_type a = gmm::vect_norm2(x);
      if (a >= radius) { 
        gmm::scale(g, radius/a);
        // gmm::rank_one_update(g, gmm::scaled(x, -radius/(a*a*a)), x);
        for (size_type i = 0; i < x.size(); ++i)
          for (size_type j = 0; j < x.size(); ++j)
            g(i,j) -= radius*x[i]*x[j] / (a*a*a);
      }
    }

    void precomp(MODEL_STATE &MS, size_type i0) {
      size_type i1 = this->mesh_fem_positions[num_fem];
      gmm::resize(RLN, gmm::mat_nrows(BN));
      gmm::resize(RLT, gmm::mat_nrows(BT));
      SUBU = gmm::sub_interval(i0 + i1, mf_u->nb_dof());
      SUBN = gmm::sub_interval(i0 + sub_problem.nb_dof(), gmm::mat_nrows(BN));
      SUBT = gmm::sub_interval(i0 + sub_problem.nb_dof() + gmm::mat_nrows(BN),
                               gmm::mat_nrows(BT));
      gmm::add(gmm::sub_vector(MS.state(), SUBN), gmm::scaled(gap, r*alpha), RLN);
      if (gmm::vect_size(WN) > 0)
        gmm::mult_add(BN, gmm::scaled(WN, -r*alpha), RLN);
      gmm::mult_add(BN, gmm::scaled(gmm::sub_vector(MS.state(), SUBU),
                                    -r*alpha), RLN);
      if (gmm::mat_nrows(AUG_M) > 0)
        gmm::mult_add(AUG_M, gmm::scaled(gmm::sub_vector(MS.state(),SUBN),-r),
                      RLN);
      if (!contact_only) {
        gmm::copy(gmm::sub_vector(MS.state(), SUBT), RLT);
        if (gmm::vect_size(WT) > 0)
          gmm::mult_add(BT, gmm::scaled(WT, -r*beta), RLT);
        if (!really_stationary)
          gmm::mult_add(BT, gmm::scaled(gmm::sub_vector(MS.state(), SUBU),
                                        -r*beta), RLT);
      }
    }

    void proper_update(void) {
      mf_u = this->mesh_fems[num_fem];
      d = mf_u->linked_mesh().dim();
      gmm::resize(BN, nbc, mf_u->nb_dof());
      gmm::resize(BT, nbc*(d-1), mf_u->nb_dof());
      gmm::resize(gap, nbc); gmm::resize(friction_coef, nbc);
      gmm::resize(threshold, nbc);
      // gmm::resize(WT, mf_u->nb_dof()); gmm::resize(WN, mf_u->nb_dof());
      this->proper_additional_dof = gmm::mat_nrows(BN)
        + (contact_only ? 0 : gmm::mat_nrows(BT));
      this->proper_mixed_variables.clear();
      this->proper_mixed_variables.add(sub_problem.nb_dof(),
                                       this->proper_additional_dof);
    }

  public :
    
    inline size_type nb_contact_nodes(void) const
    { return gmm::mat_nrows(BN); }
    
    virtual void do_compute_tangent_matrix(MODEL_STATE &MS, size_type i0,
                                           size_type) {
      precomp(MS, i0);
      
      RT_MATRIX BBN(gmm::mat_nrows(BN), gmm::mat_ncols(BN));
      RT_MATRIX MM(nb_contact_nodes(), nb_contact_nodes());
      gmm::copy(gmm::scaled(BN, -alpha), BBN);
      if (gmm::mat_nrows(AUG_M) > 0)
        gmm::copy(gmm::scaled(AUG_M, -value_type(1)), MM);
      for (size_type i=0; i < nb_contact_nodes(); ++i) {
	if (gmm::mat_nrows(CH_M) > 0) {
	  if (!CH_M(i, 0)) {
	    gmm::clear(BBN[i]);
	    if (gmm::mat_nrows(AUG_M) > 0) gmm::clear(MM[i]);
	    MM(i, i) = -value_type(1)/r;
	  }
	}
	else if (RLN[i] >= value_type(0)) {
	  gmm::clear(BBN[i]);
	  if (gmm::mat_nrows(AUG_M) > 0) gmm::clear(MM[i]);
	  MM(i, i) = -value_type(1)/r;
	}
      }
      gmm::copy(BBN, gmm::sub_matrix(MS.tangent_matrix(), SUBN, SUBU));
      gmm::copy(MM, gmm::sub_matrix(MS.tangent_matrix(), SUBN));

//       gmm::copy(gmm::scaled(BN, -alpha),
//              gmm::sub_matrix(MS.tangent_matrix(), SUBN, SUBU));
//       gmm::clear(gmm::sub_matrix(MS.tangent_matrix(), SUBN));
//       if (gmm::mat_nrows(AUG_M) > 0)
//      gmm::copy(gmm::scaled(AUG_M, -value_type(1)),
//                gmm::sub_matrix(MS.tangent_matrix(), SUBN));
//       for (size_type i=0; i < nb_contact_nodes(); ++i) {
//      if (RLN[i] >= value_type(0)) {
//        gmm::clear(gmm::sub_matrix(MS.tangent_matrix(),
//                                   gmm::sub_interval(SUBN.first()+i,1),
//                                   SUBU));
//        if (gmm::mat_nrows(AUG_M) > 0)
//          gmm::clear(gmm::sub_matrix(MS.tangent_matrix(),
//                                gmm::sub_interval(SUBN.first()+i,1), SUBN));
//        MS.tangent_matrix()(SUBN.first()+i, SUBN.first()+i)=-value_type(1)/r;
//      }
//       }
      
      if (!contact_only) {
        base_matrix pg(d-1, d-1);
        base_vector vg(d-1);

        RT_MATRIX BBT(gmm::mat_nrows(BT), gmm::mat_ncols(BT));
        gmm::dense_matrix<value_type> BTi(d-1,  gmm::mat_ncols(BT));
        
        for (size_type i=0; i < nb_contact_nodes(); ++i) {
          gmm::sub_interval SUBI(i*(d-1), d-1);
          gmm::sub_interval SUBJ(SUBT.first()+i*(d-1),(d-1));
          gmm::sub_interval SUBJJ(i*(d-1),(d-1));
          value_type th = Tresca_version ? threshold[i]
            : - (MS.state())[SUBN.first()+i] * friction_coef[i];
	  std::vector<double> rlt_CH(1);

	  if (mat_nrows(CH_M) > 0) {
	    if (!CH_M(i, 0)) th = 0.0;
	    else th = 1.0;
	    if (!CH_M(i, 1)) rlt_CH[0] = -2.0;
	    else if (!CH_M(i, 2)) rlt_CH[0] = 2.0;
	    else rlt_CH[0] = 0.0;
	    ball_projection_grad(rlt_CH, th, pg);
          } else
	    ball_projection_grad(gmm::sub_vector(RLT, SUBI), th, pg);
          if (!really_stationary)
            gmm::mult(gmm::scaled(pg, -beta), 
                      gmm::sub_matrix(BT, SUBI,
                                      gmm::sub_interval(0, gmm::mat_ncols(BT))),
                      BTi);
          gmm::copy(BTi, gmm::sub_matrix(BBT, SUBJJ, SUBU));

          if (!Tresca_version) {
	    if (mat_nrows(CH_M) > 0)
	      ball_projection_grad_r(rlt_CH, th, vg);
	    else
	      ball_projection_grad_r(gmm::sub_vector(RLT, SUBI), th, vg);
            for (size_type k = 0; k < d-1; ++k)
              MS.tangent_matrix()(SUBT.first()+i*(d-1)+k, SUBN.first()+i)
                = - friction_coef[i] * vg[k] / r;
          }
          for (size_type j = 0; j < d-1; ++j) pg(j,j) -= value_type(1);
          gmm::copy(gmm::scaled(pg,value_type(1)/r), 
                    gmm::sub_matrix(MS.tangent_matrix(), SUBJ));
        }
        T_MATRIX BBBT(gmm::mat_nrows(BT), gmm::mat_ncols(BT));
        gmm::copy(BBT, BBBT);
        gmm::copy(BBBT, gmm::sub_matrix(MS.tangent_matrix(), SUBT, SUBU));
      }

//       if (!contact_only) {
//      base_matrix pg(d-1, d-1);
//      base_vector vg(d-1);
        
//      for (size_type i=0; i < nb_contact_nodes(); ++i) {
//        gmm::sub_interval SUBI(i*(d-1), d-1);
//        gmm::sub_interval SUBJ(SUBT.first()+i*(d-1),(d-1));
//        value_type th = Tresca_version ? threshold[i]
//          : - (MS.state())[SUBN.first()+i] * friction_coef[i];
          
//        ball_projection_grad(gmm::sub_vector(RLT, SUBI), th, pg);
//        if (!really_stationary)
//          gmm::mult(gmm::scaled(pg, -beta), 
//                    gmm::sub_matrix(BT, SUBI,
//                                    gmm::sub_interval(0,gmm::mat_ncols(BT))),
//                    gmm::sub_matrix(MS.tangent_matrix(), SUBJ, SUBU));
          
//        if (!Tresca_version) {
//          ball_projection_grad_r(gmm::sub_vector(RLT, SUBI), th, vg);
//          for (size_type k = 0; k < d-1; ++k)
//            MS.tangent_matrix()(SUBT.first()+i*(d-1)+k, SUBN.first()+i)
//              = - friction_coef[i] * vg[k] / r;
//        }
//        for (size_type j = 0; j < d-1; ++j) pg(j,j) -= value_type(1);
//        gmm::copy(gmm::scaled(pg,value_type(1)/r), 
//                  gmm::sub_matrix(MS.tangent_matrix(), SUBJ));
//      }
//       }

      
      if (symmetrized) {
        T_MATRIX tmp(mf_u->nb_dof(), mf_u->nb_dof());
        
        gmm::resize(tmp, mf_u->nb_dof(), gmm::mat_nrows(BN));
        gmm::copy(gmm::transposed(gmm::sub_matrix(MS.tangent_matrix(),
                                                  SUBN, SUBU)), tmp);
        gmm::copy(tmp, gmm::sub_matrix(MS.tangent_matrix(), SUBU, SUBN));
        gmm::resize(tmp, mf_u->nb_dof(), mf_u->nb_dof());
        gmm::mult(gmm::transposed(gmm::scaled(BN,-r*alpha)),
                  gmm::sub_matrix(MS.tangent_matrix(), SUBN, SUBU), tmp);
        gmm::add(tmp, gmm::sub_matrix(MS.tangent_matrix(), SUBU));
        
        if (!contact_only) {
          gmm::mult(gmm::transposed(gmm::scaled(BT,-r*beta)),
                    gmm::sub_matrix(MS.tangent_matrix(), SUBT, SUBU), tmp);
          gmm::add(tmp, gmm::sub_matrix(MS.tangent_matrix(), SUBU));
          gmm::resize(tmp, mf_u->nb_dof(), gmm::mat_nrows(BT));
          gmm::copy(gmm::transposed(gmm::sub_matrix(MS.tangent_matrix(),
                                                    SUBT, SUBU)), tmp);
          gmm::copy(tmp, gmm::sub_matrix(MS.tangent_matrix(), SUBU, SUBT));
        }
      }
      else {
        gmm::copy(gmm::scaled(gmm::transposed(BN), value_type(-1)),
                  gmm::sub_matrix(MS.tangent_matrix(), SUBU, SUBN));
        if (!contact_only)
          gmm::copy(gmm::scaled(gmm::transposed(BT), value_type(-1)), 
                    gmm::sub_matrix(MS.tangent_matrix(), SUBU, SUBT));
      }
    }
    
    virtual void do_compute_residual(MODEL_STATE &MS, size_type i0,size_type) {
      precomp(MS, i0);
      value_type c1(1);
      
      for (size_type i=0; i < nb_contact_nodes(); ++i) {
        RLN[i] = std::min(value_type(0), RLN[i]);
        if (!contact_only)
          ball_projection(gmm::sub_vector(RLT, gmm::sub_interval(i*(d-1),d-1)),
                          Tresca_version ? threshold[i]
                          : -friction_coef[i]*(MS.state())[SUBN.first()+i]);
      }
      
      if (symmetrized) {
        gmm::mult_add(gmm::transposed(BN), gmm::scaled(RLN, -c1),
                      gmm::sub_vector(MS.residual(), SUBU));
        if (!contact_only)
          gmm::mult_add(gmm::transposed(BT), gmm::scaled(RLT, -c1),
                        gmm::sub_vector(MS.residual(), SUBU));
      } else {
        gmm::mult_add(gmm::transposed(BN),
                      gmm::scaled(gmm::sub_vector(MS.state(), SUBN),-c1),
                      gmm::sub_vector(MS.residual(), SUBU));
        if (!contact_only)
          gmm::mult_add(gmm::transposed(BT),
                        gmm::scaled(gmm::sub_vector(MS.state(), SUBT),-c1),
                        gmm::sub_vector(MS.residual(), SUBU));
      }
      
      /* residual on LN */
      gmm::add(gmm::scaled(gmm::sub_vector(MS.state(), SUBN), -c1/r),
               gmm::scaled(RLN, c1/r), gmm::sub_vector(MS.residual(), SUBN));
      // gmm::scale(gmm::sub_vector(MS.residual(), SUBN), c1 / r);

      /* residual on LT */
      if (!contact_only) {
        gmm::add(gmm::scaled(gmm::sub_vector(MS.state(),SUBT), -c1/r),
                 gmm::scaled(RLT, c1/r), gmm::sub_vector(MS.residual(), SUBT));
        // gmm::scale(gmm::sub_vector(MS.residual(), SUBT), c1 / r);
      }
    }

    void init(void) {
      this->add_sub_brick(sub_problem);
      this->proper_is_linear_ = this->proper_is_coercive_ = false;
      this->proper_is_symmetric_ = symmetrized && contact_only;
      r = value_type(1);
      beta = value_type(1);
      alpha = value_type(1);
      this->force_update();
    }

    void set_stationary(bool b) { really_stationary = b; }
    void set_beta(value_type a) { beta = a; }
    value_type get_beta(void) const { return beta; }
    void set_alpha(value_type a) { alpha = a; }
    value_type get_alpha(void) const { return alpha; }
    template<typename MAT> void set_augmented_matrix(const MAT &M) {
      gmm::resize(AUG_M, gmm::mat_nrows(M), gmm::mat_ncols(M));
      gmm::copy(M, AUG_M);
    }

    void clear_character_matrix(void) { resize(CH_M, 0, 0); }
    template<typename MAT> void set_character_matrix(const MAT &M) {
      gmm::resize(CH_M, gmm::mat_nrows(M), gmm::mat_ncols(M));
      gmm::copy(M, CH_M);
    }

    void set_r(value_type r_) { r = r_; }
    value_type get_r(void) const { return r; }
    template <class VEC> void set_WT(const VEC &WT_)
    { gmm::resize(WT, gmm::vect_size(WT_)); gmm::copy(WT_, WT); }
    template <class VEC> void set_WN(const VEC &WN_)
    { gmm::resize(WN, gmm::vect_size(WN_)); gmm::copy(WN_, WN); }

    VECTOR &get_gap(void) { return gap; }
    const VECTOR &get_gap(void) const { return gap; }

    VECTOR &get_friction_coef(void) { return friction_coef; }
    const VECTOR &get_friction_coef(void) const { return friction_coef; }

    SUBVECTOR get_LN(MODEL_STATE &MS) {
      SUBN = gmm::sub_interval(this->first_index() + sub_problem.nb_dof(),
                               gmm::mat_nrows(BN));
      return gmm::sub_vector(MS.state(), SUBN);
    }

    SUBVECTOR get_LT(MODEL_STATE &MS) {
      SUBT = gmm::sub_interval(this->first_index() + sub_problem.nb_dof()
                               + gmm::mat_nrows(BN),  gmm::mat_nrows(BT));
      return gmm::sub_vector(MS.state(), SUBT);
    }

    /* contact and friction */
    template <class MAT, class VEC> mdbrick_Coulomb_friction
    (mdbrick_abstract<MODEL_STATE> &problem, const MAT &BN_, const VEC &gap_,
     scalar_type FC_, const MAT &BT_, size_type num_fem_=0)
      : sub_problem(problem), num_fem(num_fem_) {
      contact_only = Tresca_version = symmetrized = really_stationary = false;
      nbc = gmm::mat_nrows(BN_);
      init();
      gmm::copy(BN_, BN); gmm::copy(BT_, BT); gmm::copy(gap_, gap);
      std::fill(friction_coef.begin(), friction_coef.end(), FC_);
    }

    /* contact only.        */
    template <class MAT, class VEC> mdbrick_Coulomb_friction
    (mdbrick_abstract<MODEL_STATE> &problem, const MAT &BN_, const VEC &gap_,
     size_type num_fem_=0) : sub_problem(problem), num_fem(num_fem_) {
      contact_only = true;
      Tresca_version = symmetrized = really_stationary = false;
      nbc = gmm::mat_nrows(BN_);
      init();
      gmm::copy(BN_, BN); gmm::copy(gap_, gap);
    }

  };


}  /* end of namespace getfem.                                             */


#endif /* GETFEM_COULOMB_FRICTION_H__ */
