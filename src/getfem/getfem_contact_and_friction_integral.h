/* -*- c++ -*- (enables emacs c++ mode) */
/*===========================================================================

 Copyright (C) 2011-2015 Yves Renard, Konstantinos Poulios.

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

/** @file getfem_contact_and_friction_integral.h
    @author Yves Renard <Yves.Renard@insa-lyon.fr>
    @author Konstantinos Poulios <logari81@googlemail.com>
    @date November, 2011.
    @brief Unilateral contact and Coulomb friction condition brick.
 */
#ifndef GETFEM_CONTACT_AND_FRICTION_INTEGRAL_H__
#define GETFEM_CONTACT_AND_FRICTION_INTEGRAL_H__

#include "getfem_models.h"
#include "getfem_assembling_tensors.h"

namespace getfem {


  /** Add a frictionless contact condition with a rigid obstacle
      to the model, which is defined in an integral way. It is the direct
      approximation of an augmented  Lagrangian formulation (see Getfem user
      documentation) defined at the continuous level. The advantage should be
      a better scalability: the number of Newton iterations should be more or
      less independent of the mesh size.
      The condition is applied on the variable `varname_u`
      on the boundary corresponding to `region`. The rigid obstacle should
      be described with the data `dataname_obstacle` being a signed distance to
      the obstacle (interpolated on a finite element method).
      `multname_n` should be a fem variable representing the contact stress.
      An inf-sup condition between `multname_n` and `varname_u` is required.
      The augmentation parameter `dataname_r` should be chosen in a
      range of acceptable values.
      Possible values for `option` is 1 for the non-symmetric Alart-Curnier
      augmented Lagrangian method, 2 for the symmetric one, 3 for the
      non-symmetric Alart-Curnier method with an additional augmentation.
      The default value is 1.
  */
  size_type add_integral_contact_with_rigid_obstacle_brick
  (model &md, const mesh_im &mim, const std::string &varname_u,
   const std::string &multname_n, const std::string &dataname_obs,
   const std::string &dataname_r, size_type region, int option = 1);

  /** Add a contact with friction condition with a rigid obstacle
      to the model, which is defined in an integral way. It is the direct
      approximation of an augmented  Lagrangian formulation (see Getfem user
      documentation) defined at the continuous level. The advantage should be
      a better scalability: the number of Newton iterations should be more or
      less independent of the mesh size.
      The condition is applied on the variable `varname_u`
      on the boundary corresponding to `region`. The rigid obstacle should
      be described with the data `dataname_obstacle` being a signed distance
      to the obstacle (interpolated on a finite element method).
      `multname` should be a fem variable representing the contact and
      friction stress.
      An inf-sup condition between `multname` and `varname_u` is required.
      The augmentation parameter `dataname_r` should be chosen in a
      range of acceptable values.
      The parameter `dataname_friction_coeffs` contains the Coulomb friction
      coefficient and optionally an adhesional shear stress threshold and the
      tresca limit shear stress. For constant coefficients its size is from
      1 to 3. For coefficients described on a finite element method, this
      vector contains a number of single values, value pairs or triplets
      equal to the number of the corresponding mesh_fem's basic dofs.
      Possible values for `option` is 1 for the non-symmetric Alart-Curnier
      augmented Lagrangian method, 2 for the symmetric one, 3 for the
      non-symmetric Alart-Curnier method with an additional augmentation
      and 4 for a new unsymmetric method. The default value is 1.
      (Option 4 ignores any adhesional stress and tresca limit coefficients.)
      `dataname_alpha` and `dataname_wt` are optional parameters to solve
      evolutionary friction problems. `dataname_gamma` and `dataname_vt`
      represent optional data for adding a parameter-dependent sliding
      velocity to the friction condition.
  */
  size_type add_integral_contact_with_rigid_obstacle_brick
  (model &md, const mesh_im &mim, const std::string &varname_u,
   const std::string &multname, const std::string &dataname_obs,
   const std::string &dataname_r, const std::string &dataname_friction_coeffs,
   size_type region, int option = 1, const std::string &dataname_alpha = "",
   const std::string &dataname_wt = "",
   const std::string &dataname_gamma = "",
   const std::string &dataname_vt = "");


  /** Add a penalized contact frictionless condition with a rigid obstacle
      to the model.
      The condition is applied on the variable `varname_u`
      on the boundary corresponding to `region`. The rigid obstacle should
      be described with the data `dataname_obstacle` being a signed distance to
      the obstacle (interpolated on a finite element method).
      The penalization parameter `dataname_r` should be chosen
      large enough to prescribe an approximate non-penetration condition
      but not too large not to deteriorate too much the conditionning of
      the tangent system. `dataname_n` is an optional parameter used if option
      is 2. In that case, the penalization term is shifted by lambda_n (this
      allows the use of an Uzawa algorithm on the corresponding augmented
      Lagrangian formulation)
  */
  size_type add_penalized_contact_with_rigid_obstacle_brick
  (model &md, const mesh_im &mim, const std::string &varname_u,
   const std::string &dataname_obs, const std::string &dataname_r,
   size_type region, int option = 1, const std::string &dataname_lambda_n = "");

  /** Add a penalized contact condition with Coulomb friction with a
      rigid obstacle to the model.
      The condition is applied on the variable `varname_u`
      on the boundary corresponding to `region`. The rigid obstacle should
      be described with the data `dataname_obstacle` being a signed distance to
      the obstacle (interpolated on a finite element method).
      The parameter `dataname_friction_coeffs` contains the Coulomb friction
      coefficient and optionally an adhesional shear stress threshold and the
      tresca limit shear stress. For constant coefficients its size is from
      1 to 3. For coefficients described on a finite element method, this
      vector contains a number of single values, value pairs or triplets
      equal to the number of the corresponding mesh_fem's basic dofs.
      The penalization parameter `dataname_r` should be chosen
      large enough to prescribe approximate non-penetration and friction
      conditions but not too large not to deteriorate too much the
      conditionning of the tangent system.
      `dataname_lambda` is an optional parameter used if option
      is 2. In that case, the penalization term is shifted by lambda (this
      allows the use of an Uzawa algorithm on the corresponding augmented
      Lagrangian formulation)
      `dataname_alpha` and `dataname_wt` are optional parameters to solve
      evolutionary friction problems.
  */
  size_type add_penalized_contact_with_rigid_obstacle_brick
  (model &md, const mesh_im &mim, const std::string &varname_u,
   const std::string &dataname_obs, const std::string &dataname_r,
   const std::string &dataname_friction_coeffs,
   size_type region, int option = 1, const std::string &dataname_lambda = "",
   const std::string &dataname_alpha = "",
   const std::string &dataname_wt = "");


  /** Add a frictionless contact condition between nonmatching meshes
      to the model. This brick adds a contact which is defined
      in an integral way. It is the direct approximation of an augmented
      Lagrangian formulation (see Getfem user documentation) defined at the
      continuous level. The advantage should be a better scalability:
      the number of Newton iterations should be more or less independent
      of the mesh size.
      The condition is applied on the variables `varname_u1` and `varname_u2`
      on the boundaries corresponding to `region1` and `region2`.
      `multname_n` should be a fem variable representing the contact stress.
      An inf-sup condition between `multname_n` and `varname_u1` and
      `varname_u2` is required.
      The augmentation parameter `dataname_r` should be chosen in a
      range of acceptable values.
      Possible values for `option` is 1 for the non-symmetric Alart-Curnier
      augmented Lagrangian method, 2 for the symmetric one, 3 for the
      non-symmetric Alart-Curnier method with an additional augmentation.
      The default value is 1.
  */
  size_type add_integral_contact_between_nonmatching_meshes_brick
  (model &md, const mesh_im &mim, const std::string &varname_u1,
   const std::string &varname_u2, const std::string &multname_n,
   const std::string &dataname_r,
   size_type region1, size_type region2, int option = 1);

  /** Add a contact with friction condition between nonmatching meshes
      to the model. This brick adds a contact which is defined
      in an integral way. It is the direct approximation of an augmented
      Lagrangian formulation (see Getfem user documentation) defined at the
      continuous level. The advantage should be a better scalability:
      the number of Newton iterations should be more or less independent
      of the mesh size.
      The condition is applied on the variables `varname_u1` and `varname_u2`
      on the boundaries corresponding to `region1` and `region2`.
      `multname` should be a fem variable representing the contact and
      friction stress.
      An inf-sup condition between `multname` and `varname_u1` and
      `varname_u2` is required.
      The augmentation parameter `dataname_r` should be chosen in a
      range of acceptable values.
      The parameter `dataname_friction_coeffs` contains the Coulomb friction
      coefficient and optionally an adhesional shear stress threshold and the
      tresca limit shear stress. For constant coefficients its size is from
      1 to 3. For coefficients described on a finite element method on the
      same mesh as `varname_u1`, this vector contains a number of single
      values, value pairs or triplets equal to the number of the
      corresponding mesh_fem's basic dofs.
      Possible values for `option` is 1 for the non-symmetric Alart-Curnier
      augmented Lagrangian method, 2 for the symmetric one, 3 for the
      non-symmetric Alart-Curnier method with an additional augmentation
      and 4 for a new unsymmetric method. The default value is 1.
      (Option 4 ignores any adhesional stress and tresca limit coefficients.)
      `dataname_alpha`, `dataname_wt1` and `dataname_wt2` are optional
      parameters to solve evolutionary friction problems.
  */
  size_type add_integral_contact_between_nonmatching_meshes_brick
  (model &md, const mesh_im &mim, const std::string &varname_u1,
   const std::string &varname_u2, const std::string &multname,
   const std::string &dataname_r, const std::string &dataname_friction_coeffs,
   size_type region1, size_type region2, int option = 1,
   const std::string &dataname_alpha = "",
   const std::string &dataname_wt1 = "",
   const std::string &dataname_wt2 = "");


  /** Add a penalized contact frictionless condition between nonmatching
      meshes to the model.
      The condition is applied on the variables `varname_u1` and  `varname_u2`
      on the boundaries corresponding to `region1` and `region2`.
      The penalization parameter `dataname_r` should be chosen
      large enough to prescribe an approximate non-penetration condition
      but not too large not to deteriorate too much the conditionning of
      the tangent system. `dataname_n` is an optional parameter used if option
      is 2. In that case, the penalization term is shifted by lambda_n (this
      allows the use of an Uzawa algorithm on the corresponding augmented
      Lagrangian formulation)
  */
  size_type add_penalized_contact_between_nonmatching_meshes_brick
  (model &md, const mesh_im &mim, const std::string &varname_u1,
   const std::string &varname_u2, const std::string &dataname_r,
   size_type region1, size_type region2,
   int option = 1, const std::string &dataname_lambda_n = "");


  /** Add a penalized contact condition with Coulomb friction between
      nonmatching meshes to the model.
      The condition is applied on the variables `varname_u1` and  `varname_u2`
      on the boundaries corresponding to `region1` and `region2`.
      The penalization parameter `dataname_r` should be chosen
      large enough to prescribe approximate non-penetration and friction
      conditions but not too large not to deteriorate too much the
      conditionning of the tangent system.
      The parameter `dataname_friction_coeffs` contains the Coulomb friction
      coefficient and optionally an adhesional shear stress threshold and the
      tresca limit shear stress. For constant coefficients its size is from
      1 to 3. For coefficients described on a finite element method on the
      same mesh as `varname_u1`, this vector contains a number of single
      values, value pairs or triplets equal to the number of the
      corresponding mesh_fem's basic dofs.
      `dataname_lambda` is an optional parameter used if option
      is 2. In that case, the penalization term is shifted by lambda (this
      allows the use of an Uzawa algorithm on the corresponding augmented
      Lagrangian formulation)
      `dataname_alpha`, `dataname_wt1` and `dataname_wt2` are optional parameters
      to solve evolutionary friction problems.
  */
  size_type add_penalized_contact_between_nonmatching_meshes_brick
  (model &md, const mesh_im &mim, const std::string &varname_u1,
   const std::string &varname_u2, const std::string &dataname_r,
   const std::string &dataname_friction_coeffs,
   size_type region1, size_type region2, int option = 1,
   const std::string &dataname_lambda = "",
   const std::string &dataname_alpha = "",
   const std::string &dataname_wt1 = "",
   const std::string &dataname_wt2 = "");


  enum contact_nonlinear_term_version {  RHS_L_V1,
                                         RHS_L_V2,
                                         K_LL_V1,
                                         K_LL_V2,
                                         UZAWA_PROJ,
                                         CONTACT_FLAG,
                                         CONTACT_PRESSURE,

                                         RHS_U_V1,
                                         RHS_U_V2,
                                         RHS_U_V4,
                                         RHS_U_V5,
                                         RHS_U_FRICT_V1,
                                         RHS_U_FRICT_V4,
                                         RHS_U_FRICT_V6,
                                         RHS_U_FRICT_V7,
                                         RHS_U_FRICT_V8,
                                         RHS_U_FRICT_V5,
                                         RHS_L_FRICT_V1,
                                         RHS_L_FRICT_V2,
                                         RHS_L_FRICT_V4,
                                         K_UL_V1,
                                         K_UL_V2,
                                         K_UL_V3,
                                         UZAWA_PROJ_FRICT,
                                         UZAWA_PROJ_FRICT_SAXCE,

                                         K_UU_V1,
                                         K_UU_V2,
                                         K_UL_FRICT_V1, // negative EYE
                                         K_UL_FRICT_V2,
                                         K_UL_FRICT_V3,
                                         K_UL_FRICT_V4,
                                         K_UL_FRICT_V5,
                                         K_UL_FRICT_V7,
                                         K_UL_FRICT_V8,
                                         K_LL_FRICT_V1,
                                         K_LL_FRICT_V2,
                                         K_LL_FRICT_V3,
                                         K_LL_FRICT_V4,
                                         K_UU_FRICT_V1,
                                         K_UU_FRICT_V2,
                                         K_UU_FRICT_V3,
                                         K_UU_FRICT_V4,
                                         K_UU_FRICT_V5
  };

  class contact_nonlinear_term : public nonlinear_elem_term {

  protected:
    // the following variables are used in the compute method and their values
    // have to be calculated during the preceding calls to the prepare method
    // at the current interpolation context
    base_small_vector lnt, lt; // multiplier lambda and its tangential component lambda_t
    scalar_type ln;            // normal component lambda_n of the multiplier
    base_small_vector zt;      // tangential relative displacement
    scalar_type un;            // normal relative displacement (positive when the first
                               // elastic body surface moves outwards)
    base_small_vector no;      // surface normal, pointing outwards with respect
                               // to the (first) elastic body
    scalar_type g;             // gap value
    scalar_type f_coeff;       // coefficient of friction value
    scalar_type tau_adh;       // adhesional shear resistance of the interface
    scalar_type tresca_lim;    // tresca shear limit for the interface

    // these variables are used as temporary storage and they will usually contain
    // garbage from previous calculations
    base_small_vector aux1, auxN, V; // helper vectors of size 1, N and N respectively
    base_matrix GP;                  // helper matrix of size NxN

    void adjust_tensor_size(void);

  public:
    dim_type N;
    size_type option;
    scalar_type r;
    bool contact_only;
    scalar_type alpha;

    bgeot::multi_index sizes_;

    contact_nonlinear_term(dim_type N_, size_type option_, scalar_type r_,
                           bool contact_only_ = true,
                           scalar_type alpha_ = scalar_type(1)) :
      tau_adh(0), tresca_lim(gmm::default_max(scalar_type())),
      N(N_), option(option_), r(r_), contact_only(contact_only_), alpha(alpha_) {

      adjust_tensor_size();
    }

    const bgeot::multi_index &sizes(size_type) const { return sizes_; }

    virtual void friction_law(scalar_type p, scalar_type &tau);
    virtual void friction_law(scalar_type p, scalar_type &tau, scalar_type &tau_grad);

    virtual void compute(fem_interpolation_context&, bgeot::base_tensor &t);
    virtual void prepare(fem_interpolation_context& /*ctx*/, size_type /*nb*/)
    { GMM_ASSERT1(false, "the prepare method has to be reimplemented in "
                         "a derived class"); }

  };


  class contact_rigid_obstacle_nonlinear_term : public contact_nonlinear_term {

    // temporary variables to be used inside the prepare method
    base_small_vector vt; // of size N
    base_vector coeff; // of variable size
    base_matrix grad;  // of size 1xN

  public:
    // class specific objects to take into account inside the prepare method
    const mesh_fem &mf_u;       // mandatory
    const mesh_fem &mf_obs;     // mandatory
    const mesh_fem *pmf_lambda; // optional for terms involving lagrange multipliers
    const mesh_fem *pmf_coeff;  // optional for terms involving fem described coefficient of friction
    base_vector U, obs, lambda, friction_coeff, tau_adhesion, tresca_limit, WT, VT;
    scalar_type gamma;

    template <typename VECT1>
    contact_rigid_obstacle_nonlinear_term
    (size_type option_, scalar_type r_,
     const mesh_fem &mf_u_, const VECT1 &U_,
     const mesh_fem &mf_obs_, const VECT1 &obs_,
     const mesh_fem *pmf_lambda_ = 0, const VECT1 *lambda_ = 0,
     const mesh_fem *pmf_coeff_ = 0, const VECT1 *f_coeffs_ = 0,
     scalar_type alpha_ = scalar_type(1), const VECT1 *WT_ = 0,
     scalar_type gamma_ = scalar_type(1), const VECT1 *VT_ = 0
    )
      : contact_nonlinear_term(mf_u_.linked_mesh().dim(), option_, r_,
                               (f_coeffs_ == 0), alpha_
                              ),
        mf_u(mf_u_), mf_obs(mf_obs_),
        pmf_lambda(pmf_lambda_), pmf_coeff(pmf_coeff_),
        U(mf_u.nb_basic_dof()), obs(mf_obs.nb_basic_dof()),
        lambda(0), friction_coeff(0), tau_adhesion(0), tresca_limit(0),
        WT(0), VT(0), gamma(gamma_)
    {

      mf_u.extend_vector(U_, U);
      mf_obs.extend_vector(obs_, obs);

      if (pmf_lambda) {
        lambda.resize(pmf_lambda->nb_basic_dof());
        pmf_lambda->extend_vector(*lambda_, lambda);
      }

      if (!contact_only) {
        if (!pmf_coeff) {
          f_coeff = (*f_coeffs_)[0];
          if (gmm::vect_size(*f_coeffs_) > 1) tau_adh = (*f_coeffs_)[1];
          if (gmm::vect_size(*f_coeffs_) > 2) tresca_lim = (*f_coeffs_)[2];
        }
        else {
          size_type ncoeffs = gmm::vect_size(*f_coeffs_)/pmf_coeff->nb_dof();
          GMM_ASSERT1(ncoeffs >= 1 && ncoeffs <= 3, "Wrong vector dimension for friction coefficients");
          gmm::resize(friction_coeff, pmf_coeff->nb_basic_dof());
          pmf_coeff->extend_vector(gmm::sub_vector(*f_coeffs_, gmm::sub_slice(0,pmf_coeff->nb_dof(),ncoeffs)),
                                   friction_coeff);
          if (ncoeffs > 1) {
            gmm::resize(tau_adhesion, pmf_coeff->nb_basic_dof());
            pmf_coeff->extend_vector(gmm::sub_vector(*f_coeffs_, gmm::sub_slice(1,pmf_coeff->nb_dof(),ncoeffs)),
                                     tau_adhesion);
          }
          if (ncoeffs > 2) {
            gmm::resize(tresca_limit, pmf_coeff->nb_basic_dof());
            pmf_coeff->extend_vector(gmm::sub_vector(*f_coeffs_, gmm::sub_slice(2,pmf_coeff->nb_dof(),ncoeffs)),
                                     tresca_limit);
          }
        }

        if (WT_ && gmm::vect_size(*WT_)) {
          WT.resize(mf_u.nb_basic_dof());
          mf_u.extend_vector(*WT_, WT);
        }

        if (VT_ && gmm::vect_size(*VT_)) {
          VT.resize(mf_u.nb_basic_dof());
          mf_u.extend_vector(*VT_, VT);
        }
      }

      vt.resize(N);
      gmm::resize(grad, 1, N);

      GMM_ASSERT1(mf_u.get_qdim() == N, "wrong qdim for the mesh_fem");
    }

    // this methode prepares all necessary data for the compute method
    // of the base class
    virtual void prepare(fem_interpolation_context& ctx, size_type nb);

  };


  class contact_nonmatching_meshes_nonlinear_term : public contact_nonlinear_term {

    // temporary variables to be used inside the prepare method
    base_vector coeff; // of variable size
    base_matrix grad;  // of size 1xN

  public:
    const mesh_fem &mf_u1;      // displacements on the non-mortar side
    const mesh_fem &mf_u2;      // displacements of the mortar side projected on the non-mortar side
    const mesh_fem *pmf_lambda; // Lagrange multipliers defined on the non-mortar side
    const mesh_fem *pmf_coeff;  // coefficient of friction defined on the non-mortar side
    base_vector U1, U2, lambda, friction_coeff, tau_adhesion, tresca_limit, WT1, WT2;

    template <typename VECT1>
    contact_nonmatching_meshes_nonlinear_term
    (size_type option_, scalar_type r_,
     const mesh_fem &mf_u1_, const VECT1 &U1_,
     const mesh_fem &mf_u2_, const VECT1 &U2_,
     const mesh_fem *pmf_lambda_ = 0, const VECT1 *lambda_ = 0,
     const mesh_fem *pmf_coeff_ = 0, const VECT1 *f_coeffs_ = 0,
     scalar_type alpha_ = scalar_type(1),
     const VECT1 *WT1_ = 0, const VECT1 *WT2_ = 0
    )
      : contact_nonlinear_term(mf_u1_.linked_mesh().dim(), option_, r_,
                               (f_coeffs_ == 0), alpha_
                              ),
        mf_u1(mf_u1_), mf_u2(mf_u2_),
        pmf_lambda(pmf_lambda_), pmf_coeff(pmf_coeff_),
        U1(mf_u1.nb_basic_dof()), U2(mf_u2.nb_basic_dof()),
        lambda(0), friction_coeff(0), tau_adhesion(0), tresca_limit(0),
        WT1(0), WT2(0)
    {

      GMM_ASSERT1(mf_u2.linked_mesh().dim() == N,
                  "incompatible mesh dimensions for the given mesh_fem's");

      mf_u1.extend_vector(U1_, U1);
      mf_u2.extend_vector(U2_, U2);

      if (pmf_lambda) {
        lambda.resize(pmf_lambda->nb_basic_dof());
        pmf_lambda->extend_vector(*lambda_, lambda);
      }

      if (!contact_only) {
        if (!pmf_coeff) {
          f_coeff = (*f_coeffs_)[0];
          if (gmm::vect_size(*f_coeffs_) > 1) tau_adh = (*f_coeffs_)[1];
          if (gmm::vect_size(*f_coeffs_) > 2) tresca_lim = (*f_coeffs_)[2];
        }
        else {
          size_type ncoeffs = gmm::vect_size(*f_coeffs_)/pmf_coeff->nb_dof();
          GMM_ASSERT1(ncoeffs >= 1 && ncoeffs <= 3, "Wrong vector dimension for friction coefficients");
          gmm::resize(friction_coeff, pmf_coeff->nb_basic_dof());
          pmf_coeff->extend_vector(gmm::sub_vector(*f_coeffs_, gmm::sub_slice(0,pmf_coeff->nb_dof(),ncoeffs)),
                                   friction_coeff);
          if (ncoeffs > 1) {
            gmm::resize(tau_adhesion, pmf_coeff->nb_basic_dof());
            pmf_coeff->extend_vector(gmm::sub_vector(*f_coeffs_, gmm::sub_slice(1,pmf_coeff->nb_dof(),ncoeffs)),
                                     tau_adhesion);
          }
          if (ncoeffs > 2) {
            gmm::resize(tresca_limit, pmf_coeff->nb_basic_dof());
            pmf_coeff->extend_vector(gmm::sub_vector(*f_coeffs_, gmm::sub_slice(2,pmf_coeff->nb_dof(),ncoeffs)),
                                     tresca_limit);
          }
        }
        if (WT1_ && WT2_ && gmm::vect_size(*WT1_) && gmm::vect_size(*WT2_)) {
          WT1.resize(mf_u1.nb_basic_dof());
          mf_u1.extend_vector(*WT1_, WT1);
          WT2.resize(mf_u2.nb_basic_dof());
          mf_u2.extend_vector(*WT2_, WT2);
        }
      }

      gmm::resize(grad, 1, N);

      GMM_ASSERT1(mf_u1.get_qdim() == N, "wrong qdim for the 1st mesh_fem");
      GMM_ASSERT1(mf_u2.get_qdim() == N, "wrong qdim for the 2nd mesh_fem");
    }

    // this method prepares all necessary data for the compute method
    // of the base class
    virtual void prepare(fem_interpolation_context& ctx, size_type nb);

  };


  /** Specific assembly procedure for the use of an Uzawa algorithm to solve
      contact with rigid obstacle problems with friction.
  */
  template<typename VECT1>
  void asm_integral_contact_Uzawa_proj
  (VECT1 &R, const mesh_im &mim,
   const getfem::mesh_fem &mf_u, const VECT1 &U,
   const getfem::mesh_fem &mf_obs, const VECT1 &obs,
   const getfem::mesh_fem &mf_lambda, const VECT1 &lambda,
   const getfem::mesh_fem *pmf_coeff, const VECT1 &f_coeff, const VECT1 *WT,
   scalar_type r, scalar_type alpha, const mesh_region &rg, int option = 1) {

    contact_rigid_obstacle_nonlinear_term
      nterm1((option == 1) ? UZAWA_PROJ_FRICT : UZAWA_PROJ_FRICT_SAXCE, r,
             mf_u, U, mf_obs, obs, &mf_lambda, &lambda,
             pmf_coeff, &f_coeff, alpha, WT);

    getfem::generic_assembly assem;
    if (pmf_coeff) // variable coefficient of friction
      assem.set("V(#3)+=comp(NonLin$1(#1,#1,#2,#3,#4).vBase(#3))(i,:,i); ");
    else // constant coefficient of friction
      assem.set("V(#3)+=comp(NonLin$1(#1,#1,#2,#3).vBase(#3))(i,:,i); ");
    assem.push_mi(mim);
    assem.push_mf(mf_u);
    assem.push_mf(mf_obs);
    assem.push_mf(mf_lambda);
    if (pmf_coeff)
      assem.push_mf(*pmf_coeff);
    assem.push_nonlinear_term(&nterm1);
    assem.push_vec(R);
    assem.assembly(rg);
  }

  /** Specific assembly procedure for the use of an Uzawa algorithm to solve
      contact problems.
  */
  template<typename VECT1>
  void asm_integral_contact_Uzawa_proj
  (VECT1 &R, const mesh_im &mim,
   const getfem::mesh_fem &mf_u, const VECT1 &U,
   const getfem::mesh_fem &mf_obs, const VECT1 &obs,
   const getfem::mesh_fem &mf_lambda, const VECT1 &lambda,
   scalar_type r, const mesh_region &rg) {

    contact_rigid_obstacle_nonlinear_term
      nterm1(UZAWA_PROJ, r, mf_u, U, mf_obs, obs, &mf_lambda, &lambda);

    getfem::generic_assembly assem;
    assem.set("V(#3)+=comp(NonLin$1(#1,#1,#2,#3).Base(#3))(i,:); ");
    assem.push_mi(mim);
    assem.push_mf(mf_u);
    assem.push_mf(mf_obs);
    assem.push_mf(mf_lambda);
    assem.push_nonlinear_term(&nterm1);
    assem.push_vec(R);
    assem.assembly(rg);
  }


  /** Specific assembly procedure for the use of an Uzawa algorithm to solve
      contact problems.
  */
  template<typename VEC>
  void asm_level_set_normal_source_term
  (VEC &R, const mesh_im &mim,
   const getfem::mesh_fem &mf_u,
   const getfem::mesh_fem &mf_obs, const VEC &obs,
   const getfem::mesh_fem &mf_lambda, const VEC &lambda,
   const mesh_region &rg) {

    bool contact_only = (mf_lambda.get_qdim() == 1);

    VEC U;
    gmm::resize(U, mf_u.nb_dof());
    scalar_type dummy_r(0.);
    VEC dummy_f_coeff; gmm::resize(dummy_f_coeff,1);
    contact_rigid_obstacle_nonlinear_term
      nterm(RHS_U_V1, dummy_r, mf_u, U, mf_obs, obs, &mf_lambda, &lambda,
            0, contact_only ? 0 : &dummy_f_coeff);

    getfem::generic_assembly assem;
    assem.set("V(#1)+=comp(NonLin$1(#1,#1,#2,#3).vBase(#1))(i,:,i); ");
    assem.push_mi(mim);
    assem.push_mf(mf_u);
    assem.push_mf(mf_obs);
    assem.push_mf(mf_lambda);
    assem.push_nonlinear_term(&nterm);
    assem.push_vec(R);
    assem.assembly(rg);
  }

  template<typename VEC>
  scalar_type asm_level_set_contact_area
  (const mesh_im &mim,
   const getfem::mesh_fem &mf_u, const VEC &U,
   const getfem::mesh_fem &mf_obs, const VEC &obs,
   const mesh_region &rg, scalar_type threshold_factor=0.0,
   const getfem::mesh_fem *mf_lambda=0, const VEC *lambda=0,
   scalar_type threshold_pressure_factor=0.0) {

    if (!rg.get_parent_mesh())
      rg.from_mesh(mim.linked_mesh());
    getfem::mesh_fem mf_ca(mf_u.linked_mesh());
    mf_ca.set_classical_finite_element(rg.index(),1);

    VEC mesh_size(mf_ca.nb_dof());
    VEC mesh_size2(mf_ca.nb_dof());
    { // assemble an estimator of the mesh size
      getfem::generic_assembly assem_mesh_size;
      assem_mesh_size.set("V(#1)+=comp(Base(#1))");
      assem_mesh_size.push_mi(mim);
      assem_mesh_size.push_mf(mf_ca);
      assem_mesh_size.push_vec(mesh_size2);
      assem_mesh_size.assembly(rg);
      for (dal::bv_visitor_c dof(mf_ca.basic_dof_on_region(rg));
           !dof.finished(); ++dof)
        mesh_size[dof] = sqrt(mesh_size2[dof]);
    }

    VEC threshold(mf_ca.nb_dof());
    if (mf_lambda && lambda) {
      VEC pressure(mf_ca.nb_dof());
      VEC dummy_f_coeff(1);
      bool contact_only = (mf_lambda->get_qdim() == 1);
      contact_rigid_obstacle_nonlinear_term
        nterm_pressure(CONTACT_PRESSURE, 0., mf_u, U, mf_obs, obs,
                       mf_lambda, lambda, 0, contact_only ? 0 : &dummy_f_coeff);

      getfem::generic_assembly assem_pressure;
      assem_pressure.set("V(#4)+=comp(NonLin(#1,#1,#2,#3).Base(#4))(i,:)");
      assem_pressure.push_mi(mim);
      assem_pressure.push_mf(mf_u);
      assem_pressure.push_mf(mf_obs);
      assem_pressure.push_mf(*mf_lambda);
      assem_pressure.push_mf(mf_ca);
      assem_pressure.push_nonlinear_term(&nterm_pressure);
      assem_pressure.push_vec(pressure);
      assem_pressure.assembly(rg);
      for (dal::bv_visitor_c dof(mf_ca.basic_dof_on_region(rg));
           !dof.finished(); ++dof)
        pressure[dof] /= mesh_size2[dof];

      // in areas where pressure is clearly non-zero set a low threshold
      // in order to avoid false negative contact detection
      scalar_type threshold_pressure(threshold_pressure_factor *
                                     gmm::vect_norminf(pressure));
      gmm::copy(gmm::scaled(mesh_size, scalar_type(-1)), threshold);
      for (getfem::mr_visitor v(rg); !v.finished(); ++v) {
        size_type nbdof = mf_ca.nb_basic_dof_of_face_of_element(v.cv(), v.f());
        mesh_fem::ind_dof_face_ct::const_iterator
          itdof = mf_ca.ind_basic_dof_of_face_of_element(v.cv(), v.f()).begin();
        bool all_positive = true;
        for (size_type k=0; k < nbdof; ++k, ++itdof)
          if (pressure[*itdof] < threshold_pressure) { all_positive = false; break; }
        if (!all_positive) {
          itdof = mf_ca.ind_basic_dof_of_face_of_element(v.cv(), v.f()).begin();
          for (size_type k=0; k < nbdof; ++k, ++itdof)
            threshold[*itdof] = threshold_factor * mesh_size[*itdof];
        }
      }
    }
    else
      gmm::copy(gmm::scaled(mesh_size, threshold_factor), threshold);

    // compute the total contact area
    // remark: the CONTACT_FLAG option misuses lambda as a threshold field
    contact_rigid_obstacle_nonlinear_term
      nterm(CONTACT_FLAG, 0., mf_u, U, mf_obs, obs, &mf_ca, &threshold);

    getfem::generic_assembly assem;
    assem.set("V()+=comp(NonLin(#1,#1,#2,#3))(i)");
    assem.push_mi(mim);
    assem.push_mf(mf_u);
    assem.push_mf(mf_obs);
    assem.push_mf(mf_ca);
    assem.push_nonlinear_term(&nterm);
    std::vector<scalar_type> v(1);
    assem.push_vec(v);
    assem.assembly(rg);
    return v[0];
  }


  template<typename VEC>
  void asm_nonmatching_meshes_normal_source_term
  (VEC &R, const mesh_im &mim,
   const getfem::mesh_fem &mf_u1,
   const getfem::mesh_fem &mf_u2_proj,
   const getfem::mesh_fem &mf_lambda, const VEC &lambda,
   const mesh_region &rg) {

    bool contact_only = (mf_lambda.get_qdim() == 1);

    VEC U1, U2_proj;
    gmm::resize(U1, mf_u1.nb_dof());
    gmm::resize(U2_proj, mf_u2_proj.nb_dof());
    scalar_type dummy_r(0);
    VEC dummy_f_coeff; gmm::resize(dummy_f_coeff,1);
    contact_nonmatching_meshes_nonlinear_term
      nterm(RHS_U_V1, dummy_r, mf_u1, U1, mf_u2_proj, U2_proj, &mf_lambda, &lambda,
            0, contact_only ? 0 : &dummy_f_coeff);

    getfem::generic_assembly assem;
    assem.set("V(#1)+=comp(NonLin(#1,#1,#2,#3).vBase(#1))(i,:,i)");
    assem.push_mi(mim);
    assem.push_mf(mf_u1);
    assem.push_mf(mf_u2_proj);
    assem.push_mf(mf_lambda);
    assem.push_nonlinear_term(&nterm);
    assem.push_vec(R);
    assem.assembly(rg);
  }

  template<typename VEC>
  scalar_type asm_nonmatching_meshes_contact_area
  (const mesh_im &mim,
   const getfem::mesh_fem &mf_u1, const VEC &U1,
   const getfem::mesh_fem &mf_u2_proj, const VEC &U2_proj,
   const mesh_region &rg, scalar_type threshold_factor=0.0,
   const getfem::mesh_fem *mf_lambda=0, const VEC *lambda=0,
   scalar_type threshold_pressure_factor=0.0) {

    if (!rg.get_parent_mesh())
      rg.from_mesh(mim.linked_mesh());
    getfem::mesh_fem mf_ca(mf_u1.linked_mesh());
    mf_ca.set_classical_finite_element(rg.index(),1);

    VEC mesh_size(mf_ca.nb_dof());
    VEC mesh_size2(mf_ca.nb_dof());
    { // assemble an estimator of the mesh size
      getfem::generic_assembly assem_mesh_size;
      assem_mesh_size.set("V(#1)+=comp(Base(#1))");
      assem_mesh_size.push_mi(mim);
      assem_mesh_size.push_mf(mf_ca);
      assem_mesh_size.push_vec(mesh_size2);
      assem_mesh_size.assembly(rg);
      for (dal::bv_visitor_c dof(mf_ca.basic_dof_on_region(rg));
           !dof.finished(); ++dof)
        mesh_size[dof] = sqrt(mesh_size2[dof]);
    }

    VEC threshold(mf_ca.nb_dof());
    if (mf_lambda && lambda) {
      VEC pressure(mf_ca.nb_dof());
      VEC dummy_f_coeff(1);
      bool contact_only = (mf_lambda->get_qdim() == 1);
      contact_nonmatching_meshes_nonlinear_term
        nterm_pressure(CONTACT_PRESSURE, 0., mf_u1, U1, mf_u2_proj, U2_proj,
                       mf_lambda, lambda, 0, contact_only ? 0 : &dummy_f_coeff);

      getfem::generic_assembly assem_pressure;
      assem_pressure.set("V(#4)+=comp(NonLin(#1,#1,#2,#3).Base(#4))(i,:)");
      assem_pressure.push_mi(mim);
      assem_pressure.push_mf(mf_u1);
      assem_pressure.push_mf(mf_u2_proj);
      assem_pressure.push_mf(*mf_lambda);
      assem_pressure.push_mf(mf_ca);
      assem_pressure.push_nonlinear_term(&nterm_pressure);
      assem_pressure.push_vec(pressure);
      assem_pressure.assembly(rg);
      for (dal::bv_visitor_c dof(mf_ca.basic_dof_on_region(rg));
           !dof.finished(); ++dof)
        pressure[dof] /= mesh_size2[dof];

      // in areas where pressure is clearly non-zero set a low threshold
      // in order to avoid false negative contact detection
      scalar_type threshold_pressure(threshold_pressure_factor *
                                     gmm::vect_norminf(pressure));
      gmm::copy(gmm::scaled(mesh_size, scalar_type(-1)), threshold);
      for (getfem::mr_visitor v(rg); !v.finished(); ++v) {
        size_type nbdof = mf_ca.nb_basic_dof_of_face_of_element(v.cv(), v.f());
        mesh_fem::ind_dof_face_ct::const_iterator
          itdof = mf_ca.ind_basic_dof_of_face_of_element(v.cv(), v.f()).begin();
        bool all_positive = true;
        for (size_type k=0; k < nbdof; ++k, ++itdof)
          if (pressure[*itdof] < threshold_pressure) { all_positive = false; break; }
        if (!all_positive) {
          itdof = mf_ca.ind_basic_dof_of_face_of_element(v.cv(), v.f()).begin();
          for (size_type k=0; k < nbdof; ++k, ++itdof)
            threshold[*itdof] = threshold_factor * mesh_size[*itdof];
        }
      }
    }
    else
      gmm::copy(gmm::scaled(mesh_size, threshold_factor), threshold);


    // compute the total contact area
    // remark: the CONTACT_FLAG option misuses lambda as a threshold field
    contact_nonmatching_meshes_nonlinear_term
      nterm(CONTACT_FLAG, 0., mf_u1, U1, mf_u2_proj, U2_proj, &mf_ca, &threshold);

    getfem::generic_assembly assem;
    assem.set("V()+=comp(NonLin(#1,#1,#2,#3))(i)");
    assem.push_mi(mim);
    assem.push_mf(mf_u1);
    assem.push_mf(mf_u2_proj);
    assem.push_mf(mf_ca);
    assem.push_nonlinear_term(&nterm);
    std::vector<scalar_type> v(1);
    assem.push_vec(v);
    assem.assembly(rg);
    return v[0];
  }


  void compute_integral_contact_area_and_force
  (model &md, size_type indbrick,
   scalar_type &area, model_real_plain_vector &Forces);



  /** Adds a contact condition with or without Coulomb friction on the variable
      `varname_u` and the mesh boundary `region`.  `Neumannterm`
      is the expression of the Neumann term (obtained by the Green formula)
      described as an expression of the high-level
      generic assembly language. This term can be obtained with
      md.Neumann_term(varname, region) once all volumic bricks have
      been added to the model. The contact condition
      is prescribed with Nitsche's method. The rigid obstacle should
      be described with the data `dataname_obstacle` being a signed distance to
      the obstacle (interpolated on a finite element method).
      `gamma0name` is the Nitsche's method parameter.
      `theta` is a scalar value which can be
      positive or negative. `theta = 1` corresponds to the standard symmetric
      method which is conditionnaly coercive for  `gamma0` small.
      `theta = -1` corresponds to the skew-symmetric method which is
      inconditionnaly coercive. `theta = 0` is the simplest method
      for which the second derivative of the Neumann term is not necessary.
      The optional parameter `dataexpr_friction_coeff` is the friction
      coefficient which could be any expression of the assembly language.
      Returns the brick index in the model.
  */
  size_type add_Nitsche_contact_with_rigid_obstacle_brick
  (model &md, const mesh_im &mim, const std::string &varname_u,
   const std::string &Neumannterm,
   const std::string &expr_obs, const std::string &dataname_gamma0,
   scalar_type theta_,
   std::string dataexpr_friction_coeff,
   const std::string &dataname_alpha,
   const std::string &dataname_wt,
   size_type region);

#ifdef EXPERIMENTAL_PURPOSE_ONLY

  /** Adds a contact condition with or without Coulomb friction on the variable
      `varname_u` and the mesh boundary `region`. The contact condition
      is prescribed with Nitsche's method. The rigid obstacle should
      be described with the data `dataname_obstacle` being a signed distance to
      the obstacle (interpolated on a finite element method).
      `gamma0name` is the Nitsche's method parameter.
      `theta` is a scalar value which can be
      positive or negative. `theta = 1` corresponds to the standard symmetric
      method which is conditionnaly coercive for  `gamma0` small.
      `theta = -1` corresponds to the skew-symmetric method which is
      inconditionnaly coercive. `theta = 0` is the simplest method
      for which the second derivative of the Neumann term is not necessary.
      The optional parameter `dataname_friction_coeff` is the friction
      coefficient which could be constant or defined on a finite element
      method.
      Returns the brick index in the model.
  */
  size_type add_Nitsche_contact_with_rigid_obstacle_brick_modified_midpoint
  (model &md, const mesh_im &mim, const std::string &varname_u,
   const std::string &Neumannterm, const std::string &Neumannterm_wt,
   const std::string &obs, const std::string &dataname_gamma0,
   scalar_type theta_,
   std::string dataname_friction_coeff,
   const std::string &dataname_alpha,
   const std::string &dataname_wt,
   size_type region);

#endif



  /** Adds a contact condition with or without Coulomb friction between
 two bodies in a fictitious domain. The contact condition is applied on
 the variable `varname_u1` corresponds with the first
 and slave body with Nitsche's method and on the variable `varname_u2`
 corresponds with the second and master body with Nitsche's method.
 The contact condition is evaluated on the fictitious slave bondary.
 The first body should be described by the level-set `dataname_d1`
 and the second body should be described by the level-set
`dataname_d2`. `gamma0name` is the Nitsche's method parameter.
 `theta` is a scalar value which can be positive or negative.
 `theta = 1` corresponds to the standard symmetric method which
 is conditionnaly coercive for  `gamma0` small.
 `theta = -1` corresponds to the skew-symmetric method which is
 inconditionnaly coercive. `theta = 0` is the simplest method for
 which the second derivative of the Neumann term is not necessary. 
The optional parameter `dataname_friction_coeff` is the friction 
coefficient which could be constant or defined on a finite element method. 
CAUTION: This brick has to be added in the model
 after all the bricks corresponding to partial differential
 terms having a Neumann term. Moreover, This brick can only
 be applied to bricks declaring their Neumann terms. Returns the brick index in the model. 

  */
  size_type add_Nitsche_fictitious_domain_contact_brick
  (model &md, const mesh_im &mim, const std::string &varname_u1,
   const std::string &varname_u2, const std::string &dataname_d1,
   const std::string &dataname_d2, const std::string &dataname_gamma0,
   scalar_type theta,
   const std::string &dataname_friction_coeff,
   const std::string &dataname_alpha,
   const std::string &dataname_wt1, const std::string &dataname_wt2);


}  /* end of namespace getfem.                                             */


#endif /* GETFEM_CONTACT_AND_FRICTION_INTEGRAL_H__ */
