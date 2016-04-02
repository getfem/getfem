/* -*- c++ -*- (enables emacs c++ mode) */
/*===========================================================================

 Copyright (C) 2002-2015 Amandine Cottaz, Yves Renard

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

/**@file getfem_plasticity.h
   @author  Yves Renard <Yves.Renard@insa-lyon.fr>,
   @author  Julien Pommier <Julien.Pommier@insa-toulouse.fr>
   @author  Amandine Cottaz
   @date June 2, 2010.
   @brief Plasticty bricks.
*/
#ifndef GETFEM_PLASTICITY_H__
#define GETFEM_PLASTICITY_H__

#include "getfem_models.h"
#include "getfem_assembling_tensors.h"
#include "getfem_derivatives.h"
#include "getfem_interpolation.h"
#include "gmm/gmm_dense_qr.h"

namespace getfem {


  //=================================================================
  // Abstract contraints projection
  //=================================================================


  /** Abstract projection of a stress tensor onto a set of admissible
      stress tensors.
  */
  class abstract_constraints_projection {
  protected :
    size_type flag_hyp;

  public :
      /* if flag_proj=0 the output will be Proj(tau)
       * if flag_proj=1 the output will be gradProj(tau)
       * no others values allowed for flag_proj
       */
    virtual void do_projection(const base_matrix& tau,
                               scalar_type stress_threshold,
                               base_matrix& proj,
                               size_type flag_proj) const = 0;
    abstract_constraints_projection (size_type flag_hyp_ = 0) :
      flag_hyp(flag_hyp_) {}
    virtual ~abstract_constraints_projection () {}
  };

  typedef std::shared_ptr<const abstract_constraints_projection>
  pconstraints_projection;

  //=================================================================
  // Von Mises projection
  //=================================================================


  /** Von Mises projection */
  class VM_projection : public abstract_constraints_projection {

    /* used to compute the projection */
    template<typename MAT>
    void tau_m_Id(const MAT& tau, MAT &taumid) const {
      scalar_type trace = gmm::mat_trace(tau);
      size_type size_of_tau = gmm::mat_nrows(tau);
      gmm::copy(gmm::identity_matrix(),taumid);
      gmm::scale(taumid, trace / scalar_type(size_of_tau));
    }

    /* used to compute the projection */
    template<typename MAT>
    void tau_d(const MAT& tau, MAT &taud) const {
      tau_m_Id(tau, taud);
      gmm::scale(taud, scalar_type(-1));
      gmm::add(tau, taud);
    }


    public :

    /** the Von Mises projection computation */
    /* on input : tau matrix, on output : the projection of tau */
    virtual void do_projection(const base_matrix& tau,
                               scalar_type stress_threshold,
                               base_matrix& proj,
                               size_type flag_proj)  const {

      /* be sure that flag_proj has a correct value */
      GMM_ASSERT1(flag_proj == 0 || flag_proj ==1,
                  "wrong value for the projection flag, "
                  "must be 0 or 1 ");

      /* be sure that stress_threshold has a correct value */
      GMM_ASSERT1(stress_threshold>=0., "s is not a positive number "
                  << stress_threshold << ". You need to set "
                  << "s as a positive number");

      size_type N = gmm::mat_nrows(tau);
      size_type projsize = (flag_proj == 0) ? N : gmm::sqr(N);
      scalar_type normtaud;

      /* calculate tau_m*Id */
      base_matrix taumId(N, N);
      tau_m_Id(tau, taumId);

      // calcul du deviateur de tau, taud
      base_matrix taud(N,N);
      gmm::add(gmm::scaled(taumId, scalar_type(-1)), tau, taud);

      /* plane constraints */
      if (flag_hyp == 1) {  // To be done ...
        N /= 2;
        GMM_ASSERT1(!N, "wrong value for CALCULATION HYPOTHESIS, "
                    "must be /=1 SINCE n/=2");
        // we form the 3D tau tensor considering
        // that tau(3,j)=tau(i,3)=0
        base_matrix tau_aux(3,3); gmm::clear(tau_aux);
        gmm::copy(tau,gmm::sub_matrix
                  (tau_aux,gmm::sub_interval(0,2)));
        // we calculate tau deviator and its norms
        base_matrix taud_aux(3,3);
        tau_d(tau_aux, taud_aux);
        normtaud=gmm::mat_euclidean_norm(taud_aux);
      }
      else normtaud=gmm::mat_euclidean_norm(taud);


      /* dimension and initialization of proj matrix or
         its derivative */
      gmm::resize(proj, projsize, projsize);

      if (normtaud <= stress_threshold) {
        switch(flag_proj) {
        case 0: gmm::copy(tau, proj); break;
        case 1: gmm::copy(gmm::identity_matrix(), proj); break;
        }
      }
      else {
        switch(flag_proj) {
        case 0:
          gmm::copy(gmm::scaled(taud, stress_threshold/normtaud),
                    proj);
          gmm::add(taumId,proj);
          break;
        case 1:
          base_matrix Igrad(projsize, projsize);
          gmm::copy(gmm::identity_matrix(),Igrad);
          base_matrix Igrad2(projsize, projsize);

          // build vector[1 0 0 1  0 0 1...] to be copied in certain
          // columns of Igrad(*)Igrad
          base_vector aux(projsize);
          for (size_type i=0; i < N; ++i)
            aux[i*N + i] = scalar_type(1);

          // Copy in a selection of columns of Igrad(*)Igrad
          for (size_type i=0; i < N; ++i)
            gmm::copy(aux, gmm::mat_col(Igrad2, i*N + i));

          // Compute Id_grad
          base_matrix Id_grad(projsize, projsize);
          scalar_type rr = scalar_type(1)/scalar_type(N);
          gmm::copy(gmm::scaled(Igrad2, -rr), Id_grad);
          gmm::add(Igrad, Id_grad);


          // Compute ngrad(*)ngrad
          base_matrix ngrad2(projsize, projsize);
          // Compute the normal n
          base_matrix un(N, N);
          gmm::copy(gmm::scaled(taud, 1./normtaud),un);

          // Copy of the normal in a column vector
          // in the Fortran order
          std::copy(un.begin(), un.end(), aux.begin());

          // Loop on the columns of ngrad(*)ngrad
          for (size_type j=0; j < projsize; ++j)
            gmm::copy(gmm::scaled(aux,aux[j]),
                      gmm::mat_col(ngrad2,j));


          // Final computation of the projection gradient
          gmm::copy(gmm::identity_matrix(), proj);
          gmm::add(gmm::scaled(ngrad2, scalar_type(-1)), proj);
          base_matrix aux2(projsize, projsize);
          gmm::copy(gmm::scaled(proj, stress_threshold/normtaud),
                    aux2);
          gmm::mult(aux2,Id_grad,proj);
          gmm::add(gmm::scaled(Igrad2, rr),proj);
          break;
        }
      }
    }


    VM_projection(size_type flag_hyp_ = 0) :
      abstract_constraints_projection (flag_hyp_) {}
  };




  //=================================================================
  //
  //  Elastoplasticity Brick
  //
  //=================================================================


  /** Add a nonlinear elastoplasticity term to the model for small
      deformations and an isotropic material, with respect
      to the variable `varname`.
      Note that the constitutive lawtype of projection
      to be used is described by `ACP` which should not be
      freed while the model is used.
      `datalambda` and `datamu` describe the Lamé coeffcients
      of the studied material. Could be scalar or vector fields
      described on a finite element method.
      `datathreshold` represents the elasticity threshold
      of the material. It could be a scalar or a vector field
      described on the same finite element method as
      the Lamé coefficients.
      `datasigma` represents the stress constraints values
      supported by the material. It should be a vector field
      described on a finite element method.
      `previous_dep_name` represents the displacement at the previous time step.
      Moreover, if `varname` is described
      onto a K-th order mesh_fem, `datasigma` has to be described
      on a mesh_fem of order at least K-1.
  */
  size_type add_elastoplasticity_brick(model &md,
                                       const mesh_im &mim,
                                       const pconstraints_projection &ACP,
                                       const std::string &varname,
                                       const std::string &previous_dep_name,
                                       const std::string &datalambda,
                                       const std::string &datamu,
                                       const std::string &datathreshold,
                                       const std::string &datasigma,
                                       size_type region = size_type(-1));

  /** This function permits to compute the new stress constraints
      values supported by the material after a load or an unload.
      `varname` is the main unknown of the problem
      (the displacement),
      `previous_dep_name` represents the displacement at the previous time step,
      `ACP` is the type of projection to be used that could only be
      `Von Mises` for the moment,
      `datalambda` and `datamu` are the Lamé coefficients
      of the material,
      `datathreshold` is the elasticity threshold of the material,
      `datasigma` is the vector which will contains the new
      computed values. */
  void elastoplasticity_next_iter(model &md,
                                  const mesh_im &mim,
                                  const std::string &varname,
                                  const std::string &previous_dep_name,
                                  const pconstraints_projection &ACP,
                                  const std::string &datalambda,
                                  const std::string &datamu,
                                  const std::string &datathreshold,
                                  const std::string &datasigma);

  /** This function computes on mf_vm the Von Mises or Tresca stress
      of a field for elastoplasticity and return it into the vector VM.
      Note that `datasigma` should be the vector containing the new
      stress constraints values, i.e. after a load or an unload
      of the material. If `tresca` = 'true', the Tresca stress will
      be computed, otherwise it will be the Von Mises one.*/
  void compute_elastoplasticity_Von_Mises_or_Tresca
  (model &md,
   const std::string &datasigma,
   const mesh_fem &mf_vm,
   model_real_plain_vector &VM,
   bool tresca);

  /** This function computes on mf_pl the plastic part, that could appear
      after a load and an unload, into the vector `plast`.
      Note that `datasigma` should be the vector containing the new
      stress constraints values, i.e. after a load or an unload
      of the material. */
  void compute_plastic_part(model &md,
                            const mesh_im &mim,
                            const mesh_fem &mf_pl,
                            const std::string &varname,
                            const std::string &previous_dep_name,
                            const pconstraints_projection &ACP,
                            const std::string &datalambda,
                            const std::string &datamu,
                            const std::string &datathreshold,
                            const std::string &datasigma,
                            model_real_plain_vector &plast);


  // Finite strain elastoplasticity

  /** Add a linear function with the name specified by `name` to represent
      linear isotropoc hardening in plasticity with initial yield limit
      `sigma_y0` and hardening modulus `H`.
      A true value of the `frobenius` argument will express the hardening
      function in terms of Frobenius norms both for the strain input and
      the stress output, instead of the corresponding Von-Mises quantities.
  */
  void ga_define_linear_hardening_function
  (const std::string &name, scalar_type sigma_y0, scalar_type H, bool frobenius=true);

  /** Add a Ramberg-Osgood hardening function with the name specified by
     `name`, for reference stress and strain given by `sigma_ref` and
      `eps_ref` respectively and for a hardening exponent `n`.
      A true value of the `frobenius` argument will express the hardening
      function in terms of Frobenius norms both for the strain input and
      the stress output, instead of the corresponding Von-Mises quantities.
  */
  void ga_define_Ramberg_Osgood_hardening_function
  (const std::string &name,
   scalar_type sigma_ref, scalar_type eps_ref, scalar_type n,
   bool frobenius=false);

  /** Add a Ramberg-Osgood hardening function with the name specified by
     `name`, for reference stress `sigma_ref`, Young's modulus `E`,
      offset parameter `alpha` and hardening parameter `n`.
      A true value of the `frobenius` argument will express the hardening
      function in terms of Frobenius norms both for the strain input and
      the stress output, instead of the corresponding Von-Mises quantities.
  */
  inline void ga_define_Ramberg_Osgood_hardening_function
  (const std::string &name, scalar_type sigma_ref, scalar_type E, 
   scalar_type alpha, scalar_type n, bool frobenius=false) {
    ga_define_Ramberg_Osgood_hardening_function
      (name, sigma_ref, alpha*sigma_ref/E, n, frobenius);
  }

  /** Add a finite strain elastoplasticity brick to the model
      with respect to the displacement variable `dispname` and the
      plastic multiplier `multname`. If `pressname` is not empty,
      a mixed displacement-pressure formulation is defined.
      For the moment there is only one supported law defined through 
      `lawname` as "Simo_Miehe" and expects as `params` a vector of
      the following five parameters:
      1) an expression for the initial bulk modulus K,
      2) an expression for the initial shear modulus G,
      3) the name of a user predefined function that decribes
         the yield limit as a function of the hardening variable
         (both the yield limit and the hardening variable values are
          assumed to be Frobenius norms of appropriate stress and strain
          tensors, respectively),
      4) the name of a (scalar) fem_data or im_data field that holds the
         plastic strain at the previous time step, and
      5) the name of a fem_data or im_data field that holds all
         non-repeated components of the inverse of the plastic right
         Cauchy-Green tensor at the previous time step
         (it has to be a 4 element vector for plane strain 2D problems
         and a 6 element vector for 3D problems).
  */
  size_type add_finite_strain_elastoplasticity_brick
  (model &md, const mesh_im &mim, const std::string &dispname,
   const std::string &multname, const std::string &pressname,
   const std::string &lawname, const std::vector<std::string> &params,
   size_type region = size_type(-1));

  /** This function permits to compute the new stress constraints
      values supported by the material after a load or an unload.
      `varname` is the main unknown of the problem
      (the displacement),
      `previous_dep_name` represents the displacement at the previous time step,
      `ACP` is the type of projection to be used that could only be
      `Von Mises` for the moment,
      `datalambda` and `datamu` are the Lamé coefficients
      of the material,
      `datathreshold` is the elasticity threshold of the material,
      `datasigma` is the vector which will contains the new
      computed values. */
  void finite_strain_elastoplasticity_next_iter
  (model &md, const mesh_im &mim, const std::string &dispname,
   const std::string &multname, const std::string &pressname,
   const std::string &lawname, const std::vector<std::string> &params,
   size_type region = size_type(-1));

  /** This function computes on mf_vm the Von Mises stress with respect to
      a finite strain elastoplasticity term.
      If `assemble` = 'true', the Von-Mises stress field is assembled on
      mf_vm, otherwise for the default value `assemble` = 'false', the field
      is simply interpolated.*/
  void compute_finite_strain_elastoplasticity_Von_Mises
  (model &md, const mesh_im &mim, const std::string &dispname,
   const std::string &multname, const std::string &pressname,
   const std::string &lawname, const std::vector<std::string> &params,
   const mesh_fem &mf_vm, model_real_plain_vector &VM, bool assemble=false,
   size_type region = size_type(-1));


} /* namespace getfem */

#endif
