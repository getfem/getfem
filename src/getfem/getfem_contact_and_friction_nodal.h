/* -*- c++ -*- (enables emacs c++ mode) */
/*===========================================================================
 
 Copyright (C) 2004-2012 Yves Renard, Konstantinos Poulios.
 
 This file is a part of GETFEM++
 
 Getfem++  is  free software;  you  can  redistribute  it  and/or modify it
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

/** @file getfem_contact_and_friction_nodal.h
    @author Yves Renard <Yves.Renard@insa-lyon.fr>
    @author Konstantinos Poulios <logari81@googlemail.com>
    @date July 6, 2004.
    @brief Unilateral contact and Coulomb friction condition brick.
 */
#ifndef GETFEM_CONTACT_AND_FRICTION_NODAL_H__
#define GETFEM_CONTACT_AND_FRICTION_NODAL_H__

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
      (see Getfem user documentation). The parameter `aug_version` indicates
      the augmentation strategy : 1 for the non-symmetric Alart-Curnier
      augmented Lagrangian, 2 for the symmetric one, 3 for the unsymmetric
      method with augmented multiplier.
  */
  size_type add_basic_contact_brick
  (model &md, const std::string &varname_u, const std::string &multname_n,
   const std::string &dataname_r, CONTACT_B_MATRIX &BN,
   std::string dataname_gap = "", std::string dataname_alpha = "",
   int aug_version=1, bool Hughes_stabilized=false);


  /** Add a contact with friction condition to the model. If U is the vector
      of degrees of freedom on which the condition is applied,
      the matrix `BN` has to be such that the contact condition is defined
      by $B_N U \le gap$ and `BT` have to be such that the relative tangential
      displacement is $B_T U$. The matrix `BT` should have as many rows as
      `BN` multiplied by $d-1$ where $d$ is the domain dimension.
      The contact condition is prescribed thank to a multiplier
      `multname_n` whose dimension should be equal to the number of rows of
      `BN` and the friction condition by a mutliplier `multname_t` whose size
      should be the number of rows of `BT`.
      The parameter `dataname_friction_coeff` describes the friction
      coefficient. It could be a scalar or a vector describing the
      coefficient on each contact condition.
      The augmentation parameter
      `r` should be chosen in a range of acceptabe values
      (see Getfem user documentation). `dataname_gap` is an
      optional parameter representing the initial gap. It can be a single value
      or a vector of value. `dataname_alpha` is an optional homogenization
      parameter for the augmentation parameter
      (see Getfem user documentation). The parameter `aug_version`
      indicates the augmentation strategy : 1 for the non-symmetric
      Alart-Curnier augmented Lagrangian, 2 for the symmetric one,
      3 for the unsymmetric
      method with augmented multiplier and 4 for the unsymmetric
      method with augmented multiplier and De Saxce projection.
  */
  size_type add_basic_contact_with_friction_brick
  (model &md, const std::string &varname_u, const std::string &multname_n,
   const std::string &multname_t, const std::string &dataname_r,
   CONTACT_B_MATRIX &BN, CONTACT_B_MATRIX &BT,
   std::string dataname_friction_coeff,
   std::string dataname_gap="", std::string dataname_alpha="",
   int aug_version=1, bool Tresca_version=false, bool Hughes_stabilized=false);

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
    parameter for the augmentation parameter. The parameter `aug_version`
    indicates the augmentation strategy : 1 for the non-symmetric
      Alart-Curnier augmented Lagrangian, 2 for the symmetric one,
      3 for the unsymmetric method with augmented multiplier.
  */
   inline size_type add_Hughes_stab_basic_contact_brick
   (model &md, const std::string &varname_u, const std::string &multname_n,
    const std::string &dataname_r, CONTACT_B_MATRIX &BN, CONTACT_B_MATRIX &DN,
    std::string dataname_gap="", std::string dataname_alpha="",
    int aug_version=1) {

    size_type indbrick = add_basic_contact_brick
      (md, varname_u, multname_n, dataname_r, BN,
       dataname_gap, dataname_alpha, aug_version, true);
    gmm::resize(contact_brick_set_DN(md, indbrick),
                gmm::mat_nrows(DN), gmm::mat_ncols(DN));
    gmm::copy(DN, contact_brick_set_DN(md, indbrick));
    return indbrick;
   }

  /**  Add Hughes stabilized friction contact condition to the model (broken ?). If U is the vector
      of degrees of freedom on which the condition is applied,
      the matrix `BN` have to be such that the contact condition is defined
      by $B_N U+DN Lambda \le 0$ (where 'DN' is the masse matrix
      relative to stabilzed term) and `BT` have to be such that the relative
      tangential displacement is $B_T U$. The matrix `BT` should have as many
      rows as `BN` multiplied b $d-1$ where $d$ is the domain dimension.
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
      (see Getfem user documentation). The parameter `aug_version`
      indicates the augmentation strategy : 1 for the non-symmetric
      Alart-Curnier augmented Lagrangian, 2 for the symmetric one,
      3 for the unsymmetric
      method with augmented multiplier and 4 for the unsymmetric
      method with augmented multiplier and De Saxce projection.
   **/
  inline size_type add_Hughes_stab_with_friction_contact_brick
  (model &md, const std::string &varname_u, const std::string &multname_n,
   const std::string &multname_t, const std::string &dataname_r,
   CONTACT_B_MATRIX &BN, CONTACT_B_MATRIX &BT, CONTACT_B_MATRIX &DN,CONTACT_B_MATRIX &DT,
   std::string dataname_friction_coeff,
   std::string dataname_gap="", std::string dataname_alpha="",
   int aug_version=1, bool Tresca_version=false) {

    size_type indbrick = add_basic_contact_with_friction_brick
      (md, varname_u, multname_n, multname_t, dataname_r, BN, BT,
       dataname_friction_coeff, dataname_gap, dataname_alpha,
       aug_version, Tresca_version, true);
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
      body, see Getfem user documentation). The parameter `aug_version`
      indicates the augmentation strategy : 1 for the non-symmetric
      Alart-Curnier augmented Lagrangian, 2 for the symmetric one,
      3 for the unsymmetric
      method with augmented multiplier.
      Basically, this brick computes the matrix BN
      and the vectors gap and alpha and calls the basic contact brick.
  */
  size_type add_contact_with_rigid_obstacle_brick
  (model &md, const mesh_im &mim, const std::string &varname_u,
   const std::string &multname_n, const std::string &dataname_r,
   size_type region, const std::string &obstacle, int aug_version=1);


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
      The parameter `aug_version`
      indicates the augmentation strategy : 1 for the non-symmetric
      Alart-Curnier augmented Lagrangian, 2 for the symmetric one,
      3 for the unsymmetric
      method with augmented multiplier and 4 for the unsymmetric
      method with augmented multiplier and De Saxce projection.
      Basically, this brick computes the matrix BN
      and the vectors gap and alpha and calls the basic contact brick.
  */
  size_type add_contact_with_friction_with_rigid_obstacle_brick
  (model &md, const mesh_im &mim, const std::string &varname_u,
   const std::string &multname_n, const std::string &multname_t,
   const std::string &dataname_r, const std::string &dataname_friction_coeff,
   size_type region, const std::string &obstacle, int aug_version=1);


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
      one of `slave1` and `slave2` is set to true. The parameter `aug_version`
      indicates the augmentation strategy : 1 for the non-symmetric
      Alart-Curnier augmented Lagrangian, 2 for the symmetric one,
      3 for the unsymmetric
      method with augmented multiplier.
      Basically, this brick computes the matrix BN and the vectors gap and
      alpha and calls the basic contact brick.
  */
  size_type add_nonmatching_meshes_contact_brick
  (model &md, const mesh_im &mim1, const mesh_im &mim2,
   const std::string &varname_u1, const std::string &varname_u2,
   std::string &multname_n, const std::string &dataname_r,
   const std::vector<size_type> &rg1, const std::vector<size_type> &rg2,
   bool slave1=true, bool slave2=false, int aug_version=1);

  inline size_type add_nonmatching_meshes_contact_brick
  (model &md, const mesh_im &mim1, const mesh_im &mim2,
   const std::string &varname_u1, const std::string &varname_u2,
   std::string &multname_n, const std::string &dataname_r,
   size_type rg1, size_type rg2, bool slave1=true, bool slave2=false,
   int aug_version=0) {

    std::vector<size_type> vrg1(1,rg1);
    std::vector<size_type> vrg2(1,rg2);
    return add_nonmatching_meshes_contact_brick
      (md, mim1, mim2, varname_u1, varname_u2, multname_n, dataname_r,
       vrg1, vrg2, slave1, slave2, aug_version);
  }

  inline size_type add_nonmatching_meshes_contact_brick
  (model &md, const mesh_im &mim, const std::string &varname_u,
   std::string &multname_n, const std::string &dataname_r,
   const std::vector<size_type> &rg1, const std::vector<size_type> &rg2,
   bool slave1=true, bool slave2=false, int aug_version=1) {

    return add_nonmatching_meshes_contact_brick
      (md, mim, mim, varname_u, varname_u, multname_n, dataname_r,
       rg1, rg2, slave1, slave2, aug_version);
  }

  inline size_type add_nonmatching_meshes_contact_brick
  (model &md, const mesh_im &mim, const std::string &varname_u,
   std::string &multname_n, const std::string &dataname_r,
   size_type rg1, size_type rg2, bool slave1=true, bool slave2=false,
   int aug_version=1) {

    std::vector<size_type> vrg1(1,rg1);
    std::vector<size_type> vrg2(1,rg2);
    return add_nonmatching_meshes_contact_brick
      (md, mim, mim, varname_u, varname_u, multname_n, dataname_r,
       vrg1, vrg2, slave1, slave2, aug_version);
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
      `slave2` is set to true.  The parameter `aug_version`
      indicates the augmentation strategy : 1 for the non-symmetric
      Alart-Curnier augmented Lagrangian, 2 for the symmetric one,
      3 for the unsymmetric
      method with augmented multiplier and 4 for the unsymmetric
      method with augmented multiplier and De Saxce projection.
      Basically, this brick computes the matrices BN and BT as well the vectors
      gap and alpha and calls the basic contact brick.
  */
  size_type add_nonmatching_meshes_contact_with_friction_brick
  (model &md, const mesh_im &mim1, const mesh_im &mim2,
   const std::string &varname_u1, const std::string &varname_u2,
   std::string &multname_n, std::string &multname_t,
   const std::string &dataname_r, const std::string &dataname_friction_coeff,
   const std::vector<size_type> &rg1, const std::vector<size_type> &rg2,
   bool slave1=true, bool slave2=false, int aug_version=1);

  inline size_type add_nonmatching_meshes_contact_with_friction_brick
  (model &md, const mesh_im &mim1, const mesh_im &mim2,
   const std::string &varname_u1, const std::string &varname_u2,
   std::string &multname_n, std::string &multname_t,
   const std::string &dataname_r, const std::string &dataname_friction_coeff,
   size_type rg1, size_type rg2, bool slave1=true, bool slave2=false,
   int aug_version=1) {

    std::vector<size_type> vrg1(1,rg1);
    std::vector<size_type> vrg2(1,rg2);
    return add_nonmatching_meshes_contact_with_friction_brick
      (md, mim1, mim2, varname_u1, varname_u2, multname_n, multname_t,
       dataname_r, dataname_friction_coeff,
       vrg1, vrg2, slave1, slave2, aug_version);
  }

  inline size_type add_nonmatching_meshes_contact_with_friction_brick
  (model &md, const mesh_im &mim, const std::string &varname_u,
   std::string &multname_n, std::string &multname_t,
   const std::string &dataname_r, const std::string &dataname_friction_coeff,
   const std::vector<size_type> &rg1, const std::vector<size_type> &rg2,
   bool slave1=true, bool slave2=false, int aug_version=1) {

    return add_nonmatching_meshes_contact_with_friction_brick
      (md, mim, mim, varname_u, varname_u, multname_n, multname_t,
       dataname_r, dataname_friction_coeff,
       rg1, rg2, slave1, slave2, aug_version);
  }

  inline size_type add_nonmatching_meshes_contact_with_friction_brick
  (model &md, const mesh_im &mim, const std::string &varname_u,
   std::string &multname_n, std::string &multname_t,
   const std::string &dataname_r, const std::string &dataname_friction_coeff,
   size_type rg1, size_type rg2, bool slave1=true, bool slave2=false,
   int aug_version=1) {

    std::vector<size_type> vrg1(1,rg1);
    std::vector<size_type> vrg2(1,rg2);
    return add_nonmatching_meshes_contact_with_friction_brick
      (md, mim, mim, varname_u, varname_u, multname_n, multname_t,
       dataname_r, dataname_friction_coeff,
       vrg1, vrg2, slave1, slave2, aug_version);
  }


}  /* end of namespace getfem.                                             */


#endif /* GETFEM_CONTACT_AND_FRICTION_NODAL_H__ */
