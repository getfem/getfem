/* -*- c++ -*- (enables emacs c++ mode) */
/*===========================================================================

 Copyright (C) 2013-2013 Yves Renard, Konstantinos Poulios.

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

/** @file getfem_contact_and_friction_integral.h
    @author Yves Renard <Yves.Renard@insa-lyon.fr>
    @author Konstantinos Poulios <logari81@googlemail.com>
    @date May, 2013.
    @brief Large sliding unilateral contact and friction condition brick.
 */
#ifndef GETFEM_CONTACT_AND_FRICTION_LARGE_SLIDING_H__
#define GETFEM_CONTACT_AND_FRICTION_LARGE_SLIDING_H__

#include "getfem_contact_and_friction_common.h"

namespace getfem {


  /** Adds a large sliding contact with friction brick to the model.
      This brick is able to deal with self-contact, contact between
      several deformable bodies and contact with rigid obstacles.
      It takes a variable of type multi_contact_frame wich describe
      the contact situation (master and slave contact boundaries,
      self-contact detection or not, and a few parameter).
      For each slave boundary (and also master boundaries if self-contact
      is asked) a multiplier variable should be defined.
  */
  size_type add_integral_large_sliding_contact_brick_raytrace
  (model &md, multi_contact_frame &mcf,
   const std::string &dataname_r,
   const std::string &dataname_friction_coeff = std::string(),
   const std::string &dataname_alpha = std::string());




  // Old brick, to be adapted ...



  /** Adds a large sliding contact with friction brick to the model.
      This brick is able to deal with auto-contact, contact between
      several deformable bodies and contact with rigid obstacles.
      The condition is applied on the variable `varname_u` on the
      boundary corresponding to `region`. `dataname_r` is the augmentation
      parameter of the augmented Lagrangian. `dataname_friction_coeff`
      is the friction coefficient. `mim` is an integration method on the
      boundary. `varname_u` is the variable on which the contact condition
      will be prescribed (should be of displacement type). `multname` is
      a multiplier defined on the boundary which will represent the contact
      force. If no additional boundary or rigid
      obstacle is added, only auto-contact will be detected. Use
      `add_boundary_to_large_sliding_contact_brick` and
      `add_rigid_obstacle_to_large_sliding_contact_brick` to add contact
      boundaries and rigid obstacles.
  */
  size_type add_integral_large_sliding_contact_brick_field_extension
  (model &md, const mesh_im &mim, const std::string &varname_u,
   const std::string &multname, const std::string &dataname_r,
   const std::string &dataname_friction_coeff, size_type region);


  /** Adds a contact boundary to an existing large sliding contact brick.
      `indbrick` is the brick index.
  */
  void add_boundary_to_large_sliding_contact_brick
  (model &md, size_type indbrick, const mesh_im &mim,
   const std::string &varname_u, const std::string &multname,
   size_type region);

  /** Adds a rigid obstacle to an existing large sliding contact brick.
      `indbrick` is the brick index, `obs` is the expression of a
      function which should be closed to a signed distance to the obstacle.
  */
  void add_rigid_obstacle_to_large_sliding_contact_brick
  (model &md, size_type indbrick, const std::string &obs);


}  /* end of namespace getfem.                                             */


#endif /* GETFEM_CONTACT_AND_FRICTION_LARGE_SLIDING_H__ */
