/* -*- c++ -*- (enables emacs c++ mode) */
/*===========================================================================
 
 Copyright (C) 2004-2015 Yves Renard, Jeremie Lasry
 
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
/**@file getfem_linearized_plates.h
   @author  Yves Renard <Yves.Renard@insa-lyon.fr>
   @date November 1, 2004.
   @brief Reissner-Mindlin plate model brick.
*/

#ifndef GETFEM_LINEARIZED_PLATES_H__
#define GETFEM_LINEARIZED_PLATES_H__

#include "getfem_models.h"


namespace getfem {

  /** Add the elementary transformation corresponding to the projection
      on rotated RT0 element for two-dimensional elements to the model.
      The name is the name given to the elementary transformation.
  */
  void add_2D_rotated_RT0_projection(model &md, std::string name);


  /** Add a term corresponding to the classical Reissner-Mindlin plate
      model for which `u3` is the transverse displacement,
      `Theta` the rotation of
      fibers normal to the midplane, 'param_E' the Young Modulus,
      `param_nu` the poisson ratio,
      `param_epsilon` the plate thickness,
      `param_kappa` the shear correction factor. Note that since this brick
      uses the high level generic assembly language, the parameter can
      be regular expression of this language.
      There are three variants.
      `variant = 0` corresponds to the an
      unreduced formulation and in that case only the integration
      method `mim` is used. Practically this variant is not usable since
      it is subject to a strong locking phenomenon.
      `variant = 1` corresponds to a reduced integration where `mim` is
      used for the rotation term and `mim_reduced` for the transverse
      shear term. `variant = 2` (default) corresponds to the projection onto
      a rotated RT0 element of the transverse shear term. For the moment, this
      is adapted to quadrilateral only (because it is not sufficient to
      remove the locking phenomenon on triangle elements). Note also that if
      you use high order elements, the projection on RT0 will reduce the order
      of the approximation.
      Returns the brick index in the model.
   */
  size_type add_Mindlin_Reissner_plate_brick
  (model &md, const mesh_im &mim, const mesh_im &mim_reduced,
   const std::string &u3,
   const std::string &Theta, const std::string &param_E,
   const std::string &param_nu, const std::string &param_epsilon,
   const std::string &param_kappa, size_type variant = size_type(2), 
   size_type region = size_type(-1));

}  /* end of namespace getfem.                                             */


#endif /* GETFEM_LINEARIZED_PLATES_H__ */
