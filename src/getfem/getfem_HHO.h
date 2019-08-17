// /* -*- c++ -*- (enables emacs c++ mode) */
/*===========================================================================

 Copyright (C) 2019-2019 Yves Renard

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
/**@file getfem_HHO.h
   @author  Yves Renard <Yves.Renard@insa-lyon.fr>
   @date August 16, 2019.
   @brief Tools for Hybrid-High-Order methods.
*/

#ifndef GETFEM_HHO_H__
#define GETFEM_HHO_H__

#include "getfem_models.h"


namespace getfem {

  /** Add the elementary transformation corresponding to the reconstruction
      of a gradient for HHO methods to the model.
      The name is the name given to the elementary transformation.
  */
  void add_HHO_reconstructed_gradient(model &md, std::string name);

  
  /** Add the elementary transformation corresponding to the reconstruction
      of a symmetrized gradient for HHO methods to the model.
      The name is the name given to the elementary transformation.
  */
  void add_HHO_reconstructed_symmetrized_gradient(model &md, std::string name);

  /** Add the elementary transformation to the model corresponding to the
      reconstruction of the variable.
      The name is the name given to the elementary transformation.
  */
  void add_HHO_reconstructed_value(model &md, std::string name);

  /** Add the elementary transformation to the model corresponding to the
      reconstruction of the variable using a symmetrized gradient.
      The name is the name given to the elementary transformation.
  */
  void add_HHO_reconstructed_symmetrized_value(model &md, std::string name);

  /** Add the elementary transformation to the model corresponding to the
      HHO stabilization operator.
      The name is the name given to the elementary transformation.
  */
  void add_HHO_stabilization(model &md, std::string name);

  /** Add the elementary transformation to the model corresponding to the
      HHO stabilization operator using a symmetrized gradient.
      The name is the name given to the elementary transformation.
  */
  void add_HHO_symmetrized_stabilization(model &md, std::string name);


}  /* end of namespace getfem.                                             */


#endif /* GETFEM_HHO_H__ */
