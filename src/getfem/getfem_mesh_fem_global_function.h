/* -*- c++ -*- (enables emacs c++ mode) */
/*===========================================================================

 Copyright (C) 2004-2020 Yves Renard
 Copyright (C) 2016-2020 Konstantinos Poulios

 This file is a part of GetFEM

 GetFEM  is  free software;  you  can  redistribute  it  and/or modify it
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

/**@file getfem_mesh_fem_global_function.h
   @author  Yves Renard <Yves.Renard@insa-lyon.fr>, J. Pommier
   @date March, 2005.
   @brief Define a mesh_fem with base functions which are global functions
          given by the user.
*/
#ifndef GETFEM_MESH_FEM_GLOBAL_FUNCTION_H__
#define GETFEM_MESH_FEM_GLOBAL_FUNCTION_H__

#include "getfem_fem_global_function.h"

namespace getfem {

  /** this is a convenience class for defining a mesh_fem with base functions
      which are global functions (functions defined across more than one
      convexes of a mesh) given by the user.
  */
  class mesh_fem_global_function : public mesh_fem {
  protected :
    getfem::pfem fem_;
  public :

    void set_functions(const std::vector<pglobal_function>& f,
                       const mesh_im &mim=dummy_mesh_im());
    // size_type memsize() const;
    virtual void clear();

    mesh_fem_global_function(const mesh &me, dim_type q=1)
      : mesh_fem(me, q), fem_(0) {}
    virtual ~mesh_fem_global_function() { clear(); }
  };

  enum class bspline_boundary { FREE=0, PERIODIC=1, SYMMETRY=2};

  /** This function will generate bspline basis functions on NX uniform
      elements along a line. The dimensions of the domain correspond to
      the bounding interval of the 1d mesh linked by mf. The generated
      bspline basis functions are then set as the basis of mf.
      In case mim is provided, this integration method will be used
      to determine the support of he basis functions more precisely.
  */
  void define_uniform_bspline_basis_functions_for_mesh_fem
  (mesh_fem_global_function &mf, size_type NX, size_type order,
   bspline_boundary bcX_low=bspline_boundary::FREE,
   bspline_boundary bcX_high=bspline_boundary::FREE,
   const mesh_im &mim=dummy_mesh_im());

  inline void define_uniform_bspline_basis_functions_for_mesh_fem
  (mesh_fem_global_function &mf, size_type NX, size_type order,
   bspline_boundary bcX_low, const mesh_im &mim=dummy_mesh_im()) {
    define_uniform_bspline_basis_functions_for_mesh_fem
    (mf, NX, order, bcX_low, bcX_low, mim);
  }


  /** This function will generate bspline basis functions in an NX x NY
      rectilinear grid. The generated basis spans the entire bounding
      box of the 2d mesh linked by mf. The generated bspline basis
      functions are then set as the basis of mf.
      In case mim is provided, this integration method will be used to
      determine the support of he basis functions more precisely.
  */
  void define_uniform_bspline_basis_functions_for_mesh_fem
  (mesh_fem_global_function &mf,
   size_type NX, size_type NY, size_type order,
   bspline_boundary bcX_low=bspline_boundary::FREE,
   bspline_boundary bcY_low=bspline_boundary::FREE,
   bspline_boundary bcX_high=bspline_boundary::FREE,
   bspline_boundary bcY_high=bspline_boundary::FREE,
   const mesh_im &mim=dummy_mesh_im());

  inline void define_uniform_bspline_basis_functions_for_mesh_fem
  (mesh_fem_global_function &mf,
   size_type NX, size_type NY, size_type order,
   bspline_boundary bcX_low, bspline_boundary bcY_low,
   const mesh_im &mim=dummy_mesh_im()) {
    define_uniform_bspline_basis_functions_for_mesh_fem
    (mf, NX, NY, order, bcX_low, bcY_low, bcX_low, bcY_low, mim);
  }

}  /* end of namespace getfem.                                            */

#endif
