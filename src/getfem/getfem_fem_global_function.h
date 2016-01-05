/* -*- c++ -*- (enables emacs c++ mode) */
/*===========================================================================

 Copyright (C) 2004-2016 Yves Renard
 Copyright (C) 2016      Konstantinos Poulios

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

/**@file getfem_fem_global_function.h
   @author  Yves Renard <Yves.Renard@insa-lyon.fr>, J. Pommier
   @date March, 2005.
   @brief Define mesh_fem whose base functions are global function given by the user.
*/
#ifndef GETFEM_FEM_GLOBAL_FUNCTION_H__
#define GETFEM_FEM_GLOBAL_FUNCTION_H__

#include "getfem_fem.h"
#include "getfem_mesh_fem.h"
#include "getfem_global_function.h"

namespace getfem {

  /// fem object with global basis functions.
  class fem_global_function : public virtual_fem, public context_dependencies {
  protected :
    std::vector<pglobal_function> functions;
    const mesh &m;
    const mesh_im &mim;
    const bool has_mesh_im;

    mutable bgeot::multi_index mib,mig,mih;
    mutable std::vector<std::vector<size_type> > index_of_global_dof_;
    mutable bgeot::pstored_point_tab pspt_override;

    void init();
    virtual void update_from_context() const;
  public :
    virtual size_type nb_dof(size_type cv) const;
    virtual size_type index_of_global_dof(size_type cv, size_type i) const;
    virtual bgeot::pconvex_ref ref_convex(size_type cv) const;
    virtual const bgeot::convex<base_node> &node_convex(size_type cv) const;
    virtual bgeot::pstored_point_tab node_tab(size_type) const
    { return pspt_override; }

    void base_value(const base_node &, base_tensor &) const;
    void grad_base_value(const base_node &, base_tensor &) const;
    void hess_base_value(const base_node &, base_tensor &) const;

    void real_base_value(const fem_interpolation_context& c,
                         base_tensor &t, bool = true) const;
    void real_grad_base_value(const fem_interpolation_context& c,
                              base_tensor &t, bool = true) const;
    void real_hess_base_value(const fem_interpolation_context&,
                              base_tensor &, bool = true) const;

    fem_global_function(const std::vector<pglobal_function> &funcs,
                        const mesh &m_);
    fem_global_function(const std::vector<pglobal_function> &funcs,
                        const mesh_im &mim_);
    virtual ~fem_global_function()
    { DAL_STORED_OBJECT_DEBUG_DESTROYED(this, "Global function fem"); }
  };

  /** create a new global function FEM. 
      @param funcs is a vector containing all global basis functions.
      @param m is the mesh to be used for numerical integration in the assembly.
  */
  pfem new_fem_global_function(const std::vector<pglobal_function> &funcs,
                               const mesh &m);

  /** create a new global function FEM. 
      @param funcs is a vector containing all global basis functions.
      @param mim is the integration method to be used in the assembly.
  */
  pfem new_fem_global_function(const std::vector<pglobal_function> &funcs,
                               const mesh_im &mim);

  /** release a global function FEM */
  inline void del_fem_global_function(const pfem &pf)
  { dal::del_stored_object(pf); }


}  /* end of namespace getfem.                                            */

#endif
