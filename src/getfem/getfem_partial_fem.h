// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Copyright (C) 2006-2007 Yves Renard
//
// This file is a part of GETFEM++
//
// Getfem++ is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
// You should have received a copy of the GNU Lesser General Public
// License along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301,
// USA.
//
//========================================================================
/**@file getfem_partial_fem.h
   @author  Yves Renard <Yves.Renard@insa-lyon.fr>, Julien Pommier <Julien.Pommier@insa-toulouse.fr>
   @date June 7, 2006.
   @brief Implement a special fem which disables some dof of another fem.
*/


#ifndef GETFEM_PARTIAL_FEM_H__
#define GETFEM_PARTIAL_FEM_H__

#include "getfem_fem.h"

namespace getfem {
  
  /** @internal FEM used in parital_mesh_fem objects. */
  class partial_fem : public virtual_fem {
    pfem org_fem;
    dal::bit_vector selected_dofs;
    std::vector<unsigned> ind;
    size_type cv;
    
  public:

    partial_fem(pfem pf,const dal::bit_vector &s, size_type cv_ = 0)
      : org_fem(pf), selected_dofs(s), cv(cv_) { init(); }
    size_type index_of_global_dof(size_type, size_type j) const;
    void init();
    void valid();
    void base_value(const base_node &x, base_tensor &t) const;
    void grad_base_value(const base_node &x, base_tensor &t) const;
    void hess_base_value(const base_node &x, base_tensor &t) const;

    void real_base_value(const fem_interpolation_context& c, 
			 base_tensor &t, bool = true) const;    
    void real_grad_base_value(const fem_interpolation_context& c, 
			      base_tensor &t, bool = true) const;
    void real_hess_base_value(const fem_interpolation_context& c, 
			      base_tensor &t, bool = true) const;
  };



}  /* end of namespace getfem.                                            */

#endif
  
