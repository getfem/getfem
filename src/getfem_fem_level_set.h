// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Library : GEneric Tool for Finite Element Methods (getfem)
// File    : getfem_fem_level_set.h : definition of a finite element
//           method reprensenting a discontinous field across some level sets.
// Date    : March 09, 2005.
// Author  : Yves Renard <Yves.Renard@insa-toulouse.fr>
//           Julien Pommier <Julien.Pommier@insa-toulouse.fr>           
//
//========================================================================
//
// Copyright (C) 2004-2005 Yves Renard
//
// This file is a part of GETFEM++
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; version 2 of the License.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software Foundation,
// Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
//
//========================================================================
// To be corrected : dependencies. The mesh fem using this fem will not
//                   depend on the mesh fem arguments.



#ifndef GETFEM_FEM_LEVEL_SET_H__
#define GETFEM_FEM_LEVEL_SET_H__

#include <getfem_mesh_level_set.h>

namespace getfem {

  class fem_level_set : public virtual_fem, public context_dependencies {
    
  protected :

    const mesh_fem &mf;    // mf represents the original finite element method
                           // to be cutted.
    const mesh_level_set &mls;
  
    void update_from_context(void) const;

  public :

    void adapt(void);

    virtual size_type nb_dof(size_type cv) const;
    virtual size_type index_of_global_dof(size_type cv, size_type i) const;
    virtual bgeot::pconvex_ref ref_convex(size_type cv) const;
    virtual const bgeot::convex<base_node> &node_convex(size_type cv) const;
    virtual bgeot::pstored_point_tab node_tab(size_type) const;
    void base_value(const base_node &, base_tensor &) const;
    void grad_base_value(const base_node &, base_tensor &) const;
    void hess_base_value(const base_node &, base_tensor &) const;
    void real_base_value(const fem_interpolation_context& c, 
			 base_tensor &t) const;
    void real_grad_base_value(const fem_interpolation_context& c, 
			      base_tensor &t) const;
    void real_hess_base_value(const fem_interpolation_context&, 
			      base_tensor &) const;


    fem_level_set(const mesh_fem &mef_, const mesh_level_set &mls_);
  };
  
  pfem new_fem_level_set(const mesh_fem &mef, const mesh_level_set &mls);
  
  inline void del_fem_level_set(pfem pf) { dal::del_stored_object(pf); }

  
}  /* end of namespace getfem.                                            */

#endif
  
