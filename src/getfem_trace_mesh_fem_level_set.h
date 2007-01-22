// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Copyright (C) 2007-2007 Yves Renard
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

/**@file getfem_trace_mesh_fem_level_set.h
   @author  Yves Renard <Yves.Renard@insa-toulouse.fr>,
   @date January 19, 2007.
   @brief a subclass of getfem::mesh_fem which allows to represent a subspace
   of a fem on a level-set. Different strategies should be provided.
*/

#ifndef GETFEM_TRACE_MESH_FEM_LEVEL_SET_H__
#define GETFEM_TRACE_MESH_FEM_LEVEL_SET_H__

#include <getfem_mesh_fem.h>
#include <getfem_mesh_im.h>
#include <getfem_mesh_level_set.h>

namespace getfem {

  /** @internal FEM used in  objects. */
  class sub_space_fem : public virtual_fem {
    pfem org_fem;
    size_type cv;
    std::vector<unsigned> ind;
    base_matrix B;
    
  public:

    sub_space_fem(pfem pf, const std::vector<unsigned> &indg,
		  const base_matrix &B_, size_type cv_)
      : org_fem(pf), cv(cv_), ind(indg), B(B_) { init(); }
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


  /**
     a subclass of mesh_fem which allows to represent a subspace of a trace
     of a fem on a level-set.
  */
  class trace_mesh_fem_level_set : public mesh_fem, public boost::noncopyable {
  protected :
    const mesh_level_set &mls;
    const mesh_fem &mf;
    unsigned degree;
    mutable std::vector<pfem> build_methods;
    mutable bool is_adapted;
    void clear_build_methods();

  public :
    void update_from_context(void) const { is_adapted = false; }

    /** build the mesh_fem keeping only the dof of the original
	mesh_fem which are listed in kept_dof. */
    void adapt(void);
    void clear(void);

    void receipt(const MESH_CLEAR &);
    void receipt(const MESH_DELETE &);
    
    size_type memsize() const {
      return mesh_fem::memsize(); // + ... ;
    }
    
    trace_mesh_fem_level_set(const mesh_level_set &me, const mesh_fem &mef,
			     unsigned degree_);

    ~trace_mesh_fem_level_set() { clear_build_methods(); }
  };
  
}  /* end of namespace getfem.                                            */

#endif
  
