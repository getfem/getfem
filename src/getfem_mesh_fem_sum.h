// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Copyright (C) 2004-2006 Yves Renard
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
/**@file getfem_mesh_fem_sum.h
   @author  Yves Renard <Yves.Renard@insa-toulouse.fr>, Julien Pommier <Julien.Pommier@insa-toulouse.fr>
   @date March 18, 2005.
   @brief Implement a special mesh_fem with merges the FEMs of two (or more) mesh_fems.
*/


#ifndef GETFEM_MESH_SUM_H__
#define GETFEM_MESH_SUM_H__

#include <getfem_mesh_fem.h>

namespace getfem {
  typedef std::vector<const std::string *> dof_enrichments;
  
  /** @internal FEM used in mesh_fem_sum objects. */
  class fem_sum : public virtual_fem {
    std::vector<pfem> pfems; /* the fems to be summed */
    size_type cv;
    
  public:

    size_type index_of_global_dof(size_type cv_, size_type j) const;
    fem_sum(const std::vector<pfem> &pfs, size_type i) : pfems(pfs), cv(i) { init(); }
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


  /** @brief Implement a special mesh_fem with merges the FEMs of two (or more) mesh_fems.*/
  class mesh_fem_sum : public mesh_fem {
  protected :
    std::vector<const mesh_fem *> mfs;

    mutable std::map< std::vector<pfem>, pfem> situations;
    mutable std::vector<pfem> build_methods;
    mutable bool is_adapted;
    void clear_build_methods();

  public :
    void adapt(void);
    void update_from_context(void) const { is_adapted = false; }
    void clear(void);
    
    void receipt(const MESH_CLEAR &);
    void receipt(const MESH_DELETE &);
    // void receipt(const MESH_ADD_CONVEX &m)  // to be modified ?
    // void receipt(const MESH_SUP_CONVEX &m)  // to be modified ?
    // void receipt(const MESH_SWAP_CONVEX &m) // to be modified ?
    
    size_type memsize() const {
      return mesh_fem::memsize(); // + ... ;
    }
    
    mesh_fem_sum(mesh &me) : mesh_fem(me) { is_adapted = false; }
    void set_mesh_fems(const std::vector<const mesh_fem *> &mefs)
    { mfs = mefs; adapt(); }
    void set_mesh_fems(const mesh_fem &mf1)
    { mfs.clear(); mfs.push_back(&mf1); adapt(); }
    void set_mesh_fems(const mesh_fem &mf1, const mesh_fem &mf2)
    { mfs.clear(); mfs.push_back(&mf1); mfs.push_back(&mf2);  adapt(); }
    void set_mesh_fems(const mesh_fem &mf1, const mesh_fem &mf2, const mesh_fem &mf3)
    { mfs.clear(); mfs.push_back(&mf1); mfs.push_back(&mf2); mfs.push_back(&mf3);  adapt(); }


    ~mesh_fem_sum() { clear_build_methods(); }
  };


}  /* end of namespace getfem.                                            */

#endif
  
