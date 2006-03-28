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
// To be corrected : dependencies. The mesh fem using this fem will not
//                   depend on the mesh fem arguments.

/**@file getfem_mesh_fem_product.h
   @author  Yves Renard <Yves.Renard@insa-toulouse.fr>, Julien Pommier <Julien.Pommier@insa-toulouse.fr>
   @date April 8, 2005.
   @brief NEED A DESCRIPTION!
*/

#ifndef GETFEM_MESH_FEM_PRODUCT_H__
#define GETFEM_MESH_FEM_PRODUCT_H__

#include <getfem_mesh_fem.h>

namespace getfem {
  
  class fem_product : public virtual_fem {
    pfem pfems[2];
    size_type cv, xfem_index;
    dal::bit_vector enriched_dof1;
    
  public:

    fem_product(pfem pf1_, pfem pf2_, size_type i, size_type xfi, const dal::bit_vector &nn)
      : cv(i), xfem_index(xfi),  enriched_dof1(nn)
    { pfems[0] = pf1_; pfems[1] = pf2_; init(); }
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


  class mesh_fem_product : public mesh_fem {
  protected :
    const mesh_fem &mf1, &mf2;

    mutable std::vector<pfem> build_methods;
    mutable bool is_adapted;
    size_type xfem_index;
    dal::bit_vector enriched_dof;
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
    
    mesh_fem_product(const mesh_fem &me1, const mesh_fem &me2)
      : mesh_fem(me1.linked_mesh()), mf1(me1), mf2(me2)
    { is_adapted = false; xfem_index = reserve_xfem_index(); }
    void set_enrichment(const dal::bit_vector &nn)
    { enriched_dof = nn; adapt(); }
    
    ~mesh_fem_product() { clear_build_methods(); }
  };


}  /* end of namespace getfem.                                            */

#endif
  
