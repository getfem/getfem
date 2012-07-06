/* -*- c++ -*- (enables emacs c++ mode) */
/*===========================================================================
 
 Copyright (C) 2004-2012 Yves Renard
 
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
// To be corrected : dependencies. The mesh fem using this fem will not
//                   depend on the mesh fem arguments.

/**@file getfem_mesh_fem_product.h
   @author Yves Renard <Yves.Renard@insa-lyon.fr>,
   @author Julien Pommier <Julien.Pommier@insa-toulouse.fr>
   @date April 8, 2005.
   @brief A kind of product of two mesh_fems. Specific for Xfem enrichment.
*/

#ifndef GETFEM_MESH_FEM_PRODUCT_H__
#define GETFEM_MESH_FEM_PRODUCT_H__

#include "getfem_mesh_fem.h"

namespace getfem {
  
  class fem_product : public virtual_fem {
    pfem pfems[2];
    size_type cv, xfem_index;
    dal::bit_vector enriched_dof1;
    
  public:

    fem_product(pfem pf1_, pfem pf2_, size_type i, size_type xfi,
		const dal::bit_vector &nn)
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
    void update_from_context(void) const { is_adapted = false; touch(); }
    void clear(void);
    
    size_type memsize() const
    { return mesh_fem::memsize(); /* + ... */ }
    
    mesh_fem_product(const mesh_fem &me1, const mesh_fem &me2)
      : mesh_fem(me1.linked_mesh()), mf1(me1), mf2(me2)
    { is_adapted = false; xfem_index = reserve_xfem_index(); }
    void set_enrichment(const dal::bit_vector &nn)
    { enriched_dof = nn; adapt(); }
    
    ~mesh_fem_product() { clear_build_methods(); }
  };


}  /* end of namespace getfem.                                            */

#endif
  
