// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Library : GEneric Tool for Finite Element Methods (getfem)
// File    : getfem_mesh_fem_sum.h : definition of a finite element
//           method reprensenting a direct sum of two mesh_fem.
// Date    : March 18, 2005.
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



#ifndef GETFEM_MESH_SUM_H__
#define GETFEM_MESH_SUM_H__

#include <getfem_mesh_fem.h>

namespace getfem {
  typedef std::vector<const std::string *> dof_enrichments;
  
  class fem_sum : public virtual_fem {
    pfem bfem; /* the base FEM which is to be enriched */
    const mesh_sum &mls;
    /* dof_ls_enrichment are stored in the parent mesh_fem 
       the pointer is NULL for non enriched dofs
     */
    std::vector< const dof_ls_enrichment* > dofzones;
    dal::bit_vector ls_index; /* lists only the significant level sets */
  public:
    template <typename IT_LS_ENRICH>
    fem_sum(IT_LS_ENRICH it,pfem pf, const mesh_sum &mls_) : 
      bfem(pf), mls(mls_) {
      if (!(bfem->is_equivalent()))
	DAL_THROW(to_be_done_error,
		  "Sorry, fem_sum for non tau-equivalent "
		  "elements to be done.");
      
      dofzones.assign(it, it + bfem->nb_dof(0));
      init();
    }
    void init();
    void valid();
    void base_value(const base_node &x, base_tensor &t) const;
    void grad_base_value(const base_node &x, base_tensor &t) const;
    void hess_base_value(const base_node &x, base_tensor &t) const;

    void real_base_value(const fem_interpolation_context& c, 
			 base_tensor &t) const;    
    void real_grad_base_value(const fem_interpolation_context& c, 
			      base_tensor &t) const;
    void real_hess_base_value(const fem_interpolation_context& c, 
			      base_tensor &t) const;
    
  };


  /// Describe an adaptable integration method linked to a mesh.
  class mesh_fem_sum : public mesh_fem {
  protected :
    std::vector<const mesh_fem *> mfs;

    const mesh_sum &mls;
    const mesh_fem &mf;
    mutable std::vector<pfem> build_methods;
    mutable bool is_adapted;
    void clear_build_methods();
    void build_method_of_convex(size_type cv);

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
    
    mesh_fem_sum(const getfem_mesh &me) : mesh_fem(me) { is_adapted = false; }
    void set_mesh_fems(const std::vector<const mesh_fem *> &mefs)
    { mfs = mefs; adapt(); }
    void set_mesh_fems(const mesh_fem &mf1, const mesh_fem &mf2)
    { mfs.clear(); mfs.push_back(mf1); mfs.push_back(mf2);  adapt(); }
    void set_mesh_fems(const mesh_fem &mf1, const mesh_fem &mf2, const mesh_fem &mf3)
    { mfs.clear(); mfs.push_back(mf1); mfs.push_back(mf2); mfs.push_back(mf3);  adapt(); }


    ~mesh_fem_sum() { clear_build_methods(); }
  private:
    mesh_fem_sum(const mesh_fem_sum &);
    mesh_fem_sum & operator=(const mesh_fem_sum &);
  };


}  /* end of namespace getfem.                                            */

#endif
  
