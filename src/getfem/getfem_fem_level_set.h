// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Copyright (C) 2004-2007 Yves Renard
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

/**@file getfem_fem_level_set.h
   @author  Yves Renard <Yves.Renard@insa-lyon.fr>, Julien Pommier <Julien.Pommier@insa-toulouse.fr>
   @date March 09, 2005.
   @brief FEM associated with getfem::mesh_fem_level_set objects.

   To be fixed : dependencies. The mesh fem using this fem will not
   depend on the mesh fem arguments.
*/


#ifndef GETFEM_FEM_LEVEL_SET_H__
#define GETFEM_FEM_LEVEL_SET_H__

#include "getfem_mesh_level_set.h"


namespace getfem {
  /*
  struct zoneset_t {    
    typedef unsigned char zid_t;
    dal::bit_vector ls_idx;
    std::vector<zid_t> table;
    zid_t &operator()(const std::string &p) {
      return table[getpos(p)];
    }
    size_type getpos(const std::string &s) {
      size_type p2 = 1, pos = 0;
      for (dal::bv_visitor i(ls_idx); !i.finished(); ++i, p2 *= 2) {
	pos += (s[i] == '+') ? p2 : 0;
      }
      return pos;
    }
    void merge(const zoneset_t &z) {
      dal::bit_vector idx2 = ls_idx | z.ls_idx;
      std::vector<size_type> s1, s2 = strides_for(idx2)
    }
  };
  */

  /** FEM associated with getfem::mesh_fem_level_set objects.
  */
  class fem_level_set : public virtual_fem {
    pfem bfem; /* the base FEM which is to be enriched */
    const mesh_level_set &mls;
    size_type xfem_index;
    std::vector< const mesh_level_set::zoneset * > dofzones;
    dal::bit_vector ls_index; /* lists only the significant level sets */
    std::string common_ls_zones;
    void find_zone_id(const fem_interpolation_context &c, 
		      std::vector<unsigned> &ids) const;
  public:
    template <typename IT_LS_ENRICH>
    fem_level_set(IT_LS_ENRICH it,pfem pf, const mesh_level_set &mls_,
		  size_type xfi) : 
      bfem(pf), mls(mls_), xfem_index(xfi) {
      dofzones.assign(it, it + bfem->nb_dof(0));
      init();
    }
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
  
