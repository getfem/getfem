// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Library : GEneric Tool for Finite Element Methods (getfem)
// File    : getfem_level_set.cc : Dealing with level set representation.
//           
// Date    : January 31, 2005.
// Author  : Yves Renard <Yves.Renard@insa-toulouse.fr>
//           Julien Pommier <Julien.Pommier@insa-toulouse.fr>
//
//========================================================================
//
// Copyright (C) 1999-2005 Yves Renard, Julien Pommier
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


#include <getfem_level_set.h>

namespace getfem {

  bool level_set::mf_key::operator <(const mf_key &a) const {
    if (pmesh < a.pmesh) return true; else
      if (a.pmesh < pmesh) return false; else
	if (order < a.order) return true; else return false;
  }

  dal::shared_ptr<mesh_fem> level_set::add_mesh_fem(getfem_mesh &mesh,
						    dim_type o) {
    mf_key key(mesh, o);
    std::map<mf_key, pmesh_fem>::iterator it = mesh_fems.find(key);
    if (it == mesh_fems.end()) {
      mesh_fem *pmf = new mesh_fem(mesh);
      pmf->set_classical_finite_element(o);
      return (mesh_fems[key] = pmesh_fem(pmf));
    }
    else return it->second;
  }
  
  void level_set::sup_mesh_fem(getfem_mesh &mesh, dim_type o) {
    mf_key key(mesh, o);
    std::map<mf_key, pmesh_fem>::iterator it = mesh_fems.find(key);
    if (it != mesh_fems.end()) {
      if (it->second.use_count() <= 2) mesh_fems.erase(key);
    }
  }
 

}  /* end of namespace getfem.                                             */

