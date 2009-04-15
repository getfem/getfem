// -*- c++ -*- (enables emacs c++ mode)
//===========================================================================
//
// Copyright (C) 2004-2009 Yves Renard
//
// This file is a part of GETFEM++
//
// Getfem++  is  free software;  you  can  redistribute  it  and/or modify it
// under  the  terms  of the  GNU  Lesser General Public License as published
// by  the  Free Software Foundation;  either version 2.1 of the License,  or
// (at your option) any later version.
// This program  is  distributed  in  the  hope  that it will be useful,  but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
// License for more details.
// You  should  have received a copy of the GNU Lesser General Public License
// along  with  this program;  if not, write to the Free Software Foundation,
// Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.
//
//===========================================================================
#include "getfem/getfem_modeling.h"

namespace getfem {


  void mdbrick_abstract_common_base::parameters_set_uptodate(void) {
    for (PARAM_MAP::iterator it = parameters.begin();
	 it != parameters.end(); ++it)
      it->second->set_uptodate();
  }
  
  bool mdbrick_abstract_common_base::parameters_is_any_modified(void) const {
    for (PARAM_MAP::const_iterator it = parameters.begin();
	 it != parameters.end(); ++it)
      if (it->second->is_modified()) return true;
    return false;
  }

  void mdbrick_abstract_common_base::update_from_context(void) const {
    nb_total_dof = 0;
    nb_total_constraints = 0;
    total_mixed_variables.clear();
    is_linear_ = proper_is_linear_;
    is_symmetric_ = proper_is_symmetric_;
    is_coercive_ = proper_is_coercive_;
    mesh_fems.resize(0); mesh_ims.resize(0); mesh_fem_positions.resize(0);
    /* get information from parent bricks */
    for (size_type i = 0; i < sub_bricks.size(); ++i) {
      sub_bricks[i]->context_check();
      for (size_type j = 0; j < sub_bricks[i]->mesh_fems.size(); ++j) {
	mesh_fems.push_back(sub_bricks[i]->mesh_fems[j]);
	mesh_fems_info.push_back(sub_bricks[i]->mesh_fems_info[j]);
	mesh_fem_positions.push_back(nb_total_dof 
				     + sub_bricks[i]->mesh_fem_positions[j]);
      }
      for (size_type j = 0; j < sub_bricks[i]->mesh_ims.size(); ++j) {
	mesh_ims.push_back(sub_bricks[i]->mesh_ims[j]);
      }
      is_linear_ = is_linear_ && sub_bricks[i]->is_linear();
      is_symmetric_ = is_symmetric_ && sub_bricks[i]->is_symmetric();
      is_coercive_ = is_coercive_ && sub_bricks[i]->is_coercive();
      if (i == 0)
	total_mixed_variables |= sub_bricks[i]->total_mixed_variables;
      else {
	for (dal::bv_visitor k(sub_bricks[i]->total_mixed_variables);
	     !k.finished(); ++k)
	  total_mixed_variables.add(k+nb_total_dof);
      }
      nb_total_dof += sub_bricks[i]->nb_total_dof;
      nb_total_constraints += sub_bricks[i]->nb_total_constraints;
    }
    /* merge with information from this brick */
    for (size_type j = 0; j < proper_mesh_fems.size(); ++j) {
      mesh_fems.push_back(proper_mesh_fems[j]);
      mesh_fems_info.push_back(proper_mesh_fems_info[j]);
      mesh_fem_positions.push_back(nb_total_dof);
      nb_total_dof += proper_mesh_fems[j]->nb_dof();
    }
    for (size_type j = 0; j < proper_mesh_ims.size(); ++j)
      mesh_ims.push_back(proper_mesh_ims[j]);
    for (size_type j = 0; j < proper_boundary_cond_info.size(); ++j) {
      mesh_fems_info[proper_boundary_cond_info[j].num_fem]
	.add_boundary(proper_boundary_cond_info[j].num_bound,
		      proper_boundary_cond_info[j].bc);
    }
    /* call the customizable update procedure */
    const_cast<mdbrick_abstract_common_base *>(this)->proper_update_();
    nb_total_dof += proper_additional_dof;
    nb_total_constraints += proper_nb_constraints;
    total_mixed_variables |= proper_mixed_variables;
    // cerr << " nb_total_dof         = " << nb_total_dof  << " (" << proper_additional_dof << ")\n";
    // cerr << " nb_total_constraints = " << nb_total_constraints  << " (" << proper_nb_constraints << ")\n";
  }
  
  


}

