// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Library : GEneric Tool for Finite Element Methods (getfem)
// File    : getfem_level_set.h : Dealing with level set representation.
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

/**@file getfem_level_set.h
   @brief Define level-sets.
*/
#ifndef GETFEM_LEVEL_SET_H__
#define GETFEM_LEVEL_SET_H__

#include <dal_shared_ptr.h>
#include <getfem_mesh_fem.h>
#include <getfem_mesher.h>

namespace getfem {
  /** @brief Define a level-set.  

      In getfem, a levelset is one or two scalar functions, defined on
      a lagrange polynomial mesh_fem.
      
      The (optional) second function is a way to limit the level-set
      to handle cracks for example.
  */

  class level_set : public context_dependencies {

  protected :
    typedef dal::shared_ptr<mesh_fem> pmesh_fem;
    pmesh_fem add_mesh_fem(getfem_mesh &mesh, dim_type o);
    void sup_mesh_fem(getfem_mesh &mesh, dim_type o);
    getfem_mesh *pmesh;
    dim_type degree_;
    pmesh_fem mf;
    std::vector<scalar_type> primary_, secondary_;
    bool with_secondary;

  public :

    void update_from_context(void) const { }    
    void reinit(void);
    std::vector<scalar_type> &values(unsigned i = 0)
    { return (i == 0) ? primary_ : secondary_; }
    const std::vector<scalar_type> &values(unsigned i = 0) const
    { return (i == 0) ? primary_ : secondary_; }

    mesher_level_set mls_of_convex(size_type cv, unsigned lsnum = 0,
				   bool inverted = false) const;
    bool has_secondary(void) const { return with_secondary; }
    mesh_fem &get_mesh_fem(void) { return *mf; }
    dim_type degree() const { return degree_; }
    level_set(getfem_mesh &mesh, dim_type deg = dim_type(1),
	      bool with_secondary_ = false)
      : pmesh(&mesh), degree_(deg), mf(add_mesh_fem(mesh, deg)),
	with_secondary(with_secondary_) {
      primary_.resize(mf->nb_dof());
      secondary_.resize(mf->nb_dof());
      this->add_dependency(*mf);
    }
    
    ~level_set() { sup_mesh_fem(*pmesh, degree_); }

  };
 

}  /* end of namespace getfem.                                             */


#endif /* GETFEM_LEVEL_SET_H__  */
