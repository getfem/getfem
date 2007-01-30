// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Copyright (C) 1999-2007 Yves Renard
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

/**@file getfem_level_set.h
   @author  Yves Renard <Yves.Renard@insa-toulouse.fr>, Julien Pommier <Julien.Pommier@insa-toulouse.fr>
   @date January 31, 2005.
   @brief Define level-sets.
*/
#ifndef GETFEM_LEVEL_SET_H__
#define GETFEM_LEVEL_SET_H__

#include "getfem_mesh_fem.h"
#include "getfem_mesher.h"

namespace getfem {
  /** @brief Define a level-set.  

      In getfem, a levelset is one or two scalar functions, defined on
      a lagrange polynomial mesh_fem.
      
      The (optional) second function is a way to limit the level-set
      to handle cracks for example.
  */

  class level_set : public context_dependencies {

  protected :
    mesh *pmesh;
    dim_type degree_;
    const mesh_fem *mf;
    std::vector<scalar_type> primary_, secondary_;
    bool with_secondary;

  public :

    void simplify(scalar_type eps = 0.01);
    void update_from_context(void) const { }
    void reinit(void);
    std::vector<scalar_type> &values(unsigned i = 0)
    { return (i == 0) ? primary_ : secondary_; }
    const std::vector<scalar_type> &values(unsigned i = 0) const
    { return (i == 0) ? primary_ : secondary_; }

    mesher_level_set mls_of_convex(size_type cv, unsigned lsnum = 0,
				   bool inverted = false) const;
    bool has_secondary(void) const { return with_secondary; }
    const mesh_fem &get_mesh_fem(void) const { return *mf; }
    const mesh &linked_mesh() const { return mf->linked_mesh(); }
    dim_type degree() const { return degree_; }
    level_set(mesh &msh, dim_type deg = dim_type(1),
	      bool with_secondary_ = false);
    ~level_set();
    size_type memsize() const;
  };
 

}  /* end of namespace getfem.                                             */


#endif /* GETFEM_LEVEL_SET_H__  */
