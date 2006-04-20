// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Copyright (C) 1999-2006 Yves Renard
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


#include <getfem_level_set.h>

namespace getfem {

  void level_set::reinit(void) {
    primary_.resize(mf->nb_dof());
    if (has_secondary()) secondary_.resize(mf->nb_dof());
  }

  mesher_level_set level_set::mls_of_convex(size_type cv, unsigned lsnum,
					    bool inverted) const {
    if (!mf->linked_mesh().convex_index().is_in(cv)) 
      DAL_THROW(dal::failure_error, "convex " << cv << " is not in the level set mesh!");
    if (!mf->fem_of_element(cv)) DAL_INTERNAL_ERROR("");
    std::vector<scalar_type> coeff(mf->nb_dof_of_element(cv));
    for (size_type i = 0; i < coeff.size(); ++i)
      coeff[i] = (!inverted ? scalar_type(1) : scalar_type(-1)) * 
	values(lsnum)[mf->ind_dof_of_element(cv)[i]];
    //cout << "mls_of_convex[lsnum=" << lsnum << "] : coeff = " << coeff << "\n";
    return mesher_level_set(mf->fem_of_element(cv), coeff);
  }
 
  size_type level_set::memsize() const {
    return sizeof(*this) + 
      primary_.capacity() * sizeof(scalar_type) + 
      secondary_.capacity() * sizeof(scalar_type);
  }
}  /* end of namespace getfem.                                             */

