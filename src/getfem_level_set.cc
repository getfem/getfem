// -*- c++ -*- (enables emacs c++ mode)
//===========================================================================
//
// Copyright (C) 1999-2008 Yves Renard
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


#include "getfem/getfem_level_set.h"

namespace getfem {

  level_set::level_set(mesh &msh, dim_type deg,
		       bool with_secondary_)
    : pmesh(&msh), degree_(deg), mf(&classical_mesh_fem(msh, deg)),
      with_secondary(with_secondary_) {
    primary_.resize(mf->nb_dof());
    secondary_.resize(mf->nb_dof());
    this->add_dependency(*mf);
  }

  level_set::~level_set() { }

  void level_set::reinit(void) {
    primary_.resize(mf->nb_dof());
    if (has_secondary()) secondary_.resize(mf->nb_dof());
    touch();
  }

  mesher_level_set level_set::mls_of_convex(size_type cv, unsigned lsnum,
					    bool inverted) const {
    assert(this); assert(mf); 
    GMM_ASSERT1(mf->linked_mesh().convex_index().is_in(cv), "convex " << cv
		<< " is not in the level set mesh!");
    GMM_ASSERT1(mf->fem_of_element(cv), "Internal error");
    GMM_ASSERT1(!mf->is_reduced(), "Internal error");
    std::vector<scalar_type> coeff(mf->nb_basic_dof_of_element(cv));
    GMM_ASSERT1(values(lsnum).size() == mf->nb_dof(),
		"Inconsistent state in the levelset: nb_dof=" << 
		mf->nb_dof() << ", values(" << lsnum << ").size=" << 
		values(lsnum).size());
    for (size_type i = 0; i < coeff.size(); ++i)
      coeff[i] = (!inverted ? scalar_type(1) : scalar_type(-1)) * 
	values(lsnum)[mf->ind_basic_dof_of_element(cv)[i]];
    //cout << "mls_of_convex[lsnum=" << lsnum << "] : coeff = " << coeff << "\n";
    return mesher_level_set(mf->fem_of_element(cv), coeff);
  }
 
  size_type level_set::memsize() const {
    return sizeof(*this) + 
      primary_.capacity() * sizeof(scalar_type) + 
      secondary_.capacity() * sizeof(scalar_type);
  }

  void level_set::simplify(scalar_type eps) {
    for (dal::bv_visitor cv(mf->linked_mesh().convex_index());
	 !cv.finished(); ++cv) {
      scalar_type h = mf->linked_mesh().convex_radius_estimate(cv);
      for (size_type i = 0; i < mf->nb_basic_dof_of_element(cv); ++i) {
	size_type dof = mf->ind_basic_dof_of_element(cv)[i];
	if (gmm::abs(primary_[dof]) < h*eps) {
	  primary_[dof] = scalar_type(0);
	  // cout << "Simplify dof " << dof << " : " << mf->point_of_dof(dof) << endl;
	}
	if (has_secondary() && gmm::abs(secondary_[dof]) < h*eps)
	    secondary_[dof] = scalar_type(0);
      }

    }
    touch();
  }

  

}  /* end of namespace getfem.                                             */

