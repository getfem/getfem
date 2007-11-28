// -*- c++ -*- (enables emacs c++ mode)
//===========================================================================
//
// Copyright (C) 2002-2008 Yves Renard
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

/**@file getfem_external_data_fem.h
   @author  Yves Renard <Yves.Renard@insa-lyon.fr>
   @date June 28, 2003.
   @brief A pseudo FEM allowing to define any function on a getfem::mesh_fem. 
*/
#include "getfem_fem.h"
#include "getfem_mesh_fem.h"
#include "bgeot_geotrans_inv.h"

namespace getfem
{

  /** A pseudo FEM allowing to define any function on a getfem::mesh_fem.

      Inheriting from this FEM, it is possible to define a mesh_fem which interpolates exactly a given function. The resulting mesh_fem is scalar, and has only 1 dof.

      Overload at least real_base_value.
  */
  class external_data_fem : public virtual_fem {
    
  public :    
    void base_value(const base_node &, base_tensor &) const
    { GMM_ASSERT1(false, "Uninstantied function"); }
    void grad_base_value(const base_node &, base_tensor &) const
    { GMM_ASSERT1(false, "Uninstantied function"); }
    void hess_base_value(const base_node &, base_tensor &) const
    { GMM_ASSERT1(false, "Uninstantied function"); }

    void real_base_value(const fem_interpolation_context &, 
			 base_tensor &, bool = true) const
    { GMM_ASSERT1(false, "Uninstantied function"); } 
    
    void real_grad_base_value(const fem_interpolation_context &, 
			      base_tensor &, bool = true) const
    { GMM_ASSERT1(false, "Uninstantied function"); }
    
    void real_hess_base_value(const fem_interpolation_context &, 
			      base_tensor &, bool = true) const
    { GMM_ASSERT1(false, "Uninstantied function"); }

    external_data_fem(const bgeot::pconvex_ref cvr_, size_type dim = 1) {
      cvr = cvr_;
      is_equiv = real_element_defined = true;
      is_polycomp = is_pol = is_lag = false;
      es_degree = 5;
      ntarget_dim = dim;
      init_cvs_node();
      base_node pt(cvr->structure()->dim()); pt.fill(0.001);
      add_node(lagrange_dof(cvr->structure()->dim()), pt);
    }
  };


}  /* end of namespace getfem.                                            */
