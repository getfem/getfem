/*===========================================================================

 Copyright (C) 2007-2015 Julien Pommier.

 This file is a part of GetFEM++

 GetFEM++  is  free software;  you  can  redistribute  it  and/or modify it
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

===========================================================================*/
#include <getfemint_levelset.h>
#include <getfemint_workspace.h>
#include <getfem/getfem_generic_assembly.h>
#include <getfem/getfem_arch_config.h>

namespace getfemint {
  getfemint_levelset* getfemint_levelset::get_from(getfem::level_set *ls,
                                                   int flags) {
    getfem_object *o =
      getfemint::workspace().object(getfem_object::internal_key_type(ls));
    getfemint_levelset *gls = 0;
    if (!o) {
      const getfem::mesh &m = ls->get_mesh_fem().linked_mesh();
      getfemint_mesh *gm =
        getfemint_mesh::get_from(const_cast<getfem::mesh*>(&m),
                                 flags);
      gls = new getfemint_levelset();
      gls->ls = ls;
      gls->ikey = getfem_object::internal_key_type(ls);
      gls->set_flags(flags);
      getfemint::workspace().push_object(gls);
      getfemint::workspace().set_dependance(gls, gm);
    } else gls = dynamic_cast<getfemint_levelset*>(o);
    assert(gls);
    return gls;

  }

  void getfemint_levelset::values_from_func(unsigned idx,
                                            const std::string &s) {

    const getfem::mesh_fem &mf = levelset().get_mesh_fem();
    size_type N = mf.linked_mesh().dim();
    getfem::ga_workspace gw;
    getfem::model_real_plain_vector pt(N);
    gw.add_fixed_size_constant("X", pt);
    if (N >= 1) gw.add_macro("x", "X(1)");
    if (N >= 2) gw.add_macro("y", "X(2)");
    if (N >= 3) gw.add_macro("z", "X(3)");
    if (N >= 4) gw.add_macro("w", "X(4)");
    getfem::ga_function f(gw, s);
    
    f.compile();
    ls->values(idx).resize(mf.nb_dof());
    
    bool is_set = 0;
    for (unsigned i=0; i < mf.nb_dof(); ++i) {
      is_set = 0;
      gmm::copy(mf.point_of_basic_dof(i), pt);
      const bgeot::base_tensor &t = f.eval();
      GMM_ASSERT1(t.size() == 1, "Wrong size of expression result " << s);
      ls->values(idx)[i] = t[0];
    }
   
  }
}
