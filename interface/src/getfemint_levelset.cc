#include <getfemint_levelset.h>
#include <getfemint_workspace.h>

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

  void getfemint_levelset::values_from_poly(unsigned idx, 
					    const std::string &s) {
    const getfem::mesh_fem &mf = levelset().get_mesh_fem();
    bgeot::base_poly p = 
      bgeot::read_base_poly(mf.linked_mesh().dim(), s);
    ls->values(idx).resize(mf.nb_dof());
    for (unsigned i=0; i < mf.nb_dof(); ++i) {
      const getfem::base_node x = mf.point_of_dof(i);
      ls->values(idx)[i] = p.eval(x.begin());
    }
  }
}
