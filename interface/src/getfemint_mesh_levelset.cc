#include <getfemint_mesh_levelset.h>
#include <getfemint_workspace.h>

namespace getfemint {
  getfemint_mesh_levelset::getfemint_mesh_levelset() {
    mls = 0;
  }

  getfemint_mesh_levelset::~getfemint_mesh_levelset() {
    if (!is_static()) delete mls; 
    mls = 0;
  }

  getfemint_mesh_levelset *
  getfemint_mesh_levelset::get_from(getfem::mesh_level_set *mls,
				    int flags) {
    getfem_object *o = 
      getfemint::workspace().object(getfem_object::internal_key_type(mls));
    getfemint_mesh_levelset *gmls = 0;
    if (!o) {
      getfemint_mesh *gm = 
	getfemint_mesh::get_from(const_cast<getfem::mesh*>(&mls->linked_mesh()),
				 flags);
      gmls = new getfemint_mesh_levelset(); 
      gmls->mls = mls;
      gmls->ikey = getfem_object::internal_key_type(mls);
      gmls->set_flags(flags); 
      getfemint::workspace().push_object(gmls);
      getfemint::workspace().set_dependance(gmls, gm);
    } else gmls = dynamic_cast<getfemint_mesh_levelset*>(o);
    assert(gmls);
    return gmls;
  }
}
