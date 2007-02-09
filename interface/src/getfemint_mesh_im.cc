#include <getfemint_mesh_im.h>
#include <getfemint_workspace.h>

namespace getfemint {
  getfemint_mesh_im::getfemint_mesh_im(getfem::mesh_im *mim_, 
				       id_type idmesh) {
    linked_mesh_id_ = idmesh;
    mim = mim_;
    ikey = getfem_object::internal_key_type(mim);
  }

  getfemint_mesh_im::~getfemint_mesh_im() {
    if (!is_static()) delete mim;
    mim = 0;
  }

  getfemint_mesh_im *
  getfemint_mesh_im::get_from(getfem::mesh_im *mim, 
			      int flags) {
    getfem_object *o = 
      getfemint::workspace().object(getfem_object::internal_key_type(mim));
    getfemint_mesh_im *gmim = 0;
    if (!o) {
      getfemint_mesh *gm = 
	getfemint_mesh::get_from(const_cast<getfem::mesh*>(&mim->linked_mesh()),
				 flags);
      gmim = new getfemint_mesh_im(mim, gm->get_id());
      gmim->set_flags(flags);
      getfemint::workspace().push_object(gmim);
      getfemint::workspace().set_dependance(gmim, gm);
    } else gmim = dynamic_cast<getfemint_mesh_im*>(o);
    assert(gmim);
    return gmim;
  }

  getfemint_mesh_im *
  getfemint_mesh_im::new_from(getfemint_mesh *m) { 
    //assert(::workspace == 0);
    getfem::mesh_im *mim = new getfem::mesh_im(m->mesh());
    getfemint_mesh_im *gmim = get_from(mim);
    assert(gmim->linked_mesh_id() == m->get_id());
    return gmim;
  }
}
