#include <getfemint_mesh_fem.h>
#include <getfemint_workspace.h>

namespace getfemint {

  getfemint_mesh_fem::getfemint_mesh_fem(getfem::mesh_fem *mf_, 
					 id_type idmesh) {
    linked_mesh_id_ = idmesh;
    mf = mf_;
    ikey = getfem_object::internal_key_type(mf);
  }

  getfemint_mesh_fem::~getfemint_mesh_fem() {
    if (!is_static()) delete mf;
    mf = 0;
  }

  getfemint_mesh_fem *
  getfemint_mesh_fem::get_from(getfem::mesh_fem *mf, 
			       int flags) {
    getfem_object *o = 
      getfemint::workspace().object(getfem_object::internal_key_type(mf));
    getfemint_mesh_fem *gmf = 0;
    if (!o) {
      getfemint_mesh *gm = 
	getfemint_mesh::get_from(const_cast<getfem::mesh*>(&mf->linked_mesh()),
				 flags);
      gmf = new getfemint_mesh_fem(mf, gm->get_id());
      gmf->set_flags(flags);
      getfemint::workspace().push_object(gmf);
      getfemint::workspace().set_dependance(gmf, gm);

    } else gmf = dynamic_cast<getfemint_mesh_fem*>(o);
    assert(gmf);
    return gmf;
  }

  getfemint_mesh_fem *
  getfemint_mesh_fem::new_from(getfemint_mesh *m, 
			       unsigned q_dim) { 
    //assert(::workspace == 0);
    getfem::mesh_fem *mf = new getfem::mesh_fem(m->mesh());
    mf->set_qdim(q_dim);
    getfemint_mesh_fem *gmf = get_from(mf);
    assert(gmf->linked_mesh_id() == m->get_id());
    return gmf;
  }
}
