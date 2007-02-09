#include <getfemint_mesh.h>
#include <getfemint_workspace.h>

namespace getfemint {
  getfemint_mesh::getfemint_mesh(getfem::mesh *m_) {
    assert(workspace == 0);
    assert(workspace == 0);
    m = m_;
    ikey = getfem_object::internal_key_type(m);
  }
  
  getfemint_mesh::~getfemint_mesh() {
    if (!is_static()) {
      m->clear();
      delete m;
    }
  }
  
  getfemint_mesh *getfemint_mesh::get_from(getfem::mesh *m, int f) {
    getfem_object *o = 
      getfemint::workspace().object(getfem_object::internal_key_type(m));
    getfemint_mesh *gm = 0;
    if (!o) {
      gm = new getfemint_mesh(m); 
      gm->set_flags(f); 
      getfemint::workspace().push_object(gm);
    } else gm = dynamic_cast<getfemint_mesh*>(o);
    assert(gm);
    return gm;
  }
}
