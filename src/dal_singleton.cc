#include <dal_singleton.h>
#include <algorithm>
namespace dal {
  std::auto_ptr<singletons_manager> singletons_manager::m;

  void singletons_manager::register_new_singleton(singleton_instance_base *p) {
    if (!m.get()) m.reset(new singletons_manager());
    m->lst.push_back(p);
  }

  static int level_compare(singleton_instance_base *a, singleton_instance_base *b) {
    return a->level() < b->level();
  }

  singletons_manager::~singletons_manager() { 
    /* sort singletons in increasing levels,
       lowest levels will be destroyed first */
    std::sort(m->lst.begin(),m->lst.end(), level_compare);
    std::vector<singleton_instance_base *>::const_iterator 
      it = m->lst.begin(), ite = m->lst.end();
    for ( ; it != ite; ++it) { delete *it; }
  }
}
