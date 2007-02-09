// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Copyright (C) 2001-2006 Yves Renard.
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


#include <getfemint_pfem.h>
#include <getfemint_workspace.h>
#include <getfem/getfem_interpolated_fem.h>

namespace getfemint
{
  size_type getfemint_pfem::memsize() const {
    const getfem::interpolated_fem*p = 
      dynamic_cast<const getfem::interpolated_fem*>(&(*pf));
    if (p) return p->memsize();
    else return 0;
  }

  getfemint_pfem *getfemint_pfem::get_from(getfem::pfem pf, int flags) {
    getfem_object *o = 
      getfemint::workspace().object(getfem_object::internal_key_type(&(*pf)));
    getfemint_pfem *gfi_pf = 0;
    if (!o) {
      gfi_pf = new getfemint_pfem(pf); 
      gfi_pf->set_flags(flags);
      getfemint::workspace().push_object(gfi_pf);
    } else gfi_pf = dynamic_cast<getfemint_pfem*>(o);
    assert(gfi_pf);
    return gfi_pf;
  }


  /*
  struct static_getfemint_pfem_list : public std::vector<getfemint_pfem*> {};
  struct map_pfem_to_id : public std::map<getfem::pfem, id_type> {};

  getfemint_pfem *getfemint_pfem::get_from(getfem::pfem pf, bool is_static) {
    getfemint_pfem *p = 0;
    if (is_static) {
      static_getfemint_pfem_list &l = 
	dal::singleton<static_getfemint_pfem_list>::instance();
      for (unsigned i=0; i < l.size(); ++i) {
	if (pf == l[i]->pfem()) { p = l[i]; break; }
      }
      if (p == 0) {
	p = new getfemint_pfem(pf);
	l.push_back(p);
	p->set_static();
	p->set_id((l.size() - 1) | 0x80000000);
	return p;
      }
    } else {
      p = new getfemint_pfem(pf);
      getfemint::workspace().push_object(p);
    }
    dal::singleton<map_pfem_to_id>::instance()[pf] = p->get_id();
    return p;
  }

  getfemint_pfem *getfemint_pfem::get_from(id_type id) {
    if (id & 0x80000000) {
      static_getfemint_pfem_list &l = 
	dal::singleton<static_getfemint_pfem_list>::instance();
      return l.at(id & 0x7fffffff);
    } else {
      getfem_object *o = 
	getfemint::workspace().object(id, 
			  name_of_getfemint_class_id(FEM_CLASS_ID));
      return object_to_pfem(o);
    }
  }
  */

  /*getfemint_pfem::~getfemint_pfem() {
    map_pfem_to_id &m = dal::singleton<map_pfem_to_id>::instance();
    assert(m[pf] == get_id());
    m.erase(pf);
    if (!is_static())
      dal::del_stored_object(pf); // force the removal of the fem object
  }
  */
  /*  id_type ind_pfem(getfem::pfem pf) {
    map_pfem_to_id &m = dal::singleton<map_pfem_to_id>::instance();
    if (m.find(pf) == m.end()) {
      getfemint_pfem *p = getfemint_pfem::get_from(pf, true);
      return p->get_id();
    } else return m[pf];
  }
  */

  getfemint_pfem* object_to_pfem(getfem_object *o) {
    if (object_is_pfem(o)) return ((getfemint_pfem*)o);
    else THROW_INTERNAL_ERROR;
  }



}  /* end of namespace getfemint.                                          */
