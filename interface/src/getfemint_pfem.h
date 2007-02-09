// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Copyright (C) 2001-2006 Yves Renard, Julien Pommier.
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


#ifndef GETFEMINT_PFEM_H__
#define GETFEMINT_PFEM_H__

#include <getfemint_object.h>
#include <getfem/getfem_fem.h>


namespace getfemint
{
  class getfemint_pfem : public getfem_object {
  private:
    getfem::pfem pf;
    bool nbdof_need_convex_number_;
  public:
    getfemint_pfem(getfem::pfem pf_) {
      assert(workspace == 0);
      pf = pf_;
      nbdof_need_convex_number_ = false;
      this->ikey = getfem_object::internal_key_type(&(*pf_));
    }
    ~getfemint_pfem() {}
    id_type class_id() const { return FEM_CLASS_ID; }
    size_type memsize() const;
    getfem::pfem pfem() { return pf; }

    static getfemint_pfem *get_from(getfem::pfem pf, int flags=0);
    bool &nbdof_need_convex_number() { return nbdof_need_convex_number_; }
  };

  inline bool object_is_pfem(getfem_object *o) {
    return (o->class_id() == FEM_CLASS_ID);
  }

  getfemint_pfem* object_to_pfem(getfem_object *o);
}  /* end of namespace getfemint.                                          */

#endif /* GETFEMINT_PFEM_H__                                               */
