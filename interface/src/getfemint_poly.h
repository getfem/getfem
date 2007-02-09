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


#ifndef GETFEMINT_POLY_H__
#define GETFEMINT_POLY_H__

#include <getfemint_std.h>
#include <getfem/bgeot_config.h>
#include <getfem/bgeot_poly.h>
#include <getfemint_object.h>

namespace getfemint
{
  class getfemint_poly : public getfem_object {
  private:
    bgeot::base_poly *p;
    
  public:
    getfemint_poly() { 
      p = new bgeot::base_poly();
    }
    ~getfemint_poly() {
      delete p;
    }
    id_type class_id() const { return POLY_CLASS_ID; }

    bgeot::base_poly& poly() { return *p; }
    const bgeot::base_poly& poly() const { return *p; }
  };

  inline bool object_is_poly(getfem_object *o) {
    return (o->class_id() == POLY_CLASS_ID);
  }

  inline getfemint_poly* object_to_poly(getfem_object *o) {
    if (object_is_poly(o)) return ((getfemint_poly*)o);
    else THROW_INTERNAL_ERROR;
  }
}  /* end of namespace getfemint.                                          */

#endif /* GETFEMINT_POLY_H__                                               */
