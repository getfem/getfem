/* -*- c++ -*- (enables emacs c++ mode) */
/*===========================================================================
 
 Copyright (C) 2001-2015 Yves Renard, Julien Pommier.
 
 This file is a part of GETFEM++
 
 Getfem++  is  free software;  you  can  redistribute  it  and/or modify it
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
 
 As a special exception, you  may use  this file  as it is a part of a free
 software  library  without  restriction.  Specifically,  if   other  files
 instantiate  templates  or  use macros or inline functions from this file,
 or  you compile this  file  and  link  it  with other files  to produce an
 executable, this file  does  not  by itself cause the resulting executable
 to be covered  by the GNU Lesser General Public License.  This   exception
 does not  however  invalidate  any  other  reasons why the executable file
 might be covered by the GNU Lesser General Public License.
 
===========================================================================*/


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
