/* -*- c++ -*- (enables emacs c++ mode) */
/*===========================================================================

 Copyright (C) 2012-2015 Tomas Ligursky, Yves Renard, Konstantinos Poulios.

 This file is a part of GetFEM++

 GetFEM++  is  free software;  you  can  redistribute  it  and/or modify it
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

/**\file getfemint_cont_struct.h
   \brief getfem::cont_struct_getfem_model interface
*/

#include <getfemint_object.h>
#include <getfemint.h>
#include <getfem/getfem_continuation.h>

namespace getfemint {

  class getfemint_cont_struct : public getfem_object {
  private:
    getfem::cont_struct_getfem_model *s;
    getfemint_cont_struct(getfem::cont_struct_getfem_model *s_);
  public:
    ~getfemint_cont_struct();
    id_type class_id() const { return CONT_STRUCT_CLASS_ID; }
    size_type memsize() const {
      return s->estimated_memsize();
    }

    static getfemint_cont_struct*
    get_from(getfem::cont_struct_getfem_model *s_, int flags = 0);

    getfem::cont_struct_getfem_model &cont_struct() { return *s; }
  };

  inline bool object_is_cont_struct(getfem_object *o) {
    return o->class_id() == CONT_STRUCT_CLASS_ID;
  }

  inline getfemint_cont_struct* object_to_cont_struct(getfem_object *o) {
    if (object_is_cont_struct(o)) return (getfemint_cont_struct*)o;
    else THROW_INTERNAL_ERROR;
  }
}
