/* -*- c++ -*- (enables emacs c++ mode) */
/*===========================================================================
 
 Copyright (C) 2013-2013 Yves Renard.
 
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

/**\file getfemint_cont_struct.h
   \brief getfem::multi_contact_frame interface
*/

#include <getfemint_std.h>
#include <getfemint_object.h>
#include <getfemint_workspace.h>
#include <getfem/getfem_contact_and_friction_common.h>

namespace getfemint {

  class getfemint_multi_contact_frame : public getfem_object {
  private:
    getfem::multi_contact_frame *s;
    getfemint_multi_contact_frame(getfem::multi_contact_frame *s_) {
      assert(workspace == 0);
      s = s_;
      ikey = getfem_object::internal_key_type(s);
    }

  public:
    ~getfemint_multi_contact_frame() {}
    id_type class_id() const { return MULTI_CONTACT_FRAME_CLASS_ID; }
    size_type memsize() const {
      return sizeof(getfem::multi_contact_frame) 
        + s->ct_pairs().size()
        * (sizeof(getfem::multi_contact_frame::contact_pair)
           + sizeof(scalar_type) * (s->dim() + 1) * 3
           );
    }

    static getfemint_multi_contact_frame*
    get_from(getfem::multi_contact_frame *ps, int flags = 0) {
      getfem_object *o =
	getfemint::workspace().object(getfem_object::internal_key_type(ps));
      getfemint_multi_contact_frame *gs = NULL;
      if (!o) {
	gs = new getfemint_multi_contact_frame(ps);
	gs->set_flags(flags);
	getfemint::workspace().push_object(gs);
      } else gs = dynamic_cast<getfemint_multi_contact_frame*>(o);
      assert(gs);
      return gs;
    }

    getfem::multi_contact_frame &multi_contact_frame() { return *s; }
  };
  
  inline bool object_is_multi_contact_frame(getfem_object *o) {
    return o->class_id() == MULTI_CONTACT_FRAME_CLASS_ID;
  }

  inline getfemint_multi_contact_frame* object_to_multi_contact_frame(getfem_object *o) {
    if (object_is_multi_contact_frame(o)) return (getfemint_multi_contact_frame*)o;
    else THROW_INTERNAL_ERROR;
  }
}
