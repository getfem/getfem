/* -*- c++ -*- (enables emacs c++ mode) */
/*===========================================================================

 Copyright (C) 2009-2015 Yves Renard.

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

/**\file getfemint_models.h
   \brief getfem::models interface
*/

#include <getfemint.h>
#include <getfemint_object.h>
#include <getfem/getfem_mesher.h>

namespace getfemint {

  class getfemint_mesher_object : public getfemint::getfem_object {
  private:
    getfem::pmesher_signed_distance msd;
    getfemint_mesher_object(getfem::pmesher_signed_distance msd_) {
      assert(workspace == 0);
      msd  = msd_;
      ikey = getfem_object::internal_key_type(msd.get());
    }

  public:
    ~getfemint_mesher_object() {}
    id_type class_id() const { return MESHER_OBJECT_CLASS_ID; }
    size_type memsize() const {
      return sizeof(getfem::mesher_signed_distance);
    }
    static getfemint_mesher_object *
    get_from(getfem::pmesher_signed_distance pmsd, int flags = 0) {
      getfem_object *o =
	getfemint::workspace().object
	(getfem_object::internal_key_type(pmsd.get()));
      getfemint_mesher_object *gpgf = NULL;
      if (!o) {
	gpgf = new getfemint_mesher_object(pmsd);
	gpgf->set_flags(flags);
	getfemint::workspace().push_object(gpgf);
      } else gpgf = dynamic_cast<getfemint_mesher_object*>(o);
      assert(gpgf);
      return gpgf;
    }

    getfem::pmesher_signed_distance &mesher_object() { return msd; }

  };

  inline bool object_is_mesher_object(getfem_object *o) {
    return o->class_id() == MESHER_OBJECT_CLASS_ID;
  }

  inline getfemint_mesher_object* object_to_mesher_object(getfem_object *o) {
    if (object_is_mesher_object(o)) return (getfemint_mesher_object*)o;
    else THROW_INTERNAL_ERROR;
  }
}
