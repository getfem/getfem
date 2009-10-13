// -*- c++ -*- (enables emacs c++ mode)
//===========================================================================
//
// Copyright (C) 2009 Luis Saavedra.
//
// This file is a part of GETFEM++
//
// Getfem++  is  free software;  you  can  redistribute  it  and/or modify it
// under  the  terms  of the  GNU  Lesser General Public License as published
// by  the  Free Software Foundation;  either version 2.1 of the License,  or
// (at your option) any later version.
// This program  is  distributed  in  the  hope  that it will be useful,  but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
// License for more details.
// You  should  have received a copy of the GNU Lesser General Public License
// along  with  this program;  if not, write to the Free Software Foundation,
// Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.
//
// As a special exception, you  may use  this file  as it is a part of a free
// software  library  without  restriction.  Specifically,  if   other  files
// instantiate  templates  or  use macros or inline functions from this file,
// or  you compile this  file  and  link  it  with other files  to produce an
// executable, this file  does  not  by itself cause the resulting executable
// to be covered  by the GNU Lesser General Public License.  This   exception
// does not  however  invalidate  any  other  reasons why the executable file
// might be covered by the GNU Lesser General Public License.
//
//===========================================================================
// $Id$
/**\file getfemint_global_function.h
   \brief getfem::abstract_xy_function wrapper class.
*/

#ifndef GETFEMINT_GLOBAL_FUNCTION_H__
#define GETFEMINT_GLOBAL_FUNCTION_H__

#include <getfemint_std.h>
#include <getfemint_object.h>
#include <getfem/getfem_mesh_fem_global_function.h>

namespace getfemint
{
  class getfemint_global_function : public getfem_object {
  private:
    getfem::abstract_xy_function *pgf;
  public:
    ~getfemint_global_function();
    id_type class_id() const { return GLOBAL_FUNCTION_CLASS_ID; }
    getfem::abstract_xy_function& global_function() { return *pgf; }
    const getfem::abstract_xy_function& global_function() const { return *pgf; }

    static getfemint_global_function* get_from(getfem::abstract_xy_function *pabs,
                                               int flags = 0);

    getfemint_global_function(getfem::abstract_xy_function *pabs);
  };

  inline bool object_is_global_function(getfem_object *o) {
    return (o->class_id() == GLOBAL_FUNCTION_CLASS_ID);
  }

  inline getfemint_global_function* object_to_global_function(getfem_object *o) {
    if (object_is_global_function(o)) return ((getfemint_global_function*)o);
    else THROW_INTERNAL_ERROR;
  }
}  /* end of namespace getfemint.                                          */
#endif /* GETFEMINT_GLOBAL_FUNCTION_H__                                    */
