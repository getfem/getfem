/* -*- c++ -*- (enables emacs c++ mode) */
/*===========================================================================

 Copyright (C) 2002-2015 Julien Pommier

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

#ifndef GETFEMINT_OBJECT_H__
#define GETFEMINT_OBJECT_H__

#include <getfemint_std.h>

namespace getfemint
{

  static const id_type anonymous_workspace = id_type(-1);
  
  /* common class for all getfem objects addressable from the interface */

  enum { STATIC_OBJ = 1, CONST_OBJ = 2 };

  class getfem_object {
    friend class workspace_data;
    friend class workspace_stack;
  public:
    typedef const void* internal_key_type;
  protected:
    //  public:
    id_type workspace;
    id_type id;
    std::vector<id_type> used_by; /* list of objects which depends on this object */
    internal_key_type ikey; /* generally the pointer to the corresponding 
			       getfem object */

    
    typedef int obj_flags_t;
    obj_flags_t flags; /* if STATIC_OBJ, the linked getfem object is not
			  deleted when the getfem_object is destroyed
			  (example: pfem, mesh_fems obtained with
			  classical_mesh_fem(..) etc. */
    
    /* fonctions reservées au workspace */
    void set_workspace(id_type w) { workspace = w; }
    void set_id(id_type i) { id = i; }

  public:

    getfem_object() : ikey(0) { workspace = 0; id = 0; flags = 0; }

    /* these functions can't be pure virtual functions !?
       it breaks the linking of getfem with libgetfemint.so
    */
    virtual ~getfem_object()
    { id = id_type(-1); ikey = 0; workspace = id = 0x77777777;}
    virtual id_type class_id() const = 0;

    virtual size_type memsize() const { return 0; }

    /* 
       the clear function is called before deletion 
       the object should prepare itself to be deleted
       after or before any other object it uses.
       (for ex. the mesh_fem object may be destroyed
       before or after its linked_mesh)
       
       be careful not do destroy objects which are marked static!
    */
    virtual void clear_before_deletion() {};
   
    /*
      mark the object as a static object (infinite lifetime)
    */
    void set_flags(int v) { flags = v; workspace = anonymous_workspace; }
    bool is_static() const { return flags & STATIC_OBJ; }
    bool is_const()  const { return flags & CONST_OBJ; }

    id_type get_workspace() const { return workspace; }
    id_type get_id() const { return id; }
    bool is_anonymous() const { return workspace == anonymous_workspace; }
    const std::vector<id_type>& get_used_by() const { return used_by; }
  };
}
#endif
