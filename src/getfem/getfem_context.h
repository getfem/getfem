// -*- c++ -*- (enables emacs c++ mode)
//===========================================================================
//
// Copyright (C) 2004-2009 Yves Renard
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

/**@file getfem_context.h
   @author  Yves Renard <Yves.Renard@insa-lyon.fr>
   @date June 17, 2004.
   @brief Deal with interdependencies of objects (getfem::context_dependencies).
*/
#ifndef GETFEM_CONTEXT_H__
#define GETFEM_CONTEXT_H__

#include "getfem_config.h"
#include <list>

namespace getfem {
  /**Deal with interdependencies of objects.

   An object can be in three different states :
     NORMAL  : no change is necessary
     CHANGED : something in the context has changed and the update function
               of the object has to be called.
     INVALID : one of the dependencies desappears, the object is invalid
               and should no longer be used.
  
   add_dependency(ct) : add a dependency to the dependency list.

   touch()            : significate to the dependent objects that something
                        has change in the object. This make the dependent
                        objects to be in the CHANGED state

   context_check()    : check if the object has to be updated. if it is the
                        case it makes first a check to the dependency list
                        and call the update function of the object.
                        (the update function of the dependencies are called
                        before the update function of the current object).

   context_valid()    : says if the object has still a valid context
  
   Remarks :

   - A protection against round dependencies exists. In this case, the
   order of call of the update functions can be arbitrary

   - Detection of context changes is very fast (control the
   state). the touch operation can cover the whole tree of dependent
   object.  But this is the case only for the first touch operation
   because once a dependent object is in the CHANGED state it will not
   be considered by next touch operations.
  */
  class context_dependencies {

  protected :
    enum context_state { CONTEXT_NORMAL, CONTEXT_CHANGED, CONTEXT_INVALID };
    mutable context_state state;
    mutable bool touched;
    mutable std::vector<const context_dependencies *> dependencies;
    mutable std::vector<const context_dependencies *> dependent;
    typedef std::vector<const context_dependencies *>::iterator iterator_list;

    void sup_dependent_(const context_dependencies &cd) const;
    void sup_dependency_(const context_dependencies &cd) const;
    void invalid_context(void) const;

  public :
    
    /** this function has to be defined and should update the object when
	the context is modified. */
    virtual void update_from_context(void) const = 0;

    void change_context(void) const
    { if (state == CONTEXT_NORMAL) { state = CONTEXT_CHANGED; touch(); } }
    void add_dependency(const context_dependencies &cd);
    void sup_dependency(const context_dependencies &cd)
    { cd.sup_dependent_(*this); sup_dependency_(cd); }
    bool is_context_valid(void) const { return (state != CONTEXT_INVALID); }
    bool is_context_changed() const { return (state == CONTEXT_CHANGED); }
    /** return true if update_from_context was called */
    bool context_check(void) const;
    void touch(void) const;
    virtual ~context_dependencies();
    context_dependencies() : state(CONTEXT_NORMAL), touched(false) {}

  };

}  /* end of namespace getfem.                                             */


#endif /* GETFEM_CONTEXT_H__  */
