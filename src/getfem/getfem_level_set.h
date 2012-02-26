/* -*- c++ -*- (enables emacs c++ mode) */
/*===========================================================================
 
 Copyright (C) 1999-2012 Yves Renard
 
 This file is a part of GETFEM++
 
 Getfem++  is  free software;  you  can  redistribute  it  and/or modify it
 under  the  terms  of the  GNU  Lesser General Public License as published
 by  the  Free Software Foundation;  either version 2.1 of the License,  or
 (at your option) any later version.
 This program  is  distributed  in  the  hope  that it will be useful,  but
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 License for more details.
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

/**@file getfem_level_set.h
   @author  Yves Renard <Yves.Renard@insa-lyon.fr>, Julien Pommier <Julien.Pommier@insa-toulouse.fr>
   @date January 31, 2005.
   @brief Define level-sets.
*/
#ifndef GETFEM_LEVEL_SET_H__
#define GETFEM_LEVEL_SET_H__

#include "getfem_mesh_fem.h"
#include "getfem_mesher.h"

namespace getfem {
  /** @brief Define a level-set.  

      In getfem, a levelset is one or two scalar functions, defined on
      a lagrange polynomial mesh_fem.
      
      The (optional) second function is a way to limit the level-set
      to handle cracks for example.
  */

  class level_set : public context_dependencies {

  protected :
    mesh *pmesh;
    dim_type degree_;
    const mesh_fem *mf;
    std::vector<scalar_type> primary_, secondary_;
    bool with_secondary;
    scalar_type shift_ls;     // for the computation of a gap on a level_set.
    // shift the level set on the ref element for mesher_level_set call

    void copy_from(const level_set &ls); // WARNING :  to be updated if
                                         //    some components are added 

  public :

    void set_shift(scalar_type shift_ls_) { shift_ls = shift_ls_; }
    void simplify(scalar_type eps = 0.01);
    void update_from_context(void) const { }
    void reinit(void);
    std::vector<scalar_type> &values(unsigned i = 0)
    { return (i == 0) ? primary_ : secondary_; }
    const std::vector<scalar_type> &values(unsigned i = 0) const
    { return (i == 0) ? primary_ : secondary_; }

    mesher_level_set mls_of_convex(size_type cv, unsigned lsnum = 0,
				   bool inverted = false) const;
    bool has_secondary(void) const { return with_secondary; }
    const mesh_fem &get_mesh_fem(void) const { return *mf; }
    const mesh &linked_mesh() const { return mf->linked_mesh(); }
    dim_type degree() const { return degree_; }
    level_set(mesh &msh, dim_type deg = dim_type(1),
	      bool with_secondary_ = false);
    level_set(const level_set &ls);
    level_set &operator =(const level_set &ls);

    ~level_set();
    size_type memsize() const;
  };
 

}  /* end of namespace getfem.                                             */


#endif /* GETFEM_LEVEL_SET_H__  */
