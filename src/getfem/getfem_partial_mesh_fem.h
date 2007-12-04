// -*- c++ -*- (enables emacs c++ mode)
//===========================================================================
//
// Copyright (C) 2006-2008 Yves Renard
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
// As a special exception, you may use this file as part of a free software
// library without restriction.  Specifically, if other files instantiate
// templates or use macros or inline functions from this file, or you compile
// this file and link it with other files to produce an executable, this
// file does not by itself cause the resulting executable to be covered by
// the GNU General Public License.  This exception does not however
// invalidate any other reasons why the executable file might be covered by
// the GNU General Public License.
//
//===========================================================================

/**@file getfem_partial_mesh_fem.h
   @author  Yves Renard <Yves.Renard@insa-lyon.fr>, Julien Pommier <Julien.Pommier@insa-toulouse.fr>
   @date June 08, 2006.
   @brief a subclass of getfem::mesh_fem which allows to eliminate a number
          of dof of the original mesh_fem.

   This elimination is done via the pseudo-fem getfem::partial_fem,
   hence it is not very efficient.
*/

#ifndef GETFEM_PARTIAL_MESH_FEM_H__
#define GETFEM_PARTIAL_MESH_FEM_H__

#include "getfem_mesh_fem.h"
#include "getfem_mesh_im.h"

namespace getfem {
  /**
     a subclass of mesh_fem which allows to eliminate a number of dof
     of the original mesh_fem.
  */
  class partial_mesh_fem : public mesh_fem, public boost::noncopyable {
  protected :
    const mesh_fem &mf;
    mutable std::vector<pfem> build_methods;
    mutable bool is_adapted;
    void clear_build_methods();

  public :
    void update_from_context(void) const { is_adapted = false; }

    /** build the mesh_fem keeping only the dof of the original
	mesh_fem which are listed in kept_dof. */
    void adapt(const dal::bit_vector &kept_dof,
	       const dal::bit_vector &rejected_elt = dal::bit_vector());
    void clear(void);

    void receipt(const MESH_CLEAR &);
    void receipt(const MESH_DELETE &);
    
    size_type memsize() const {
      return mesh_fem::memsize(); // + ... ;
    }
    
    partial_mesh_fem(const mesh_fem &mef);

    ~partial_mesh_fem() { clear_build_methods(); }
  };

  /**
     @brief Return a selection of dof who contribute significantly to the
     mass-matrix that would be computed with mf and the integration method mim.

     P represents the dimension on what the integration method operates
     (default mf.linked_mesh().dim()).

     An example of use can be found in the contrib/xfem_contact/ directory.
   */
  dal::bit_vector select_dofs_from_im(const mesh_fem &mf, const mesh_im &mim,
				      unsigned P = unsigned(-1));
  
}  /* end of namespace getfem.                                            */

#endif
  
