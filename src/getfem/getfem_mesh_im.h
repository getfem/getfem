// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Copyright (C) 2005-2007 Yves Renard
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

/**@file getfem_mesh_im.h
   @author  Yves Renard <Yves.Renard@insa-lyon.fr>
   @date January 26, 2005.
   @brief Define the getfem::mesh_im class (integration of getfem::mesh_fem).
*/
#ifndef GETFEM_MESH_IM_H__
#define GETFEM_MESH_IM_H__

#include "getfem_integration.h"
#include "getfem_fem.h"
#include "getfem_mesh.h"

namespace getfem {

  /// Describe an integration method linked to a mesh.
  class mesh_im : public mesh_receiver, public context_dependencies {
  protected :
    
    dal::dynamic_array<pintegration_method> ims;
    dal::bit_vector im_convexes;
    mesh *linked_mesh_;
    bool is_valid_;
  public :
    bool is_valid() const { return is_valid_; }
    void update_from_context(void) const {}

    /** Get the set of convexes where an integration method has been assigned.
     */
    inline const dal::bit_vector &convex_index(void) const
    { return im_convexes; }
    
    /// Give a reference to the linked mesh of type mesh.
    mesh &linked_mesh(void) const { return *linked_mesh_; }
    /** Set the integration method of a convex.

	@param cv the convex number

	@param pim the integration method, typically obtained with
	@code getfem::int_method_descriptor("IM_SOMETHING(..)") 
	@endcode
     */
    void set_integration_method(size_type cv, pintegration_method pim);
    /** Set the integration method on all the convexes of indexes in bv,
     *  which is of type dal::bit_vector.
     */
    void set_integration_method(const dal::bit_vector &cvs, 
				pintegration_method pim);
    /** shortcut for
	@code 
	set_integration_method(linked_mesh().convex_index(),ppi); 
	@endcode
    */
    void set_integration_method(pintegration_method ppi);
    /** Set an approximate integration method chosen to be exact for
	polynomials of degree 'im_degree'
    */
    void set_integration_method(const dal::bit_vector &cvs, 
				dim_type im_degree);
    
    /** return the integration method associated with an element (in
	no integration is associated, the function will crash! use the
	convex_index() of the mesh_im to check that a fem is
	associated to a given convex) */
    virtual pintegration_method int_method_of_element(size_type cv) const
    { return  ims[cv]; }
    void clear(void);
    void receipt(const MESH_CLEAR &);
    void receipt(const MESH_DELETE &);
    void receipt(const MESH_ADD_CONVEX &m) { mesh_receiver::receipt(m); }
    void receipt(const MESH_SUP_CONVEX &m);
    void receipt(const MESH_SWAP_CONVEX &m);
    void receipt(const MESH_REFINE_CONVEX &m);
    
    size_type memsize() const {
      return 
	sizeof(mesh_im) +
	ims.memsize() + im_convexes.memsize();
    }
    
    mesh_im(mesh &me);
    virtual ~mesh_im();
    /** Read the mesh_im from a stream.
        @param ist the stream. */
    void read_from_file(std::istream &ist);
    /** Read the mesh_im from a file.
        @param name the file name. */
    void read_from_file(const std::string &name);
    /** Write the mesh_im to a stream. */
    void write_to_file(std::ostream &ost) const;
    /** Write the mesh_im to a file. 

	@param name the file name

	@param with_mesh if set, then the linked_mesh() will also be
	saved to the file.
    */
    void write_to_file(const std::string &name, bool with_mesh=false) const;
  private:
    mesh_im(const mesh_im &);
    mesh_im & operator=(const mesh_im &);
  };
  
}  /* end of namespace getfem.                                             */


#endif /* GETFEM_MESH_IM_H__  */
