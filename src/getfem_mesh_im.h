// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Library : GEneric Tool for Finite Element Methods (getfem)
// File    : getfem_mesh_im.h : Integration methods on convex meshes.
//           
// Date    : January 26, 2005.
// Author  : Yves Renard <Yves.Renard@insa-toulouse.fr>
//
//========================================================================
//
// Copyright (C) 2005-2005 Yves Renard
//
// This file is a part of GETFEM++
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; version 2 of the License.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software Foundation,
// Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
//
//========================================================================


#ifndef GETFEM_MESH_IM_H__
#define GETFEM_MESH_IM_H__

#include <getfem_mesh.h>
#include <getfem_integration.h>
#include <getfem_precomp.h>

namespace getfem {

  /// Describe an integration method linked to a mesh.
  class mesh_im : public getfem_mesh_receiver, public context_dependencies {
  protected :
    
    dal::dynamic_array<pintegration_method> ims;
    dal::bit_vector im_convexes;
    getfem_mesh *linked_mesh_;
    
  public :
    
    typedef base_node point_type;
    void update_from_context(void) const {}

    /** Gives in a structure dal::bit\_vector all convexes of the
     *          mesh where an integration method is defined.
     */
    inline const dal::bit_vector &convex_index(void) const
    { return im_convexes; }
    
    /// Gives a reference to the linked mesh of type getfem\_mesh.
    getfem_mesh &linked_mesh(void) const { return *linked_mesh_; }
    /** Set the integration method on the convex of index i
     */
    void set_integration_method(size_type cv, pintegration_method pim);
    /** Set the integration method on all the convexes of indexes in bv,
     *  which is of type dal::bit\_vector.
     */
    void set_integration_method(const dal::bit_vector &cvs, 
				pintegration_method pim);
    /** shortcut for
	set_integration_method(linked_mesh().convex_index(),ppi); */
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
    pintegration_method int_method_of_element(size_type cv) const
    { return  ims[cv]; }
    void clear(void);
    /* explicit calls to parent class 
       for HP aCC and mipspro CC who complain about hidden functions 
       (they're right)
    */
    void receipt(const MESH_CLEAR &);
    void receipt(const MESH_DELETE &);
    void receipt(const MESH_ADD_CONVEX &m) {getfem_mesh_receiver::receipt(m); }
    void receipt(const MESH_SUP_CONVEX &m);
    void receipt(const MESH_SWAP_CONVEX &m);
    
    size_type memsize() const {
      return 
	sizeof(mesh_im) +
	ims.memsize() + im_convexes.memsize();
    }
    
    mesh_im(getfem_mesh &me);
    virtual ~mesh_im();
    void read_from_file(std::istream &ist);
    void read_from_file(const std::string &name);
    void write_to_file(std::ostream &ost) const;
    void write_to_file(const std::string &name, bool with_mesh=false) const;
  private:
    mesh_im(const mesh_im &);
    mesh_im & operator=(const mesh_im &);
  };
  
}  /* end of namespace getfem.                                             */


#endif /* GETFEM_MESH_IM_H__  */
