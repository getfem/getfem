/* -*- c++ -*- (enables emacs c++ mode) */
/*===========================================================================

 Copyright (C) 2005-2020 Yves Renard

 This file is a part of GetFEM

 GetFEM  is  free software;  you  can  redistribute  it  and/or modify it
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
  class mesh_im : public context_dependencies, virtual public dal::static_stored_object {
  private :
    void copy_from(const mesh_im &mim);

  protected :
    bool is_lower_dim;
    dal::dynamic_array<pintegration_method> ims;
    dal::bit_vector im_convexes;
    const mesh *linked_mesh_;
    mutable gmm::uint64_type v_num_update, v_num;
    pintegration_method auto_add_elt_pim; /* im for automatic addition     */
                          /* of element option. (0 = no automatic addition)*/

  public :
    void update_from_context(void) const;
    gmm::uint64_type version_number(void) const
    { context_check(); return v_num; }

    /** Set the im for automatic addition
     *  of element option. pim=0 disables the automatic addition.
     */
    void set_auto_add(pintegration_method pim)
    { auto_add_elt_pim = pim; }

    /** Get the set of convexes where an integration method has been assigned.
     */
    inline const dal::bit_vector &convex_index(void) const
    { context_check(); return im_convexes; }

    bool is_lower_dimensional() const { return is_lower_dim; }

    /// Give a reference to the linked mesh of type mesh.
    const mesh &linked_mesh() const
    { return linked_mesh_ ? *linked_mesh_ : dummy_mesh(); }
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
        set_integration_method(linked_mesh().convex_index(),pim);
        and set_auto_add(pim)
        @endcode
    */
    void set_integration_method(pintegration_method ppi);
    /** Set an approximate integration method chosen to be exact for
        polynomials of degree 'im_degree'.
    */
    void set_integration_method(const dal::bit_vector &cvs,
                                dim_type im_degree);

    /** Set an approximate integration method chosen to be exact for
        polynomials of degree 'im_degree' on the whole mesh.
    */
    void set_integration_method(dim_type im_degree);

    /** return the integration method associated with an element (in
        no integration is associated, the function will crash! use the
        convex_index() of the mesh_im to check that a fem is
        associated to a given convex) */
    virtual pintegration_method int_method_of_element(size_type cv) const
    { return  ims[cv]; }
    void clear(void);

    size_type memsize() const {
      context_check(); 
      return sizeof(mesh_im) + ims.memsize() + im_convexes.memsize();
    }

    void init_with_mesh(const mesh &me);
    mesh_im(const mesh &me);
    mesh_im();
    virtual ~mesh_im();
    mesh_im(const mesh_im &mim);
    mesh_im &operator=(const mesh_im &mim);

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
  };

  /** Dummy mesh_im for default parameter of functions. */
  const mesh_im &dummy_mesh_im();

}  /* end of namespace getfem.                                             */


#endif /* GETFEM_MESH_IM_H__  */
