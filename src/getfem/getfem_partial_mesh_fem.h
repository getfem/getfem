/* -*- c++ -*- (enables emacs c++ mode) */
/*===========================================================================
 
 Copyright (C) 2006-2015 Yves Renard, Julien Pommier
 
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

/**@file getfem_partial_mesh_fem.h
   @author Yves Renard <Yves.Renard@insa-lyon.fr>,
   @author Julien Pommier <Julien.Pommier@insa-toulouse.fr>
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
  class partial_mesh_fem : public mesh_fem, public boost::noncopyable,
                           public dal::static_stored_object {
  protected :
    const mesh_fem &mf;
    mutable bool is_adapted;

  public :
    void update_from_context(void) const
    { mf.context_check(); is_adapted = false; }

    /** build the mesh_fem keeping only the dof of the original
        mesh_fem which are listed in kept_dof. */
    void adapt(const dal::bit_vector &kept_dof,
               const dal::bit_vector &rejected_elt = dal::bit_vector());
    void clear(void);

    pfem fem_of_element(size_type cv) const
    { return  mf.fem_of_element(cv); }

    virtual dim_type get_qdim() const { return mf.get_qdim(); }
    virtual const bgeot::multi_index &get_qdims() const
    { return mf.get_qdims(); }

    virtual void set_qdim(dim_type) {
      GMM_ASSERT1(false, "The Qdim of a partial_mesh_fem is the same as "
                  "the original fem");
    }

    virtual void set_qdim(dim_type, dim_type) {
      GMM_ASSERT1(false, "The Qdim of a partial_mesh_fem is the same as "
                  "the original fem");
    }

    virtual void set_qdim(dim_type, dim_type, dim_type, dim_type) {
      GMM_ASSERT1(false, "The Qdim of a partial_mesh_fem is the same as "
                  "the original fem");
    }

    virtual void set_qdim(const bgeot::multi_index &) {
      GMM_ASSERT1(false, "The Qdim of a partial_mesh_fem is the same as "
                  "the original fem");
    }

    ind_dof_ct ind_basic_dof_of_element(size_type cv) const
    { return  mf.ind_basic_dof_of_element(cv); }

    ind_dof_face_ct
    ind_basic_dof_of_face_of_element(size_type cv, short_type f) const
    { return  mf.ind_basic_dof_of_face_of_element(cv, f); }

    size_type nb_basic_dof_of_face_of_element(size_type cv, short_type f) const
    { return  mf.nb_basic_dof_of_face_of_element(cv, f); }

    size_type nb_basic_dof_of_element(size_type cv) const
    { return  mf.nb_basic_dof_of_element(cv); }

    base_node point_of_basic_dof(size_type cv, size_type i) const
    { return  mf.point_of_basic_dof(cv, i); }

    base_node point_of_basic_dof(size_type d) const
    { return  mf.point_of_basic_dof(d); }

    dim_type basic_dof_qdim(size_type d) const
    { return  mf.basic_dof_qdim(d); }

    size_type first_convex_of_basic_dof(size_type d) const
    { return mf.first_convex_of_basic_dof(d); }

    const mesh::ind_cv_ct &convex_to_basic_dof(size_type d) const
    { return mf.convex_to_basic_dof(d); }

    size_type nb_dof(void) const {
      context_check();
      return use_reduction ? gmm::mat_nrows(R_) : mf.nb_dof();
    }

    size_type nb_basic_dof(void) const
    { return mf.nb_basic_dof(); }

    dal::bit_vector basic_dof_on_region(const mesh_region &b) const
    { return mf.basic_dof_on_region(b); }

    // invalid function for a mesh change.
    // dal::bit_vector retrieve_kept_dofs() const;

    void read_from_file(std::istream &)
    { GMM_ASSERT1(false, "You cannot directly read this kind of mesh_fem"); }
    void write_to_file(std::ostream &ost) const;
    void write_to_file(const std::string &name, bool with_mesh=false) const;

    partial_mesh_fem(const mesh_fem &mef);
    partial_mesh_fem(const mesh_fem *mef);

  };

  typedef boost::intrusive_ptr<partial_mesh_fem> ppartial_mesh_fem;

  /**
     @brief Return a selection of dof who contribute significantly to the
     mass-matrix that would be computed with mf and the integration method mim.

     P represents the dimension on what the integration method operates
     (default mf.linked_mesh().dim()).

     An example of use can be found in the contrib/xfem_contact/ directory.

     A more efficient algorithm is now present in gmm_range_basis.h
   */
  dal::bit_vector select_dofs_from_im(const mesh_fem &mf, const mesh_im &mim,
                                      unsigned P = unsigned(-1));


}  /* end of namespace getfem.                                            */

#endif

