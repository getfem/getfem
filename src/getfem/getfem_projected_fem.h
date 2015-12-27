/* -*- c++ -*- (enables emacs c++ mode) */
/*===========================================================================

 Copyright (C) 2012-2015 Yves Renard, Konstantinos Poulios

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
/**@file getfem_projected_fem.h
   @author Konstantinos Poulios <logari81@googlemail.com>
   @date January 12, 2012.
   @brief FEM which projects a mesh_fem on a different mesh.
*/


#ifndef GETFEM_PROJECTED_FEM_H__
#define GETFEM_PROJECTED_FEM_H__

#include "getfem_fem.h"
#include "getfem_mesh_fem.h"
#include "getfem_mesh_im.h"
#include "bgeot_kdtree.h"
#include "bgeot_geotrans_inv.h"

namespace getfem {

  struct gausspt_projection_data {
    size_type cv; // convex of the source mesh_fem opposite to the gauss point
    short_type f; // convex face of the source mesh_fem opposite to the gauss point
    size_type iflags;     // flags & 1 : there is an element or not
                          // flags & 2 : base_val is stored
                          // flags & 4 : grad_val is stored
    base_node ptref;      // coords on reference element of mf_source element
    base_node normal;     // normal vector at the projected point on the mf_source element
    scalar_type gap;      // gap distance from the gauss point to the projected point
    base_tensor base_val; // optional storage of the base values
    base_tensor grad_val; // optional storage of the grad base values
    std::map<size_type,size_type> local_dof; // correspondance between dof of the
                                             // mf_source element and dof of the projected element.
    gausspt_projection_data() :
      cv(size_type(-1)), f(short_type(-1)), iflags(size_type(-1)) {}
  };

  /** FEM which interpolates the projection of a mesh_fem on a different mesh.

    Note that the memory cost of this method is extremely high!
  */
  class projected_fem : public virtual_fem, public context_dependencies {

  protected :

    struct elt_projection_data {
      short_type f;    // face number on the target mesh
      size_type nb_dof;
      std::map<size_type,gausspt_projection_data> gausspt;
      std::vector<size_type> inddof;
      pintegration_method pim; // for DEBUG
      elt_projection_data() : f(short_type(-1)), nb_dof(0), pim(0) {}
    };

    const mesh_fem &mf_source; // mf_source represents the original finite
                               // element method to be projected.
    const mesh_im &mim_target; // mesh on which mf_source is projected.
                               // Contains also the integration method.
    mesh_region rg_source;
    mesh_region rg_target;

    bool store_values;
    dal::bit_vector blocked_dofs;

    // auxiliary variables
    mutable std::map<size_type,elt_projection_data> elements;
    mutable bgeot::kdtree tree; // Tree containing the nodes of the
                                // projected mf_source dofs
    mutable std::vector<size_type> ind_dof; /* all functions using this work
                                               array should keep it full of
                                               size_type(-1) */
    mutable bgeot::geotrans_inv_convex gic;
    mutable base_tensor taux;
    mutable fem_interpolation_context fictx;
    mutable size_type fictx_cv;
    mutable base_matrix G;
    mutable bgeot::pstored_point_tab pspt_override;
    mutable bgeot::multi_index mi2, mi3;
    mutable base_node ptref;

    void build_kdtree(void) const;

    bool find_a_projected_point(base_node pt, base_node &ptr_proj,
                                size_type &cv_proj, short_type &fc_proj) const;

    virtual void update_from_context(void) const;
    inline void actualize_fictx(pfem pf, size_type cv,
                                const base_node &ptr) const;

  public :

    virtual size_type nb_dof(size_type cv) const;
    virtual size_type index_of_global_dof(size_type cv, size_type i) const;
    virtual bgeot::pconvex_ref ref_convex(size_type cv) const;
    virtual const bgeot::convex<base_node> &node_convex(size_type cv) const;
    virtual bgeot::pstored_point_tab node_tab(size_type) const
    { return pspt_override; }
    void base_value(const base_node &, base_tensor &) const;
    void grad_base_value(const base_node &, base_tensor &) const;
    void hess_base_value(const base_node &, base_tensor &) const;
    void real_base_value(const fem_interpolation_context& c,
                         base_tensor &t, bool = true) const;
    void real_grad_base_value(const fem_interpolation_context& c,
                              base_tensor &t, bool = true) const;
    void real_hess_base_value(const fem_interpolation_context&,
                              base_tensor &, bool = true) const;

    void projection_data(const fem_interpolation_context& c,
                         base_node &normal, scalar_type &gap) const;
    void projection_data(const base_node &pt,
                         base_node &normal, scalar_type &gap) const;

    /** return the list of convexes of the projected mesh_fem which
     *  contain at least one gauss point (should be all convexes)! */
    dal::bit_vector projected_convexes() const;

    /** return the min/max/mean number of gauss points in the convexes
     *  of the projected mesh_fem */
    void gauss_pts_stats(unsigned &ming, unsigned &maxg,
                         scalar_type &meang) const;
    size_type memsize() const;

    projected_fem(const mesh_fem &mf_source_, const mesh_im &mim_target_,
                  size_type rg_source_, size_type rg_target_,
                  dal::bit_vector blocked_dofs_,
                  bool store_val);
    virtual ~projected_fem()
    { DAL_STORED_OBJECT_DEBUG_DESTROYED(this, "Projected fem"); }
  };


  /** create a new projected FEM.
      @param mf_source the mesh_fem that will be projected.
      @param mim_target the integration method on the final mesh (not the mesh
             of mf_source!).
      @param blocked_dofs list of dof of mf_source which won't be projected.
      @param store_val if true, the values/gradients of interpolated base
             function are cached at each gauss point (eats much memory).
  */
  pfem new_projected_fem(const mesh_fem &mf_source, const mesh_im &mim_target,
                         size_type rg_source_ = size_type(-1),
                         size_type rg_target_ = size_type(-1),
                         dal::bit_vector blocked_dofs = dal::bit_vector(),
                         bool store_val = true);

  /** release a projected fem */
  inline void del_projected_fem(pfem pf) { dal::del_stored_object(pf); }


}  /* end of namespace getfem.                                            */

#endif
