/* -*- c++ -*- (enables emacs c++ mode) */
/*===========================================================================
 
 Copyright (C) 2004-2015 Yves Renard
 
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
/**@file getfem_interpolated_fem.h
   @author  Yves Renard <Yves.Renard@insa-lyon.fr>
   @date October 29, 2004.
   @brief FEM which interpolates a mesh_fem on a different mesh.

   To be corrected : dependencies. The mesh fem using this fem will
   not depend on the mesh fem arguments.
*/


#ifndef GETFEM_INTERPOLATED_FEM_H__
#define GETFEM_INTERPOLATED_FEM_H__

#include "getfem_fem.h"
#include "getfem_mesh_fem.h"
#include "getfem_mesh_im.h"
#include "bgeot_rtree.h"
#include "bgeot_geotrans_inv.h"

namespace getfem {

  // Object representing global transformation. To be derived.

  struct virtual_interpolated_func {
    virtual void val(const base_node&, base_node &) const
    { GMM_ASSERT1(false, "this interpolated_func has no value"); }
    virtual void grad(const base_node&, base_matrix &) const
    { GMM_ASSERT1(false, "this interpolated_func has no gradient"); }
    virtual void hess(const base_node&, base_matrix &) const
    { GMM_ASSERT1(false, "this interpolated_func has no hessian"); }
    virtual ~virtual_interpolated_func() {}
  };


  typedef const virtual_interpolated_func *pinterpolated_func;


  struct gausspt_interpolation_data {
    size_type elt; // convex of the interpolated mesh_fem under the gauss point
    size_type iflags;     // flags & 1 : there is an element or not
                          // flags & 2 : base_val is stored
                          // flags & 4 : grad_val is stored
    base_node ptref;      // coords on reference element of mf1 element
    base_tensor base_val; // optional storage of the base values
    base_tensor grad_val; // optional storage of the grad base values
    std::vector<size_type> local_dof; // correspondance between dof of the
    // mf1 element and dof of the interpolated element.
    gausspt_interpolation_data() : elt(size_type(-1)), iflags(size_type(-1)) {}
  };

  /** FEM which interpolates a mesh_fem on a different mesh.

    Note that the memory cost of this method is extremely high!
  */
  class interpolated_fem : public virtual_fem, public context_dependencies {
    
  protected :

    struct elt_interpolation_data {
      size_type nb_dof;
      std::vector<gausspt_interpolation_data> gausspt;
      std::vector<size_type> inddof;
      pintegration_method pim; // for DEBUG
      elt_interpolation_data() : nb_dof(0), pim(0) {}
    }; 

    const mesh_fem &mf;    // mf represents the original finite element method
                           // to be interpolated.
    const mesh_im &mim;    // mesh on which mf1 is interpolated. contains
                           // also the integration method.
    pinterpolated_func pif; // optional transformation

    bool store_values;
    dal::bit_vector blocked_dof;

    // auxiliary variables
    mutable std::vector<elt_interpolation_data> elements;
    mutable bgeot::rtree boxtree; // Tree containing the bounding box
                                  // of mf1 elements
    mutable std::vector<size_type> ind_dof; /* all functions using this work
					       array should keep it full of
					       size_type(-1) */
    mutable size_type cv_stored;
    mutable bgeot::rtree::pbox_set boxlst;
    mutable bgeot::geotrans_inv_convex gic;
    mutable base_tensor taux;
    mutable fem_interpolation_context fictx;
    mutable size_type fictx_cv;
    mutable base_matrix G;
    mutable std::vector<base_node> node_tab_;
    mutable bgeot::multi_index mi2, mi3;
    mutable base_node ptref;
    mutable gmm::dense_matrix<scalar_type> trans;

    void build_rtree(void) const;

    bool find_a_point(base_node pt, base_node &ptr,
		      size_type &cv) const;

    void update_from_context(void) const;
    inline void actualize_fictx(pfem pf, size_type cv,
				const base_node &ptr) const;

  public :

    virtual size_type nb_dof(size_type cv) const;
    virtual size_type index_of_global_dof(size_type cv, size_type i) const;
    virtual bgeot::pconvex_ref ref_convex(size_type cv) const;
    virtual const bgeot::convex<base_node> &node_convex(size_type cv) const;
    virtual bgeot::pstored_point_tab node_tab(size_type) const;
    void base_value(const base_node &, base_tensor &) const;
    void grad_base_value(const base_node &, base_tensor &) const;
    void hess_base_value(const base_node &, base_tensor &) const;
    void real_base_value(const fem_interpolation_context& c, 
			 base_tensor &t, bool = true) const;
    void real_grad_base_value(const fem_interpolation_context& c, 
			      base_tensor &t, bool = true) const;
    void real_hess_base_value(const fem_interpolation_context&, 
			      base_tensor &, bool = true) const;

    /** return the list of convexes of the interpolated mesh_fem which
     *  contain at least one gauss point (should be all convexes)! */
    dal::bit_vector interpolated_convexes() const;
    
    /** return the min/max/mean number of gauss points in the convexes
     *  of the interpolated mesh_fem */
    void gauss_pts_stats(unsigned &ming, unsigned &maxg,
			 scalar_type &meang) const; 
    size_type memsize() const;
  private:
    interpolated_fem(const mesh_fem &mef, const mesh_im &meim,
		     pinterpolated_func pif_ = 0,
		     dal::bit_vector blocked_dof = dal::bit_vector(),
		     bool store_val = true);
    
    friend pfem new_interpolated_fem(const mesh_fem &mef, const mesh_im &mim,
				     pinterpolated_func pif,
				     dal::bit_vector blocked_dof,
				     bool store_val);
  };
  

  /** create a new interpolated FEM. 
      @param mef the mesh_fem that will be interpolated.
      @param mim the integration method on the final mesh (not the mesh
             of mef!).
      @param pif an optional geometric transformation applied to
             mef.linked_mesh() (used for getfem::spider_fem)
      @param blocked_dof list of dof of mef which won't be interpolated.
      @param store_val if true, the values/gradients of interpolated base
             function are cached at each gauss point (eats much memory).
  */
  pfem new_interpolated_fem(const mesh_fem &mef, const mesh_im &mim,
			    pinterpolated_func pif = 0,
			    dal::bit_vector blocked_dof = dal::bit_vector(),
			    bool store_val = true);

  /** release an interpolated fem */
  inline void del_interpolated_fem(pfem pf) { dal::del_stored_object(pf); }

  
}  /* end of namespace getfem.                                            */

#endif
  
