/* -*- c++ -*- (enables emacs c++ mode) */
/*===========================================================================

 Copyright (C) 2006-2015 Yves Renard

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

/**@file bgeot_mesh.h
   @author  Yves Renard <Yves.Renard@insa-lyon.fr>
   @date February 15, 2006.
   @brief Basic mesh definition
*/

#ifndef BGEOT_MESH_H__
#define BGEOT_MESH_H__

#include "bgeot_mesh_structure.h"
#include "bgeot_geometric_trans.h"
#include "bgeot_node_tab.h"

namespace bgeot {

  /** @internal mesh structure + points
   */
  class basic_mesh :  public bgeot::mesh_structure {

  public :
    
    // typedef basic_mesh_point_comparator pt_comp;
    typedef bgeot::node_tab PT_TAB;
    typedef bgeot::mesh_structure::ind_cv_ct ind_cv_ct;
    typedef bgeot::mesh_structure::ind_set ind_set;
    typedef bgeot::mesh_structure::ind_pt_face_ct ind_pt_face_ct;
    typedef bgeot::mesh_structure::point_ct point_ct;

    typedef gmm::tab_ref_index_ref
    <PT_TAB::const_iterator, ind_cv_ct::const_iterator> ref_mesh_pt_ct;
    typedef gmm::tab_ref_index_ref
    <PT_TAB::const_iterator, ind_pt_face_ct::const_iterator>
    ref_mesh_face_pt_ct;
    
    
  protected :
    
    PT_TAB pts;

    dal::dynamic_array<bgeot::pgeometric_trans> gtab;
    dal::bit_vector trans_exists;
 
  public :

    dim_type dim(void) const { return pts.dim(); }

    bgeot::pgeometric_trans trans_of_convex(size_type ic) const {
      GMM_ASSERT2(trans_exists[ic], "internal error");
      return gtab[ic]; 
    }

    const PT_TAB &points(void) const { return pts; }

    ref_mesh_pt_ct points_of_convex(size_type ic) const {
      const ind_cv_ct &rct = ind_points_of_convex(ic);
      return ref_mesh_pt_ct(pts.begin(), rct.begin(), rct.end());
    } 

    size_type add_point(const base_node &pt,
                        const scalar_type tol=scalar_type(0)) {
      return pts.add_node(pt, tol);
    }

    template<class ITER>
    size_type add_convex(bgeot::pgeometric_trans pgt, ITER ipts) { 
      bool present;
      size_type i = mesh_structure::add_convex(pgt->structure(), ipts,
                                               &present);
      gtab[i] = pgt; trans_exists[i] = true;
      return i;
    }

    size_type add_segment(size_type a, size_type b) { 
      size_type ipt[2]; ipt[0] = a; ipt[1] = b;
      return add_convex(simplex_geotrans(1, 1), &(ipt[0]));
    }

    size_type add_triangle(size_type a, size_type b, size_type c) {
      size_type ipt[3]; ipt[0] = a; ipt[1] = b; ipt[2] = c;
      return add_convex(simplex_geotrans(2, 1), &(ipt[0]));
    }

    size_type add_tetrahedron(size_type a, size_type b,
                              size_type c, size_type d) {
      size_type ipt[4]; ipt[0] = a; ipt[1] = b; ipt[2] = c; ipt[3] = d;
      return add_convex(simplex_geotrans(3, 1), &(ipt[0]));
    }
  };

  typedef basic_mesh *pbasic_mesh;

}

#endif /* BGEOT_MESH_H__                                         */
