/*===========================================================================

 Copyright (C) 2001-2020 Yves Renard

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

===========================================================================*/


#include "getfem/getfem_interpolation.h"

namespace getfem {

  size_type mesh_trans_inv::id_of_point(size_type ipt) const {

    if (!ids.empty()) {
      const auto it = ids.find(ipt);
      if (it != ids.end())
        return it->second;
    }
    // otherwise assume that the point id is the point index
    return ipt;
  }

  void mesh_trans_inv::points_on_convex(size_type cv,
                                        std::vector<size_type> &itab) const {
    itab.resize(pts_in_cvx[cv].size()); size_type j = 0;
    for (auto it = pts_in_cvx[cv].cbegin(); it != pts_in_cvx[cv].cend(); ++it)
      itab[j++] = *it;
  }

  size_type mesh_trans_inv::point_on_convex(size_type cv, size_type i) const {
    auto it = pts_in_cvx[cv].cbegin();
    for (size_type j = 0; it != pts_in_cvx[cv].cend() && j < i; ++it, ++j) {}
    GMM_ASSERT1(it != pts_in_cvx[cv].cend(), "internal error");
    return *it;
  }

  void mesh_trans_inv::distribute(int extrapolation, mesh_region rg_source) {
    // this function fills in these two data structures:
    // pts_in_cvx[cv_i] --> set with indices of all points belonging to convex cv_i
    // ref_coords[pt_i] --> local coordinates of point pt_i in the convex it belongs to
    rg_source.from_mesh(msh);
    rg_source.error_if_not_convexes();
    bool all_convexes = (rg_source.id() == mesh_region::all_convexes().id());
    bool projection_into_element(extrapolation == 0);

    size_type nbpts = nb_points();
    size_type nbcvx = msh.nb_allocated_convex();
    ref_coords.resize(nbpts);
    pts_in_cvx.clear();
    pts_in_cvx.resize(nbcvx);

    std::vector<double> dist(nbpts);
    std::vector<size_type> cvx_of_pt(nbpts);
    dal::bit_vector remaining_pts, cv_on_bound;
    remaining_pts.add(0, nbpts);
    scalar_type mult = scalar_type(1);
    do {
      for (dal::bv_visitor j(rg_source.index()); !j.finished(); ++j) {
        if (mult > scalar_type(1) && !(cv_on_bound.is_in(j)))
          continue;
        base_node min, max; /* bound of the box enclosing the convex */
        bgeot::pgeometric_trans pgt = msh.trans_of_convex(j);
        bounding_box(min, max, msh.points_of_convex(j), pgt);
        for (size_type k=0; k < min.size(); ++k) {
          min[k] -= EPS;
          max[k] += EPS;
        }
        if (extrapolation == 2) { // find also points far outside the convex
          if (mult == scalar_type(1)) // in the 1st pass mark all convexes on boundary
            for (short_type f = 0; f < msh.nb_faces_of_convex(j); ++f) {
              size_type neighbor_cv = msh.neighbor_of_convex(j, f);
              if (!all_convexes && neighbor_cv != size_type(-1)) {
                // check if the neighbor is also contained in rg_source ...
                if (!rg_source.is_in(neighbor_cv))
                  cv_on_bound.add(j); // ... if not, treat the element as a boundary one
              } else // boundary element of the overall mesh
                cv_on_bound.add(j);
            }
          if (cv_on_bound.is_in(j)) {
            scalar_type h = scalar_type(0);
            for (size_type k=0; k < min.size(); ++k)
              h = std::max(h, max[k] - min[k]);
            for (size_type k=0; k < min.size(); ++k) {
              min[k] -= mult*h;
              max[k] += mult*h;
            }
          }
        }
        bgeot::kdtree_tab_type boxpts;
        points_in_box(boxpts, min, max);
        if (boxpts.size() > 0)
          gic.init(msh.points_of_convex(j), pgt);
        // check and treat all points found in the bounding box around convex j
        for (const bgeot::index_node_pair &ind_node : boxpts) {
          size_type ind = ind_node.i;
          if (remaining_pts[ind] || dist[ind] > 0) {
            base_node pt_ref;
            bool converged;
            bool gicisin = gic.invert(ind_node.n, pt_ref, converged, EPS,
                                      projection_into_element);
            bool toadd = extrapolation || gicisin;
            double isin = pgt->convex_ref()->is_in(pt_ref);

            if (toadd && !(remaining_pts[ind])) {
              if (isin < dist[ind])
                pts_in_cvx[cvx_of_pt[ind]].erase(ind);
              else
                toadd = false;
            }
            if (toadd) {
//               if (mult > 1.5) {
//                 cout << "adding " << ind << "ref_coord = " << pt_ref
//                      << " cv = " << j << " mult = " << mult << endl;
//               }
              pts_in_cvx[j].insert(ind); // output
              ref_coords[ind] = pt_ref;  // output
              dist[ind] = isin;          // ephemeral
              cvx_of_pt[ind] = j;        // ephemeral
              remaining_pts.sup(ind);    // ephemeral
            }
          }
        }
      }
      mult *= scalar_type(2);
    } while (remaining_pts.card() > 0 && extrapolation == 2);
  }
}  /* end of namespace getfem.                                             */

