/*===========================================================================
 
 Copyright (C) 2001-2015 Yves Renard
 
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
 
===========================================================================*/


#include "getfem/getfem_interpolation.h"

namespace getfem {

  size_type mesh_trans_inv::id_of_point(size_type ipt) const {

    if (!ids.empty()) {
      map_iterator it=ids.find(ipt);
      if (it != ids.end())
        return it->second;
    }
    // otherwise assume that the point id is the point index
    return ipt;
  }

  void mesh_trans_inv::points_on_convex(size_type cv,
                                        std::vector<size_type> &itab) const {
    itab.resize(pts_cvx[cv].size()); size_type j = 0;
    for (set_iterator it = pts_cvx[cv].begin(); it != pts_cvx[cv].end(); ++it)
      itab[j++] = *it;
  }

  size_type mesh_trans_inv::point_on_convex(size_type cv, size_type i) const {
    set_iterator it = pts_cvx[cv].begin();
    for (size_type j = 0; it != pts_cvx[cv].end() && j < i; ++it, ++j);
    GMM_ASSERT1(it != pts_cvx[cv].end(), "internal error");    
    return *it;
  }

  void mesh_trans_inv::distribute(int extrapolation, mesh_region rg_source) {

    rg_source.from_mesh(msh);
    rg_source.error_if_not_convexes();
    bool all_convexes = (rg_source.id() == mesh_region::all_convexes().id());

    size_type nbpts = nb_points();
    size_type nbcvx = msh.convex_index().last_true() + 1;
    ref_coords.resize(nbpts);
    std::vector<double> dist(nbpts);
    std::vector<size_type> cvx_pts(nbpts);
    pts_cvx.clear(); pts_cvx.resize(nbcvx);
    base_node min, max, pt_ref; /* bound of the box enclosing the convex */
    bgeot::kdtree_tab_type boxpts;
    dal::bit_vector npt, cv_on_bound;
    npt.add(0, nbpts);
    scalar_type mult = scalar_type(1);

    do {
      for (dal::bv_visitor j(rg_source.index()); !j.finished(); ++j) {
        if (mult > scalar_type(1) && !(cv_on_bound.is_in(j))) continue;
        bgeot::pgeometric_trans pgt = msh.trans_of_convex(j);
        bounding_box(min, max, msh.points_of_convex(j), pgt);
        for (size_type k=0; k < min.size(); ++k) { min[k]-=EPS; max[k]+=EPS; }
        if (extrapolation == 2) {
          if (mult == scalar_type(1))
            for (short_type f = 0; f < msh.nb_faces_of_convex(j); ++f) {
              size_type neighbour_cv = msh.neighbour_of_convex(j, f);
              if (!all_convexes && neighbour_cv != size_type(-1)) {
                // check if the neighbour is also contained in rg_source ...
                if (!rg_source.is_in(neighbour_cv)) 
                  cv_on_bound.add(j); // ... if not, treat the element as a boundary one
              }
              else // boundary element of the overall mesh
                cv_on_bound.add(j);
            }
          if (cv_on_bound.is_in(j)) {
            scalar_type h = scalar_type(0);
            for (size_type k=0; k < min.size(); ++k)
              h = std::max(h, max[k] - min[k]);
            for (size_type k=0; k < min.size(); ++k)
              { min[k]-=mult*h; max[k]+=mult*h; }
          }
        }
        points_in_box(boxpts, min, max);

        if (boxpts.size() > 0) gic.init(msh.points_of_convex(j), pgt);
        
        for (size_type l = 0; l < boxpts.size(); ++l) {
          size_type ind = boxpts[l].i;
          if (npt[ind] || dist[ind] > 0) {
            bool converged;
            bool gicisin = gic.invert(boxpts[l].n, pt_ref, converged, EPS);
            bool toadd = extrapolation || gicisin;
            double isin = pgt->convex_ref()->is_in(pt_ref);
            
            if (toadd && !(npt[ind])) {
              if (isin < dist[ind]) pts_cvx[cvx_pts[ind]].erase(ind);
              else toadd = false;
            }
            if (toadd) {
//               if (mult > 1.5) {
//                 cout << "adding " << ind << "ref_coord = " << pt_ref
//                      << " cv = " << j << " mult = " << mult << endl;
//               }
              ref_coords[ind] = pt_ref;
              dist[ind] = isin; cvx_pts[ind] = j;
              pts_cvx[j].insert(ind);
              npt.sup(ind);
            }
          }
        }
      }
      mult *= scalar_type(2);
    } while (npt.card() > 0 && extrapolation == 2);
  }
}  /* end of namespace getfem.                                             */

