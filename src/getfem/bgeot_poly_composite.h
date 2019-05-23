/* -*- c++ -*- (enables emacs c++ mode) */
/*===========================================================================

 Copyright (C) 2002-2017 Yves Renard

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

/**@file bgeot_poly_composite.h
   @author  Yves Renard <Yves.Renard@insa-lyon.fr>
   @date August 26, 2002.
   @brief Handle composite polynomials.

   Composite polynomials are used in hierarchical FEM, composite geometric
   transformations and composite fems.
*/

#ifndef BGEOT_POLY_COMPOSITE_H__
#define BGEOT_POLY_COMPOSITE_H__

#include "bgeot_poly.h"
#include "bgeot_mesh.h"
#include "bgeot_rtree.h"

// TODO : Use of rtree instead of dal::dynamic_tree_sorted<base_node,
//        imbricated_box_less>


namespace bgeot {

  /// A comparison function for bgeot::base_node
  struct imbricated_box_less
  {
    mutable int exp_max, exp_min;
    mutable scalar_type c_max;
    unsigned base;

    /// comparaison function
    int operator()(const base_node &x, const base_node &y) const;

    imbricated_box_less(unsigned ba = 10, int emi = -15, int ema = -2) {
      base = ba; exp_max = ema; exp_min = emi;
      c_max = pow(double(base), double(-exp_max));
    }
  };



  struct mesh_precomposite {

    typedef dal::dynamic_tree_sorted<base_node, imbricated_box_less> PTAB;

    const basic_mesh *msh;
    PTAB vertexes;
    rtree box_tree;
    std::map<size_type, std::vector<size_type>> box_to_convexes_map;
    std::vector<base_matrix> gtrans;
    std::vector<scalar_type> det;
    std::vector<base_node> orgs;
    
    const basic_mesh &linked_mesh(void) const { return *msh; }
    size_type nb_convex(void) const { return gtrans.size(); }
    dim_type dim(void) const { return msh->dim(); }
    pgeometric_trans trans_of_convex(size_type ic) const
    { return msh->trans_of_convex(ic); }
    void initialise(const basic_mesh &m);
    
    mesh_precomposite(const basic_mesh &m);
    mesh_precomposite(void) : msh(0) {}
  };

  typedef const mesh_precomposite *pmesh_precomposite;

  class polynomial_composite {

  protected :
    const mesh_precomposite *mp;
    std::map<size_type, dal::pstatic_stored_object_key> polytab;
    bool local_coordinate;  // are the polynomials described on the local
                            // coordinates of each sub-element or on global coordinates.
    base_poly default_poly;

  public :

    template <class ITER> scalar_type eval(const ITER &it, size_type l = -1) const;
    void derivative(short_type k);
    void set_poly_of_subelt(size_type l, const base_poly &poly);
    const base_poly &poly_of_subelt(size_type l) const;
    size_type nb_subelt() const { return polytab.size(); }

    polynomial_composite(bool lc = true) : local_coordinate(lc) {}
    polynomial_composite(const mesh_precomposite &m, bool lc = true);

  };

  inline std::ostream &operator <<
  (std::ostream &o, const polynomial_composite& P) {
    o << "poly_composite [";
    for (size_type i = 0; i < P.nb_subelt(); ++i) {
      if (i != 0) o << ", " << P.poly_of_subelt(i);
    }
    o << "]";
    return o;
  }

  template <class ITER>
  void mult_diff_transposed(
    const base_matrix &M, const ITER &it, const base_node &p1, base_node &p2) {
    for (dim_type d = 0; d < p2.size(); ++d) {
      p2[d] = 0;
      auto &col = mat_col(M, d);
      for (dim_type i = 0; i < p1.size(); ++i) p2[d] += col[i] * (*(it + i) - p1[i]);
    }
  }

  template <class ITER>
  scalar_type polynomial_composite::eval(const ITER &it, size_type l) const {

    if (l != size_type(-1)) {
      if (!local_coordinate) return poly_of_subelt(l).eval(it);
      base_node p(mp->dim());
      std::copy(it, it + mp->dim(), p.begin());
      mult_diff_transposed(mp->gtrans[l], it, mp->orgs[l], p);
      return poly_of_subelt(l).eval(p.begin());
    }

    auto dim = mp->dim();
    base_node p0(dim);
    std::copy(it, it + mp->dim(), p0.begin());

    auto &box_tree = mp->box_tree;
    rtree::pbox_set box_list;
    box_tree.find_boxes_at_point(p0, box_list);
    while (box_list.size())
    {
      auto pmax_box = *box_list.begin();

      if (box_list.size() > 1)
      {
        auto max_rate = -1.0;
        for (auto &&pbox : box_list)
        {
          auto box_rate = 1.0;
          for (size_type i = 0; i < dim; ++i)
          {
            auto h = pbox->max->at(i) - pbox->min->at(i);
            if (h > 0)
            {
              auto rate = std::min(pbox->max->at(i) - p0[i], p0[i] - pbox->min->at(i)) / h;
              box_rate = std::min(rate, box_rate);
            }
          }
          if (box_rate > max_rate)
          {
            pmax_box = pbox;
            max_rate = box_rate;
          }
        }
      }

      for (auto cv : mp->box_to_convexes_map.at(pmax_box->id)) {
        mult_diff_transposed(mp->gtrans[cv], it, mp->orgs[cv], p0);
        if (mp->trans_of_convex(cv)->convex_ref()->is_in(p0) < 1E-10) {
          return to_scalar(poly_of_subelt(cv).eval(local_coordinate ? p0.begin() : it));
        }
      }
      if (box_list.size() == 1) break;
      box_list.erase(pmax_box);
    }
    GMM_ASSERT1(
      false, "Element not found in composite polynomial: " << base_node(*it, *it + mp->dim()));
  }

  void structured_mesh_for_convex(pconvex_ref cvr, short_type k,
				  pbasic_mesh &pm, pmesh_precomposite &pmp,
				  bool force_simplexification=false);

  /** simplexify a convex_ref.
      @param cvr the convex_ref.
      @param k the refinement level.
      @return a pointer to a statically allocated mesh. Do no free it!
  */
  const basic_mesh *
  refined_simplex_mesh_for_convex(pconvex_ref cvr, short_type k);

  /** simplexify the faces of a convex_ref

      @param cvr the convex_ref.

      @param k the refinement level.

      @return vector of pointers to a statically allocated
      mesh_structure objects. Do no free them! The point numbers in
      the mesh_structure refer to the points of the mesh given by
      refined_simplex_mesh_for_convex.
  */      
  const std::vector<std::unique_ptr<mesh_structure>>&
  refined_simplex_mesh_for_convex_faces(pconvex_ref cvr, short_type k);
}  /* end of namespace bgeot.                                            */


#endif
