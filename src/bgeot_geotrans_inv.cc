/*===========================================================================

 Copyright (C) 2000-2026 Yves Renard

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
 along  with  this program. If not, see https://www.gnu.org/licenses/.

===========================================================================*/

#include "getfem/bgeot_geotrans_inv.h"
#include "getfem/bgeot_mesh_structure.h"
#include "getfem/bgeot_torus.h"
#include "gmm/gmm_solver_bfgs.h"


namespace bgeot
{
  void project_into_convex(base_node &xref, const pgeometric_trans pgt) {

    auto pgt_torus = std::dynamic_pointer_cast<const torus_geom_trans>(pgt);
    const pgeometric_trans
      orig_pgt = pgt_torus ? pgt_torus->get_original_transformation()
                           : pgt;
    pconvex_structure cvs = orig_pgt->basic_structure();
    dim_type n = cvs->dim();

    if (cvs == pyramid_QK_structure(1)) {
      // reference pyramid: 0 <= x2 <= 1 and |x0|, |x1| <= 1 - x2
      xref[2] = std::min(std::max(xref[2], 0.0), 1.0);
      scalar_type h = 1.0 - xref[2];
      for (dim_type k = 0; k < 2; ++k)
        xref[k] = std::min(std::max(xref[k], -h), h);
      return;
    }

    for (auto &coord : xref) {
      if (coord < 0.0) coord = 0.0;
      if (coord > 1.0) coord = 1.0;
    }

    if (cvs == simplex_structure(n)) {
      scalar_type sum_coordinates = 0.0;
      for (const auto &coord : xref) sum_coordinates += coord;

      if (sum_coordinates > 1.0) gmm::scale(xref, 1.0 / sum_coordinates);
    }
    else if (n >= 3 && cvs == prism_P1_structure(n)) {
      // simplex in the first n-1 coordinates
      scalar_type sum_coordinates = 0.0;
      for (dim_type k = 0; k < n-1; ++k) sum_coordinates += xref[k];
      if (sum_coordinates > 1.0)
        for (dim_type k = 0; k < n-1; ++k) xref[k] /= sum_coordinates;
    }
  }

  bool geotrans_inv_convex::invert(const base_node& x, base_node& xref,
                                   scalar_type IN_EPS,
                                   bool project_into_element) {
    bool converged = true;
    return invert(x, xref, converged, IN_EPS, project_into_element);
  }

  bool geotrans_inv_convex::invert(const base_node& x, base_node& xref,
                                   bool &converged,
                                   scalar_type IN_EPS,
                                   bool project_into_element) {
    assert(pgt);
    xref.resize(pgt->structure()->dim());
    converged = true;
    if (pgt->is_linear())
      return invert_lin(x, xref, IN_EPS);
    else
      return invert_nonlin(x, xref, IN_EPS, converged, false,
                           project_into_element);
  }

  /* inversion for linear geometric transformations */
  bool geotrans_inv_convex::invert_lin(const base_node& x, base_node& xref,
                                       scalar_type IN_EPS) {
    // the point relative to the first node, as Grel: the residual does not
    // depend on the element position
    base_node xrel(x); gmm::add(gmm::scaled(mat_col(G, 0), -1.0), xrel);
    mult(transposed(B), xrel, xref);
    base_node y = pgt->transform(xref, Grel);
    add(gmm::scaled(xrel, -1.0), y);

    return (pgt->convex_ref()->is_in(xref) < IN_EPS) &&
           (gmm::vect_norm2(y) < IN_EPS);
  }

  bool geotrans_inv_convex::update_B() {
    if (P != N) {
      pgt->compute_K_matrix(G, pc, K);
      gmm::mult(gmm::transposed(K), K, CS);
      double det = bgeot::lu_inverse(&(*(CS.begin())), P, false);
      if (det == 0)
        return false;
      gmm::mult(K, CS, B);
    } else {
      // L'inversion peut être optimisée par le non calcul global de B
      // et la resolution d'un système linéaire.
      base_matrix KT(K.nrows(), K.ncols());
      pgt->compute_K_matrix(G, pc, KT);
      gmm::copy(gmm::transposed(KT), K);
      gmm::copy(K, B);
      double det = bgeot::lu_inverse(&(*(K.begin())), P, false);
      if (det == 0)
        return false;
      B.swap(K);
    }
    return true;
  }

  class geotrans_inv_convex_bfgs {
    geotrans_inv_convex &gic;
    base_node x;
  public:
    geotrans_inv_convex_bfgs(geotrans_inv_convex &gic_,
                             const base_node &x_) : gic(gic_), x(x_) {}
    scalar_type operator()(const base_node& xref) const {
      base_node r = gic.pgt->transform(xref, gic.G) - x;
      return gmm::vect_norm2_sqr(r)/2.;
    }
    void operator()(const base_node& xref, base_small_vector& gr) const {
      gic.pgt->poly_vector_grad(xref, gic.pc);
      gic.update_B();
      base_node r = gic.pgt->transform(xref, gic.G) - x;
      gr.resize(xref.size());
      gmm::mult(gmm::transposed(gic.K), r, gr);
    }
  };

  void geotrans_inv_convex::update_linearization() {

    const convex_ind_ct &dir_pt_ind = pgt->structure()->ind_dir_points();
    const stored_point_tab &ref_nodes = pgt->geometric_nodes();

    has_linearized_approx = true;

    auto n_points = dir_pt_ind.size();
    auto N_ref = ref_nodes.begin()->size();

    std::vector<base_node> dir_pts, dir_pts_ref;
    for (auto i : dir_pt_ind) {
      dir_pts.push_back(base_node(N));
      gmm::copy(mat_col(G, i), dir_pts.back());
      dir_pts_ref.push_back(ref_nodes[i]);
    }

    base_matrix K_lin(N, n_points - 1),
                B_transp_lin(n_points - 1, N),
                K_ref_lin(N_ref, n_points - 1);

    P_lin = dir_pts[0];
    P_ref_lin = dir_pts_ref[0];

    for (size_type i = 1; i < n_points; ++i) {
      add(dir_pts[i], gmm::scaled(P_lin, -1.0), mat_col(K_lin, i - 1));
      add(dir_pts_ref[i], gmm::scaled(P_ref_lin, -1.0),
          mat_col(K_ref_lin, i - 1));
    }

    if (K_lin.nrows() == K_lin.ncols()) {
      lu_inverse(K_lin);
      gmm::copy(K_lin, B_transp_lin);
    }
    else {
      base_matrix temp(n_points - 1, n_points - 1);
      mult(transposed(K_lin), K_lin, temp);
      lu_inverse(temp);
      mult(temp, transposed(K_lin), B_transp_lin);
    }

    K_ref_B_transp_lin.base_resize(N_ref, N);
    mult(K_ref_lin, B_transp_lin, K_ref_B_transp_lin);
  }


  /* inversion for non-linear geometric transformations
     (Newton on Grad(pgt)(y - pgt(x)) = 0 )
  */
  bool geotrans_inv_convex::invert_nonlin(const base_node& x,
                                          base_node& xref,
                                          scalar_type IN_EPS,
                                          bool &converged,
                                          bool /* throw_except */,
                                          bool project_into_element) {
    converged = false;
    base_node x0_ref(P), x0rel(N), diff(N), xref_prev(P);
    // the point relative to the first node, as Grel
    base_node xrel(x); gmm::add(gmm::scaled(mat_col(G, 0), -1.0), xrel);
    scalar_type maxnormG = gmm::mat_maxnorm(Grel); // element size
    { // find initial guess
      x0_ref = pgt->geometric_nodes()[0];
      scalar_type res = gmm::vect_dist2(mat_col(Grel, 0), xrel);
      for (size_type j = 1; j < pgt->nb_points(); ++j) {
        scalar_type res0 = gmm::vect_dist2(mat_col(Grel, j), xrel);
        if (res0 < res) {
          res = res0;
          x0_ref = pgt->geometric_nodes()[j];
        }
      }
      if (res < maxnormG * IN_EPS/100.) {
        gmm::copy(x0_ref, xref);
        if (project_into_element) project_into_convex(xref, pgt);
        x0rel = pgt->transform(xref, Grel);
        add(x0rel, gmm::scaled(xrel, -1.0), diff);
        res = gmm::vect_norm2(diff);
        if (res < maxnormG * IN_EPS/100. &&
            pgt->convex_ref()->is_in(xref) < IN_EPS) {
          converged = true;
          return true;
        }
      }

      scalar_type res0 = std::numeric_limits<scalar_type>::max();
      if (has_linearized_approx) {

        add(x, gmm::scaled(P_lin, -1.0), diff);
        mult(K_ref_B_transp_lin, diff, xref);
        gmm::add(P_ref_lin, xref);

        if (project_into_element) project_into_convex(xref, pgt);
        res0 = gmm::vect_dist2(pgt->transform(xref, Grel), xrel);
      }

      if (res < res0) gmm::copy(x0_ref, xref);
      if (res < IN_EPS)      // TODO: fix this too intrusive hack
        xref *= 0.999888783; // For pyramid element to avoid the singularity
    }
    
    add(pgt->transform(xref, Grel), gmm::scaled(xrel, -1.0), diff);
    scalar_type res = gmm::vect_norm2(diff);
    scalar_type res0 = std::numeric_limits<scalar_type>::max();
    scalar_type factor = 1.0;

    while (res >= maxnormG * IN_EPS/100.) {
      if ((gmm::abs(res - res0) < IN_EPS/100.) || (factor < IN_EPS)) { // stalled?
        // relaxed convergence criterion depending on the size and position
        // of the real element
        converged = (res < maxnormG * IN_EPS/100.);
        if (!converged) { // stalled at rounding level, e.g. near a pyramid
          // apex: accept if the Newton correction is below the tolerance
          pgt->poly_vector_grad(xref, pc);
          if (update_B()) {
            mult(transposed(B), diff, x0_ref);
            converged = (gmm::vect_norm2(x0_ref) < IN_EPS);
            if (converged && P != N) { // the residual normal to the element
              base_node rn(N);         // is not reduced by the correction
              gmm::mult(K, gmm::scaled(x0_ref, -1.0), diff, rn);
              converged = (gmm::vect_norm2(rn) < maxnormG * IN_EPS/100.);
            }
          }
        }
        return converged && (pgt->convex_ref()->is_in(xref) < IN_EPS);
      }
      if (res > res0) { // back to the previous iterate, with a shorter step
        gmm::copy(xref_prev, xref);
        x0rel = pgt->transform(xref, Grel);
        add(x0rel, gmm::scaled(xrel, -1.0), diff);
        factor *= 0.5;
      } else {
        if (factor < 1.0-IN_EPS) factor *= 2.0;
        res0 = res;
      }
      pgt->poly_vector_grad(xref, pc);
      if (!update_B()) {
        converged = false;
        return false;
      }
      mult(transposed(B), diff, x0_ref);
      gmm::copy(xref, xref_prev);
      add(gmm::scaled(x0_ref, -factor), xref);
      if (project_into_element) project_into_convex(xref, pgt);
      x0rel = pgt->transform(xref, Grel);
      add(x0rel, gmm::scaled(xrel, -1.0), diff);
      res = gmm::vect_norm2(diff);
    }
    converged = true;
    return (pgt->convex_ref()->is_in(xref) < IN_EPS);
  }

}  /* end of namespace bgeot.                                             */
