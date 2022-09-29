/*===========================================================================

 Copyright (C) 2000-2020 Yves Renard

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

#include "getfem/bgeot_geotrans_inv.h"
#include "getfem/bgeot_mesh_structure.h"
#include "getfem/bgeot_torus.h"
#include "gmm/gmm_solver_bfgs.h"


namespace bgeot
{
  void project_into_convex(base_node &x, const pgeometric_trans pgt) {

    for (auto &coord : x) {
      if (coord < 0.0) coord = 0.0;
      if (coord > 1.0) coord = 1.0;
    }


    auto pgt_torus = std::dynamic_pointer_cast<const torus_geom_trans>(pgt);
    const pgeometric_trans
      orig_pgt = pgt_torus ? pgt_torus->get_original_transformation()
                           : pgt;

    auto pbasic_convex_ref = basic_convex_ref(orig_pgt->convex_ref());
    auto nb_simplices = pbasic_convex_ref->simplexified_convex()->nb_convex();

    if (nb_simplices == 1) { // simplex
      auto sum_coordinates = 0.0;
      for (const auto &coord : x) sum_coordinates += coord;

      if (sum_coordinates > 1.0) gmm::scale(x, 1.0 / sum_coordinates);
    }
    else if (pgt->dim() == 3 && nb_simplices != 4) { // prism
      auto sum_coordinates = x[0] + x[1];
      if (sum_coordinates > 1.0) {
        x[0] /= sum_coordinates;
        x[1] /= sum_coordinates;
      }
    }
  }

  bool geotrans_inv_convex::invert(const base_node& n, base_node& n_ref,
                                   scalar_type IN_EPS,
                                   bool project_into_element) {
    bool converged = true;
    return invert(n, n_ref, converged, IN_EPS, project_into_element);
  }

  bool geotrans_inv_convex::invert(const base_node& n, base_node& n_ref,
                                   bool &converged,
                                   scalar_type IN_EPS,
                                   bool project_into_element) {
    assert(pgt);
    n_ref.resize(pgt->structure()->dim());
    converged = true;
    if (pgt->is_linear())
      return invert_lin(n, n_ref, IN_EPS);
    else
      return invert_nonlin(n, n_ref, IN_EPS, converged, false,
                           project_into_element);
  }

  /* inversion for linear geometric transformations */
  bool geotrans_inv_convex::invert_lin(const base_node& n, base_node& n_ref,
                                       scalar_type IN_EPS) {
    base_node y(n); for (size_type i=0; i < N; ++i) y[i] -= G(i,0);
    mult(transposed(B), y, n_ref);
    y = pgt->transform(n_ref, G);
    add(gmm::scaled(n, -1.0), y);

    return (pgt->convex_ref()->is_in(n_ref) < IN_EPS) &&
           (gmm::vect_norm2(y) < IN_EPS);
  }

  void geotrans_inv_convex::update_B() {
    if (P != N) {
      pgt->compute_K_matrix(G, pc, K);
      gmm::mult(gmm::transposed(K), K, CS);
      bgeot::lu_inverse(&(*(CS.begin())), P);
      gmm::mult(K, CS, B);
    }
    else {
      // L'inversion peut être optimisée par le non calcul global de B
      // et la resolution d'un système linéaire.
      base_matrix KT(K.nrows(), K.ncols());
      pgt->compute_K_matrix(G, pc, KT);
      gmm::copy(gmm::transposed(KT), K);
      gmm::copy(K, B);
      bgeot::lu_inverse(&(*(K.begin())), P); B.swap(K);
    }
  }

  class geotrans_inv_convex_bfgs {
    geotrans_inv_convex &gic;
    base_node xreal;
  public:
    geotrans_inv_convex_bfgs(geotrans_inv_convex &gic_,
                             const base_node &xr) : gic(gic_), xreal(xr) {}
    scalar_type operator()(const base_node& x) const {
      base_node r = gic.pgt->transform(x, gic.G) - xreal;
      return gmm::vect_norm2_sqr(r)/2.;
    }
    void operator()(const base_node& x, base_small_vector& gr) const {
      gic.pgt->poly_vector_grad(x, gic.pc);
      gic.update_B();
      base_node r = gic.pgt->transform(x, gic.G) - xreal;
      gr.resize(x.size());
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
  bool geotrans_inv_convex::invert_nonlin(const base_node& xreal,
                                          base_node& x, scalar_type IN_EPS,
                                          bool &converged,
                                          bool /* throw_except */,
                                          bool project_into_element) {
    converged = true;
    base_node x0_ref(P), diff(N);

    { // find initial guess
      x0_ref = pgt->geometric_nodes()[0];
      scalar_type res = gmm::vect_dist2(mat_col(G, 0), xreal);
      for (size_type j = 1; j < pgt->nb_points(); ++j) {
        scalar_type res0 = gmm::vect_dist2(mat_col(G, j), xreal);
        if (res0 < res) {
          res = res0;
          x0_ref = pgt->geometric_nodes()[j];
        }
      }

      scalar_type res0 = std::numeric_limits<scalar_type>::max();
      if (has_linearized_approx) {

        add(xreal, gmm::scaled(P_lin, -1.0), diff);
        mult(K_ref_B_transp_lin, diff, x);
        gmm::add(P_ref_lin, x);

        if (project_into_element) project_into_convex(x, pgt);
        res0 = gmm::vect_dist2(pgt->transform(x, G), xreal);
      }

      if (res < res0) gmm::copy(x0_ref, x);
      if (res < IN_EPS)
        x *= 0.999888783; // For pyramid element to avoid the singularity
    }
    
    add(pgt->transform(x, G), gmm::scaled(xreal, -1.0), diff);
    scalar_type res = gmm::vect_norm2(diff);
    scalar_type res0 = std::numeric_limits<scalar_type>::max();
    scalar_type factor = 1.0;

    base_node x0_real(N);
    while (res > IN_EPS/100.) {
      if ((gmm::abs(res - res0) < IN_EPS/100.) || (factor < IN_EPS)) {
        // relaxed convergence criterion depending on the size and position
        // of the real element
        converged = (res < gmm::mat_maxnorm(G) * IN_EPS/100.);
        return (pgt->convex_ref()->is_in(x) < IN_EPS) && (res < IN_EPS);
      }
      if (res > res0) {
        add(gmm::scaled(x0_ref, factor), x);
        x0_real = pgt->transform(x, G);
        add(x0_real, gmm::scaled(xreal, -1.0), diff);
        factor *= 0.5;
      }
      else {
        if (factor < 1.0-IN_EPS) factor *= 2.0;
        res0 = res;
      }
      pgt->poly_vector_grad(x, pc);
      update_B();
      mult(transposed(B), diff, x0_ref);
      add(gmm::scaled(x0_ref, -factor), x);
      if (project_into_element) project_into_convex(x, pgt);
      x0_real = pgt->transform(x, G);
      add(x0_real, gmm::scaled(xreal, -1.0), diff);
      res = gmm::vect_norm2(diff);
    }
    return (pgt->convex_ref()->is_in(x) < IN_EPS) && (res < IN_EPS);
  }

}  /* end of namespace bgeot.                                             */
