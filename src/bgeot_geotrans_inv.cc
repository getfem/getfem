/*===========================================================================

 Copyright (C) 2000-2017 Yves Renard

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

===========================================================================*/

#include "getfem/bgeot_geotrans_inv.h"
#include "getfem/bgeot_mesh_structure.h"
#include "getfem/bgeot_torus.h"
#include "gmm/gmm_solver_bfgs.h"

using namespace gmm;
using namespace std;

namespace bgeot
{ 
  bool geotrans_inv_convex::invert(const base_node& n, base_node& n_ref,
				   scalar_type IN_EPS) {
    assert(pgt);
    n_ref.resize(pgt->structure()->dim());
    bool converged = true;
    if (pgt->is_linear()) {
      return invert_lin(n, n_ref,IN_EPS);
    } else {
      return invert_nonlin(n, n_ref,IN_EPS,converged,true);
    }
  }

  bool geotrans_inv_convex::invert(const base_node& n, base_node& n_ref, 
				   bool &converged, 
				   scalar_type IN_EPS) {
    assert(pgt);
    n_ref.resize(pgt->structure()->dim());
    converged = true;
    if (pgt->is_linear()) {
      return invert_lin(n, n_ref,IN_EPS);
    } else return invert_nonlin(n, n_ref,IN_EPS,converged, false);
  }

  bool point_in_convex(const geometric_trans &geoTrans,
		       const base_node &x,
		       scalar_type res,
		       scalar_type IN_EPS) {
    // Test un peu sevère peut-être en ce qui concerne res.
    return (geoTrans.convex_ref()->is_in(x) < IN_EPS) && (res < IN_EPS);
  }

  /* inversion for linear geometric transformations */
  bool geotrans_inv_convex::invert_lin(const base_node& n, base_node& n_ref,
				       scalar_type IN_EPS) {
    base_node y(n); for (size_type i=0; i < N; ++i) y[i] -= G(i,0);
    mult(transposed(B), y, n_ref);
	y = pgt->transform(n_ref, G);
    add(scaled(n, -1.0), y);

    return point_in_convex(*pgt, n_ref, vect_norm2(y), IN_EPS);
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

  nonlinear_storage_struct::linearised_structure::linearised_structure
  (const convex_ind_ct &direct_points_indices,
   const stored_point_tab &reference_nodes,
   const std::vector<base_node> &real_nodes) :
    diff(real_nodes.begin()->size()), diff_ref(direct_points_indices.size()-1) {
    auto n_points = direct_points_indices.size();
    std::vector<base_node> direct_points;
    std::vector<base_node> direct_points_ref;
    direct_points.reserve(n_points);
    direct_points_ref.reserve(n_points);
    
    for (auto i : direct_points_indices) {
      direct_points.push_back(real_nodes[i]);
      direct_points_ref.push_back(reference_nodes[i]);
    }
    
    auto N = direct_points.begin()->size();
    auto N_ref = direct_points_ref.begin()->size();
    base_matrix K_linear(N, n_points - 1);
    B_linear.base_resize(N, n_points - 1);
    K_ref_linear.base_resize(N_ref, n_points - 1);
    P_linear = direct_points[0];
    P_ref_linear = direct_points_ref[0];
    
    for (size_type i = 1; i < n_points; ++i) {
      add(direct_points[i], scaled(P_linear, -1.0), mat_col(K_linear, i - 1));
      add(direct_points_ref[i], scaled(P_ref_linear, -1.0),
	  mat_col(K_ref_linear, i - 1));
    }
    
    if (K_linear.nrows() == K_linear.ncols()) {
      lu_inverse(K_linear);
      copy(transposed(K_linear), B_linear);
    }
    else {
      base_matrix temp(n_points - 1, n_points - 1);
      mult(transposed(K_linear), K_linear, temp);
      lu_inverse(temp);
      mult(K_linear, temp, B_linear);
    }
  }

  void nonlinear_storage_struct::linearised_structure::invert
  (const base_node &x_real, base_node &x_ref, scalar_type /* IN_EPS */) const {
    add(x_real, scaled(P_linear, -1.0), diff);
    mult(transposed(B_linear), diff, diff_ref);
    mult(K_ref_linear, diff_ref, x_ref);
    add(P_ref_linear, x_ref);
  }

  void project_into_convex(base_node &x, const geometric_trans *pgeo_trans,
			   bool project) {
    if (project == false) return;
    
    auto dim = pgeo_trans->dim();
    
    for (auto d = 0; d < dim; ++d) {
      if (x[d] < 0.0) x[d] = 0.0;
      if (x[d] > 1.0) x[d] = 1.0;
    }

    auto poriginal_trans = pgeo_trans;
    
    if (auto ptorus_trans = dynamic_cast<const torus_geom_trans*>(pgeo_trans)) {
      poriginal_trans = ptorus_trans->get_original_transformation().get();
    }
    
    auto pbasic_convex_ref = basic_convex_ref(poriginal_trans->convex_ref());
    auto nb_simplices = pbasic_convex_ref->simplexified_convex()->nb_convex();
    
    if (nb_simplices == 1) { // simplex
      auto sum_coordinates = 0.0;
      
      for (auto d = 0; d < dim; ++d) sum_coordinates += x[d];
      
      if (sum_coordinates > 1.0) scale(x, 1.0 / sum_coordinates);
    }
    else if ((dim == 3) && (nb_simplices != 4)) { // prism
      auto sum_coordinates = x[0] + x[1];
      
      if (sum_coordinates > 1.0) {
	x[0] /= sum_coordinates;
	x[1] /= sum_coordinates;
      }
    }
  }
  
  void find_initial_guess(base_node &x,
			  nonlinear_storage_struct &storage,
			  const base_node &xreal,
			  const base_matrix &G,
			  const geometric_trans *pgt,
			  scalar_type IN_EPS) {
    storage.x_ref = pgt->geometric_nodes()[0];
    
    auto res = vect_dist2(mat_col(G, 0), xreal);
    double res0;
    
    for (size_type j = 1; j < pgt->nb_points(); ++j) { 
      res0 = vect_dist2(mat_col(G, j), xreal);
      if (res > res0) {
	res = res0;
	storage.x_ref = pgt->geometric_nodes()[j];
      }
    }
    
    res0 = std::numeric_limits<scalar_type>::max();
    
    if (storage.plinearised_structure != nullptr) {
      storage.plinearised_structure->invert(xreal, x, IN_EPS);
      project_into_convex(x, pgt, storage.project_into_element);
      res0 = vect_dist2(pgt->transform(x, G), xreal);
    }
    
    if (res < res0) copy(storage.x_ref, x);
  }
  

  /* inversion for non-linear geometric transformations 
     (Newton on Grad(pgt)(y - pgt(x)) = 0 )
  */
  bool geotrans_inv_convex::invert_nonlin(const base_node& xreal,
	       			  base_node& x, scalar_type IN_EPS,
				  bool &converged, bool /* throw_except */) {
    converged = true;
    find_initial_guess(x, nonlinear_storage, xreal, G, pgt.get(), IN_EPS);
    add(pgt->transform(x, G), scaled(xreal, -1.0), nonlinear_storage.diff);
    auto res = vect_norm2(nonlinear_storage.diff);
    auto res0 = std::numeric_limits<scalar_type>::max();
    double factor = 1.0;
    auto cnt = 0;

    while (res > IN_EPS) {
      if ((abs(res - res0) < IN_EPS) || (factor < IN_EPS)) {
        converged = false;
        return point_in_convex(*pgt, x, res, IN_EPS);
      }

      if (res > res0) {
        add(scaled(nonlinear_storage.x_ref, factor), x);
        nonlinear_storage.x_real = pgt->transform(x, G);
        add(nonlinear_storage.x_real, scaled(xreal, -1.0),
	    nonlinear_storage.diff);
        factor *= 0.5;
      }
      else {
        if (factor < 1.0) factor *= 2.0;
        res0 = res;
      }

      pgt->poly_vector_grad(x, pc);
      update_B();
      mult(transposed(B), nonlinear_storage.diff, nonlinear_storage.x_ref);
      add(scaled(nonlinear_storage.x_ref, -1.0 * factor), x);
      project_into_convex(x, pgt.get(), nonlinear_storage.project_into_element);
      nonlinear_storage.x_real = pgt->transform(x, G);
      add(nonlinear_storage.x_real, scaled(xreal, -1.0),
	  nonlinear_storage.diff);
      res = vect_norm2(nonlinear_storage.diff);
      ++cnt;
    }

    return point_in_convex(*pgt, x, res, IN_EPS);
  }

}  /* end of namespace bgeot.                                             */
