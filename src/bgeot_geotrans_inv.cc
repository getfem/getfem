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


  /* inversion for linear geometric transformations */
  bool geotrans_inv_convex::invert_lin(const base_node& n, base_node& n_ref, scalar_type IN_EPS) {
    base_node y(n); for (size_type i=0; i < N; ++i) y[i] -= G(i,0);
    gmm::mult(gmm::transposed(B), y, n_ref);
    if (pgt->convex_ref()->is_in(n_ref) < IN_EPS) {
      if (P == N) return true;
      else {
	gmm::mult(K,gmm::scaled(n_ref,-1.0),y,y);
	//        y -= K * n_ref;
        if (gmm::vect_norm2(y) < IN_EPS) return true;
      }
    }
    return false;
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

string element_type_of_pgt(const string &pgt_name) {
  return pgt_name.substr(3, pgt_name.find("(") - 3);
}

pair<string, string> get_pgt_names(const string &pgt_name) {
  auto pos_start = pgt_name.find("(") + 1;
  auto pos_end = pgt_name.find(")") + 1;
  auto pgt_name_1 = pgt_name.substr(pos_start, pos_end - pos_start);

  pos_start = pos_end + 1;
  pos_end = pgt_name.find(")", pos_start) + 1;
  auto pgt_name_2 = pgt_name.substr(pos_start, pos_end - pos_start);

  return {pgt_name_1, pgt_name_2};
}

dim_type get_dim_of_pgt(const string &pgt_name)
{
  auto element_type = element_type_of_pgt(pgt_name);
  if ((element_type == "PK") || (element_type == "QK") || (element_type == "PRISM")
      || (element_type == "Q2_INCOMPLETE")) {
    auto pos_start = pgt_name.find("(") + 1;
    auto end_symbol = (element_type == "Q2_INCOMPLETE") ? ")" : ",";
    return stoi(pgt_name.substr(pos_start, pgt_name.find(end_symbol) - pos_start));
  }
  else if (element_type == "PRODUCT") {
    auto pgt_names = get_pgt_names(pgt_name);
    return get_dim_of_pgt(pgt_names.first) + get_dim_of_pgt(pgt_names.second);
  }
  else
  {
    GMM_THROW_DEFAULT("Could not determine type of " + pgt_name);
    return -1;
  }
}

string create_linear_pgt_name(const string &original_pgt) {
  string linear_pgt;
  auto element_type = element_type_of_pgt(original_pgt);

  if ((element_type == "PK") || (element_type == "QK") || (element_type == "PRISM"))
  {
    linear_pgt = original_pgt;

    auto start_pos = linear_pgt.find(",") + 1;
    linear_pgt.replace(start_pos, linear_pgt.find(")") - start_pos, "1");
  }
  else if (element_type == "Q2_INCOMPLETE") {
    linear_pgt = "GT_QK(" + std::to_string(get_dim_of_pgt(original_pgt)) + ",1)";
  }
  else if (element_type == "PRODUCT") {
    auto pgt_names = get_pgt_names(original_pgt);
    auto linear_pgt1 = create_linear_pgt_name(pgt_names.first);
    auto linear_pgt2 = create_linear_pgt_name(pgt_names.second);
    linear_pgt = "GT_PRODUCT(" + linear_pgt1 + "," + linear_pgt2 + ")";
  }
  else GMM_THROW_DEFAULT("Could not determine type of " + original_pgt);

  return linear_pgt;
}

pgeometric_trans create_linear_pgt(pgeometric_trans poriginal_gt) {
  auto is_axisymmetric = is_torus_geom_trans(poriginal_gt);
  auto name_original_pgt = name_of_geometric_trans(poriginal_gt);
  auto plinear_pgt = geometric_trans_descriptor(create_linear_pgt_name(name_original_pgt));

  return is_axisymmetric ? std::make_shared<const torus_geom_trans>(move(plinear_pgt))
                         : plinear_pgt;
}

base_vector get_degrees_of_pgt(const string &pgt_name, dim_type dim) {
  auto element_type = element_type_of_pgt(pgt_name);
  base_vector degrees(dim);

  if ((element_type == "PK") || (element_type == "QK") || (element_type == "PRISM")) {
    auto pos = pgt_name.find(",") + 1;
    fill(degrees, stoi(pgt_name.substr(pos, pgt_name.find(")") - pos)));
  }
  else if (element_type == "Q2_INCOMPLETE") fill(degrees, 2);
  else if (element_type == "PRODUCT") {
    auto pgt_names = get_pgt_names(pgt_name);
    auto dim1 = get_dim_of_pgt(pgt_names.first);
    auto degrees1 = sub_vector(degrees, sub_interval(0, dim1));
    copy(get_degrees_of_pgt(pgt_names.first, dim1), degrees1);
    auto dim2 = get_dim_of_pgt(pgt_names.second);
    auto degrees2 = sub_vector(degrees, sub_interval(dim1, dim2));
    copy(get_degrees_of_pgt(pgt_names.second, dim2), degrees2);
  }
  else GMM_THROW_DEFAULT("Could not determine type of " + pgt_name);

  return degrees;
}

vector<size_type> get_linear_nodes_indices(
  const string &pgt_name, const base_vector &degrees, dim_type dim)
{
  auto element_type = element_type_of_pgt(pgt_name);
  vector<size_type> indices;

  if (element_type == "PK") {
    indices.push_back(0);
    indices.push_back(degrees[0]);
    if (dim > 1) indices.push_back((degrees[0] + 1) * (degrees[1] + 2) / 2 - 1);
    if (dim > 2) {
      indices.push_back(
        (degrees[0] + 1) * (degrees[1] + 2) * ((2 * degrees[2] + 3) / 3 + 1) / 4 - 1);
    }
  }
  else if (element_type == "QK") {
    indices.push_back(0);
    indices.push_back(degrees[0]);
    if (dim > 1) {
      indices.push_back(degrees[0] * (degrees[1] + 1));
      indices.push_back((degrees[0] + 1) * (degrees[1] + 1) - 1);
    }
    if (dim > 2) {
      auto nb_nodes_plane = (degrees[0] + 1) * (degrees[1] + 1);
      indices.push_back(nb_nodes_plane * degrees[2]);
      indices.push_back(nb_nodes_plane * degrees[2] + degrees[0]);
      indices.push_back(nb_nodes_plane * degrees[2] + degrees[0] * (degrees[1] + 1));
      indices.push_back(nb_nodes_plane * (degrees[2] + 1) - 1);
    }
  }
  else if (element_type == "PRISM") {
    indices.push_back(0);
    indices.push_back(degrees[0]);
    indices.push_back((degrees[0] + 1) * (degrees[1] + 2) / 2 - 1);
    indices.push_back((degrees[0] + 1) * (degrees[1] + 2) / 2);
    indices.push_back((degrees[0] + 1) * (degrees[1] + 2) * (degrees[2] + 1) / 2 - 1);
  }
  else if (element_type == "Q2_INCOMPLETE") {
    indices.push_back(0);
    indices.push_back(2);
    if (dim > 1) {
      indices.push_back(5);
      indices.push_back(7);
    }
    if (dim > 2) {
      indices.push_back(12);
      indices.push_back(14);
      indices.push_back(17);
      indices.push_back(19);
    }
  }
  else if (element_type == "PRODUCT") {
    auto pgt_name1 = get_pgt_names(pgt_name).first;
    auto dim1 = get_dim_of_pgt(pgt_name1);
    auto indices_plane = get_linear_nodes_indices(
                           pgt_name1, get_degrees_of_pgt(pgt_name1, dim1), dim1);

    for (auto i : indices_plane) indices.push_back(i);

    for (auto d = dim1; d < dim; ++d) {
      auto indices_old = indices;

      for (auto i : indices_old) indices.push_back(degrees[d] * indices_old.size() + i);
    }
  }

  return indices;
}

void project_into_convex(base_node &x, pgeometric_trans &pgeo_trans) {
  auto dim = pgeo_trans->dim();

  for (auto d = 0; d < dim; ++d) {
    if (x[d] < 0.0) x[d] = 0.0;
    if (x[d] > 1.0) x[d] = 1.0;
  }

  auto poriginal_trans = pgeo_trans.get();

  if (auto ptorus_trans = dynamic_cast<const torus_geom_trans*>(pgeo_trans.get())) {
    poriginal_trans = ptorus_trans->get_original_transformation().get();
  }

  auto nb_simplices = poriginal_trans->convex_ref()->simplexified_convex()->nb_convex();

  if (nb_simplices == 1) { //simplex
    auto sum_coordinates = 0.0;

    for (auto d = 0; d < dim; ++d) sum_coordinates += x[d];

    if (sum_coordinates > 1.0) scale(x, 1.0 / sum_coordinates);
  }
  else if ((dim == 3) && (nb_simplices != 4)) { //prism
    auto sum_coordinates = x[0] + x[1];

    if (sum_coordinates > 1.0) {
      x[0] /= sum_coordinates;
      x[1] /= sum_coordinates;
    }
  }
}

vector<size_type> get_linear_nodes_indices(pgeometric_trans pgt) {
  auto pgt_name = name_of_geometric_trans(pgt);
  auto original_dim = is_torus_geom_trans(pgt) ? 2 : pgt->dim();
  return get_linear_nodes_indices(
    pgt_name, get_degrees_of_pgt(pgt_name, pgt->dim()), original_dim);
}

  /* inversion for non-linear geometric transformations 
     (Newton on Grad(pgt)(y - pgt(x)) = 0 )
  */
  bool geotrans_inv_convex::invert_nonlin(const base_node& xreal,
	       			  base_node& x, scalar_type IN_EPS,
				  bool &converged, bool throw_except) {
    converged = true;
    /* find an initial guess */
    nonlinear_storage.plinear_inversion->invert(xreal, x);
    project_into_convex(x, nonlinear_storage.plinear_inversion->pgt);
    nonlinear_storage.x_real = pgt->transform(x, G);
    add(nonlinear_storage.x_real, scaled(xreal, -1.0), nonlinear_storage.diff);

    auto res = vect_norm2(nonlinear_storage.diff);
    double res0 = 1e10;
    auto cnt = 0;

    while (res > IN_EPS) {
      if (abs(res - res0) < IN_EPS) {
        converged = false;
        return false;
      }

      pgt->poly_vector_grad(x, pc);
      update_B();
      mult(transposed(B), nonlinear_storage.diff, nonlinear_storage.x_ref);
      add(scaled(nonlinear_storage.x_ref, -1.0), x);
      project_into_convex(x, nonlinear_storage.plinear_inversion->pgt);
      nonlinear_storage.x_real = pgt->transform(x, G);
      add(nonlinear_storage.x_real, scaled(xreal, -1.0), nonlinear_storage.diff);
      res0 = res;
      res = vect_norm2(nonlinear_storage.diff);
      ++cnt;
    }

    return true;
  }

}  /* end of namespace bgeot.                                             */
